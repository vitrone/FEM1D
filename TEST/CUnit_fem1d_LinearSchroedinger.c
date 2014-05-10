/*============================================================================+/
 | File: Cunit_fem1d_LinearSchroedinger.c 
 | Description: Test 
 |
 | Boundary Condition: Dirichlet
 | Initial Conditions: Compactly supported
 |
/+============================================================================*/
#include <pthread.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>

/* MKL */ 
#include "mkl.h"
#include "mkl_pardiso.h"

#define NDEBUG
#define MATLIB_NTRACE_DATA

#include "legendre.h"
#include "fem1d.h"
#include "assert.h"

/* CUnit modules */
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>
#include <CUnit/TestDB.h>

static const matlib_real TOL = 1e-9;
/*============================================================================*/
int init_suite(void)
{
      return 0;
}

int clean_suite(void)
{
      return 0;
}

/*============================================================================*/

void fem1d_zm_sparse_GSM
/* Global Stiffness Matrix */ 
(
    matlib_index     N,
    matlib_complex   coeff,
    matlib_zm_sparse M
)
{
    debug_enter("%s", "");
    matlib_index i;

    M.elem_p[0] =  coeff/2.0 + M.elem_p[0];
    M.elem_p[1] = -coeff/2.0 + M.elem_p[1];
    for(i=1; i<N; i++)
    {
        M.elem_p[M.rowIn[i]]   =  coeff     + M.elem_p[M.rowIn[i]  ];
        M.elem_p[M.rowIn[i]+1] = -coeff/2.0 + M.elem_p[M.rowIn[i]+1];
    }
    M.elem_p[M.rowIn[i]] = coeff/2.0 + M.elem_p[M.rowIn[i]];
    
    for(i=N+1; i<M.lenc; i++)
    {
        M.elem_p[M.rowIn[i]] = coeff + M.elem_p[M.rowIn[i]];
    }
    debug_exit("%s", "");
}

void fem1d_zm_nsparse_GSM
/* Global Stiffness Matrix */ 
(
    matlib_index      N,
    matlib_complex    coeff,
    matlib_zm_nsparse M
)
{
    debug_enter("%s", "");
    matlib_complex* ptr;
    matlib_index i, k;
    for( k=0; k<M.nsparse; k++)
    {
        ptr = M.elem_p[k];
        ptr[0] =  coeff/2.0 + ptr[0];
        ptr[1] = -coeff/2.0 + ptr[1];
        for(i=1; i<N; i++)
        {
            ptr[M.rowIn[i]]   =  coeff     + ptr[M.rowIn[i]  ];
            ptr[M.rowIn[i]+1] = -coeff/2.0 + ptr[M.rowIn[i]+1];
        }
        ptr[M.rowIn[i]] = coeff/2.0 + ptr[M.rowIn[i]];
        
        for(i=N+1; i<M.lenc; i++)
        {
            ptr[M.rowIn[i]] = coeff + ptr[M.rowIn[i]];
        }
    }
    debug_exit("%s", "");
}

/*============================================================================+/
| Example 1.:
| Schroedinger Equation (IVP):
| iu_t + u_xx + phi(x) u= 0
| t=0: u(x,t) = u_0(x) -> Initial data
| This example deal with time-independent real potential.
|
| Weak formulation in time discrete form: For any test functions 'psi'
| i/(rho*J^2)(v_xi, psi_xi)+(phi_1*v,psi) = (u,psi) where u -> 2 * v - u
| phi_1 = 1 + (-i/rho) phi(x)
| rho = 2/dt
|
/+============================================================================*/
matlib_real solve_Schroedinger_IVP1
(
    matlib_index p,
    matlib_index nr_LGL,
    matlib_index N,
    matlib_real domain[2],
    matlib_real dt,
    matlib_index Nt,
    void (*func_p)(matlib_xv, matlib_real, matlib_zv),
    void (*potential_p)(matlib_complex [], matlib_xv, matlib_zv)
)
{

    debug_enter( "polynomial degree: %d, nr. of LGL points: %d", p, nr_LGL );
    matlib_index i, j;
    matlib_index P = nr_LGL-1;

    matlib_real x_l = domain[0];
    matlib_real x_r = domain[1];

    matlib_xv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_xm FM, IM, Q;
    matlib_create_xm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_xm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    legendre_LGLdataIM( xi, IM);
    fem1d_quadM( quadW, IM, &Q);

    /* generate the grid */ 
    matlib_xv x;
    fem1d_ref2mesh (xi, N, x_l, x_r, &x);
    debug_body("length of x: %d", x.len);

    matlib_real ta[2] = {0, Nt*dt};

    /* Provide the initial condition */
    matlib_zv u0, u_exact, u_computed;
    matlib_create_zv( x.len,         &u0, MATLIB_COL_VECT);
    matlib_create_zv( x.len,    &u_exact, MATLIB_COL_VECT);
    matlib_create_zv( x.len, &u_computed, MATLIB_COL_VECT);
    
    (*func_p)(x, ta[0], u0);
    DEBUG_PRINT_ZV(u0, "%s :", "Initial condition");
    
    /* real potential */ 
    matlib_zv phi;
    matlib_create_zv( x.len, &phi, MATLIB_COL_VECT);

    matlib_real irho = dt/2.0;
    matlib_real J    = 0.5*(x_r-x_l)/N;

    matlib_complex s_coeff = 0.0 + I*irho/(J*J); /* Stiffness Matrix Coeff.*/ 
    matlib_complex m_coeff[2] = {1.0, -I*irho};  /* Mass Matrix Coeff.*/
    debug_body("coeff: %0.16f%+0.16fi", s_coeff);
    debug_body("coeff: %0.16f%+0.16fi", m_coeff);

    (*potential_p)(m_coeff, x, phi);

    matlib_zm_sparse M;
    fem1d_zm_sparse_GMM(p, Q, phi, &M);

    
    /* Global Stiffness-Mass Matrix */ 
    fem1d_zm_sparse_GSM(N, s_coeff, M);

    BEGIN_DTRACE
        debug_print("dimension of the sparse square matrix: %d", M.lenc);
        for(i=0; i<M.lenc; i++)
        {
            for(j=M.rowIn[i]; j<M.rowIn[i+1]; j++)
            {
                debug_print( "M(%d,%d): % 0.16f %+0.16fi", 
                             i, M.colIn[j], M.elem_p[j]);
            }
        }
    END_DTRACE

    /* Solve the equation while marching in time 
     *
     * Linear system M * V_vb = Pvb 
     * */
    
    matlib_zv U_tmp, V_tmp, V_vb, Pvb;
    matlib_real dim = N*(p+1);
    matlib_create_zv( M.lenc,   &Pvb, MATLIB_COL_VECT);
    matlib_create_zv( M.lenc,  &V_vb, MATLIB_COL_VECT);
    matlib_create_zv(    dim, &V_tmp, MATLIB_COL_VECT);
    matlib_create_zv(    dim, &U_tmp, MATLIB_COL_VECT);

    fem1d_ZFLT( N, FM, u0, U_tmp);

    matlib_real t = 0, norm_actual;
    matlib_xv e_relative;
    matlib_create_xv( Nt, &e_relative, MATLIB_COL_VECT);

    pardiso_solver_t data = { .nsparse  = 1, 
                              .mnum     = 1, 
                              .sol_enum = PARDISO_LHS, 
                              .mtype    = PARDISO_COMPLEX_SYM,
                              .smat_p   = (void*)&M,
                              .rhs_p    = (void*)&Pvb,
                              .sol_p    = (void*)&V_vb};

    data.phase_enum = PARDISO_INIT;
    matlib_pardiso(&data);

    data.phase_enum = PARDISO_ANALYSIS_AND_FACTOR;
    matlib_pardiso(&data);
    for (i=1; i<Nt+1; i++)
    {
        debug_body("begin iteration: %d", i);
        fem1d_ZPrjL2F(p, U_tmp, Pvb);

        data.phase_enum = PARDISO_SOLVE_AND_REFINE;
        matlib_pardiso(&data);
        
        fem1d_ZF2L(p, V_vb, V_tmp);

        /* 2.0 * V_tmp -U_tmp --> U_tmp*/ 
        matlib_zaxpby(2.0, V_tmp, -1.0, U_tmp );
        t += dt;

        (*func_p)(x, t, u_exact);
        fem1d_ZFLT(N, FM, u_exact, V_tmp);

        /* Error analysis */ 
        norm_actual = fem1d_ZNorm2(p, N, V_tmp);
        matlib_zaxpy(-1.0, U_tmp, V_tmp );
        e_relative.elem_p[i-1] = fem1d_ZNorm2(p, N, V_tmp)/norm_actual;
    }
    

    (*func_p)(x, ta[1], u_exact);
    fem1d_ZILT(N, IM, U_tmp, u_computed);

    /* Write the final evolution data to a file */ 
    matlib_xzvwrite_csv( "initial_state.dat", 1, &x, 1, &u0);

    matlib_zv tmp_mat[2] = {u_exact, u_computed }; 
    matlib_xzvwrite_csv( "final_state.dat", 1, &x, 2, tmp_mat);
    
    matlib_xvwrite_csv("e_relative.dat", 1, &e_relative);
    DEBUG_PRINT_XV(e_relative, "%s :", "relative error");

    data.phase_enum = PARDISO_FREE;
    matlib_pardiso(&data);

    matlib_free(u0.elem_p);
    matlib_free(u_exact.elem_p);
    matlib_free(u_computed.elem_p);
    matlib_free(phi.elem_p);
    matlib_free(U_tmp.elem_p);
    matlib_free(V_tmp.elem_p);
    matlib_free(V_vb.elem_p);
    matlib_free(Pvb.elem_p);

    /* Free the sparse matrix */ 
    matlib_free(M.elem_p);
    matlib_free(M.colIn);
    matlib_free(M.rowIn);

    matlib_real r = e_relative.elem_p[Nt-1];

    debug_exit("Relative error: % 0.16g", r);
    return(r);
}

/* Gaussian Wave Packet in Constant Potential */ 
/*============================================================================*/
void solution_Gaussian_WP
(
    matlib_xv   x,
    matlib_real t,
    matlib_zv   u
)
{
    debug_enter("%s", "");
    matlib_complex coeff = I/(-4.0*t+I);
    matlib_real *xptr;
    for (xptr = x.elem_p; xptr<(x.elem_p+x.len); xptr++)
    {
        *(u.elem_p) = csqrt(coeff)*cexp(-*xptr**xptr*coeff+I*t);
        u.elem_p++;
    }
    debug_exit("%s", "");
}

void constant_zpotential
( 
          matlib_complex m_coeff[2], 
    const matlib_xv      x, 
          matlib_zv      y
)
{
    matlib_index i;
    if(y.len >= x.len)
    {
        for (i=0; i<x.len; i++)
        {
            *(y.elem_p) = m_coeff[0] + m_coeff[1];
            y.elem_p++;
        }
    }
}

void test_solve_Schroedinger_IVP1a(void)
{
    matlib_index p = 4;
    matlib_index nr_LGL = 2*p+1;

    matlib_index N  = 800;
    matlib_index Nt = 1000;

    matlib_real dt = 0.5e-3;
    matlib_real domain[2] = {-15.0, 15.0 };

    matlib_real e_relative;

    struct timespec tb, te;
    matlib_real dn;
    START_TIMMING(tb);

    e_relative = solve_Schroedinger_IVP1( p, nr_LGL, N, domain, dt, Nt,
                                          solution_Gaussian_WP,
                                          constant_zpotential);
    GET_DURATION(tb, te, dn);
    debug_print("Execution Time [msec]: %0.6f", dn);
    CU_ASSERT_TRUE(e_relative<1e-6);

}

/* Hermite-Gaussian WP in Harmonic Potential */ 
/*============================================================================*/
void harmonic_zpotential
( 
          matlib_complex m_coeff[2], 
    const matlib_xv x, 
          matlib_zv y
)
{
    matlib_index i;
    if(y.len >= x.len)
    {
        for (i=0; i<x.len; i++)
        {
            *(y.elem_p) = m_coeff[0]
                          + m_coeff[1]*(-*(x.elem_p+i)**(x.elem_p+i));
            y.elem_p++;
        }
    }
}

void solution_HermiteGaussian_WP
(
    matlib_xv   x,
    matlib_real t,
    matlib_zv   u
)
/* 
 * iu_t+u_xx-x^2u=0
 * u = psi_0(x)exp(-it)+psi_1(x)exp(-i3*t)
 * psi_0 = (1/pi)^1/4 * exp(-x^2/2)
 * psi_1 = (1/pi)^1/4 * sqrt(2) x exp(-x^2/2)
 *
 *
 * */
{
    debug_enter("%s", "");
    
    matlib_real coeff[2] = {1.0/pow(M_PI, 0.25), sqrt(2.0)/pow(M_PI, 0.25)};
    matlib_complex e[2];
    matlib_real *xptr;
    for (xptr = x.elem_p; xptr<(x.elem_p+x.len); xptr++)
    {
        e[0] = cexp(-I*t); 
        e[1] = cexp(-I*3*t);
        *(u.elem_p) = (e[0]*coeff[0]+e[1]*coeff[1]**xptr)*exp(-0.5**xptr**xptr);
        u.elem_p++;
    }
    debug_exit("%s", "");
}

void test_solve_Schroedinger_IVP1b(void)
{
    matlib_index p = 4;
    matlib_index nr_LGL = 2*p+1;
    matlib_index N = 800;
    matlib_index Nt = 2000;
    matlib_real dt = 0.5e-3;
    matlib_real domain[2] = {-15.0, 15.0 };
    matlib_real e_relative;

    struct timespec tb, te;
    matlib_real dn;
    START_TIMMING(tb);

    e_relative = solve_Schroedinger_IVP1( p, nr_LGL, N, domain, dt, Nt,
                                          solution_HermiteGaussian_WP,
                                          harmonic_zpotential);
    GET_DURATION(tb, te, dn);
    debug_print("Execution Time [msec]: %0.6f", dn);
    CU_ASSERT_TRUE(e_relative<1e-6);

}

/*============================================================================+/
| Example 2.:
| Schroedinger Equation (IVP):
| iu_t + u_xx + phi(x,t) u= 0
| t=0: u(x,t) = u_0(x) -> Initial data
| This example deals with time-dependent real potential.
|   phi = x cos(mu*t).
/+============================================================================*/
void fem1d_linespace
(
    matlib_real  t0, 
    matlib_real  dt, 
    matlib_index Nt, 
    matlib_xv    *t
)
{
    matlib_create_xv( Nt+1, t, MATLIB_COL_VECT);

    t->elem_p[0] = t0;

    for(matlib_index i=1; i<Nt+1; i++)
    {
        t->elem_p[i] = t->elem_p[i-1] + dt;
    }
}

matlib_real solve_Schroedinger_IVP2
(
    matlib_index p,
    matlib_index nr_LGL,
    matlib_index N,
    matlib_real  domain[2],
    matlib_real  dt,
    matlib_index Nt,
    void* params,
    void  (*fun_p)(void*, matlib_xv, matlib_real, matlib_zv),
    void  (*pot_p)(void*, matlib_complex*, matlib_xv, matlib_xv, matlib_zm)
)
{

    debug_enter( "polynomial degree: %d, nr. of LGL points: %d", p, nr_LGL );
    matlib_index i, j;
    matlib_index P = nr_LGL-1;

    matlib_real x_l = domain[0];
    matlib_real x_r = domain[1];

    matlib_xv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_xm FM, IM, Q;
    matlib_create_xm(    p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_xm( xi.len,    p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    legendre_LGLdataIM( xi, IM);
    fem1d_quadM( quadW, IM, &Q);

    /* generate the grid */ 
    matlib_real t0 = 0;
    matlib_xv x, t;
    fem1d_ref2mesh (xi, N, x_l, x_r, &x);
    fem1d_linespace(t0, dt, Nt, &t);
    debug_body("length of x: %d", x.len);
    debug_body("length of t: %d", t.len);

    /* time discretization */ 
    matlib_zv u0, u_exact, u_computed;
    matlib_create_zv( x.len,         &u0, MATLIB_COL_VECT);
    matlib_create_zv( x.len,    &u_exact, MATLIB_COL_VECT);
    matlib_create_zv( x.len, &u_computed, MATLIB_COL_VECT);
    
    /* Provide the initial condition */
    (*fun_p)(params, x, t0, u0);
    

    matlib_real irho = dt/2.0;
    matlib_real J    = 0.5*(x_r-x_l)/N;

    matlib_complex s_coeff    = 0.0 + I*irho/(J*J);
    matlib_complex m_coeff[2] = {1.0, -I*irho};
    debug_body("s_coeff: %0.16f%+0.16fi"   , s_coeff);
    debug_body("m_coeff[0]: %0.16f%+0.16fi", m_coeff[0]);
    debug_body("m_coeff[1]: %0.16f%+0.16fi", m_coeff[1]);

    /* Assemble the global mass matrix */ 
    matlib_zm phi, q;
    matlib_zm_nsparse M;
    matlib_index nsparse = 20;
    fem1d_zm_nsparse_GMM(p, N, nsparse, Q, &phi,   &q, &M, FEM1D_GMM_INIT);
    fem1d_zm_nsparse_GMM(p, N, nsparse, Q, NULL, NULL, &M, FEM1D_GET_SPARSITY_ONLY);
    
    /* Setup the sparse linear system */
    /* Solve the equation while marching in time 
     *
     * Linear system M * V_vb = Pvb 
     * */
    
    matlib_zv U_tmp, V_tmp, V_vb, U_vb, Pvb;
    matlib_real dim = N*(p+1);
    matlib_create_zv( M.lenc,   &Pvb, MATLIB_COL_VECT);
    matlib_create_zv( M.lenc,  &V_vb, MATLIB_COL_VECT);
    matlib_create_zv(    dim, &V_tmp, MATLIB_COL_VECT);
    matlib_create_zv(    dim, &U_tmp, MATLIB_COL_VECT);

    fem1d_ZFLT( N, FM, u0, U_tmp);


    matlib_real norm_actual;
    matlib_xv e_relative;
    matlib_create_xv( Nt, &e_relative, MATLIB_COL_VECT);

    pardiso_solver_t data = { .nsparse  = nsparse, 
                              .mnum     = 1, 
                              .mtype    = PARDISO_COMPLEX_SYM,
                              .sol_enum = PARDISO_LHS, 
                              .smat_p   = (void*)&M,
                              .rhs_p    = (void*)&Pvb,
                              .sol_p    = (void*)&V_vb};

    debug_body("%s", "Solver data initialized");

    BEGIN_DTRACE
    matlib_complex *ptr;
        for(matlib_index k = 0; k< nsparse; k++)
        {
            debug_print("dimension of the sparse square matrix: %d", M.lenc);
            ptr = M.elem_p[k];
            for(i=0; i<M.lenc; i++)
            {
                for(j=M.rowIn[i]; j<M.rowIn[i+1]; j++)
                {
                    debug_print("M(%d,%d): % 0.16f %+0.16fi", i, M.colIn[j], ptr[j]);
                }
            }
        }
    END_DTRACE

    data.phase_enum = PARDISO_INIT;
    matlib_pardiso(&data);

    matlib_index Nt_ = Nt/nsparse;
    matlib_xv t_tmp  = {.len = nsparse}; 

    for (i=0; i<Nt_; i++)
    {
        t_tmp.elem_p = t.elem_p + i*nsparse;
        (*pot_p)( params, m_coeff, x, t_tmp, phi);
        fem1d_zm_nsparse_GMM(p, N, nsparse, Q, &phi, &q, &M, FEM1D_GET_NZE_ONLY);
        fem1d_zm_nsparse_GSM(N, s_coeff, M);

        for(j=0; j<nsparse; j++)
        {
            debug_body("begin iteration: %d", i*nsparse+j);
            fem1d_ZPrjL2F(p, U_tmp, Pvb);

            data.mnum = j+1;
            data.phase_enum = PARDISO_ANALYSIS_AND_FACTOR;
            matlib_pardiso(&data);

            data.phase_enum = PARDISO_SOLVE_AND_REFINE;
            matlib_pardiso(&data);
            
            fem1d_ZF2L(p, V_vb, V_tmp);

            /* 2.0 * V_tmp -U_tmp --> U_tmp*/ 
            matlib_zaxpby(2.0, V_tmp, -1.0, U_tmp );

            (*fun_p)(params, x, t_tmp.elem_p[j+1], u_exact);
            fem1d_ZFLT(N, FM, u_exact, V_tmp);

            /* Error analysis */ 
            norm_actual = fem1d_ZNorm2(p, N, V_tmp);
            matlib_zaxpy(-1.0, U_tmp, V_tmp );
            e_relative.elem_p[j+i*nsparse] = fem1d_ZNorm2(p, N, V_tmp)/norm_actual;
        }
    }

    (*fun_p)(params, x, t.elem_p[t.len-1], u_exact);
    fem1d_ZILT(N, IM, U_tmp, u_computed);

    /* Write the final evolution data to a file */ 
    matlib_xzvwrite_csv( "initial_state.dat", 1, &x, 1, &u0);

    matlib_zv tmp_mat[2] = {u_exact, u_computed }; 
    matlib_xzvwrite_csv( "final_state.dat", 1, &x, 2, tmp_mat);
    
    matlib_xvwrite_csv("e_relative.dat", 1, &e_relative);
    DEBUG_PRINT_XV(e_relative, "%s :", "relative error");


    data.phase_enum = PARDISO_FREE;
    matlib_pardiso(&data);
    fem1d_zm_nsparse_GMM(p, N, nsparse, Q, &phi, &q, &M, FEM1D_GMM_FREE);

    matlib_free(u0.elem_p);
    matlib_free(u_exact.elem_p);
    matlib_free(u_computed.elem_p);
    matlib_free(U_tmp.elem_p);
    matlib_free(V_tmp.elem_p);
    matlib_free(V_vb.elem_p);
    matlib_free(Pvb.elem_p);


    matlib_real r = e_relative.elem_p[Nt-1];

    debug_exit("Relative error: % 0.16g", r);
    return(r);
}

/* Gaussian WP  in time-dependent linear potential */
/*============================================================================*/
void solution_Gaussian_WP_linpot
(
    void*        params_,
    matlib_xv    x,
    matlib_real  t,
    matlib_zv    u
)
/* 
 * 
 * beta   = g_0 *sin(mu*t)/mu 
 * nu     = (2*g_0*cos(mu*t))/mu^2 - (2*g_0)/mu^2;
 * ibeta2 = (1/mu^3)(mu*(g_0^2*t)/2 - (g_0^2*sin(2*mu*t))/4;
 * Xi     = beta*x-ibeta2;
 *
 * */ 
{
    debug_enter("%s", "");

    matlib_real* params = (matlib_real*)params_;

    matlib_complex coeff = 1.0/(4.0*I*t+1);
    matlib_real *xptr;

    /* parameters 
     * param[0] = g_0, 
     * param[1] = mu
     * */ 
    matlib_real Xi, x1;
    
    matlib_real theta     = params[1]*t;
    matlib_real beta      = params[0]*sin(theta)/params[1];
    matlib_real nu        = 2.0*params[0]*(cos(theta)-1.0)/(params[1]*params[1]);
    matlib_real ibeta2    = params[0]*params[0]*
                            (theta/2.0 - sin(2*theta)/4.0)/(params[1]*params[1]*params[1]);

    for (xptr = x.elem_p; xptr<(x.elem_p+x.len); xptr++)
    {
        Xi = beta**xptr - ibeta2;
        x1 = *xptr + nu;
        *(u.elem_p) = csqrt(coeff)*cexp(-x1*x1*coeff+I*Xi);
        u.elem_p++;
    }
    debug_exit("%s", "");
}

void linear_timedependent_zpotential
( 
    void*          params_,
    matlib_complex m_coeff[2], 
    matlib_xv      x, 
    matlib_xv      t, 
    matlib_zm      y
)
{
    debug_enter("%s", "");
    matlib_index i, j;

    matlib_real* params = (matlib_real*) params_;

    for(j=0; j<y.lenr; j++)
    {
        for (i=0; i<x.len; i++)
        {
            *(y.elem_p) =  m_coeff[0]+
                           0.5*m_coeff[1]**(x.elem_p+i)*params[0]*
                           (cos(params[1]*(t.elem_p[j]))+cos(params[1]*(t.elem_p[j+1])));
            y.elem_p++;
        }
    }
    debug_exit("%s", "");
}

void test_solve_Schroedinger_IVP2a(void)
{
    matlib_index p = 4;
    matlib_index nr_LGL = 2*p+1;
    matlib_index N  = 400;
    matlib_index Nt = 2000;
    matlib_real dt = 0.5e-3;
    matlib_real domain[2] = {-15.0, 15.0 };
    matlib_real params[2] = {1.0, 2*M_PI};
    matlib_real e_relative;

    struct timespec tb, te;
    matlib_real dn;
    START_TIMMING(tb);


    e_relative = solve_Schroedinger_IVP2( p, nr_LGL, N, domain, dt, Nt, (void*)params,
                                          solution_Gaussian_WP_linpot,
                                          linear_timedependent_zpotential);
    GET_DURATION(tb, te, dn);
    debug_print("Execution Time [msec]: %0.6f", dn);
    CU_ASSERT_TRUE(e_relative<3e-6);

}
/*============================================================================*/
void solution_Gaussian_WP_timepot
(
    void*       params_,
    matlib_xv   x,
    matlib_real t,
    matlib_zv   u
)
/* 
 * 
 * phi = g_0 *cos(mu*t)
 * Phi = g_0*sin(mu*t)/mu
 *
 * */ 
{
    debug_enter("%s", "");
    matlib_complex coeff = 1.0/(4.0*I*t+1);
    matlib_real *xptr;
    /* parameters 
     * param[0] = g_0, 
     * param[1] = mu
     *
     * */ 
    matlib_real* params = (matlib_real*)params_;
    matlib_real theta = params[1]*t;
    matlib_real Phi   = params[0]*sin(theta)/params[1];

    for (xptr = x.elem_p; xptr<(x.elem_p+x.len); xptr++)
    {
        *(u.elem_p) = csqrt(coeff)*cexp(-*xptr**xptr*coeff+I*Phi);
        u.elem_p++;
    }
    debug_exit("%s", "");
}

void timedependent_zpotential
( 
    void*          params_,
    matlib_complex m_coeff[2], 
    matlib_xv      x, 
    matlib_xv      t, 
    matlib_zm      y
)
{
    debug_enter("%s", "");
    matlib_index i, j;

    matlib_real* params = (matlib_real*) params_;

    for(j=0; j<y.lenr; j++)
    {
        for (i=0; i<x.len; i++)
        {
            *(y.elem_p) =  m_coeff[0]+
                           0.5*m_coeff[1]*params[0]*
                           (cos(params[1]*(t.elem_p[j]))+cos(params[1]*(t.elem_p[j+1])));
            y.elem_p++;
        }
    }
    debug_exit("%s", "");
}

void test_solve_Schroedinger_IVP2b(void)
{
    matlib_index p = 4;
    matlib_index nr_LGL = 2*p+1;
    matlib_index N = 400;
    matlib_index Nt = 1000;
    matlib_real dt = 0.5e-3;
    matlib_real domain[2] = {-15.0, 15.0 };
    matlib_real params[2] = {1.0, 2*M_PI};
    matlib_real e_relative;

    struct timespec tb, te;
    matlib_real dn;
    START_TIMMING(tb);

    e_relative = solve_Schroedinger_IVP2( p, nr_LGL, N, domain, dt, Nt, 
                                          (void*)params,
                                          solution_Gaussian_WP_timepot,
                                          timedependent_zpotential);
    GET_DURATION(tb, te, dn);
    debug_print("Execution Time [msec]: %0.6f", dn);
    CU_ASSERT_TRUE(e_relative<3e-6);

}



/*============================================================================*/
void _Schroedinger_IVP3_assemble_GSMM
(
    matlib_index     N,
    matlib_complex   s_coeff,
    matlib_complex   m_coeff[2],
    matlib_zm_sparse M1,
    matlib_zm_sparse M2,
    matlib_zm_sparse GSMM
)
{
    matlib_index i, j;

    for(i=0; i<GSMM.lenc; i++)
    {
        for(j=GSMM.rowIn[i]; j<GSMM.rowIn[i+1]; j++)
        {
            GSMM.elem_p[j] = m_coeff[0]*M1.elem_p[j] + m_coeff[1]*M2.elem_p[j];
        }
    }
    GSMM.elem_p[0] =  s_coeff/2.0 + GSMM.elem_p[0];
    GSMM.elem_p[1] = -s_coeff/2.0 + GSMM.elem_p[1];
    for(i=1; i<N; i++)
    {
        GSMM.elem_p[GSMM.rowIn[i]]   =  s_coeff     + GSMM.elem_p[GSMM.rowIn[i]  ];
        GSMM.elem_p[GSMM.rowIn[i]+1] = -s_coeff/2.0 + GSMM.elem_p[GSMM.rowIn[i]+1];
    }
    GSMM.elem_p[GSMM.rowIn[i]] = s_coeff/2.0 + GSMM.elem_p[GSMM.rowIn[i]];
    
    for(i=N+1; i<GSMM.lenc; i++)
    {
        GSMM.elem_p[GSMM.rowIn[i]] = s_coeff + GSMM.elem_p[GSMM.rowIn[i]];
    }
}

matlib_real solve_Schroedinger_IVP3
(
    matlib_index p,
    matlib_index nr_LGL,
    matlib_index N,
    matlib_real domain[2],
    matlib_real dt,
    matlib_index Nt,
    void* params_,
    void  (*fun_p)(void*, matlib_xv, matlib_real, matlib_zv)
)
{

    debug_enter( "polynomial degree: %d, nr. of LGL points: %d", p, nr_LGL );

    matlib_real* params = (matlib_real*)params_;
    matlib_index i, j;
    matlib_index P = nr_LGL-1;

    matlib_real x_l = domain[0];
    matlib_real x_r = domain[1];

    matlib_xv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_xm FM, IM, Q;
    matlib_create_xm(    p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_xm( xi.len,    p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    legendre_LGLdataIM( xi, IM);
    fem1d_quadM( quadW, IM, &Q);

    /* generate the grid */ 
    matlib_real t0 = 0;
    matlib_xv x, t;
    fem1d_ref2mesh (xi, N, x_l, x_r, &x);
    fem1d_linespace(t0, dt, Nt, &t);
    debug_body("length of x: %d", x.len);
    debug_body("length of t: %d", t.len);

    /* Provide the initial condition */
    matlib_zv u0, u_exact, u_computed;
    matlib_create_zv( x.len,         &u0, MATLIB_COL_VECT);
    matlib_create_zv( x.len,    &u_exact, MATLIB_COL_VECT);
    matlib_create_zv( x.len, &u_computed, MATLIB_COL_VECT);
    
    (*fun_p)(params, x, t0, u0);
    //DEBUG_PRINT_ZV(u0, "%s :", "Initial condition");
    
    matlib_real irho = dt/2.0;
    matlib_real J    = 0.5*(x_r-x_l)/N;

    matlib_complex m_coeff[2], s_coeff = I*irho/(J*J);

    /* Assemble the global mass matrix */ 
    matlib_zv phi;
    matlib_create_zv( x.len, &phi, MATLIB_COL_VECT);

    matlib_zm_sparse M, M1;
    for(i=0; i<phi.len; i++)
    {
        phi.elem_p[i] = 1.0;
    }
    fem1d_zm_sparse_GMM(p, Q, phi, &M);
    
    for(i=0; i<phi.len; i++)
    {
        phi.elem_p[i] = x.elem_p[i];
    }
    fem1d_zm_sparse_GMM(p, Q, phi, &M1);

    matlib_complex GSMMa[M.rowIn[M.lenc]];
    matlib_zm_sparse GSMM = {.lenc   = M.lenc,
                             .lenr   = M.lenr,
                             .rowIn  = M.rowIn,
                             .colIn  = M.colIn,
                             .elem_p = GSMMa };
    /* Setup the sparse linear system */
    /* Solve the equation while marching in time 
     *
     * Linear system M * V_vb = Pvb 
     * */
    
    matlib_zv U_tmp, V_tmp, V_vb, U_vb, Pvb;
    matlib_real dim = N*(p+1);
    matlib_create_zv( M.lenc,   &Pvb, MATLIB_COL_VECT);
    matlib_create_zv( M.lenc,  &V_vb, MATLIB_COL_VECT);
    matlib_create_zv(    dim, &V_tmp, MATLIB_COL_VECT);
    matlib_create_zv(    dim, &U_tmp, MATLIB_COL_VECT);

    fem1d_ZFLT( N, FM, u0, U_tmp);


    matlib_real norm_actual;
    matlib_xv e_relative;
    matlib_create_xv( Nt, &e_relative, MATLIB_COL_VECT);

    pardiso_solver_t data = { .nsparse  = 1, 
                              .mnum     = 1, 
                              .mtype    = PARDISO_COMPLEX_SYM,
                              .sol_enum = PARDISO_LHS, 
                              .smat_p   = (void*) &GSMM,
                              .rhs_p    = (void*) &Pvb,
                              .sol_p    = (void*) &V_vb};

    debug_body("%s", "Solver data initialized");

    BEGIN_DTRACE
        debug_print("dimension of the sparse square matrix: %d", M.lenc);
        for(i=0; i<M.lenc; i++)
        {
            for(j=M.rowIn[i]; j<M.rowIn[i+1]; j++)
            {
                debug_print( "M(%d,%d): % 0.16f %+0.16fi", 
                             i, M.colIn[j], M.elem_p[j]);
            }
        }
    END_DTRACE

    data.phase_enum = PARDISO_INIT;
    matlib_pardiso(&data);

    for (i=0; i<Nt; i++)
    {

        debug_body("begin iteration: %d", i);
        m_coeff[0] = 1.0;
        m_coeff[1] = -I*irho*0.5*params[0]*(cos(params[1]*t.elem_p[i])
                     + cos(params[1]*t.elem_p[i+1]));
        debug_body("m_coeff[0]: %0.16f%+0.16fi", m_coeff[0]);
        debug_body("m_coeff[1]: %0.16f%+0.16fi", m_coeff[1]);

        _Schroedinger_IVP3_assemble_GSMM(N, s_coeff, m_coeff, M, M1, GSMM);

        fem1d_ZPrjL2F(p, U_tmp, Pvb);

        data.phase_enum = PARDISO_ANALYSIS_AND_FACTOR;
        matlib_pardiso(&data);

        data.phase_enum = PARDISO_SOLVE_AND_REFINE;
        matlib_pardiso(&data);
        
        fem1d_ZF2L(p, V_vb, V_tmp);

        /* 2.0 * V_tmp -U_tmp --> U_tmp*/ 
        matlib_zaxpby(2.0, V_tmp, -1.0, U_tmp );

        (*fun_p)(params, x, t.elem_p[i+1], u_exact);
        fem1d_ZFLT(N, FM, u_exact, V_tmp);

        /* Error analysis */ 
        norm_actual = fem1d_ZNorm2(p, N, V_tmp);
        matlib_zaxpy(-1.0, U_tmp, V_tmp );
        e_relative.elem_p[i] = fem1d_ZNorm2(p, N, V_tmp)/norm_actual;
    }

    (*fun_p)(params, x, t.elem_p[t.len-1], u_exact);
    fem1d_ZILT(N, IM, U_tmp, u_computed);

    /* Write the final evolution data to a file */ 
    matlib_xzvwrite_csv( "initial_state.dat", 1, &x, 1, &u0);

    matlib_zv tmp_mat[2] = {u_exact, u_computed }; 
    matlib_xzvwrite_csv( "final_state.dat", 1, &x, 2, tmp_mat);
    
    matlib_xvwrite_csv("e_relative.dat", 1, &e_relative);
    DEBUG_PRINT_XV(e_relative, "%s :", "relative error");

    data.phase_enum = PARDISO_FREE;
    matlib_pardiso(&data);

    matlib_free(u0.elem_p);
    matlib_free(u_exact.elem_p);
    matlib_free(u_computed.elem_p);
    matlib_free(phi.elem_p);
    matlib_free(U_tmp.elem_p);
    matlib_free(V_tmp.elem_p);
    matlib_free(V_vb.elem_p);
    matlib_free(Pvb.elem_p);

    /* Free the sparse matrix */ 
    matlib_free(M.elem_p);
    matlib_free(M.colIn);
    matlib_free(M.rowIn);

    matlib_real r = e_relative.elem_p[Nt-1];
    debug_exit("Relative error: % 0.16g", r);
    return(r);
}

/*============================================================================*/
void test_solve_Schroedinger_IVP3a(void)
{
    matlib_index p = 4;
    matlib_index nr_LGL = 2*p+1;
    matlib_index N = 400;
    matlib_index Nt = 2000;
    matlib_real dt = 0.5e-3;
    matlib_real domain[2] = {-15.0, 15.0 };
    matlib_real params[2] = {1.0, 2*M_PI};
    matlib_real e_relative;


    struct timespec tb, te;
    matlib_real dn;
    START_TIMMING(tb);

    e_relative = solve_Schroedinger_IVP3( p, nr_LGL, N, domain, dt, Nt, 
                                          (void*)params,
                                          solution_Gaussian_WP_linpot);
    GET_DURATION(tb, te, dn);
    debug_print("Execution Time [msec]: %0.6f", dn);
    CU_ASSERT_TRUE(e_relative<3.0e-6);

}

/*============================================================================+/
 | Test runner
 |
 |
 +============================================================================*/

int main(void)
{
    CU_pSuite pSuite = NULL;

    /* initialize the CUnit test registry */
    if (CUE_SUCCESS != CU_initialize_registry())
    {
        return CU_get_error();
    }

    /* Create a test array */
    CU_TestInfo test_array[] = 
    {
        { "Constant potential"              , test_solve_Schroedinger_IVP1a},
        { "Harmonic potential"              , test_solve_Schroedinger_IVP1b},
        { "Linear time-dependent potential" , test_solve_Schroedinger_IVP2a},
        { "Time-dependent potential"        , test_solve_Schroedinger_IVP2b},
        { "Linear time-dependent potential2", test_solve_Schroedinger_IVP3a},
        CU_TEST_INFO_NULL,
    };

    /* Create the test suite */ 
    CU_SuiteInfo suites[] = 
    {
        { "Linear Schroedinger Equation", init_suite, clean_suite, test_array },
        CU_SUITE_INFO_NULL,
    }; 

    /* Register test suites */ 
    CU_ErrorCode CU_error = CU_register_suites(suites); 
    if (CU_error != CUE_SUCCESS) 
    {
        debug_body("%s", CU_get_error_msg());
        CU_cleanup_registry();
        return CU_get_error();
    }

   /* Run all tests using the CUnit Basic interface */
   CU_basic_set_mode(CU_BRM_VERBOSE);
   CU_basic_run_tests();
   CU_cleanup_registry();
   return CU_get_error();
}

