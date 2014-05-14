/*============================================================================+/
 | File: Cunit_fem1d_NonlinearSchroedinger.c 
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
| Nonlinear Schrodinger Eqation (IVP):
| iu_t + u_xx + chi |u|^2u= 0
| t=0: u(x,t) = u_0(x) -> Initial data
| 
/+============================================================================*/
matlib_real solve_NLSEquation_IVP
(
    matlib_index p,
    matlib_index nr_LGL,
    matlib_index N,
    matlib_real  domain[2],
    matlib_real  dt,
    matlib_index Nt,
    matlib_real  Chi,
    void  (*func_p)(matlib_real*, matlib_xv, matlib_real, matlib_zv)
)
/* 
 * */ 
{
    struct timespec tb0, te0;
    matlib_real dn0;
    START_TIMMING(tb0);

    struct timespec tb, te;
    matlib_real dn;
    START_TIMMING(tb);

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

    /* time discretization */ 

    matlib_real ta[2] = {0, Nt*dt};

    matlib_zv u0, phi, u_exact, u_computed;
    matlib_create_zv( x.len,         &u0, MATLIB_COL_VECT);
    matlib_create_zv( x.len,        &phi, MATLIB_COL_VECT);
    matlib_create_zv( x.len,    &u_exact, MATLIB_COL_VECT);
    matlib_create_zv( x.len, &u_computed, MATLIB_COL_VECT);
    
    /* Provide the initial condition: u0 */
    matlib_real params[2] = {1.0, 2*M_PI};
    (*func_p)(params, x, ta[0], u0);
    //DEBUG_PRINT_ZV(u0, "%s :", "Initial condition");

    matlib_real irho = dt/2.0;
    matlib_real J    = 0.5*(x_r-x_l)/N;

    matlib_complex coeff = 0.0 + I*irho/(J*J);
    matlib_complex m_coeff[2] = {1.0, -I*irho};
    debug_body("coeff: %0.16f%+0.16fi", coeff);


    /* Assemble the global mass matrix */ 
    matlib_zm q;
    matlib_zm_sparse M;
    for(i=0; i<phi.len; i++)
    {
        phi.elem_p[i] = 1.0;
    }
    fem1d_zm_sparse_GMM(p, Q, phi, &M);

    /* Global Stiffness-Mass Matrix */ 
    fem1d_zm_sparse_GSM(N, coeff, M);

    BEGIN_DTRACE
        debug_print("dimension of the sparse square matrix: %d", M.lenc);
        for(i=0; i<M.lenc; i++)
        {
            for(j=M.rowIn[i]; j<M.rowIn[i+1]; j++)
            {
                debug_print("M(%d,%d): % 0.16f %+0.16fi", i, M.colIn[j], M.elem_p[j]);
            }
        }
    END_DTRACE
    
    /* Setup the sparse linear system */
    /* Solve the equation while marching in time 
     *
     * Linear system M * V_vb = Pvb 
     * */
    
    matlib_zv U_tmp, V_tmp, V_vb, U_vb, Pvb, u_NL, U_NL, PNL_vb;
    matlib_real dim = N*(p+1);
    matlib_create_zv( M.lenc,    &Pvb, MATLIB_COL_VECT);
    matlib_create_zv( M.lenc,   &V_vb, MATLIB_COL_VECT);
    matlib_create_zv(  x.len,   &u_NL, MATLIB_COL_VECT);
    matlib_create_zv(    dim,   &U_NL, MATLIB_COL_VECT);
    matlib_create_zv( M.lenc, &PNL_vb, MATLIB_COL_VECT);
    matlib_create_zv(    dim,  &V_tmp, MATLIB_COL_VECT);
    matlib_create_zv(    dim,  &U_tmp, MATLIB_COL_VECT);

    fem1d_ZFLT( N, FM, u0, U_tmp);
    for(i=0; i<V_tmp.len; i++)
    {
        V_tmp.elem_p[i] = 0;
    }

    matlib_real t = 0, norm_actual;
    matlib_xv e_relative;
    matlib_create_xv( Nt, &e_relative, MATLIB_COL_VECT);

    pardiso_solver_t data = { .nsparse  = 1, 
                              .mnum     = 1, 
                              .mtype    = PARDISO_COMPLEX_SYM,
                              .sol_enum = PARDISO_RHS, 
                              .smat_p   = (void*)&M,
                              .rhs_p    = (void*)&PNL_vb,
                              .sol_p    = (void*)&V_vb};

    debug_body("%s", "Solver data initialized");

    data.phase_enum = PARDISO_INIT;
    matlib_pardiso(&data);

    matlib_real res;
    matlib_index Nth = 50, iter;

    GET_DURATION(tb, te, dn);
    debug_print("Time of initialize[msec]: %0.6f", dn);

    data.phase_enum = PARDISO_ANALYSIS_AND_FACTOR;
    matlib_pardiso(&data);

    GET_DURATION(te, tb, dn);
    debug_print( "Time for factorization[msec]: %0.6f", dn);
    for (i=0; i<Nt; i++)
    {
        debug_body("time-step: %d", i);
        fem1d_ZPrjL2F(p, U_tmp, Pvb);


        /* The non-linear term is I*irho*chi*|v|^2v on the RHS 
         * which is computed using the solution obtained in the 
         * previous iteration 
         * */ 
        iter = 0;
        res  = 1.0;
        START_TIMMING(tb);
        while((res>TOL) && (iter<Nth))
        {
            fem1d_ZILT(N, IM, V_tmp, u_NL);
            for (j=0; j<u_NL.len; j++)
            {
                u_NL.elem_p[j] = I*irho*Chi*u_NL.elem_p[j]*u_NL.elem_p[j]*conj(u_NL.elem_p[j]);
            }
            fem1d_ZFLT( N, FM, u_NL, U_NL);
            fem1d_ZPrjL2F(p, U_NL, PNL_vb);
            matlib_zaxpy(1.0, Pvb, PNL_vb);
            //DEBUG_PRINT_ZV(PNL_vb, "%s :", "projection");

            //GET_DURATION(tb, te, dn);
            //debug_print( "Time for other crap at %d-th time-step[msec]: %0.6f", i, dn);
            data.phase_enum = PARDISO_SOLVE_AND_REFINE;
            data.mnum = i+1;

            matlib_pardiso(&data);
            
            matlib_zaxpby(1.0, V_tmp, 0, U_NL );

            fem1d_ZF2L(p, PNL_vb, V_tmp);
            matlib_zaxpy(-1.0, V_tmp, U_NL );
            res = sqrt(J) * fem1d_ZNorm2(p, N, U_NL);
            iter++;
            debug_body("time-step: %d, iteration: %d, residue: %0.16f", i, iter, res);
        }
        GET_DURATION(tb, te, dn);
        //debug_print( "Time for %d-th time-step[msec]: %0.6f", i, dn);
        if((iter==Nth) || (isnan(res)))
        {
            term_exec("Iteration reached threshold with res: %0.16f", res);
        }

        /* 2.0 * V_tmp -U_tmp --> U_tmp*/ 
        matlib_zaxpby(2.0, V_tmp, -1.0, U_tmp );
        t += dt;

        (*func_p)(params, x, t, u_exact);
        fem1d_ZFLT(N, FM, u_exact, V_tmp);

        /* Error analysis */ 
        norm_actual = fem1d_ZNorm2(p, N, V_tmp);
        matlib_zaxpy(-1.0, U_tmp, V_tmp );
        e_relative.elem_p[i] = fem1d_ZNorm2(p, N, V_tmp)/norm_actual;
        debug_body("relative error: %0.16f", e_relative.elem_p[i]);
    }

    GET_DURATION(tb0, te0, dn0);
    debug_print("Time for solution [sec]: %0.6f", dn0/1e3);

    assert(fabs(ta[1]-t)<TOL);
    data.phase_enum = PARDISO_FREE;
    matlib_pardiso(&data);

    (*func_p)(params, x, t, u_exact);
    fem1d_ZILT(N, IM, U_tmp, u_computed);

    matlib_xzvwrite_csv("u_initial.dat", 1, &x, 1, &u0);
    matlib_xvwrite_csv("e_relative.dat", 1, &e_relative);

    matlib_zv tmp_mat[2] = {u_exact, u_computed }; 
    matlib_xzvwrite_csv("final_state.dat", 1, &x, 2, tmp_mat);

    DEBUG_PRINT_XV(e_relative, "%s :", "relative error");

    matlib_real r = e_relative.elem_p[Nt-1];
    debug_exit("Relative error: % 0.16g", r);
    return(r);
    
}
/* Bright 1-Soliton Solution */ 
/*============================================================================*/
void solution_BrightSolitonNLS
( 
    matlib_real *params,
    matlib_xv   x, 
    matlib_real t, 
    matlib_zv   u
)
/*                  e^(i*eta^2*t)   
 * u(x,t) = eta ---------------------
 *                  cosh[eta*(x)]
 *
 * */ 
{
    debug_enter("%s", "");
    
    matlib_real *xptr;
    matlib_real eta = 1.0;
    matlib_complex e = cexp(I*eta*eta*t);

    for (xptr = x.elem_p; xptr<(x.elem_p+x.len); xptr++)
    {
        *(u.elem_p) = eta*e/cosh(eta*(*xptr));
        u.elem_p++;
    }
    debug_exit("%s", "");
}

void test_solve_NLSEquation_IVPa(void)
{
    matlib_index p = 4;
    matlib_index nr_LGL = 2*p+1;
    matlib_index N = 400;
    matlib_index Nt = 4000;
    matlib_real dt = 0.5e-3;
    matlib_real domain[2] = {-15.0, 15.0 };
    matlib_real e_relative, Chi = 2.0;
    
    struct timespec tb, te;
    matlib_real dn;
    START_TIMMING(tb);

    e_relative = solve_NLSEquation_IVP( p, nr_LGL, N, domain, dt, Nt, Chi,
                                        solution_BrightSolitonNLS);

    GET_DURATION(tb, te, dn);
    debug_print("Execution Time [msec]: %0.6f", dn);

    CU_ASSERT_TRUE(e_relative<1e-6);

}

/* Bright 2-Soliton Solution */ 
/*============================================================================*/
void solution_BrightSolitonNLS2
( 
    matlib_real *params,
    matlib_xv   x, 
    matlib_real t, 
    matlib_zv   u
)
{
    debug_enter("%s", "");
    
    matlib_real *xptr;
    matlib_real eta = 0.5;
    matlib_complex e = cexp(I*eta*eta*t);

    for (xptr = x.elem_p; xptr<(x.elem_p+x.len); xptr++)
    {
        *(u.elem_p) = 4.0*eta*(cosh(3.0*eta**xptr)+3.0*cpow(e,8.0)*cosh(eta**xptr))*e
                      /(cosh(4.0*eta**xptr)+4.0*cosh(2.0*eta**xptr)+3.0*cos(8.0*eta*eta*t));
        u.elem_p++;
    }
    debug_exit("%s", "");
}

void test_solve_NLSEquation_IVPb(void)
{
    matlib_index p = 4;
    matlib_index nr_LGL = 2*p+1;
    matlib_index N = 1000;
    matlib_index Nt = 4000;
    matlib_real dt = 2.5e-4;
    matlib_real domain[2] = {-25.0, 25.0 };
    matlib_real e_relative, Chi = 2.0;
    
    struct timespec tb, te;
    matlib_real dn;
    START_TIMMING(tb);

    e_relative = solve_NLSEquation_IVP( p, nr_LGL, N, domain, dt, Nt, Chi,
                                        solution_BrightSolitonNLS2);

    GET_DURATION(tb, te, dn);
    debug_print("Execution Time [msec]: %0.6f", dn);

    CU_ASSERT_TRUE(e_relative<3e-6);

}
/*============================================================================+/
| Example 3.:
| Gross-Pitaevvskii Eqation (IVP):
| iu_t + u_xx + phi(x,t) u + chi |u|^2u= 0
| t=0: u(x,t) = u_0(x) -> Initial data
| This example deals with time-dependent real potential.
| chi=2, phi(x,t) = x cos(mu*t)
| 
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


matlib_real solve_GPEquation_IVP
(
    matlib_index p,
    matlib_index nr_LGL,
    matlib_index N,
    matlib_real  domain[2],
    matlib_real  dt,
    matlib_index Nt,
    matlib_real  Chi,
    void* params,    
    void  (*fun_p)(void*, matlib_xv, matlib_real, matlib_zv),
    void  (*pot_p)(void*, matlib_complex*, matlib_xv, matlib_xv, matlib_zm)
)
/* 
 * */ 
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
    
    matlib_zv U_tmp, V_tmp, V_vb, U_vb, Pvb, u_NL, U_NL, PNL_vb;
    matlib_real dim = N*(p+1);
    matlib_create_zv( M.lenc,    &Pvb, MATLIB_COL_VECT);
    matlib_create_zv( M.lenc,   &V_vb, MATLIB_COL_VECT);
    matlib_create_zv(  x.len,   &u_NL, MATLIB_COL_VECT);
    matlib_create_zv(    dim,   &U_NL, MATLIB_COL_VECT);
    matlib_create_zv( M.lenc, &PNL_vb, MATLIB_COL_VECT);
    matlib_create_zv(    dim,  &V_tmp, MATLIB_COL_VECT);
    matlib_create_zv(    dim,  &U_tmp, MATLIB_COL_VECT);

    fem1d_ZFLT( N, FM, u0, U_tmp);
    for(i=0; i<V_tmp.len; i++)
    {
        V_tmp.elem_p[i] = 0;
    }


    matlib_real norm_actual;
    matlib_xv e_relative;
    matlib_create_xv( Nt, &e_relative, MATLIB_COL_VECT);

    pardiso_solver_t data = { .nsparse  = nsparse, 
                              .mnum     = 1, 
                              .mtype    = PARDISO_COMPLEX_SYM,
                              .sol_enum = PARDISO_LHS, 
                              .smat_p   = (void*)&M,
                              .rhs_p    = (void*)&PNL_vb,
                              .sol_p    = (void*)&V_vb};

    debug_body("%s", "Solver data initialized");

    data.phase_enum = PARDISO_INIT;
    matlib_pardiso(&data);

    matlib_real res;
    matlib_index Nth = 50, iter;
    matlib_index Nt_ = Nt/nsparse;
    debug_body("Nt_: %d", Nt_);
    matlib_xv t_tmp  = {.len = nsparse}; 

    for(matlib_index k=0; k<Nt_; k++ )
    {
        t_tmp.elem_p = t.elem_p + k*nsparse;
        (*pot_p)( params, m_coeff, x, t_tmp, phi);
        fem1d_zm_nsparse_GMM(p, N, nsparse, Q, &phi, &q, &M, FEM1D_GET_NZE_ONLY);
        fem1d_zm_nsparse_GSM(N, s_coeff, M);

        for (i=0; i<nsparse; i++)
        {
            debug_body("time-step: %d", i+k*nsparse);
            fem1d_ZPrjL2F(p, U_tmp, Pvb);

            data.phase_enum = PARDISO_ANALYSIS_AND_FACTOR;
            data.mnum = i+1;
            matlib_pardiso(&data);

            /* The non-linear term is I*irho*chi*|v|^2v on the RHS 
             * which is computed using the solution obtained in the 
             * previous iteration 
             * */ 
            iter = 0;
            res  = 1.0;
            while((res>TOL) && (iter<Nth))
            {
                fem1d_ZILT(N, IM, V_tmp, u_NL);
                for (j=0; j<u_NL.len; j++)
                {
                    u_NL.elem_p[j] = I*irho*Chi*u_NL.elem_p[j]*u_NL.elem_p[j]*conj(u_NL.elem_p[j]);
                }
                fem1d_ZFLT( N, FM, u_NL, U_NL);
                fem1d_ZPrjL2F(p, U_NL, PNL_vb);
                matlib_zaxpy(1.0, Pvb, PNL_vb);
                //DEBUG_PRINT_ZV(PNL_vb, "%s :", "projection");

                data.phase_enum = PARDISO_SOLVE_AND_REFINE;
                matlib_pardiso(&data);
                
                matlib_zaxpby(1.0, V_tmp, 0, U_NL );

                fem1d_ZF2L(p, V_vb, V_tmp);
                matlib_zaxpy(-1.0, V_tmp, U_NL );
                res = sqrt(J) * fem1d_ZNorm2(p, N, U_NL);
                iter++;
                debug_body("time-step: %d, iteration: %d, residue: %0.16f", i+k*nsparse, iter, res);
            }
            if((iter==Nth) || (isnan(res)))
            {
                term_exec("Iteration reached threshold with res: %0.16f", res);
            }

            /* 2.0 * V_tmp -U_tmp --> U_tmp*/ 
            matlib_zaxpby(2.0, V_tmp, -1.0, U_tmp );

            (*fun_p)(params, x, t_tmp.elem_p[i+1], u_exact);
            fem1d_ZFLT(N, FM, u_exact, V_tmp);

            /* Error analysis */ 
            norm_actual = fem1d_ZNorm2(p, N, V_tmp);
            matlib_zaxpy(-1.0, U_tmp, V_tmp );
            e_relative.elem_p[k*nsparse+i] = fem1d_ZNorm2(p, N, V_tmp)/norm_actual;
            debug_body("relative error: %0.16f", e_relative.elem_p[k*nsparse+i]);
        }
    }
    data.phase_enum = PARDISO_FREE;
    matlib_pardiso(&data);
    fem1d_zm_nsparse_GMM(p, N, nsparse, Q, &phi, &q, &M, FEM1D_GMM_FREE);

    (*fun_p)(params, x, t.elem_p[t.len-1], u_exact);
    fem1d_ZILT(N, IM, U_tmp, u_computed);

    matlib_xzvwrite_csv("u_initial.dat", 1, &x, 1, &u0);
    matlib_xvwrite_csv("e_relative.dat", 1, &e_relative);

    matlib_zv tmp_mat[2] = {u_exact, u_computed }; 
    matlib_xzvwrite_csv("final_state.dat", 1, &x, 2, tmp_mat);

    matlib_real r = e_relative.elem_p[Nt-1];
    debug_exit("Relative error: % 0.16g", r);
    return(r);
    
}
/* Bright 1-Soliton Soltuion */
/*============================================================================*/
void solution_BrightSolitonGPE
( 
    void*       params_,
    matlib_xv   x, 
    matlib_real t, 
    matlib_zv   u
)
/*                  e^(i*eta^2*t)   
 * u(x,t) = eta ---------------------e^(i*Xi)
 *               cosh[eta*(x-nu)]
 *
 * beta   = g_0 *sin(mu*t)/mu 
 * nu     = (2*g_0*cos(mu*t))/mu^2 - (2*g_0)/mu^2;
 * ibeta2 = (1/mu^3)(mu*(g_0^2*t)/2 - (g_0^2*sin(2*mu*t))/4;
 * Xi     = beta*x-ibeta2;
 *
 * */ 
{

    debug_enter("%s", "");
    
    matlib_real* params = (matlib_real*) params_;
    matlib_real *xptr;
    /* parameters 
     * param[0] = g_0, 
     * param[1] = mu
     * */ 
    matlib_real eta = 1, Xi;
    
    matlib_complex e   = cexp(I*eta*eta*t);
    matlib_real theta  = params[1]*t;
    matlib_real beta   = params[0]*sin(theta)/params[1];
    matlib_real nu     = 2.0*params[0]*(cos(theta)-1.0)/(params[1]*params[1]);
    matlib_real ibeta2 = params[0]*params[0]*
                         (theta/2.0 - sin(2*theta)/4.0)/(params[1]*params[1]*params[1]);

    for (xptr = x.elem_p; xptr<(x.elem_p+x.len); xptr++)
    {
        Xi = beta**xptr - ibeta2;
        *(u.elem_p) = eta*e*cexp(I*Xi)/cosh(eta*(*xptr+nu));
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
            *(y.elem_p) =  m_coeff[0]
                           +0.5*m_coeff[1]**(x.elem_p+i)*params[0]*
                           (cos(params[1]*(t.elem_p[j]))+cos(params[1]*(t.elem_p[j+1])));
            y.elem_p++;
        }
    }
    debug_exit("%s", "");
}

void test_solve_GPEquation_IVP(void)
{
    matlib_index p = 4;
    matlib_index nr_LGL = 2*p+1;
    matlib_index N = 1000;
    matlib_index Nt = 4000;
    matlib_real dt = 0.5e-3;
    matlib_real domain[2] = {-15.0, 15.0 };
    matlib_real e_relative, Chi = 2.0;

    matlib_real params[2] = {1.0, 2*M_PI};
    
    struct timespec tb, te;
    matlib_real dn;
    START_TIMMING(tb);
 
    e_relative = solve_GPEquation_IVP( p, nr_LGL, N, domain, dt, Nt, Chi, 
                                       (void*)params,
                                       solution_BrightSolitonGPE,
                                       linear_timedependent_zpotential);
    GET_DURATION(tb, te, dn);
    debug_print("Execution Time [msec]: %0.6f", dn);

    CU_ASSERT_TRUE(e_relative<1e-6);

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
        { "NLS Equation 1-Soliton", test_solve_NLSEquation_IVPa},
        { "NLS Equation 2-Soliton", test_solve_NLSEquation_IVPb},
        { "GP Equation, linear time-dependent potential", test_solve_GPEquation_IVP},
        CU_TEST_INFO_NULL,
    };

    /* Create the test suite */ 
    CU_SuiteInfo suites[] = 
    {
        { "Nonlinear Schroedinger Equation", init_suite, clean_suite, test_array },
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

