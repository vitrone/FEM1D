/*============================================================================+/
 | File: Cunit_SchroedingerEquation.c 
 | Description: Test 
 |
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
#include "omp.h"

//#define NDEBUG

#include "legendre.h"
#include "fem1d.h"
#include "pfem1d.h"
#include "pthpool.h"
#include "assert.h"

/* CUnit modules */
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>
#include <CUnit/TestDB.h>

static const double TOL = 1e-9;
#define NUM_THREADS 5
/*============================================================================*/
int init_suite(void)
{
      return 0;
}
int clean_suite(void)
{
      return 0;
}

/*============================================================================+/
| Example 1.:
| Schroedinger Equation (IVP):
| iu_t + u_xx + phi(x) u= 0
| t=0: u(x,t) = u_0(x) -> Initial data
| This example deal with time-independent real potential.
/+============================================================================*/
void* solution_Gaussian(void* mp)
{
    debug_enter("%s", "");
    matlib_index i;
    
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_dv x = *((matlib_dv*) (ptr->shared_data[0]));
    matlib_zv u = *((matlib_zv*) (ptr->shared_data[1]));
    double t    = *((double*)    (ptr->shared_data[2]));

    matlib_complex coeff = I/(-4.0*t+I);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    

    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);
    
    (u.elem_p) += (start_end_index[0]);

    for(i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        *(u.elem_p) = csqrt(coeff)*cexp(-x.elem_p[i]*x.elem_p[i]*coeff+I*t);
        u.elem_p++;
    }
    debug_exit("%s", "");
}
void* constant_zpotential(void* mp)
{
    debug_enter("%s", "");
    matlib_index i;
    
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_dv x = *((matlib_dv*) (ptr->shared_data[0]));
    matlib_zv u = *((matlib_zv*) (ptr->shared_data[1]));
    double* time_coeff = (double*)(ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    

    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);
    
    (u.elem_p) += (start_end_index[0]);

    for(i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        *(u.elem_p) = time_coeff[0]+time_coeff[1];
        u.elem_p++;
    }
}

void harmonic_zpotential
( 
          matlib_complex time_coeff[2], 
    const matlib_dv x, 
          matlib_zv y
)
{
    matlib_index i;
    if(y.len >= x.len)
    {
        for (i=0; i<x.len; i++)
        {
            *(y.elem_p) =  time_coeff[0]-time_coeff[1]**(x.elem_p+i)**(x.elem_p+i);
            y.elem_p++;
        }
    }
}

void solution_Hermite
(
    matlib_dv x,
    double    t,
    matlib_zv u
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
    
    double coeff[2] = {1.0/pow(M_PI, 0.25), sqrt(2.0)/pow(M_PI, 0.25)};
    matlib_complex e[2];
    double *xptr;
    for (xptr = x.elem_p; xptr<(x.elem_p+x.len); xptr++)
    {
        e[0] = cexp(-I*t); 
        e[1] = cexp(-I*3*t);
        *(u.elem_p) = (e[0]*coeff[0]+e[1]*coeff[1]**xptr)*exp(-0.5**xptr**xptr);
        u.elem_p++;
    }
    debug_exit("%s", "");
}
/*============================================================================*/
double solve_Schroedinger_IVP
(
    matlib_index p,
    matlib_index nr_LGL,
    matlib_index N,
    double domain[2],
    double dt,
    matlib_index Nt
)
/* 
 * */ 
{

    debug_enter( "polynomial degree: %d, nr. of LGL points: %d", p, nr_LGL );
    struct timespec tb, te;
    double dn;
    
    START_TIMMING(tb);

    matlib_index num_threads = 4;
    pthpool_data_t mp[num_threads];
    
    pthpool_create_threads(num_threads, mp);

    matlib_index i, j;
    matlib_index P = nr_LGL-1;

    double x_l = domain[0];
    double x_r = domain[1];

    matlib_dv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_dm FM, IM, Q;
    matlib_create_dm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_dm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    legendre_LGLdataIM( xi, IM);
    fem1d_quadM( quadW, IM, &Q);

    /* generate the grid */ 
    matlib_dv x;
    fem1d_ref2mesh (xi, N, x_l, x_r, &x);
    debug_body("length of x: %d", x.len);

    /* time discretization */ 

    double ta[2] = {0, Nt*dt};

    /* Provide the initial condition */
    matlib_zv u0;
    matlib_create_zv( x.len, &u0, MATLIB_COL_VECT);

    matlib_index Np[2] = {x.len/num_threads, x.len};
    void* shared_data[3] = { (void*) &x,
                             (void*) &u0,
                             (void*) &ta[0]};
 
    pthpool_func( Np, shared_data, solution_Gaussian, num_threads, mp);
    
    //DEBUG_PRINT_ZV(u0, "%s :", "Initial condition");
    
    /* Assemble the global mass matrix */ 
    matlib_zv phi;
    matlib_create_zv( x.len, &phi, MATLIB_COL_VECT);

    double irho = dt/2.0;
    double J    = 0.5*(x_r-x_l)/N;
    matlib_complex coeff = 0.0 + I*irho/(J*J);
    matlib_complex time_coeff[2] = {1.0, -I*irho};
    debug_body("coeff: %0.16f%+0.16fi", coeff);

    shared_data[1] = (void*) &phi;
    shared_data[2] = (void*) time_coeff;
    pthpool_func( Np, shared_data, constant_zpotential, num_threads, mp);

    matlib_zsparsem M;
    fem1d_ZSparseMGMM(p, Q, phi, &M);

    
    /* Setup the sparse linear system */

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

    //BEGIN_DEBUG
    //    debug_print("dimension of the sparse square matrix: %d", M.lenc);
    //    for(i=0; i<M.lenc; i++)
    //    {
    //        for(j=M.rowIn[i]; j<M.rowIn[i+1]; j++)
    //        {
    //            debug_print("M(%d,%d): % 0.16f %+0.16fi", i, M.colIn[j], M.elem_p[j]);
    //        }
    //    }
    //END_DEBUG

    /* Solve the equation while marching in time 
     *
     * Linear system M * V_vb = Pvb 
     * */
    
    matlib_zv U_tmp, V_tmp, V_vb, U_vb, Pvb;
    double dim = N*(p+1);
    matlib_create_zv( M.lenc,   &Pvb, MATLIB_COL_VECT);
    matlib_create_zv( M.lenc,  &V_vb, MATLIB_COL_VECT);
    matlib_create_zv(    dim, &V_tmp, MATLIB_COL_VECT);
    matlib_create_zv(    dim, &U_tmp, MATLIB_COL_VECT);

    fem1d_ZFLT( N, FM, u0, U_tmp);

    /* ONLY FOR ERROR: Compute the exact solution and make a forward LT */ 
    matlib_dv xii, quadWi, x1;
    legendre_LGLdataLT1( p, TOL, &xii, &quadWi);
    fem1d_ref2mesh (xii, N, x_l, x_r, &x1);
    matlib_dm IMi;
    matlib_create_dm( xii.len, p+1, &IMi, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    
    legendre_LGLdataIM( xii, IMi);
    matlib_zv u_exact, u_computed;
    matlib_create_zv( x1.len,    &u_exact, MATLIB_COL_VECT);
    matlib_create_zv( x1.len, &u_computed, MATLIB_COL_VECT);

    double t = 0, norm_actual;
    matlib_dv e_relative;
    matlib_create_dv( Nt, &e_relative, MATLIB_COL_VECT);

    pardiso_solver_t data = { .nsparse = 1, 
                                 .mnum    = 1, 
                                 .mtype   = PARDISO_COMPLEX_SYM,
                                 .sol_enum = PARDISO_LHS,
                                 .smat_p  = (void*)&M,
                                 .rhs_p   = (void*)&Pvb,
                                 .sol_p   = (void*)&V_vb};

    data.phase_enum = PARDISO_INIT;
    matlib_pardiso(&data);

    data.phase_enum = PARDISO_ANALYSIS_AND_FACTOR;
    matlib_pardiso(&data);
    GET_DURATION(tb, te, dn);
    debug_print("Init time[msec]: %0.4f", dn);

    Np[0] = x1.len/num_threads;
    Np[1] = x1.len;
    shared_data[0] = (void*) &x1;
    shared_data[1] = (void*) &u_exact;
    START_TIMMING(tb);
    for (i=1; i<Nt+1; i++)
    {
        //debug_print("begin iteration: %d", i);
        pfem1d_ZPrjL2F(p, U_tmp, Pvb, num_threads, mp);

        data.phase_enum = PARDISO_SOLVE_AND_REFINE;
        matlib_pardiso(&data);
        
        pfem1d_ZF2L(p, V_vb, V_tmp, num_threads, mp);

        //for(j=0; j<U_tmp.len; j++)
        //{
        //    U_tmp.elem_p[j] = 2.0 * V_tmp.elem_p[j]-U_tmp.elem_p[j];
        //}

        /* 2.0 * V_tmp -U_tmp --> U_tmp*/ 
        matlib_zaxpby(2.0, V_tmp, -1.0, U_tmp );
        t += dt;

        shared_data[2] = (void*) &t;
        pthpool_func( Np, shared_data, solution_Gaussian, num_threads, mp);
        pfem1d_ZFLT(N, FM, u_exact, U_tmp, num_threads, mp);

        /* Error analysis */ 
        norm_actual = matlib_znrm2(u_exact);
        matlib_zaxpy(-1.0, u_exact, u_computed);
        e_relative.elem_p[i-1] = matlib_znrm2(u_computed)/norm_actual;
    }
    data.phase_enum = PARDISO_FREE;
    matlib_pardiso(&data);

    shared_data[2] = (void*) &ta[1];
    pthpool_func( Np, shared_data, solution_Gaussian, num_threads, mp);
    fem1d_ZILT(N, IMi, U_tmp, u_computed);

    GET_DURATION(tb, te, dn);
    debug_print("Solution done time[sec]:%0.4f", dn/1e3);
    DEBUG_PRINT_DV(e_relative, "%s :", "relative error");

    double r = e_relative.elem_p[Nt-1];
    debug_exit("Relative error: % 0.16g", r);
    debug_body("%s", "signal threads to exit!");
    pthpool_destroy_threads(num_threads, mp);

    return(r);
    
    
}
void test_solve_Schroedinger_IVP_1a(void)
{
    matlib_index p = 4;
    matlib_index nr_LGL = 2*p+1;
    matlib_index N = 800;
    matlib_index Nt = 1000;
    double dt = 0.5e-3;
    double domain[2] = {-15.0, 15.0 };
    double e_relative;


    e_relative = solve_Schroedinger_IVP( p, nr_LGL, N, domain, dt, Nt);
    CU_ASSERT_TRUE(e_relative<1e-3);

}
//void test_solve_Schroedinger_IVP_1b(void)
//{
//    matlib_index p = 4;
//    matlib_index nr_LGL = 2*p+1;
//    matlib_index N = 800;
//    matlib_index Nt = 2000;
//    double dt = 0.5e-3;
//    double domain[2] = {-15.0, 15.0 };
//    double e_relative;
//
//
//    e_relative = solve_Schroedinger_IVP( p, nr_LGL, N, domain, dt, Nt,
//                                         solution_Hermite,
//                                         harmonic_zpotential);
//    CU_ASSERT_TRUE(e_relative<1e-3);
//
//}

/*============================================================================+/
| Example 2.:
| Schroedinger Equation (IVP):
| iu_t + u_xx + phi(x,t) u= 0
| t=0: u(x,t) = u_0(x) -> Initial data
| This example deals with time-dependent real potential.
|   2a.: Let phi = phit(t)*phix(x).
|   2b.: Let phi = phi(x,t).
/+============================================================================*/
void linear_timedependent_zpotential
( 
    matlib_complex time_coeff[2], 
    matlib_dv x, 
    double    t0, 
    double    dt, 
    matlib_zm y
)
{
    debug_enter("%s", "");
    matlib_index i, j;
    for(j=0; j<y.lenr; j++)
    {
        for (i=0; i<x.len; i++)
        {
            *(y.elem_p) =  time_coeff[0]-
                           0.5*time_coeff[1]**(x.elem_p+i)*
                           (cos(2*M_PI*(t0))+cos(2*M_PI*(t0+dt)));
            y.elem_p++;
        }
        t0 += dt;
    }
    debug_exit("%s", "");
}
void fem1d_ZNSparseMGSM
(
    matlib_index            N,
    matlib_complex   coeff,
    matlib_znsparsem M
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
double solve_Schroedinger_IVP2
(
    matlib_index p,
    matlib_index nr_LGL,
    matlib_index N,
    double domain[2],
    double dt,
    matlib_index Nt,
    void  (*func_p)(matlib_dv, double, matlib_zv),
    void  (*potential_p)(matlib_complex*, matlib_dv, double, double, matlib_zm)
)
/* 
 * */ 
{

    debug_enter( "polynomial degree: %d, nr. of LGL points: %d", p, nr_LGL );
    matlib_index i, j;
    matlib_index P = nr_LGL-1;

    double x_l = domain[0];
    double x_r = domain[1];

    matlib_dv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_dm FM, IM, Q;
    matlib_create_dm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_dm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    legendre_LGLdataIM( xi, IM);
    fem1d_quadM( quadW, IM, &Q);

    /* generate the grid */ 
    matlib_dv x;
    fem1d_ref2mesh (xi, N, x_l, x_r, &x);
    debug_body("length of x: %d", x.len);

    /* time discretization */ 

    double ta[2] = {0, Nt*dt};

    /* Provide the initial condition */
    matlib_zv u0;
    matlib_create_zv( x.len, &u0, MATLIB_COL_VECT);
    
    (*func_p)(x, ta[0], u0);
    //DEBUG_PRINT_ZV(u0, "%s :", "Initial condition");
    
    /* Assemble the global mass matrix */ 

    double irho = dt/2.0;
    double J   = 0.5*(x_r-x_l)/N;
    matlib_complex coeff = 0.0 + I*irho/(J*J);
    matlib_complex time_coeff[2] = {1.0, -I*irho};
    debug_body("coeff: %0.16f%+0.16fi", coeff);


    matlib_zm phi, q;
    matlib_znsparsem M;
    matlib_index nsparse = Nt;
    fem1d_ZNSparseMGMM(p, N, nsparse, Q, &phi, &q, &M, FEM1D_GMM_INIT);
    (*potential_p)(time_coeff, x, 0, dt, phi);
    fem1d_ZNSparseMGMM(p, N, nsparse, Q, &phi, &q, &M, FEM1D_GET_SPARSITY_NZE);
    fem1d_ZNSparseMGSM(N, coeff, M);
    
    /* Setup the sparse linear system */


    /* Solve the equation while marching in time 
     *
     * Linear system M * V_vb = Pvb 
     * */
    
    matlib_zv U_tmp, V_tmp, V_vb, U_vb, Pvb;
    double dim = N*(p+1);
    matlib_create_zv( M.lenc,   &Pvb, MATLIB_COL_VECT);
    matlib_create_zv( M.lenc,  &V_vb, MATLIB_COL_VECT);
    matlib_create_zv(    dim, &V_tmp, MATLIB_COL_VECT);
    matlib_create_zv(    dim, &U_tmp, MATLIB_COL_VECT);

    fem1d_ZFLT( N, FM, u0, U_tmp);

    /* ONLY FOR ERROR: Compute the exact solution and make a forward LT */ 
    matlib_dv xii, quadWi, x1;
    legendre_LGLdataLT1( p, TOL, &xii, &quadWi);
    fem1d_ref2mesh (xii, N, x_l, x_r, &x1);
    matlib_dm IMi;
    matlib_create_dm( xii.len, p+1, &IMi, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    
    legendre_LGLdataIM( xii, IMi);
    matlib_zv u_exact, u_computed;
    matlib_create_zv( x1.len,    &u_exact, MATLIB_COL_VECT);
    matlib_create_zv( x1.len, &u_computed, MATLIB_COL_VECT);

    double t = 0, norm_actual;
    matlib_dv e_relative;
    matlib_create_dv( Nt, &e_relative, MATLIB_COL_VECT);

    pardiso_solver_t data = { .nsparse = Nt, 
                                 .mnum    = 1, 
                                 .mtype   = PARDISO_COMPLEX_SYM,
                                 .sol_enum = PARDISO_LHS,
                                 .smat_p  = (void*)&M,
                                 .rhs_p   = (void*)&Pvb,
                                 .sol_p   = (void*)&V_vb};
    debug_body("%s", "SolverDATA initialized");
    matlib_complex *ptr;
    BEGIN_DEBUG
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
    END_DEBUG

    data.phase_enum = PARDISO_INIT;
    matlib_pardiso(&data);

    for (i=1; i<Nt+1; i++)
    {
        debug_print("begin iteration: %d", i);
        fem1d_ZPrjL2F(p, U_tmp, Pvb);

        data.phase_enum = PARDISO_ANALYSIS_AND_FACTOR;
        matlib_pardiso(&data);

        data.phase_enum = PARDISO_SOLVE_AND_REFINE;
        data.mnum = i;
        matlib_pardiso(&data);
        
        fem1d_ZF2L(p, V_vb, V_tmp);

        //for(j=0; j<U_tmp.len; j++)
        //{
        //    U_tmp.elem_p[j] = 2.0 * V_tmp.elem_p[j]-U_tmp.elem_p[j];
        //}

        /* 2.0 * V_tmp -U_tmp --> U_tmp*/ 
        matlib_zaxpby(2.0, V_tmp, -1.0, U_tmp );
        t += dt;

        //(*func_p)(x1, t, u_exact);
        //fem1d_ZILT(N, IMi, U_tmp, u_computed);

        /* Error analysis */ 
        //norm_actual = matlib_znrm2(u_exact);
        //matlib_zaxpy(-1.0, u_exact, u_computed);
        //e_relative.elem_p[i-1] = matlib_znrm2(u_computed)/norm_actual;
    }
    data.phase_enum = PARDISO_FREE;
    matlib_pardiso(&data);
    fem1d_ZNSparseMGMM(p, N, nsparse, Q, &phi, &q, &M, FEM1D_GMM_FREE);

    return(0);
    
    
}

//void test_solve_Schroedinger_IVP_2a(void)
//{
//    matlib_index p = 4;
//    matlib_index nr_LGL = 2*p+1;
//    matlib_index N = 800;
//    matlib_index Nt = 4;
//    double dt = 0.5e-3;
//    double domain[2] = {-15.0, 15.0 };
//    double e_relative;
//
//
//    e_relative = solve_Schroedinger_IVP2( p, nr_LGL, N, domain, dt, Nt,
//                                         solution_Gaussian,
//                                         linear_timedependent_zpotential);
//    CU_ASSERT_TRUE(e_relative<1e-3);
//
//}

/*============================================================================+/
| Example 3.:
| Nonlinear Schrodinger Eqation (IVP):
| iu_t + u_xx + chi |u|^2u= 0
| t=0: u(x,t) = u_0(x) -> Initial data
| 
/+============================================================================*/


void solution_BrightSolitonNLS
( 
    double    *params,
    matlib_dv x, 
    double    t, 
    matlib_zv u
)
/*                  e^(i*eta^2*t)   
 * u(x,t) = eta ---------------------
 *                  cosh[eta*(x)]
 *
 * */ 
{

    debug_enter("%s", "");
    
    double *xptr;
    double eta = 1.0;
    matlib_complex e = cexp(I*eta*eta*t);

    for (xptr = x.elem_p; xptr<(x.elem_p+x.len); xptr++)
    {
        *(u.elem_p) = eta*e/cosh(eta*(*xptr));
        u.elem_p++;
    }
    debug_exit("%s", "");
}



double solve_NLSEquation_IVP
(
    matlib_index p,
    matlib_index nr_LGL,
    matlib_index N,
    double domain[2],
    double dt,
    matlib_index Nt,
    void  (*func_p)(double*, matlib_dv, double, matlib_zv)
)
/* 
 * */ 
{
    clock_t start = clock(), duration;
    debug_enter( "polynomial degree: %d, nr. of LGL points: %d", p, nr_LGL );
    matlib_index i, j;
    matlib_index P = nr_LGL-1;

    double x_l = domain[0];
    double x_r = domain[1];

    matlib_dv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_dm FM, IM, Q;

    matlib_create_dm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_dm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    legendre_LGLdataIM( xi, IM);
    fem1d_quadM( quadW, IM, &Q);

    /* generate the grid */ 
    matlib_dv x;
    fem1d_ref2mesh (xi, N, x_l, x_r, &x);
    debug_body("length of x: %d", x.len);

    /* time discretization */ 

    double ta[2] = {0, Nt*dt};

    /* Provide the initial condition */
    matlib_zv u0, phi;
    matlib_create_zv( x.len, &u0, MATLIB_COL_VECT);
    matlib_create_zv( x.len, &phi, MATLIB_COL_VECT);
    
    double params[2] = {1.0, 2*M_PI};
    (*func_p)(params, x, ta[0], u0);
    //DEBUG_PRINT_ZV(u0, "%s :", "Initial condition");
    
    /* Assemble the global mass matrix */ 

    double irho = dt/2.0;
    double J   = 0.5*(x_r-x_l)/N;
    matlib_complex coeff = 0.0 + I*irho/(J*J);
    matlib_complex time_coeff[2] = {1.0, -I*irho};
    debug_body("coeff: %0.16f%+0.16fi", coeff);


    matlib_zm q;
    matlib_zsparsem M;
    for(i=0; i<phi.len; i++)
    {
        phi.elem_p[i] = 1.0;
    
    }
    fem1d_ZSparseMGMM(p, Q, phi, &M);
    
    /* Setup the sparse linear system */

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

    BEGIN_DEBUG
        debug_print("dimension of the sparse square matrix: %d", M.lenc);
        for(i=0; i<M.lenc; i++)
        {
            for(j=M.rowIn[i]; j<M.rowIn[i+1]; j++)
            {
                debug_print("M(%d,%d): % 0.16f %+0.16fi", i, M.colIn[j], M.elem_p[j]);
            }
        }
    END_DEBUG
    
    /* Setup the sparse linear system */


    /* Solve the equation while marching in time 
     *
     * Linear system M * V_vb = Pvb 
     * */
    
    matlib_zv U_tmp, V_tmp, V_vb, U_vb, Pvb, u_NL, U_NL, PNL_vb;
    double dim = N*(p+1);
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

    /* ONLY FOR ERROR: Compute the exact solution and make a forward LT */ 
    matlib_dv xii, quadWi, x1;
    legendre_LGLdataLT1( p, TOL, &xii, &quadWi);
    fem1d_ref2mesh (xii, N, x_l, x_r, &x1);
    matlib_dm IMi;
    matlib_create_dm( xii.len, p+1, &IMi, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    
    legendre_LGLdataIM( xii, IMi);
    matlib_zv u_exact, u_computed;
    matlib_create_zv( x1.len,    &u_exact, MATLIB_COL_VECT);
    matlib_create_zv( x1.len, &u_computed, MATLIB_COL_VECT);

    double t = 0, norm_actual;
    matlib_dv e_relative;
    matlib_create_dv( Nt, &e_relative, MATLIB_COL_VECT);

    pardiso_solver_t data = { .nsparse = 1, 
                                 .mnum    = 1, 
                                 .mtype   = PARDISO_COMPLEX_SYM,
                                 .sol_enum = PARDISO_LHS,
                                 .smat_p  = (void*)&M,
                                 .rhs_p   = (void*)&PNL_vb,
                                 .sol_p   = (void*)&V_vb};
    debug_body("%s", "SolverDATA initialized");

    data.phase_enum = PARDISO_INIT;
    matlib_pardiso(&data);

    double chi = 2.0, res;
    matlib_index Nth = 50, iter;

    duration = clock()-start;
    debug_print("Time of initialize[msec]: %0.6f", duration*1000.0/(double)CLOCKS_PER_SEC);

    start = clock();
    clock_t start_, duration_;
    data.phase_enum = PARDISO_ANALYSIS_AND_FACTOR;
    start_ = clock();
    matlib_pardiso(&data);
    duration_ = clock()-start_;
    debug_print( "Time for factorization[msec]: %0.6f", 
                 duration_*1000.0/(double)CLOCKS_PER_SEC);
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
        while((res>TOL) && (iter<Nth))
        {
            start_ = clock();
            fem1d_ZILT(N, IM, V_tmp, u_NL);
            for (j=0; j<u_NL.len; j++)
            {
                u_NL.elem_p[j] = I*irho*chi*u_NL.elem_p[j]*u_NL.elem_p[j]*conj(u_NL.elem_p[j]);
            }
            fem1d_ZFLT( N, FM, u_NL, U_NL);
            fem1d_ZPrjL2F(p, U_NL, PNL_vb);
            matlib_zaxpy(1.0, Pvb, PNL_vb);
            //DEBUG_PRINT_ZV(PNL_vb, "%s :", "projection");

            duration_ = clock()-start_;
            debug_print( "Time for other crap at %d-th time-step[msec]: %0.6f", 
                     i, duration_*1000.0/(double)CLOCKS_PER_SEC);
            data.phase_enum = PARDISO_SOLVE_AND_REFINE;
            data.mnum = i+1;
            start_ = clock();
            matlib_pardiso(&data);
            duration_ = clock()-start_;
            debug_print( "Time for %d-th time-step[msec]: %0.6f", 
                     i, duration_*1000.0/(double)CLOCKS_PER_SEC);
            
            matlib_zaxpby(1.0, V_tmp, 0, U_NL );

            fem1d_ZF2L(p, V_vb, V_tmp);
            matlib_zaxpy(-1.0, V_tmp, U_NL );
            res = sqrt(J) * fem1d_ZNorm2(p, N, U_NL);
            iter++;
            debug_body("time-step: %d, iteration: %d, residue: %0.16f", i, iter, res);
        }
        if((iter==Nth) || (isnan(res)))
        {
            term_exec("Iteration reached threshold with res: %0.16f", res);
        }

        /* 2.0 * V_tmp -U_tmp --> U_tmp*/ 
        matlib_zaxpby(2.0, V_tmp, -1.0, U_tmp );
        t += dt;

        (*func_p)(params, x1, t, u_exact);
        fem1d_ZILT(N, IMi, U_tmp, u_computed);

        /* Error analysis */ 
        norm_actual = matlib_znrm2(u_exact);
        matlib_zaxpy(-1.0, u_exact, u_computed);
        e_relative.elem_p[i] = matlib_znrm2(u_computed)/norm_actual;
        debug_body("relative error: %0.16f", e_relative.elem_p[i]);
    }
    assert(fabs(ta[1]-t)<TOL);
    data.phase_enum = PARDISO_FREE;
    matlib_pardiso(&data);
    duration = clock()-start;
    debug_print("Time for solution [sec]: %0.6f", duration/(double)CLOCKS_PER_SEC);

    //(*func_p)(params, x1, ta[2], u_exact);
    //fem1d_ZILT(N, IMi, U_tmp, u_computed);

    DEBUG_PRINT_DV(e_relative, "%s :", "relative error");

    double r = e_relative.elem_p[Nt-1];
    debug_exit("Relative error: % 0.16g", r);
    return(r);
    
}

void test_solve_NLSEquation_IVP(void)
{
    matlib_index p = 4;
    matlib_index nr_LGL = 2*p+1;
    matlib_index N = 2000;
    matlib_index Nt = 2000;
    double dt = 0.5e-3;
    double domain[2] = {-15.0, 15.0 };
    double e_relative;
    
    e_relative = solve_NLSEquation_IVP( p, nr_LGL, N, domain, dt, Nt,
                                       solution_BrightSolitonNLS);

    CU_ASSERT_TRUE(true);

}
/*============================================================================*/

typedef struct
{
    void** param; /* array of pointers */ 
} FEM1D_PTHREAD_DATA;

void* test_initialize_var(void* mp)
{
    FEM1D_PTHREAD_DATA *ptr = (FEM1D_PTHREAD_DATA*) mp;

    matlib_index     *p       = (matlib_index*)(ptr->param[0]);
    matlib_dv *xi_p    = (matlib_dv*)ptr->param[1];
    matlib_dv *quadW_p = (matlib_dv*)ptr->param[2];
    matlib_dm *IM_p    = (matlib_dm*)ptr->param[3];
    matlib_dm *Q_p     = (matlib_dm*)ptr->param[4];

    matlib_create_dm( xi_p->len, *p+1, IM_p, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataIM( *xi_p, *IM_p);
    fem1d_quadM( *quadW_p, *IM_p, Q_p);
    int status = 0;
    pthread_exit((void*)&status);
}
void* test_eval_nonlinear_term(void* mp)
{
    FEM1D_PTHREAD_DATA *ptr = (FEM1D_PTHREAD_DATA*) mp;

    matlib_index     *p       = (matlib_index*)(ptr->param[0]);
    matlib_dv *xi_p    = (matlib_dv*)ptr->param[1];
    matlib_dv *quadW_p = (matlib_dv*)ptr->param[2];
    matlib_dm *IM_p    = (matlib_dm*)ptr->param[3];
    matlib_dm *Q_p     = (matlib_dm*)ptr->param[4];

    matlib_create_dm( xi_p->len, *p+1, IM_p, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataIM( *xi_p, *IM_p);
    fem1d_quadM( *quadW_p, *IM_p, Q_p);
    int status = 0;
    pthread_exit((void*)&status);
}

double solve_NLSEquation_IVP_Pthread
(
    matlib_index p,
    matlib_index nr_LGL,
    matlib_index N,
    double domain[2],
    double dt,
    matlib_index Nt,
    void  (*func_p)(double*, matlib_dv, double, matlib_zv)
)
/* 
 * */ 
{
    pthread_t mythread[NUM_THREADS];
    pthread_attr_t attr;
    /* initialize and set thread detached attribute */
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    
    int pthread_r;

    clock_t start = clock(), duration;
    debug_enter( "polynomial degree: %d, nr. of LGL points: %d", p, nr_LGL );
    matlib_index i, j;
    matlib_index P = nr_LGL-1;

    double x_l = domain[0];
    double x_r = domain[1];

    matlib_dv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_dm FM, IM, Q;

    void* param[5] = { (void*) &p, 
                       (void*) &xi, 
                       (void*) &quadW, 
                       (void*) &IM,
                       (void*) &Q     };

    FEM1D_PTHREAD_DATA ptr = {.param = param};
    
    pthread_r = pthread_create(&mythread[0], &attr, test_initialize_var, (void *) &ptr); 
    if(pthread_r)
    {
        term_exec("pthread failed to create; return value: %d", pthread_r);
    }


    matlib_create_dm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    //matlib_create_dm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    //legendre_LGLdataIM( xi, IM);
    //fem1d_quadM( quadW, IM, &Q);

    /* generate the grid */ 
    matlib_dv x;
    fem1d_ref2mesh (xi, N, x_l, x_r, &x);
    debug_body("length of x: %d", x.len);

    /* time discretization */ 

    double ta[2] = {0, Nt*dt};

    /* Provide the initial condition */
    matlib_zv u0, phi;
    matlib_create_zv( x.len, &u0, MATLIB_COL_VECT);
    matlib_create_zv( x.len, &phi, MATLIB_COL_VECT);
    
    double params[2] = {1.0, 2*M_PI};
    (*func_p)(params, x, ta[0], u0);
    //DEBUG_PRINT_ZV(u0, "%s :", "Initial condition");
    

    /* ONLY FOR ERROR ANALYSIS: Compute the exact solution and make a forward LT */ 
    matlib_dv xii, quadWi, x1;
    legendre_LGLdataLT1( p, TOL, &xii, &quadWi);
    fem1d_ref2mesh (xii, N, x_l, x_r, &x1);
    matlib_dm IMi;
    matlib_create_dm( xii.len, p+1, &IMi, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    
    legendre_LGLdataIM( xii, IMi);
    matlib_zv u_exact, u_computed;
    matlib_create_zv( x1.len,    &u_exact, MATLIB_COL_VECT);
    matlib_create_zv( x1.len, &u_computed, MATLIB_COL_VECT);

    /* Assemble the global mass matrix */ 

    double irho = dt/2.0;
    double J   = 0.5*(x_r-x_l)/N;
    matlib_complex coeff = 0.0 + I*irho/(J*J);
    matlib_complex time_coeff[2] = {1.0, -I*irho};
    debug_body("coeff: %0.16f%+0.16fi", coeff);


    matlib_zsparsem M;
    for(i=0; i<phi.len; i++)
    {
        phi.elem_p[i] = 1.0;
    
    }

    matlib_zv U_tmp, V_tmp, u_NL, U_NL;
    double dim = N*(p+1);
    matlib_create_zv( x.len,  &u_NL, MATLIB_COL_VECT);
    matlib_create_zv(   dim,  &U_NL, MATLIB_COL_VECT);
    matlib_create_zv(   dim, &V_tmp, MATLIB_COL_VECT);
    matlib_create_zv(   dim, &U_tmp, MATLIB_COL_VECT);

    matlib_dv e_relative;
    matlib_create_dv( Nt, &e_relative, MATLIB_COL_VECT);

    fem1d_ZFLT( N, FM, u0, U_tmp);
    for(i=0; i<V_tmp.len; i++)
    {
        V_tmp.elem_p[i] = 0;
    }
    
    /* join the thread */ 
    void* status;
    pthread_r = pthread_join(mythread[0], &status); 

    fem1d_ZSparseMGMM(p, Q, phi, &M);
    
    /* Setup the sparse linear system */

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

    BEGIN_DEBUG
        debug_print("dimension of the sparse square matrix: %d", M.lenc);
        for(i=0; i<M.lenc; i++)
        {
            for(j=M.rowIn[i]; j<M.rowIn[i+1]; j++)
            {
                debug_print("M(%d,%d): % 0.16f %+0.16fi", i, M.colIn[j], M.elem_p[j]);
            }
        }
    END_DEBUG
    
    /* Setup the sparse linear system */


    /* Solve the equation while marching in time 
     *
     * Linear system M * V_vb = Pvb 
     * */
    
    matlib_zv V_vb, U_vb, Pvb, PNL_vb;
    matlib_create_zv( M.lenc,    &Pvb, MATLIB_COL_VECT);
    matlib_create_zv( M.lenc,   &V_vb, MATLIB_COL_VECT);
    matlib_create_zv( M.lenc, &PNL_vb, MATLIB_COL_VECT);

    double t = 0, norm_actual;

    pardiso_solver_t data = { .nsparse = 1, 
                                 .mnum    = 1, 
                                 .mtype   = PARDISO_COMPLEX_SYM,
                                 .sol_enum = PARDISO_LHS,
                                 .smat_p  = (void*)&M,
                                 .rhs_p   = (void*)&PNL_vb,
                                 .sol_p   = (void*)&V_vb};
    debug_body("%s", "SolverDATA initialized");

    data.phase_enum = PARDISO_INIT;
    matlib_pardiso(&data);

    double chi = 2.0, res;
    matlib_index Nth = 50, iter;
    duration = clock()-start;
    debug_print("Time of initialize[msec]: %0.6f", duration*1000.0/(double)CLOCKS_PER_SEC);


    start = clock();
    clock_t start_, duration_;
    data.phase_enum = PARDISO_ANALYSIS_AND_FACTOR;
    start_ = clock();
    matlib_pardiso(&data);
    duration_ = clock()-start_;
    debug_print( "Time for factorization[msec]: %0.6f", 
                 duration_*1000.0/(double)CLOCKS_PER_SEC);
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
        while((res>TOL) && (iter<Nth))
        {
            start_ = clock();
            fem1d_ZILT(N, IM, V_tmp, u_NL);
            for (j=0; j<u_NL.len; j++)
            {
                u_NL.elem_p[j] = I*irho*chi*u_NL.elem_p[j]*u_NL.elem_p[j]*conj(u_NL.elem_p[j]);
            }
            fem1d_ZFLT( N, FM, u_NL, U_NL);
            fem1d_ZPrjL2F(p, U_NL, PNL_vb);
            matlib_zaxpy(1.0, Pvb, PNL_vb);
            //DEBUG_PRINT_ZV(PNL_vb, "%s :", "projection");

            duration_ = clock()-start_;
            debug_print( "Time for other crap at %d-th time-step[msec]: %0.6f", 
                     i, duration_*1000.0/(double)CLOCKS_PER_SEC);
            data.phase_enum = PARDISO_SOLVE_AND_REFINE;
            data.mnum = i+1;
            start_ = clock();
            matlib_pardiso(&data);
            duration_ = clock()-start_;
            debug_print( "Time for %d-th time-step[msec]: %0.6f", 
                     i, duration_*1000.0/(double)CLOCKS_PER_SEC);
            
            matlib_zaxpby(1.0, V_tmp, 0, U_NL );

            fem1d_ZF2L(p, V_vb, V_tmp);
            matlib_zaxpy(-1.0, V_tmp, U_NL );
            res = sqrt(J) * fem1d_ZNorm2(p, N, U_NL);
            iter++;
            debug_body("time-step: %d, iteration: %d, residue: %0.16f", i, iter, res);
        }
        if((iter==Nth) || (isnan(res)))
        {
            term_exec("Iteration reached threshold with res: %0.16f", res);
        }

        /* 2.0 * V_tmp -U_tmp --> U_tmp*/ 
        matlib_zaxpby(2.0, V_tmp, -1.0, U_tmp );
        t += dt;

        (*func_p)(params, x1, t, u_exact);
        fem1d_ZILT(N, IMi, U_tmp, u_computed);

        /* Error analysis */ 
        norm_actual = matlib_znrm2(u_exact);
        matlib_zaxpy(-1.0, u_exact, u_computed);
        e_relative.elem_p[i] = matlib_znrm2(u_computed)/norm_actual;
        debug_body("relative error: %0.16f", e_relative.elem_p[i]);
    }
    assert(fabs(ta[1]-t)<TOL);
    data.phase_enum = PARDISO_FREE;
    matlib_pardiso(&data);
    duration = clock()-start;
    debug_print("Time for solution [sec]: %0.6f", duration/(double)CLOCKS_PER_SEC);

    //(*func_p)(params, x1, ta[2], u_exact);
    //fem1d_ZILT(N, IMi, U_tmp, u_computed);

    DEBUG_PRINT_DV(e_relative, "%s :", "relative error");

    double r = e_relative.elem_p[Nt-1];

    
    debug_exit("Relative error: % 0.16g", r);
    
    return(r);
    
}

void test_solve_NLSEquation_IVP_Pthread(void)
{
    matlib_index p = 4;
    matlib_index nr_LGL = 2*p+1;
    matlib_index N = 1200;
    matlib_index Nt = 2000;
    double dt = 0.5e-3;
    double domain[2] = {-15.0, 15.0 };
    double e_relative;
    
    e_relative = solve_NLSEquation_IVP_Pthread( p, nr_LGL, N, domain, dt, Nt,
                                       solution_BrightSolitonNLS);

    CU_ASSERT_TRUE(true);

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
void solution_BrightSolitonGPE
( 
    double    *params,
    matlib_dv x, 
    double    t, 
    matlib_zv u
)
/*                  e^(i*eta^2*t)   
 * u(x,t) = eta ---------------------e^(i*Xi)
 *               cosh[eta*(x-nu)]
 *
 * beta   = g_0 *sin(mu*t)/mu 
 * nu     = (2*g_0*cos(mu*t))/mu^2 - (2*g_0)/mu^2;
 * ibeta2 = (mu*(g_0^2*t)/2 - (g_0^2*sin(2*mu*t))/4;
 * Xi     = beta*x-ibeta2;
 *
 * */ 
{

    debug_enter("%s", "");
    
    double *xptr;
    /* parameters 
     * param[0] = g_0, 
     * param[1] = mu
     * */ 
    double eta = 1, Xi;
    
    matlib_complex e = cexp(I*eta*eta*t);
    double theta     = params[1]*t;
    double beta      = params[0]*sin(theta)/params[1];
    double nu        = 2.0*params[0]*(cos(theta)-1.0)/(params[1]*params[1]);
    double ibeta2    = params[0]*params[0]*(theta/2.0 - sin(2*theta)/4.0);

    for (xptr = x.elem_p; xptr<(x.elem_p+x.len); xptr++)
    {
        Xi = beta**xptr - ibeta2;
        *(u.elem_p) = eta*e*exp(I*Xi)/cosh(eta*(*xptr+nu));
        u.elem_p++;
    }
    debug_exit("%s", "");
}



double solve_GPEquation_IVP
(
    matlib_index p,
    matlib_index nr_LGL,
    matlib_index N,
    double domain[2],
    double dt,
    matlib_index Nt,
    void  (*func_p)(double*, matlib_dv, double, matlib_zv),
    void  (*potential_p)(matlib_complex*, matlib_dv, double, double, matlib_zm)
)
/* 
 * */ 
{

    debug_enter( "polynomial degree: %d, nr. of LGL points: %d", p, nr_LGL );
    matlib_index i, j;
    matlib_index P = nr_LGL-1;

    double x_l = domain[0];
    double x_r = domain[1];

    matlib_dv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_dm FM, IM, Q;
    matlib_create_dm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_dm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    legendre_LGLdataIM( xi, IM);
    fem1d_quadM( quadW, IM, &Q);

    /* generate the grid */ 
    matlib_dv x;
    fem1d_ref2mesh (xi, N, x_l, x_r, &x);
    debug_body("length of x: %d", x.len);

    /* time discretization */ 

    double ta[2] = {0, Nt*dt};

    /* Provide the initial condition */
    matlib_zv u0;
    matlib_create_zv( x.len, &u0, MATLIB_COL_VECT);
    
    double params[2] = {1.0, 2*M_PI};
    (*func_p)(params, x, ta[0], u0);
    //DEBUG_PRINT_ZV(u0, "%s :", "Initial condition");
    
    /* Assemble the global mass matrix */ 

    double irho = dt/2.0;
    double J   = 0.5*(x_r-x_l)/N;
    matlib_complex coeff = 0.0 + I*irho/(J*J);
    matlib_complex time_coeff[2] = {1.0, -I*irho};
    debug_body("coeff: %0.16f%+0.16fi", coeff);


    matlib_zm phi, q;
    matlib_znsparsem M;
    matlib_index nsparse = 20;
    fem1d_ZNSparseMGMM(p, N, nsparse, Q, &phi, &q, &M, FEM1D_GMM_INIT);
    
    /* Setup the sparse linear system */


    /* Solve the equation while marching in time 
     *
     * Linear system M * V_vb = Pvb 
     * */
    
    matlib_zv U_tmp, V_tmp, V_vb, U_vb, Pvb, u_NL, U_NL, PNL_vb;
    double dim = N*(p+1);
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

    /* ONLY FOR ERROR: Compute the exact solution and make a forward LT */ 
    matlib_dv xii, quadWi, x1;
    legendre_LGLdataLT1( p, TOL, &xii, &quadWi);
    fem1d_ref2mesh (xii, N, x_l, x_r, &x1);
    matlib_dm IMi;
    matlib_create_dm( xii.len, p+1, &IMi, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    
    legendre_LGLdataIM( xii, IMi);
    matlib_zv u_exact, u_computed;
    matlib_create_zv( x1.len,    &u_exact, MATLIB_COL_VECT);
    matlib_create_zv( x1.len, &u_computed, MATLIB_COL_VECT);

    double t = 0, norm_actual;
    matlib_dv e_relative;
    matlib_create_dv( Nt, &e_relative, MATLIB_COL_VECT);

    pardiso_solver_t data = { .nsparse = nsparse, 
                                 .mnum    = 1, 
                                 .mtype   = PARDISO_COMPLEX_SYM,
                                 .sol_enum = PARDISO_LHS,
                                 .smat_p  = (void*)&M,
                                 .rhs_p   = (void*)&PNL_vb,
                                 .sol_p   = (void*)&V_vb};
    debug_body("%s", "SolverDATA initialized");

    data.phase_enum = PARDISO_INIT;
    matlib_pardiso(&data);

    double chi = 2.0, res;
    matlib_index Nth = 50, iter;
    matlib_index Nt_ = Nt/nsparse;
    debug_body("Nt_: %d", Nt_);

    for(matlib_index k=0; k<Nt_; k++ )
    {
        debug_print("Compute the Global Mass-Stiffness Matrix: %d", k);
        (*potential_p)(time_coeff, x, t, dt, phi);
        fem1d_ZNSparseMGMM(p, N, nsparse, Q, &phi, &q, &M, FEM1D_GET_SPARSITY_NZE);
        fem1d_ZNSparseMGSM(N, coeff, M);

        for (i=0; i<nsparse; i++)
        {
            debug_print("time-step: %d", i+k*nsparse);
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
                    u_NL.elem_p[j] = I*irho*chi*u_NL.elem_p[j]*u_NL.elem_p[j]*conj(u_NL.elem_p[j]);
                }
                fem1d_ZFLT( N, FM, u_NL, U_NL);
                fem1d_ZPrjL2F(p, U_NL, PNL_vb);
                matlib_zaxpy(1.0, Pvb, PNL_vb);
                //DEBUG_PRINT_ZV(PNL_vb, "%s :", "projection");

                data.phase_enum = PARDISO_SOLVE_AND_REFINE;
                data.mnum = i+1;
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
            t += dt;

            (*func_p)(params, x1, t, u_exact);
            fem1d_ZILT(N, IMi, U_tmp, u_computed);

            /* Error analysis */ 
            norm_actual = matlib_znrm2(u_exact);
            matlib_zaxpy(-1.0, u_exact, u_computed);
            e_relative.elem_p[k*nsparse+i] = matlib_znrm2(u_computed)/norm_actual;
            debug_body("relative error: %0.16f", e_relative.elem_p[k*nsparse+i]);
        }
    }
    assert(fabs(ta[1]-t)<TOL);
    data.phase_enum = PARDISO_FREE;
    matlib_pardiso(&data);
    fem1d_ZNSparseMGMM(p, N, nsparse, Q, &phi, &q, &M, FEM1D_GMM_FREE);

    (*func_p)(params, x1, ta[2], u_exact);
    fem1d_ZILT(N, IMi, U_tmp, u_computed);

    DEBUG_PRINT_DV(e_relative, "%s :", "relative error");

    double r = e_relative.elem_p[Nt-1];
    debug_exit("Relative error: % 0.16g", r);
    return(r);
    
}

void test_solve_GPEquation_IVP(void)
{
    matlib_index p = 4;
    matlib_index nr_LGL = 2*p+1;
    matlib_index N = 1000;
    matlib_index Nt = 2000;
    double dt = 0.5e-3;
    double domain[2] = {-15.0, 15.0 };
    double e_relative;
    
    e_relative = solve_GPEquation_IVP( p, nr_LGL, N, domain, dt, Nt,
                                       solution_BrightSolitonGPE,
                                       linear_timedependent_zpotential);

    CU_ASSERT_TRUE(true);

}





/*============================================================================+/
 | Test runner
 |
 |
 +============================================================================*/

int main(void)
{
    mkl_domain_set_num_threads(1, MKL_BLAS);    
    CU_pSuite pSuite = NULL;

    /* initialize the CUnit test registry */
    if (CUE_SUCCESS != CU_initialize_registry())
    {
        return CU_get_error();
    }

    /* Create a test array */
    CU_TestInfo test_array[] = 
    {
        { "Schroedinger Equation, const. potential"  , test_solve_Schroedinger_IVP_1a},
        //{ "Schroedinger Equation, harmonic potential", test_solve_Schroedinger_IVP_1b},
        //{ "Schroedinger Equation, linear time-dependent potential", test_solve_Schroedinger_IVP_2a},
        //{ "GP Equation, linear time-dependent potential", test_solve_GPEquation_IVP},
        //{ "Pthreded NLS Equation, linear time-dependent potential", test_solve_NLSEquation_IVP_Pthread},
        //{ "NLS Equation", test_solve_NLSEquation_IVP},
        CU_TEST_INFO_NULL,
    };

    /* Create the test suite */ 
    CU_SuiteInfo suites[] = 
    {
        { "Linear solvers", init_suite, clean_suite, test_array },
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

