/*============================================================================+/
 | File: Cunit_matlib_solver.c 
 | Description: Test 
 |
 |
/+============================================================================*/
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>

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

/*============================================================================*/

static const matlib_real TOL = 1e-9;

int init_suite(void)
{
      return 0;
}
int clean_suite(void)
{
      return 0;
}

/*============================================================================*/
void Gaussian_xfunc
( 
    const matlib_xv x, 
          matlib_xv y
)
{
    matlib_index i;
    if(y.len >= x.len)
    {
        for (i=0; i<x.len; i++)
        {
            *(y.elem_p) = exp(-0.5*(*(x.elem_p+i)**(x.elem_p+i)));
            y.elem_p++;
        }
    }
}

void Gaussian_zfunc
( 
    const matlib_xv x, 
          matlib_zv z
)
{
    matlib_index i;
    if(z.len >= x.len)
    {
        for (i=0; i<x.len; i++)
        {
            *(z.elem_p) = cexp(-(0.5+I*0.8)*(*(x.elem_p+i)**(x.elem_p+i)));
            z.elem_p++;
        }
    }
}

/*============================================================================*/

matlib_real test_matlib_xsolver_general
(
    matlib_index p,
    matlib_index nr_LGL,
    matlib_index N,
    matlib_real domain[2],
    void  (*func_p)(matlib_xv, matlib_xv), 
    void  (*potential_p)(matlib_xv, matlib_xv)
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

    matlib_xv u, U;
    matlib_create_xv( x.len,    &u, MATLIB_COL_VECT);
    matlib_create_xv( N*(p+1),  &U, MATLIB_COL_VECT);

    (*func_p)(x, u);
    fem1d_XFLT( N, FM, u, U);


    /* Aseemble the global mass matrix */ 
    matlib_xv phi, Phi;
    matlib_create_xv(   x.len, &phi, MATLIB_COL_VECT);
    matlib_create_xv( N*(p+1), &Phi, MATLIB_COL_VECT);

    (*potential_p)(x, phi);
    
    matlib_xm_sparse M;
    fem1d_xm_sparse_GMM(p, Q, phi, &M);

    /* Transform it to FEM-basis: vertex and bubble basis functions */  
    matlib_xv Uvb;
    matlib_create_xv( M.lenc, &Uvb, MATLIB_COL_VECT);
    fem1d_XL2F(p, U, Uvb);
    /* 
     * Finding Uvb by solving the linear system 
     * M * Uvb1 = P_vb where Uvb1 is the unknown and 
     * P_vb is the projection of phi(x)u(x).
     *
     * */ 
    matlib_xv Uvb1, P_vb;
    matlib_create_xv( M.lenc, &Uvb1, MATLIB_COL_VECT);
    matlib_create_xv( M.lenc, &P_vb, MATLIB_COL_VECT);
    for(i=0; i<x.len; i++)
    {
        phi.elem_p[i] = u.elem_p[i] * phi.elem_p[i];
    }

    /* compute projection from the LP representation of phi(x)*u(x)*/
    fem1d_XFLT( N, FM, phi, Phi);
    fem1d_XPrjL2F(p, Phi, P_vb);

    /* Initialize solver data */ 
    pardiso_solver_t data = { .nsparse  = 1, 
                              .mnum     = 1, 
                              .mtype    = PARDISO_REAL_SYM_PDEF,
                              .sol_enum = PARDISO_LHS, 
                              .smat_p   = (void*)&M,
                              .rhs_p    = (void*)&P_vb,
                              .sol_p    = (void*)&Uvb1};

    debug_body("%s", "Solver data initialized");
    data.phase_enum = PARDISO_INIT;
    matlib_pardiso(&data);

    data.phase_enum = PARDISO_ANALYSIS_AND_FACTOR;
    matlib_pardiso(&data);

    data.phase_enum = PARDISO_SOLVE_AND_REFINE;
    matlib_pardiso(&data);

    /* Do the error analysis */ 
    BEGIN_DTRACE
        for(i=0;i<Uvb.len; i++)
        {
            debug_print( "[%d]-> Uvb: % 0.16f, Uvb1: % 0.16f", 
                         i, *(Uvb.elem_p+i), *(Uvb1.elem_p+i));
        }
    END_DTRACE

    matlib_real norm_actual = matlib_xnrm2(Uvb);
    matlib_xaxpy(-1.0, Uvb1, Uvb);
    matlib_real e_relative = matlib_xnrm2(Uvb)/norm_actual;
    debug_exit("Relative error: % 0.16g", e_relative);

    data.phase_enum = PARDISO_FREE;
    matlib_pardiso(&data);

    return(e_relative);
}

void constant_xpotential( const matlib_xv x, matlib_xv y)
{
    matlib_index i;
    if(y.len >= x.len)
    {
        for (i=0; i<x.len; i++)
        {
            *(y.elem_p) = 1.0;
            y.elem_p++;
        }
    }
}

void harmonic_xpotential( const matlib_xv x, matlib_xv y)
{
    matlib_index i;
    if(y.len >= x.len)
    {
        for (i=0; i<x.len; i++)
        {
            *(y.elem_p) = 0.5*(*(x.elem_p+i)**(x.elem_p+i));
            y.elem_p++;
        }
    }
}

void test_matlib_xsolver(void)
{
    matlib_index p, nr_LGL;
    matlib_index N = 1000;
    matlib_real domain[2] = {-5.0, 5.0};
    matlib_real e_relative;
    p = 2;
    nr_LGL = 4*p+1;
    e_relative = test_matlib_xsolver_general( p, nr_LGL, N, domain,
                                              Gaussian_xfunc, 
                                              constant_xpotential);
    CU_ASSERT_TRUE(e_relative<TOL);

    N = 800;
    for( p=3; p<13; p++)
    {
        nr_LGL = 3*p+1;
        e_relative = test_matlib_xsolver_general( p, nr_LGL, N, domain,
                                                  Gaussian_xfunc, 
                                                  constant_xpotential);
        CU_ASSERT_TRUE(e_relative<TOL);

        e_relative = test_matlib_xsolver_general( p, nr_LGL, N, domain,
                                                  Gaussian_xfunc, 
                                                  harmonic_xpotential);
        CU_ASSERT_TRUE(e_relative<TOL);
    } 
}
/*============================================================================*/

matlib_real test_matlib_zsolver_general
(
    matlib_index p,
    matlib_index nr_LGL,
    matlib_index N,
    matlib_real domain[2],
    void  (*func_p)(matlib_xv, matlib_zv), 
    void  (*potential_p)(matlib_xv, matlib_zv)
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

    matlib_zv u, U;
    matlib_create_zv( x.len,    &u, MATLIB_COL_VECT);
    matlib_create_zv( N*(p+1),  &U, MATLIB_COL_VECT);

    (*func_p)(x, u);
    fem1d_ZFLT( N, FM, u, U);


    /* Aseemble the global mass matrix */ 
    matlib_zv phi, Phi;
    matlib_create_zv(   x.len, &phi, MATLIB_COL_VECT);
    matlib_create_zv( N*(p+1), &Phi, MATLIB_COL_VECT);

    (*potential_p)(x, phi);
    
    matlib_zm_sparse M;
    fem1d_zm_sparse_GMM(p, Q, phi, &M);

    /* Transform it to FEM-basis: vertex and bubble basis functions */  
    matlib_zv Uvb;
    matlib_create_zv( M.lenc, &Uvb, MATLIB_COL_VECT);
    fem1d_ZL2F(p, U, Uvb);
    /* 
     * Finding Uvb by solving the linear system 
     * M * Uvb1 = P_vb where Uvb1 is the unknown and 
     * P_vb is the projection of phi(x)u(x).
     *
     * */ 
    matlib_zv Uvb1, P_vb;
    matlib_create_zv( M.lenc, &Uvb1, MATLIB_COL_VECT);
    matlib_create_zv( M.lenc, &P_vb, MATLIB_COL_VECT);
    for(i=0; i<x.len; i++)
    {
        phi.elem_p[i] = u.elem_p[i] * phi.elem_p[i];
    }

    /* compute projection from the LP representation of phi(x)*u(x)*/
    fem1d_ZFLT( N, FM, phi, Phi);
    fem1d_ZPrjL2F(p, Phi, P_vb);

    /* Initialize solver data */ 
    pardiso_solver_t data = { .nsparse  = 1, 
                              .mnum     = 1, 
                              .mtype    = PARDISO_COMPLEX_SYM,
                              .sol_enum = PARDISO_LHS, 
                              .smat_p   = (void*)&M,
                              .rhs_p    = (void*)&P_vb,
                              .sol_p    = (void*)&Uvb1};

    debug_body("%s", "Solver data initialized");
    data.phase_enum = PARDISO_INIT;
    matlib_pardiso(&data);

    data.phase_enum = PARDISO_ANALYSIS_AND_FACTOR;
    matlib_pardiso(&data);

    data.phase_enum = PARDISO_SOLVE_AND_REFINE;
    matlib_pardiso(&data);

    /* Do the error analysis */ 
    BEGIN_DTRACE
        for(i=0;i<Uvb.len; i++)
        {
            debug_print( "[%d]-> Uvb: % 0.16f%+0.16f, Uvb1: % 0.16f%+0.16f", 
                         i, *(Uvb.elem_p+i), *(Uvb1.elem_p+i));
        }
    END_DTRACE

    matlib_real norm_actual = matlib_znrm2(Uvb);
    matlib_zaxpy(-1.0, Uvb1, Uvb);
    matlib_real e_relative = matlib_znrm2(Uvb)/norm_actual;
    debug_exit("Relative error: % 0.16g", e_relative);

    data.phase_enum = PARDISO_FREE;
    matlib_pardiso(&data);

    return(e_relative);
}

void constant_zpotential( const matlib_xv x, matlib_zv y)
{
    matlib_index i;
    if(y.len >= x.len)
    {
        for (i=0; i<x.len; i++)
        {
            *(y.elem_p) = 0.5+I*0.5;
            y.elem_p++;
        }
    }
}

void harmonic_zpotential( const matlib_xv x, matlib_zv y)
{
    matlib_index i;
    if(y.len >= x.len)
    {
        for (i=0; i<x.len; i++)
        {
            *(y.elem_p) = (0.25+I*0.25)*(*(x.elem_p+i)**(x.elem_p+i));
            y.elem_p++;
        }
    }
}
void test_matlib_zsolver(void)
{
    matlib_index p, nr_LGL;
    matlib_index N = 2000;
    matlib_real domain[2] = {-5.0, 5.0};
    matlib_real e_relative;

    p = 2;
    nr_LGL = 2*p+1;
    e_relative = test_matlib_zsolver_general( p, nr_LGL, N, domain,
                                              Gaussian_zfunc, 
                                              constant_zpotential);
    CU_ASSERT_TRUE(e_relative<TOL);

    e_relative = test_matlib_zsolver_general( p, nr_LGL, N, domain,
                                              Gaussian_zfunc, 
                                              harmonic_zpotential);
    CU_ASSERT_TRUE(e_relative<TOL);
    
    N = 1000;
    for(matlib_index p=3; p<13; p++)
    {
        nr_LGL = 3*p+1;
        e_relative = test_matlib_zsolver_general( p, nr_LGL, N, domain,
                                                  Gaussian_zfunc, 
                                                  constant_zpotential);
        CU_ASSERT_TRUE(e_relative<TOL);

        e_relative = test_matlib_zsolver_general( p, nr_LGL, N, domain,
                                                  Gaussian_zfunc, 
                                                  harmonic_zpotential);
        CU_ASSERT_TRUE(e_relative<TOL);

    } 
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
        { "Solve real linear system"   , test_matlib_xsolver},
        { "Solve complex linear system", test_matlib_zsolver},
        CU_TEST_INFO_NULL,
    };

    /* Create the test suite */ 
    CU_SuiteInfo suites[] = 
    {
        { "PARDISO Linear solver", init_suite, clean_suite, NULL, NULL, test_array },
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

