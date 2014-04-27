/*============================================================================+/
 | File: Cunit_solver.c 
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

//#define NDEBUG

#include "legendre.h"
#include "fem1d.h"
#include "assert.h"

/* CUnit modules */
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>
#include <CUnit/TestDB.h>

/*============================================================================*/

static const double TOL = 1e-9;

int init_suite(void)
{
      return 0;
}
int clean_suite(void)
{
      return 0;
}

/*============================================================================*/
void poly_func
( 
    const matlib_dv x, 
          matlib_dv y
)
{
    int deg = 10;
    matlib_index i;
    if(y.len >= x.len)
    {
        for (i=0; i<x.len; i++)
        {
            *(y.elem_p) = pow(*(x.elem_p+i), deg);
            y.elem_p++;
        }
    }
}
/*============================================================================*/

void Gaussian_func
( 
    const matlib_dv x, 
          matlib_dv y
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

void Gaussian_func_t
( 
    const matlib_dv x, 
    const matlib_dv t, 
          matlib_dm y
)
{
    matlib_index i, j;
    bool dim_OK = (y.lenc == x.len) && (t.len == y.lenr);

    if(dim_OK)
    {
        for(j=0;j<t.len; j++)
        {
            for (i=0; i<x.len; i++)
            {
                *(y.elem_p) = exp(-0.5*(*(x.elem_p+i)**(x.elem_p+i))/(1+*(t.elem_p+j)));
                y.elem_p++;
            }
        }
    }
}

void Gaussian_zfunc
( 
    const matlib_dv x, 
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

void Gaussian_zfunc_t
( 
    const matlib_dv x, 
    const matlib_dv t, 
          matlib_zm y
)
{
    matlib_index i, j;
    bool dim_OK = (y.lenc == x.len) && (t.len == y.lenr);

    if(dim_OK)
    {
        for(j=0;j<t.len; j++)
        {
            for (i=0; i<x.len; i++)
            {
                *(y.elem_p) = cexp(-(0.5+0.8*I)*(*(x.elem_p+i)**(x.elem_p+i))/(1+*(t.elem_p+j)));
                y.elem_p++;
            }
        }
    }
}
/*============================================================================+/
 | Linear solver based on paradiso
/+============================================================================*/
void test_solver_dsym
(
    matlib_dsparsem A,
    matlib_dv       b,
    matlib_dv       x
)
{
    debug_enter("%s", "");
    matlib_index i;
    int mtype = -2; /* Matrix type: real symmetric */ 
    int nrhs  =  1; /* Number of right hand sides  */ 

    /* SETUP PARADISO CONTROL PARAMETERS */
    void* ptr[64];
    int iparam[64];
    for (i = 0; i < 64; i++)
    {
        iparam[i] = 0;
        ptr[i] = 0;
    }
    iparam[0] = 1; /* Don't use default values */ 
    iparam[1] = 2; /* Fill-in reducing odering for input matrix */ 
    iparam[3] = 0; /* Preconditioning */ 
    iparam[4] = 0; /*  */ 
    iparam[5] = 0; /* Write solution into x */ 
    iparam[7] = 2; /* Maximum number of iterative refinement steps, output 
                      reported in iparam[6] */ 
    iparam[9] = 13; /* Perturbing pivot elements */ 
    iparam[10] = 0; /* Disable scaling  */ 
    iparam[12] = 0; /*  */ 
    iparam[17] = -1; /* Enable reporting of nnz  */ 
    iparam[18] =  1; /*  */ 
    iparam[34] =  1; /* Zero-based indexing  */ 

    int maxfct, mnum, phase, error, msglvl;
    maxfct = 1; /* Maximum number of numerical factorizations. */
    mnum   = 1; /* Which factorization to use. */
    msglvl = 0; /* Print statistical information in file */
    error  = 0; /* Initialize error flag */
    
    /* Reordering and Symbolic Factorization. This step also allocates all 
     * memory that is necessary for the factorization. 
     * */
    phase = 11; /* Analysis */
    double ddummy;
    int idummy;
    debug_body("%s", "Start testing PARADISO");

    PARDISO ( ptr, &maxfct, &mnum, &mtype, &phase,
              &A.lenc, A.elem_p, A.rowIn, A.colIn, NULL, 
              &nrhs, iparam, &msglvl, NULL, NULL, &error);
    if (error != 0)
    {
      term_exec("ERROR during symbolic factorization: %d", error);
    }
    debug_body("%s", "Analysis completed");
    
    /* Numerical factorization. 
     * */
    phase = 22;
    PARDISO ( ptr, &maxfct, &mnum, &mtype, &phase,
              &A.lenc, A.elem_p, A.rowIn, A.colIn, NULL, 
              &nrhs, iparam, &msglvl, NULL, NULL, &error);
    if (error != 0)
    {
      term_exec("ERROR during symbolic factorization: %d", error);
    }
    debug_body("%s", "Factorization completed");

    /* Back substitution and iterative refinement.
     * */ 
    phase = 33;
    PARDISO ( ptr, &maxfct, &mnum, &mtype, &phase,
              &A.lenc, A.elem_p, A.rowIn, A.colIn, NULL, 
              &nrhs, iparam, &msglvl, b.elem_p, x.elem_p, &error);
    if (error != 0)
    {
      term_exec("ERROR during back substitution: %d", error);
    }
    /* Termination and release of memory 
     * */ 
    phase = -1;
    PARDISO ( ptr, &maxfct, &mnum, &mtype, &phase,
              &A.lenc, NULL, A.rowIn, A.colIn, NULL, 
              NULL, iparam, &msglvl, NULL, NULL, &error);
    
    debug_exit("%s", "");

}

/*============================================================================*/

void test_solver_zsym
(
    matlib_zsparsem A,
    matlib_zv       b,
    matlib_zv       x
)
{
    debug_enter("%s", "");
    matlib_index i;
    int mtype = 6; /* Matrix type: real symmetric */ 
    int nrhs  = 1; /* Number of right hand sides  */ 

    /* SETUP PARADISO CONTROL PARAMETERS */
    void* ptr[64];
    int iparam[64];
    for (i = 0; i < 64; i++)
    {
        iparam[i] = 0;
        ptr[i] = 0;
    }
    iparam[0] = 1; /* Don't use default values */ 
    iparam[1] = 2; /* Fill-in reducing odering for input matrix */ 
    iparam[3] = 0; /* Preconditioning */ 
    iparam[4] = 0; /*  */ 
    iparam[5] = 0; /* Write solution into x */ 
    iparam[7] = 2; /* Maximum number of iterative refinement steps, output 
                      reported in iparam[6] */ 
    iparam[9]  = 13; /* Perturbing pivot elements */ 
    iparam[10] = 0;  /* Disable scaling  */ 
    iparam[12] = 0;  /*  */ 
    iparam[17] = -1; /* Enable reporting of nnz  */ 
    iparam[18] =  1; /*  */ 
    iparam[34] =  1; /* Zero-based indexing  */ 

    int maxfct, mnum, phase, error, msglvl;
    maxfct = 1; /* Maximum number of numerical factorizations. */
    mnum   = 1; /* Which factorization to use. */
    msglvl = 0; /* Print statistical information in file */
    error  = 0; /* Initialize error flag */
    
    /* Reordering and Symbolic Factorization. This step also allocates all 
     * memory that is necessary for the factorization. 
     * */
    phase = 11; /* Analysis */
    double ddummy;
    int idummy;
    debug_body("%s", "Start testing PARADISO");

    PARDISO ( ptr, &maxfct, &mnum, &mtype, &phase,
              &A.lenc, A.elem_p, A.rowIn, A.colIn, &idummy, 
              &nrhs, iparam, &msglvl, &ddummy, &ddummy, &error);
    if (error != 0)
    {
      term_exec("ERROR during symbolic factorization: %d", error);
    }
    debug_body("%s", "Analysis completed");
    
    /* Numerical factorization. 
     * */
    phase = 22;
    PARDISO ( ptr, &maxfct, &mnum, &mtype, &phase,
              &A.lenc, A.elem_p, A.rowIn, A.colIn, &idummy, 
              &nrhs, iparam, &msglvl, &ddummy, &ddummy, &error);
    if (error != 0)
    {
      term_exec("ERROR during symbolic factorization: %d", error);
    }
    debug_body("%s", "Factorization completed");

    /* Back substitution and iterative refinement.
     * */ 
    phase = 33;
    PARDISO ( ptr, &maxfct, &mnum, &mtype, &phase,
              &A.lenc, A.elem_p, A.rowIn, A.colIn, &idummy, 
              &nrhs, iparam, &msglvl, b.elem_p, x.elem_p, &error);
    if (error != 0)
    {
      term_exec("ERROR during back substitution: %d", error);
    }
    /* Termination and release of memory 
     * */ 
    phase = -1;
    PARDISO ( ptr, &maxfct, &mnum, &mtype, &phase,
              &A.lenc, A.elem_p, A.rowIn, A.colIn, &idummy, 
              &nrhs, iparam, &msglvl, &ddummy, &ddummy, &error);
    
    debug_exit("%s", "");

}

/*============================================================================*/
double test_fem1d_DGMM_general
(
    matlib_index p,
    matlib_index nr_LGL,
    matlib_index N,
    double domain[2],
    void  (*func_p)(matlib_dv, matlib_dv), 
    void  (*potential_p)(matlib_dv, matlib_dv)
)
/* 
 * potential function: phi(x)
 * Compute mass matrix for phi(x)u(x)
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

    matlib_dv u, U;
    matlib_create_dv( x.len,    &u, MATLIB_COL_VECT);
    matlib_create_dv( N*(p+1),  &U, MATLIB_COL_VECT);

    (*func_p)(x, u);
    fem1d_DFLT( N, FM, u, U);

    matlib_index dim = N*p+1, nnz = N*p*(p+3)/2+1;
    matlib_index row[dim+1], col[nnz];
    double ugpmm[nnz];

    /* Transform it to FEM-basis: vertex and bubble basis functions */  
    matlib_dv Uvb;
    matlib_create_dv( dim, &Uvb, MATLIB_COL_VECT);

    fem1d_DL2F(p, U, Uvb);


    DEBUG_PRINT_DV(U,  "%s: ", "Legendre Transform");
    DEBUG_PRINT_DV(Uvb, "%s: ", "coeff. fem basis func.");

    /* Aseemble the global mass matrix */ 
    matlib_dv Y, y;
    matlib_create_dv( x.len,   &y, MATLIB_COL_VECT);
    matlib_create_dv( N*(p+1), &Y, MATLIB_COL_VECT);
   
    (*potential_p)(x,y);

    debug_body("nnz: %d", nnz);

    matlib_dv q;
    matlib_create_dv( Q.lenc*N, &q, MATLIB_COL_VECT);
    fem1d_DFLT( N, Q, y, q);
    DEBUG_PRINT_DV(q, "%s: ", "inner product");

    fem1d_DCSRGMM(p, N, q, row, col, ugpmm);

    debug_body("nnz: %d, nnz actual: %d", row[dim], nnz);

    matlib_dsparsem M = { .lenc = dim, 
                          .lenr = dim, 
                          .rowIn = row, 
                          .colIn = col, 
                          .elem_p = ugpmm};

    BEGIN_DEBUG
        debug_print("dimension of the sparse square matrix: %d", dim);
        for(i=0; i<dim; i++)
        {
            for(j=row[i]; j<row[i+1]; j++)
            {
                debug_print("UGPMM(%d,%d): % 0.16f", i, col[j], ugpmm[j]);
            }
        }
    END_DEBUG
    
    /* Compute projection onto FEM-basis using the GMM */ 
    //matlib_dv Pvb;
    //matlib_create_dv( dim, &Pvb, MATLIB_COL_VECT);
    //matlib_dcsrsymv(MATLIB_UPPER, M, Uvb, Pvb);
       
    /* compute projection from the LP representation of phi(x)*u(x)*/
    matlib_dv Pvb1;
    matlib_create_dv( dim, &Pvb1, MATLIB_COL_VECT);


    for(i=0; i<x.len; i++)
    {
        y.elem_p[i] = u.elem_p[i] * y.elem_p[i];
    }
    fem1d_DFLT( N, FM, y, Y);
    fem1d_DPrjL2F(p, Y, Pvb1);

    matlib_dv Uvb1;
    matlib_create_dv( dim, &Uvb1, MATLIB_COL_VECT);

    test_solver_dsym(M, Pvb1, Uvb1);

    BEGIN_DEBUG
        for(i=0;i<dim; i++)
        {
            debug_print( "[%d]-> Uvb: % 0.16f, Uvb1: % 0.16f", 
                         i, *(Uvb.elem_p+i), *(Uvb1.elem_p+i));
        }
    END_DEBUG

    double norm_actual = matlib_dnrm2(Uvb);
    
    matlib_daxpy(-1.0, Uvb1, Uvb);

    double e_relative = matlib_dnrm2(Uvb)/norm_actual;

    debug_exit("Relative error: % 0.16g", e_relative);

    return(e_relative);
}
void constant_dpotential( const matlib_dv x, matlib_dv y)
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
void harmonic_dpotential( const matlib_dv x, matlib_dv y)
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
void sinusoidal_dpotential( const matlib_dv x, matlib_dv y)
{
    matlib_index i;
    if(y.len >= x.len)
    {
        for (i=0; i<x.len; i++)
        {
            *(y.elem_p) = sin(2*M_PI**(x.elem_p+i));
            y.elem_p++;
        }
    }

}
void test_fem1d_DGMM1(void)
{
    matlib_index p, nr_LGL;
    matlib_index N = 1000;
    double domain[2] = {-5.0, 5.0};
    double e_relative;
    p = 2;
    nr_LGL = 4*p+1;
    e_relative = test_fem1d_DGMM_general( p, nr_LGL, N, domain,
                                          Gaussian_func, 
                                          constant_dpotential);
    CU_ASSERT_TRUE(e_relative<TOL);

    N = 800;
    for( p=3; p<13; p++)
    {
        nr_LGL = 3*p+1;
        e_relative = test_fem1d_DGMM_general( p, nr_LGL, N, domain,
                                              Gaussian_func, 
                                              constant_dpotential);
        CU_ASSERT_TRUE(e_relative<TOL);

        e_relative = test_fem1d_DGMM_general( p, nr_LGL, N, domain,
                                              Gaussian_func, 
                                              harmonic_dpotential);
        CU_ASSERT_TRUE(e_relative<TOL);
    } 
}
/*============================================================================*/
double test_fem1d_ZGMM_general
(
    matlib_index p,
    matlib_index nr_LGL,
    matlib_index N,
    double domain[2],
    void  (*func_p)(matlib_dv, matlib_zv), 
    void  (*potential_p)(matlib_dv, matlib_zv)
)
/* 
 * potential function: phi(x)
 * Compute mass matrix for phi(x)u(x)
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

    matlib_zv u, U;
    matlib_create_zv( x.len,   &u, MATLIB_COL_VECT);
    matlib_create_zv( N*(p+1), &U, MATLIB_COL_VECT);

    (*func_p)(x, u);
    fem1d_ZFLT( N, FM, u, U);

    matlib_index dim = N*p+1, nnz = N*p*(p+3)/2+1;
    matlib_index row[dim+1], col[nnz];
    matlib_complex ugpmm[nnz];

    /* Transform it to FEM-basis: vertex and bubble basis functions */  
    matlib_zv Uvb;
    matlib_create_zv( dim, &Uvb, MATLIB_COL_VECT);

    fem1d_ZL2F(p, U, Uvb);


    DEBUG_PRINT_ZV(U,  "%s: ", "Legendre Transform");
    DEBUG_PRINT_ZV(Uvb, "%s: ", "coeff. fem basis func.");

    /* Aseemble the global mass matrix */ 
    matlib_zv Y, y;
    matlib_create_zv( x.len,   &y, MATLIB_COL_VECT);
    matlib_create_zv( N*(p+1), &Y, MATLIB_COL_VECT);
   
    (*potential_p)(x,y);

    debug_body("nnz: %d", nnz);

    matlib_zv q;
    matlib_create_zv( Q.lenc*N, &q, MATLIB_COL_VECT);
    fem1d_ZFLT( N, Q, y, q);
    DEBUG_PRINT_ZV(q, "%s: ", "inner product");

    fem1d_ZCSRGMM(p, N, q, row, col, ugpmm);

    debug_body("nnz: %d, nnz actual: %d", row[dim], nnz);

    matlib_zsparsem M = { .lenc = dim, 
                          .lenr = dim, 
                          .rowIn = row, 
                          .colIn = col, 
                          .elem_p = ugpmm};

    BEGIN_DEBUG
        debug_print("dimension of the sparse square matrix: %d", dim);
        for(i=0; i<dim; i++)
        {
            for(j=row[i]; j<row[i+1]; j++)
            {
                debug_print("UGPMM(%d,%d): % 0.16f %+0.16fi", i, col[j], ugpmm[j]);
            }
        }
    END_DEBUG
    
    /* Compute projection onto FEM-basis using the GMM */ 
    //matlib_zv Pvb;
    //matlib_create_zv( dim, &Pvb, MATLIB_COL_VECT);
    //matlib_zcsrsymv(MATLIB_UPPER, M, Uvb, Pvb);
       
    /* compute projection from the LP representation of phi(x)*u(x)*/
    matlib_zv Pvb1;
    matlib_create_zv( dim, &Pvb1, MATLIB_COL_VECT);

    for(i=0; i<x.len; i++)
    {
        y.elem_p[i] = u.elem_p[i] * y.elem_p[i];
    
    }
    fem1d_ZFLT( N, FM, y, Y);
    DEBUG_PRINT_ZV(Y, "%s: ", "Transfomed Vector");
    fem1d_ZPrjL2F(p, Y, Pvb1);

    matlib_zv Uvb1;
    matlib_create_zv( dim, &Uvb1, MATLIB_COL_VECT);

    test_solver_zsym(M, Pvb1, Uvb1);

    BEGIN_DEBUG
        for(i=0;i<dim; i++)
        {
            debug_print( "[%d]-> Pvb-Pvb1: % 0.16f%+0.16fi", 
                         i, *(Uvb.elem_p+i)-*(Uvb1.elem_p+i));
        }
    END_DEBUG

    double norm_actual = matlib_znrm2(Uvb1);
    
    matlib_zaxpy(-1.0, Uvb1, Uvb);

    double e_relative = matlib_znrm2(Uvb)/norm_actual;

    debug_exit("Relative error: % 0.16g", e_relative);

    return(e_relative);
}
void constant_zpotential( const matlib_dv x, matlib_zv y)
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
void harmonic_zpotential( const matlib_dv x, matlib_zv y)
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

void test_fem1d_ZGMM1(void)
{
    matlib_index p, nr_LGL;
    matlib_index N = 2000;
    double domain[2] = {-5.0, 5.0};
    double e_relative;

    p = 2;
    nr_LGL = 2*p+1;
    e_relative = test_fem1d_ZGMM_general( p, nr_LGL, N, domain,
                                          Gaussian_zfunc, 
                                          constant_zpotential);
    CU_ASSERT_TRUE(e_relative<TOL);

    e_relative = test_fem1d_ZGMM_general( p, nr_LGL, N, domain,
                                          Gaussian_zfunc, 
                                          harmonic_zpotential);
    CU_ASSERT_TRUE(e_relative<TOL);
    
    N = 1000;
    for(matlib_index p=3; p<13; p++)
    {
        nr_LGL = 3*p+1;
        e_relative = test_fem1d_ZGMM_general( p, nr_LGL, N, domain,
                                              Gaussian_zfunc, 
                                              constant_zpotential);
        CU_ASSERT_TRUE(e_relative<TOL);

        e_relative = test_fem1d_ZGMM_general( p, nr_LGL, N, domain,
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
        { "Global mass matrix for Gaussian real"   , test_fem1d_DGMM1    },
        { "Global mass matrix for Gaussian complex", test_fem1d_ZGMM1    },
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

