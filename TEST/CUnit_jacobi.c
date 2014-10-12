/*============================================================================
 | File: CUnit_jacobi.c
 | Description: Unittest for Legendre Transform Library
 |              Depends of Matrix Library for Linear Algebra.
 |
 |
 |============================================================================*/


/*============================================================================*/
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include "mkl.h"

#define NDEBUG
#include "jacobi.h"

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
matlib_real jacobi_calcjp2
(
    matlib_index n, 
    matlib_real a, 
    matlib_real b, 
    matlib_real x
)
/* 
 * Calculate Jacobi polynomials: $J^{(a,b)}_{n-1}(x),J^{(a,b)}_{n}(x)$ 
 * $a,b>-1$                                                             
 *
 * */ 
{
    matlib_index i;
    matlib_real tmp, tmp1, JP[2];
    matlib_real C, D, E;
    
    *(JP+0) = 1;
    /* 
     * C[0] = 0.5*(a-b);
     * C[1] = 0.5*(a+b+2); 
     *
     * */
     *(JP+1) = 0.5*(a-b)+0.5*(a+b+2)*x;
    if(n>1) 
    {
        for( i=2; i<(n+1); i++)
        {
            
            tmp1 = 2*i*(i+a+b)*(2*i-2+a+b);
            C = (2*i-1+a+b)*(a*a-b*b)/tmp1;
            D = (2*i-2+a+b)*(2*i-1+a+b)*(2*i+a+b)/tmp1;
            E = 2*(i-1+a)*(i-1+b)*(2*i+a+b)/tmp1;
            tmp = *(JP+1);
            *(JP+1) = (C+D*x)**(JP+1)-E**(JP+0);
            *(JP+0) = tmp;
        }
    }
    return(JP[1]);
}

/* Test the zeros of Jacobi polynomial J^(a,b)_n(x), a=b=-0.5
 * Identity:  Chebyshev polynimial T_n(x)
 * T_n(x) = J^(a,b)_n(x)/J^(a,b)_n(1)
 *
 * */ 
void test_jacobi_zeros_general(matlib_index n)
{
    matlib_real a = -0.4;
    matlib_real b = -0.3;
    matlib_real T;

    matlib_xv x;
    matlib_create_xv(n, &x, MATLIB_COL_VECT);

    jacobi_find_zeros(n, a, b, x.elem_p, TOL);
    for (matlib_index i=0; i< x.len; i++)
    {
        T = jacobi_calcjp2(n, a, b, x.elem_p[i]);
        CU_ASSERT_TRUE(fabs(T)<TOL);
    }
}
void test_jacobi_zeros(void)
{
    matlib_index n_max = 10;
    for(matlib_index n=1; n<n_max; n++)
    {
        test_jacobi_zeros_general(n);
    }
}

/*============================================================================
 | Test runner
 |
 |
 |
 |
 |
 +============================================================================*/

int main()
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
        { "Roots of Jacobi poly.", test_jacobi_zeros            },
        CU_TEST_INFO_NULL,
    };
    /* Create the test suite */ 
    CU_SuiteInfo suites[] = 
    {
        { "Jacobi polynomials", init_suite, clean_suite, NULL, NULL, test_array },
        CU_SUITE_INFO_NULL,
    }; 

    /* Register test suites */ 
    CU_ErrorCode CU_error = CU_register_suites(suites); 
    if (CU_error != CUE_SUCCESS) 
    {
        CU_cleanup_registry();
        return CU_get_error();
    }

   /* Run all tests using the CUnit Basic interface */
   CU_basic_set_mode(CU_BRM_VERBOSE);
   CU_basic_run_tests();
   CU_cleanup_registry();
   return CU_get_error();
}

