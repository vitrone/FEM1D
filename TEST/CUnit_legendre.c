/*============================================================================
 | File: CUnit_legendre.c
 | Description: Unittest for Legendre Transform Library
 |              Depends of Matrix Library for Linear Algebra.
 |
 |
 |
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
#define MATLIB_NTRACE_DATA
#include "legendre.h"

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
static void poly_func
( 
    const matlib_xv x, 
          matlib_xv y
)
{
    int deg = 5;
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

static void Gaussian_func
( 
    matlib_xv x, 
    matlib_xv y
)
{
    matlib_index i;
    if(y.len >= x.len)
    {
        for (i=0; i<x.len; i++)
        {
            *(y.elem_p) = exp(-16*(*(x.elem_p+i)**(x.elem_p+i)));
            y.elem_p++;
        }
    }
}
static void LP10
( 
    matlib_xv x, 
    matlib_xv y
)
{
    /* lp = (1.0/256.0)*(46189.0*x^10-109395.0*x^8+90090.0*x^6-30030.0*x^4+3465.0*x^2-63)
     * */;
    matlib_real C_0  = (1.0/256.0)*(63)      ;
    matlib_real C_2  = (1.0/256.0)*(3465.0)  ;
    matlib_real C_4  = (1.0/256.0)*(30030.0) ;
    matlib_real C_6  = (1.0/256.0)*(90090.0) ;
    matlib_real C_8  = (1.0/256.0)*(109395.0);
    matlib_real C_10 = (1.0/256.0)*(46189.0) ;
    matlib_real x_2, x_4, x_6, x_8, x_10;

    matlib_index i;
    if(y.len >= x.len)
    {
        for (i=0; i<x.len; i++)
        {
            x_2  = *x.elem_p * *x.elem_p;
            x_4  = x_2 * x_2;
            x_6  = x_2 * x_4;
            x_8  = x_6 * x_2;
            x_10 = x_8 * x_2;

            *(y.elem_p) = C_10*x_10-C_8*x_8+C_6*x_6-C_4*x_4+C_2*x_2-C_0;
            y.elem_p++;
            x.elem_p++;
        }
    }
}

static void LP09
( 
    matlib_xv x, 
    matlib_xv y
)
{
    /* lp = (1/128)*(12155*x**9-25740*x**7+18018*x**5-4620*x**3+315*x)
     * */;
    
    matlib_real C_1  = (1/128.0)*(315.0)  ;
    matlib_real C_3  = (1/128.0)*(4620.0) ;
    matlib_real C_5  = (1/128.0)*(18018.0);
    matlib_real C_7  = (1/128.0)*(25740.0);
    matlib_real C_9  = (1/128.0)*(12155.0);
    matlib_real x_2, x_3, x_5, x_7, x_9;

    matlib_index i;
    if(y.len >= x.len)
    {
        for (i=0; i<x.len; i++)
        {
            x_2  = *x.elem_p * *x.elem_p;
            x_3  = x_2 * *x.elem_p;
            x_5  = x_2 * x_3;
            x_7  = x_2 * x_5;
            x_9  = x_2 * x_7;

            *(y.elem_p) = C_9*x_9 - C_7*x_7 + C_5*x_5 - C_3*x_3 + C_1**x.elem_p;
            y.elem_p++;
            x.elem_p++;
        }
    }
}

/*============================================================================*/
static void test_Gauss_points(void)
{
    matlib_index p = 10;
    debug_enter( "degree of polynomial: %d", p);
    matlib_xv gzeros, y; /* must be matlib_freed */ 
    matlib_create_xv(p, &gzeros, MATLIB_COL_VECT);
    matlib_create_xv(p, &y, MATLIB_COL_VECT);

    find_Gauss_points( p, TOL, gzeros.elem_p);
    LP10(gzeros,y);

    /* Verification */ 
    matlib_real *ptr;
    for(ptr= y.elem_p; ptr < y.elem_p+y.len; ptr++)
    {
        CU_ASSERT_TRUE(fabs(*ptr)<TOL);
        debug_body( "|LP10(gzeros)| = %0.16f", fabs(*ptr));
    }
    matlib_free(gzeros.elem_p);
    matlib_free(y.elem_p);
    debug_exit( "%s", "");
}

static void test_LGL_points(void)
{
    matlib_index p      = 10;
    matlib_index nr_LGL = p+1;

    debug_enter( "nr. LGL-points: %d", nr_LGL);
    
    matlib_xv gzeros, zeros, quadW, y, yy; /* must be matlib_freed */

    matlib_create_xv(     p, &gzeros, MATLIB_COL_VECT);
    matlib_create_xv(nr_LGL,  &zeros, MATLIB_COL_VECT);
    matlib_create_xv(nr_LGL,  &quadW, MATLIB_COL_VECT);
    matlib_create_xv(nr_LGL,      &y, MATLIB_COL_VECT);
    matlib_create_xv(nr_LGL,     &yy, MATLIB_COL_VECT);

    find_Gauss_points( p, TOL, gzeros.elem_p);
    find_LGL_points( p, TOL, zeros.elem_p, quadW.elem_p, gzeros.elem_p);

    LP09(zeros,y);
    LP10(zeros,yy);

    /* Verification */ 
    matlib_real *ptr, *ptr1 = zeros.elem_p, *ptr2 = yy.elem_p ;
    for(ptr= y.elem_p; ptr < y.elem_p+y.len; ptr++)
    {
        *ptr2 = p**ptr - p**ptr1**ptr2;
        CU_ASSERT_TRUE(fabs(*ptr2)<TOL);
        debug_body( "|LGL11(zeros)| = %0.16f", fabs(*ptr2));
        ptr2++;
        ptr1++;
    }
    /* Quadrature test */
    matlib_index max_deg;

    matlib_index i;
    matlib_real sum;
    for (i=1; i<2*p-1; i++)
    {
        max_deg = 2*p-i;
        sum = 0;
        ptr1 = zeros.elem_p;
        for (ptr = quadW.elem_p; ptr<quadW.elem_p+quadW.len; ptr++)
        {
            sum += (pow(*ptr1, max_deg)**ptr );
            ptr1++;
        }
        if((max_deg % 2)==0)
        {
            CU_ASSERT_TRUE(fabs(sum-2.0/(max_deg+1))<TOL);
            debug_body( "Num. integral of x^%d = %0.16f, actual: %0.16f", 
                        max_deg, sum, 2.0/(max_deg+1));
        }
        else
        {
            CU_ASSERT_TRUE(fabs(sum)<TOL);
            debug_body( "Num. integral of x^%d = %0.16f", max_deg, sum);
        }
    }

    matlib_free(gzeros.elem_p);
    matlib_free(zeros.elem_p);
    matlib_free(quadW.elem_p);
    matlib_free(y.elem_p);
    matlib_free(yy.elem_p);
    debug_exit( "%s", "");
}
/*============================================================================
 | Decomposition of monomials into Legendre polynomials 
 | x^2 = 1/3[P_0(x)+2P_2(x)]	
 | x^3 = 1/5[3P_1(x)+2P_3(x)]	
 | x^4 = 1/(35)[7P_0(x)+20P_2(x)+8P_4(x)]	
 | x^5 = 1/(63)[27P_1(x)+28P_3(x)+8P_5(x)]	
 | x^6 = 1/(231)[33P_0(x)+110P_2(x)+72P_4(x)+16P_6(x)]
 +============================================================================*/
static void monomial_func
( 
    const int       p,
    const matlib_xv x, 
          matlib_xv y
)
{
    matlib_index i;
    if(y.len >= x.len)
    {
        for (i=0; i<x.len; i++)
        {
            *(y.elem_p) = pow(*(x.elem_p+i), p);
            y.elem_p++;
        }
    }
}

static void test_forward_transform(void)
{
    debug_enter( "%s", "");
    /* Define Legendre transforms of monomials 
     * */
    matlib_real U_2[] = {   1.0/3,       0,     2.0/3                                      };
    matlib_real U_3[] = {       0,   3.0/5,         0,   2.0/5                             };
    matlib_real U_4[] = {  7.0/35,       0,   20.0/35,       0,   8.0/35                   };
    matlib_real U_5[] = {       0, 27.0/63,         0, 28.0/63,        0, 8.0/63           };
    matlib_real U_6[] = {33.0/231,       0, 110.0/231,       0, 72.0/231,      0, 16.0/231 };
    
    int i, j, nr_MLT = 5, p, nr_LGL;
    /* An array of pointers */ 
    matlib_real *MLT_array[5] = {U_2, U_3, U_4, U_5, U_6};
    
    matlib_xv xi, quadW, u, U, U_exact;
    matlib_xm FM;
    matlib_real e_relative, norm_exact;

    for (i=0; i<nr_MLT; i++)
    {
        p = i+2;
        legendre_LGLdataLT1( p, TOL, &xi, &quadW);
        
        matlib_create_xm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);
        legendre_LGLdataFM( xi, FM);

        nr_LGL = p+1;
        matlib_create_xv( nr_LGL, &u,       MATLIB_COL_VECT);
        matlib_create_xv( nr_LGL, &U,       MATLIB_COL_VECT);
        matlib_create_xv( nr_LGL, &U_exact, MATLIB_COL_VECT);
        
        monomial_func(p, xi, u);
        /* U: LT of u */ 
        matlib_xgemv( 1.0, FM, u, 0, U);
        U_exact.len    = nr_LGL;
        U_exact.elem_p = *(MLT_array+i);

        debug_body("Monomial degree: %d", p);
        BEGIN_DTRACE
            for(j=0; (j<nr_LGL); j++)
            {
                debug_print( "U_exact[%d]: %0.16f, U_computed[%d]: %0.16f", 
                             j, *(U_exact.elem_p+j), j, *(U.elem_p+j) );
            }
        END_DTRACE

        norm_exact = matlib_xnrm2(U_exact);
        matlib_xaxpy(-1.0, U_exact, U);
        e_relative = matlib_xnrm2(U)/norm_exact;
        debug_body("Monomial degree: %d, Relative error: %0.16f", p, e_relative);
        CU_ASSERT_TRUE(e_relative<TOL);

        matlib_free( FM.elem_p);
        matlib_free( xi.elem_p);
        matlib_free( u.elem_p);
        matlib_free( U.elem_p);
        matlib_free( quadW.elem_p);
    }
    debug_exit( "%s", "");
}



/*============================================================================*/

static matlib_real test_legendre_interpolation_general
(
    matlib_index p,
    matlib_index nr_samples,
    void  (*func_p)(matlib_xv, matlib_xv)
)
{
    debug_enter( "degree of polynomial: %d, nr. LGL-points for sampling: %d",
                 p, nr_samples );
    matlib_index i;

    matlib_index nr_LGL = p+1;

    matlib_xm FM;
    matlib_xv xi, quadW, u, U;

    legendre_LGLdataLT1( p, TOL, &xi, &quadW);
    
    matlib_create_xm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);
    legendre_LGLdataFM( xi, FM);

    matlib_create_xv( nr_LGL, &u, MATLIB_COL_VECT);
    matlib_create_xv( nr_LGL, &U, MATLIB_COL_VECT);

                          
    (*func_p)(xi, u);

    /* U: LT of u */ 
    matlib_xgemv( 1.0, FM, u, 0, U);

    DEBUG_PRINT_XV(U, "%s: ", "transfomred matrix");

    /* Verification =========================*/
    matlib_index P = nr_samples - 1;

    matlib_xm IMi;
    matlib_xv xii, u_interpol, u_actual;

    matlib_create_xv( nr_samples, &xii,        MATLIB_COL_VECT);
    matlib_create_xv( nr_samples, &u_interpol, MATLIB_COL_VECT);
    matlib_create_xv( nr_samples, &u_actual,   MATLIB_COL_VECT);

    /* generate the grid on the reference interval */ 
    matlib_real dxii = 2.0/P;
    debug_body("dxi=%0.16f", dxii);

    *(xii.elem_p) = -1.0;
    for (i=1; i<P; i++)
    {
        *(xii.elem_p+i) = *(xii.elem_p+i-1) + dxii;
    }
    *(xii.elem_p+P) = 1.0;
    
    matlib_create_xm( xii.len, p+1, &IMi, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);
    legendre_LGLdataIM(xii, IMi);

    matlib_xgemv( 1.0, IMi, U, 0, u_interpol);

    (*func_p)(xii, u_actual);

    /* Compute the error */
    BEGIN_DTRACE
        for (i=0; (i<nr_samples); i++)
        {
            debug_print( "x=%0.3f, interpolated y:%0.16f, actual y= %0.16f",
                        xii.elem_p[i], 
                        u_interpol.elem_p[i], 
                        u_actual.elem_p[i]);
        }
    END_DTRACE
    matlib_real norm_actual = matlib_xnrm2(u_actual);
    matlib_xaxpy(-1.0, u_actual, u_interpol);
    matlib_real e_relative = matlib_xnrm2(u_interpol)/norm_actual;
    
    /* Free allocated memory*/ 
    matlib_free( (void*) FM.elem_p);
    matlib_free( (void*) xi.elem_p);
    matlib_free( (void*) u.elem_p);
    matlib_free( (void*) U.elem_p);
    matlib_free( (void*) quadW.elem_p);
    matlib_free( (void*) IMi.elem_p);
    matlib_free( (void*) xii.elem_p);
    matlib_free( (void*) u_actual.elem_p);
    matlib_free( (void*) u_interpol.elem_p);
    debug_exit("Relative interpolation error: %0.16g\n", e_relative);

    return(e_relative);
}
/*============================================================================*/

static void test_legendre_interpolation1(void)
{

    matlib_index p = 10;
    matlib_index nr_samples = 5;
    matlib_real e_relative;
    e_relative = test_legendre_interpolation_general( p, 
                                                      nr_samples,
                                                      poly_func);

    CU_ASSERT_TRUE(e_relative<TOL);

}
static void test_legendre_interpolation2(void)
{

    matlib_index p = 40;
    matlib_index nr_samples = 10;
    matlib_real e_relative;
    e_relative = test_legendre_interpolation_general( p, 
                                                      nr_samples,
                                                      Gaussian_func);

    CU_ASSERT_TRUE(e_relative<TOL);

}
static void test_FM_and_IM(void)
{

    matlib_index p  = 30;
    matlib_index ni = 20;
    matlib_index i, P;

    matlib_xm FM, IM, FMi, IMi;
    matlib_xv xi, quadW, u, U, Ui;

    matlib_real e_relative, norm_U;

    for (i=0; i<ni; i++)
    {
        P = (i+1)*p;
        debug_body("nr. of samples : %d", P+1);
        legendre_LGLdataLT2( P, TOL, &xi, &quadW, &FM, &IM);
        matlib_create_xv( xi.len, &u, MATLIB_COL_VECT);
        matlib_create_xv( P+1,    &U, MATLIB_COL_VECT);
        
        Gaussian_func(xi, u);
        /* U: LT of u */ 
        matlib_xgemv( 1.0, FM, u, 0, U);

        matlib_create_xm( p+1, xi.len, &FMi, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);
        legendre_LGLdataFM( xi, FMi);

        matlib_create_xv( p+1,    &Ui, MATLIB_COL_VECT);

        /* U: LT of u */ 
        matlib_xgemv( 1.0, FMi, u, 0, Ui);
        U.len = Ui.len;
        norm_U = matlib_xnrm2(U);
        matlib_xaxpy(-1.0, U, Ui);
        e_relative = matlib_xnrm2(Ui)/norm_U;
        debug_body("relative error: % 0.16g", e_relative);
        CU_ASSERT_TRUE(e_relative==0);

        matlib_free( FM.elem_p);
        matlib_free( FMi.elem_p);
        matlib_free( u.elem_p);
        matlib_free( U.elem_p);
        matlib_free( Ui.elem_p);
        matlib_free( quadW.elem_p);
        matlib_free( xi.elem_p);
    }
    debug_exit("%s", "");

}
static void test_parsevals_theorem(void)
{

    matlib_index p  = 100;
    matlib_index ni = 20;
    matlib_index i, k, P;

    matlib_xm FM, IM;
    matlib_xv xi, quadW, u, U, coeff;
    
    matlib_real e_relative[ni], snorm_spectral, snorm_actual = sqrt(M_PI/32.0);

    for (i=0; i<ni; i++)
    {
        P = (i+1)*p;
        debug_body("nr. of samples : %d", P+1);
        legendre_LGLdataLT1( P, TOL, &xi, &quadW);
        matlib_create_xm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);
        legendre_LGLdataFM( xi, FM);

        matlib_create_xv( xi.len,     &u, MATLIB_COL_VECT);
        matlib_create_xv(    p+1,     &U, MATLIB_COL_VECT);
        matlib_create_xv(    p+1, &coeff, MATLIB_COL_VECT);
                          
        Gaussian_func(xi, u);

        /* U: LT of u */ 
        matlib_xgemv( 1.0, FM, u, 0, U);
        for(k=0; k<p+1; k++ )
        {
            U.elem_p[k]  = U.elem_p[k] * U.elem_p[k];
            coeff.elem_p[k] = 2.0/(2*k+1.0);
        }

        snorm_spectral = matlib_xdot( U, coeff);
        e_relative[i] = fabs(snorm_spectral-snorm_actual)/snorm_actual;
        debug_body("relative error: % 0.16g", e_relative[i]);
        CU_ASSERT_TRUE(e_relative[i]<TOL);

        matlib_free( FM.elem_p);
        matlib_free( u.elem_p);
        matlib_free( U.elem_p);
        matlib_free( quadW.elem_p);
        matlib_free( xi.elem_p);
    }
    for(i=0; i<ni-1; i++)
    {
        CU_ASSERT_TRUE(e_relative[i]>e_relative[ni-1]);
    }
    debug_exit("%s", "");
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
        { "Gauss-points"                , test_Gauss_points            },
        { "LGL-points"                  , test_LGL_points              },
        { "forward transform"           , test_forward_transform       },
        { "interpolation for polynomial", test_legendre_interpolation1 },
        { "interpolation for Gaussian"  , test_legendre_interpolation2 },
        //{ "Aliasing effect", test_aliasing_effect1 },
        { "Parseval's Theorem", test_parsevals_theorem },
        { "FM and IM ", test_FM_and_IM },
        CU_TEST_INFO_NULL,
    };
    /* Create the test suite */ 
    CU_SuiteInfo suites[] = 
    {
        { "Legendre Transforms", init_suite, clean_suite, test_array },
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

