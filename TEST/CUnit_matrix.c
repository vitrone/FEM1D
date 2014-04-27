#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include "mkl.h"


#define DEBUG
#include "basic.h"
#include "matlib.h"
#include "debug.h"
#include "ehandler.h"

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
/* In C99 complex numbers are defined as an adjacent pair of memory locations.
 * However, there are some syntactical differences which require type-casting:
 * (1) Interpreting a complex array as N-by-2 double matrix in rwo major format.
 * (2) Interpreting a double matrix in row major format as an array of complex
 *     numbers.
 *
 * */ 

/* The following type is fully compatible with C99 complex.*/ 
typedef double matlib_c[2];

void test_complex(void)
{
    matlib_c xa[] = { {1.0, 2.0 }, 
                      {3.0, 4.5 }};

    debug_body("xa[0] = %0.16f, %0.16f", (*(xa+0))[0], (*(xa+0))[1]);
    debug_body("xa[1] = %0.16f, %0.16f", (*(xa+1))[0], (*(xa+1))[1]);
    matlib_dm xm = { .lenc   = 2,
                     .lenr  = 2,
                     .order  = MATLIB_ROW_MAJOR,
                     .elem_p = &(*xa)[0]};
    
    DEBUG_PRINT_DM(xm, "%s", "");
    matlib_zv xcm = { .len    = 2,
                      .type   = MATLIB_ROW_MAJOR,
                      .elem_p = (complex*)xa};
    
    DEBUG_PRINT_ZV(xcm, "%s", "");
    /* convert such xa to complex type */
    matlib_c ya[2];
    *((complex*)ya[0]) = csqrt(*((complex*)xa[0]));
    *((complex*)ya[1]) = csqrt(*((complex*)xa[1]));
    debug_body("ya[0]: %0.16f %+0.16fi", *((complex*)ya[0]));
    debug_body("ya[1]: %0.16f %+0.16fi", *((complex*)ya[1]));

    matlib_dm ym = { .lenc   = 2,
                     .lenr  = 2,
                     .order  = MATLIB_ROW_MAJOR,
                     .elem_p = &(*ya)[0]};
    DEBUG_PRINT_DM(ym, "%s", "");
    matlib_zv ycv = { .len    = 2,
                      .type   = MATLIB_ROW_MAJOR,
                      .elem_p = (complex*)ya};
    
    DEBUG_PRINT_ZV(ycv, "%s", "");
    
    /* Verifying if C99 complex can be converted to double matrix.*/ 

    matlib_complex za[4] = { -1.0 + I * 2.0, 
                              1.6 + I * 2.1};

    matlib_zv z = {  .len   = 2, 
                     .type  = MATLIB_COL_VECT,
                     .elem_p = za};
    
    DEBUG_PRINT_ZV(z, "%s", "");

    debug_body("za[0] = %0.16f %+0.16fi", *((double*)(za+0)), *((double*)(za+0)+1));
    debug_body("za[1] = %0.16f %+0.16fi", *((double*)(za+1)), *((double*)(za+1)+1));

    matlib_dm zm = { .lenc  = 2,
                     .lenr  = 2,
                     .order   = MATLIB_ROW_MAJOR,
                     .elem_p = ((double*)za)};

    DEBUG_PRINT_DM(zm, "%s", "");

    CU_ASSERT_TRUE(true);

}

/*============================================================================*/
void test_zgemv(void)
{
    debug_enter("%s", "");
    matlib_complex Ma[3][4] = { { 3.0 + I*1.0, 1.0+ I*1.0, 3.0+ I*1.0, 2.0+ I*1.0},
                                { 1.0 + I*1.0, 5.0+ I*1.0, 9.0+ I*1.0, 4.0+ I*1.0},
                                { 2.0 + I*1.0, 6.0+ I*1.0, 5.0+ I*1.0, 7.0+ I*1.0}
                          };

    matlib_zm M = { .lenc   = 3, 
                    .lenr   = 4, 
                    .order  = MATLIB_ROW_MAJOR,
                    .op     = MATLIB_NO_TRANS,
                    .elem_p = &Ma[0][0]};
    DEBUG_PRINT_ZM(M, "%s", "");

    matlib_complex xa[] = {-1.0 + I * 2.0, 
                           -1.0 + I * 2.0, 
                            1.0 + I * 2.0, 
                            1.0 + I * 2.0 };

    matlib_zv x = { .len    = 4, 
                    .type   = MATLIB_COL_VECT, 
                    .elem_p = xa};
    DEBUG_PRINT_ZV(x, "%s", "");

    matlib_complex ya[] = { 0, 0, 0 };

    matlib_zv y = { .len    = 3, 
                    .type   = MATLIB_COL_VECT, 
                    .elem_p = ya};

    matlib_zgemv( 1.0, M, x, 0, y);
    matlib_complex ya_actual[] = { Ma[0][0] * xa[0] + Ma[0][1] * xa[1] + Ma[0][2] * xa[2]+ Ma[0][3] * xa[3], 
                                   Ma[1][0] * xa[0] + Ma[1][1] * xa[1] + Ma[1][2] * xa[2]+ Ma[1][3] * xa[3],
                                   Ma[2][0] * xa[0] + Ma[2][1] * xa[1] + Ma[2][2] * xa[2]+ Ma[2][3] * xa[3]};

    matlib_zv y_actual = { .len    = 3, 
                           .type   = MATLIB_COL_VECT, 
                           .elem_p = ya_actual};

    DEBUG_PRINT_ZV(y, "%s", "");
    DEBUG_PRINT_ZV(y_actual, "%s", "");
    double norm_actual = matlib_znrm2(y_actual);
    matlib_zaxpy(-1.0, y_actual, y);
    double e_relative = matlib_znrm2(y)/norm_actual;

    debug_exit("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);


}
/*============================================================================*/

void test_dgemv(void)
{
    debug_enter("%s", "");

    
    double Ma[3][4] = { { 3.0, 1.0, 3.0, 2.0},
                        { 1.0, 5.0, 9.0, 4.0},
                        { 2.0, 6.0, 5.0, 7.0}
                  };
    matlib_dm M = { .lenc   = 3, 
                    .lenr   = 4, 
                    .order  = MATLIB_ROW_MAJOR,
                    .op     = MATLIB_NO_TRANS,
                    .elem_p = &Ma[0][0]};

    DEBUG_PRINT_DM(M, "%s", "");
    
    double xa[] = {-1.0, -1.0, 1.0, 1.0};
    matlib_dv x = { .len    = 4, 
                    .type   = MATLIB_COL_VECT, 
                    .elem_p = xa};
    
    DEBUG_PRINT_DV(x, "%s", "");
    
    double ya[] = { 0, 0, 0 };

    matlib_dv y = { .len    = 3, 
                    .type   = MATLIB_COL_VECT, 
                    .elem_p = ya};


    matlib_dgemv( 1, M, x, 0, y);
    double ya_actual[] = { Ma[0][0] * xa[0] + Ma[0][1] * xa[1] + Ma[0][2] * xa[2]+ Ma[0][3] * xa[3], 
                           Ma[1][0] * xa[0] + Ma[1][1] * xa[1] + Ma[1][2] * xa[2]+ Ma[1][3] * xa[3],
                           Ma[2][0] * xa[0] + Ma[2][1] * xa[1] + Ma[2][2] * xa[2]+ Ma[2][3] * xa[3]};

    matlib_dv y_actual = { .len    = 3, 
                           .type   = MATLIB_COL_VECT, 
                           .elem_p = ya_actual};

    DEBUG_PRINT_DV(y, "%s", "");
    DEBUG_PRINT_DV(y_actual, "%s", "");
    double norm_actual = matlib_dnrm2(y_actual);
    matlib_daxpy(-1.0, y_actual, y);
    double e_relative = matlib_dnrm2(y)/norm_actual;

    debug_exit("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);
    
}

void test_dgemm(void)
{
    debug_enter("%s", "");

    
    double Ma[3][4] = { { 3.0, 1.0, 3.0, 2.0},
                        { 1.0, 5.0, 9.0, 4.0},
                        { 2.0, 6.0, 5.0, 7.0}
                  };
    matlib_dm M = { .lenc   = 3, 
                    .lenr   = 4, 
                    .order  = MATLIB_ROW_MAJOR,
                    .elem_p = &Ma[0][0]};

    DEBUG_PRINT_DM(M, "%s", "");
    
    double Na[4][3] = { {  2.0,  1.0,  3.0 },
                        {  1.0,  3.0, -7.0 },
                        {  2.0, -1.0,  5.0 },
                        { -2.0,  6.0,  8.0 }
                  };

    matlib_dm N = { .lenc   = 4, 
                    .lenr   = 3, 
                    .order  = MATLIB_ROW_MAJOR,
                    .elem_p = &Na[0][0]};
    
    DEBUG_PRINT_DM(N, "%s", "");
    double Pa[3][3];
    matlib_dm P = { .lenc   = 3, 
                    .lenr   = 3, 
                    .order  = MATLIB_ROW_MAJOR,
                    .elem_p = &Pa[0][0]};
    

    matlib_dgemm( 1, M, N, 0, P);
    DEBUG_PRINT_DM(P, "%s", "");

    double Pa_actual[3][3] = { { 9.0,    15.0,    33.0},
                               {17.0,    31.0,    45.0},
                               { 6.0,    57.0,    45.0}
                             };

    matlib_dm P_actual = { .lenc   = 3, 
                           .lenr   = 3, 
                           .order  = MATLIB_ROW_MAJOR,
                           .elem_p = &Pa_actual[0][0]};


    DEBUG_PRINT_DM(P_actual, "%s", "");
    matlib_dv u = MK_VM(P);
    matlib_dv u_actual = MK_VM(P_actual);

    double norm_actual = matlib_dnrm2(u_actual);
    matlib_daxpy(-1.0, u_actual, u);
    double e_relative = matlib_dnrm2(u)/norm_actual;

    debug_exit("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);
    
}
void test_zgemm(void)
{
    debug_enter("%s", "");

    
    matlib_complex Ma[3][4] = { { 3.0 + I*2.0, 1.0 + I*2.0, 3.0 + I*2.0, 2.0 + I*2.0},
                                { 1.0 + I*2.0, 5.0 + I*2.0, 9.0 + I*2.0, 4.0 + I*2.0},
                                { 2.0 + I*2.0, 6.0 + I*2.0, 5.0 + I*2.0, 7.0 + I*2.0}
                  };
    matlib_zm M = { .lenc   = 3, 
                    .lenr   = 4, 
                    .order  = MATLIB_ROW_MAJOR,
                    .elem_p = &Ma[0][0]};

    DEBUG_PRINT_ZM(M, "%s", "");
    
    matlib_complex Na[4][3] = { {  2.0 + I*3.0,  1.0 + I*3.0,  3.0 + I*3.0 },
                                {  1.0 + I*3.0,  3.0 + I*3.0, -7.0 + I*3.0 },
                                {  2.0 + I*3.0, -1.0 + I*3.0,  5.0 + I*3.0 },
                                { -2.0 + I*3.0,  6.0 + I*3.0,  8.0 + I*3.0 }
                  };

    matlib_zm N = { .lenc   = 4, 
                    .lenr   = 3, 
                    .order  = MATLIB_ROW_MAJOR,
                    .elem_p = &Na[0][0]};
    
    DEBUG_PRINT_ZM(N, "%s", "");
    matlib_complex Pa[3][3];
    matlib_zm P = { .lenc   = 3, 
                    .lenr   = 3, 
                    .order  = MATLIB_ROW_MAJOR,
                    .elem_p = &Pa[0][0]};
    

    matlib_zgemm( 1, M, N, 0, P);
    DEBUG_PRINT_ZM(P, "%s", "");

    matlib_complex Pa_actual[3][3] = { {-15.0000 + I*33.0000,  -9.0000 + I*45.0000,   9.0000 + I*45.0000},
                                       {-7.0000  + I*63.0000,   7.0000 + I*75.0000,  21.0000 + I*75.0000},
                                       {-18.0000 + I*66.0000,  33.0000 + I*78.0000,  21.0000 + I*78.0000}};

    matlib_zm P_actual = { .lenc   = 3, 
                           .lenr   = 3, 
                           .order  = MATLIB_ROW_MAJOR,
                           .elem_p = &Pa_actual[0][0]};


    DEBUG_PRINT_ZM(P_actual, "%s", "");
    
    matlib_zv u = MK_VM(P);
    matlib_zv u_actual = MK_VM(P);

    double norm_actual = matlib_znrm2(u_actual);
    matlib_zaxpy(-1.0, u_actual, u);
    double e_relative = matlib_znrm2(u)/norm_actual;

    debug_exit("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);
    
}

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
        { "Test double matrix-vector multiplication"  , test_dgemv },
        { "Test complex matrix-vector multiplication" , test_zgemv },
        { "Test double matrix-matrix multiplication"  , test_dgemm },
        { "Test complex matrix-matrix multiplication" , test_zgemm },
        { "Test complex matrices" , test_complex },
        CU_TEST_INFO_NULL,
    };

    /* Create the test suite */ 
    CU_SuiteInfo suites[] = 
    {
        { "Matrix Library", init_suite, clean_suite, test_array },
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
