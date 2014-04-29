#include <pthread.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>

/* MKL */ 
#include "mkl.h"
#include "omp.h"

#define NDEBUG

#include "legendre.h"
#include "fem1d.h"
#include "pfem1d.h"
#include "pthpool.h"
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

void serial_Gaussian
(
    matlib_xv x,
    matlib_zv u
)
{
    debug_enter("%s", "");
    matlib_real *xptr;
    for (xptr = x.elem_p; xptr<(x.elem_p+x.len); xptr++)
    {
        *(u.elem_p) = cexp(-(1.0+1.0*I)**xptr**xptr);
        u.elem_p++;
    }
    debug_exit("%s", "");
}


void* thfunc_Gaussian(void* mp)
{

    debug_enter("%s", "");
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;
    matlib_xv x = *((matlib_xv*) (ptr->shared_data[0]));
    matlib_zv u = *((matlib_zv*) (ptr->shared_data[1]));

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    

    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);
    
    (u.elem_p) += (start_end_index[0]);

    for(i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        *(u.elem_p) = cexp(-(1.0+1.0*I)*x.elem_p[i]*x.elem_p[i]);
        u.elem_p++;
    }
    debug_exit("%s", "");
}

void test_pfunc(void)
{

    struct timespec tb, te;
    matlib_real dn;
    
    START_TIMMING(tb);

    matlib_index num_threads = 4;
    pthpool_data_t mp[num_threads];
    
    pthpool_create_threads(num_threads, mp);

    /* define the domain */ 
    matlib_real x_l = -5.0;
    matlib_real x_r =  5.0;

    matlib_index p = 4;
    matlib_xv xi, quadW;
    legendre_LGLdataLT1( p, TOL, &xi, &quadW);

    matlib_index N = 2000;
    matlib_xv x;
    fem1d_ref2mesh (xi, N, x_l, x_r, &x);
    debug_body("length of x: %d", x.len);
    
    matlib_zv u_serial, u_parallel;
    matlib_create_zv( x.len, &u_serial, MATLIB_COL_VECT);
    matlib_create_zv( x.len, &u_parallel, MATLIB_COL_VECT);

    GET_DURATION(tb, te, dn);
    debug_print("Initialization time[msec]: %0.4f", dn);
    
    START_TIMMING(tb);
    serial_Gaussian(x, u_serial);
    GET_DURATION(tb, te, dn);
    debug_print("Serial time[msec]: %0.4f", dn);

    void* param = {0};
    matlib_index Np[2] = {x.len/num_threads, x.len};
    void* shared_data[3] = { (void*) &x,
                             (void*) &u_parallel,
                             (void*) param};
    START_TIMMING(tb);
    pthpool_func( Np, shared_data, thfunc_Gaussian, num_threads, mp);
    GET_DURATION(tb, te, dn);
    debug_print("parallel time[msec]: %0.4f", dn);

    matlib_real norm_actual = matlib_znrm2(u_serial);
    matlib_zaxpy(-1.0, u_serial, u_parallel);
    matlib_real e_relative = matlib_znrm2(u_parallel)/norm_actual;
    debug_body("Relative error: % 0.16g", e_relative);
    
    CU_ASSERT_TRUE(e_relative<TOL);

    matlib_free(x.elem_p);
    matlib_free(u_serial.elem_p);
    matlib_free(u_parallel.elem_p);
    
    debug_body("%s", "signal threads to exit!");
    pthpool_destroy_threads(num_threads, mp);

}
/*============================================================================+/
 | Test runner
 |
 |
 +============================================================================*/

int main(void)
{

    //mkl_set_num_threads(1);
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
        { "Parallel Gaussian"      , test_pfunc },
        CU_TEST_INFO_NULL,
    };

    /* Create the test suite */ 
    CU_SuiteInfo suites[] = 
    {
        { "Parallel Evaluation of Expressions", init_suite, clean_suite, test_array },
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

