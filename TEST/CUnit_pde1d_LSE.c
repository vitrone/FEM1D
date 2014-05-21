/*============================================================================+/
 | File: Cunit_pde1d_LSE.c 
 | Description: Test 
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
#include "omp.h"
#include "mkl_pardiso.h"

//#define NDEBUG
//#define MATLIB_NTRACE_DATA

#include "pde1d_solver.h"
#include "debug.h"
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


/*============================================================================*/

void test_pde1d_LSE_solve_IVP_evol(void)
{
    debug_enter("%s", "");

    pde1d_LSE_data_t  input;
    pde1d_LSE_solver_t data;

    pde1d_LSE_set_defaultsIVP(&input);
    input.sol_mode = PDE1D_LSE_EVOLVE_ONLY;
    
    pde1d_LSE_init_solverIVP(&input, &data);
    
    pde1d_LSE_set_potential( &input, PDE1D_LSE_STATIC, 
                             (void*)pde1d_LSE_constant_potential);
    
    matlib_complex A_0 = 1.0 + I*0.0;
    matlib_complex a = 0.5 + I*0.5;
    matlib_real c = 1;
    matlib_complex phi_0 = 1.0 + I*0.0;

    void* params[4] = { (void*)&A_0, 
                        (void*)&a, 
                        (void*)&c, 
                        (void*)&phi_0};

    input.params = params;
    /* Provide the initial condition 
     * */ 
    pde1d_LSE_Gaussian_WP_constant_potential( params, 
                                              input.x, 
                                              (input.t.elem_p)[0],
                                              input.u_init);
    
    pde1d_LSE_solve_IVP(&input, &data);

    matlib_zm u_evol;
    matlib_create_zm( (input.x).len, 
                      (input.t).len, 
                      &(u_evol), 
                      MATLIB_COL_MAJOR, MATLIB_NO_TRANS);

    fem1d_ZILT2(input.N, data.IM, input.U_evol, u_evol);

    /* Write the final evolution data to a file 
     * */ 
    matlib_zmwrite_csv( "evol_data.dat", u_evol);
    /* Get the analytic evolution data
     * */ 
    pde1d_analytic_evol( params, pde1d_LSE_Gaussian_WP_constant_potential,
                         input.x, input.t, u_evol);
    matlib_zmwrite_csv( "evol_data_analytic.dat", u_evol);


    matlib_xvwrite_csv("x_grid.dat", 1, &input.x);
    matlib_xvwrite_csv("t_grid.dat", 1, &input.t);

    CU_ASSERT_TRUE(true);
    pde1d_LSE_destroy_solverIVP(&input, &data);

    debug_exit("%s", "");
}


void test_pde1d_LSE_solve_IVP_evol2(void)
{
    debug_enter("%s", "");

    pde1d_LSE_data_t  input;
    pde1d_LSE_solver_t data;

    pde1d_LSE_set_defaultsIVP(&input);
    input.sol_mode = PDE1D_LSE_EVOLVE_ONLY;

    pde1d_LSE_init_solverIVP(&input, &data);
    pde1d_LSE_set_potential( &input, PDE1D_LSE_STATIC, 
                             (void*)pde1d_LSE_harmonic_potential);

    void** params = NULL;
    input.params = NULL;
    /* Provide the initial condition 
     * */ 
    pde1d_LSE_HermiteGaussian_WP_harmonic_potential( params, 
                                              input.x, 
                                              (input.t.elem_p)[0],
                                              input.u_init);
    
    pde1d_LSE_solve_IVP(&input, &data);

    matlib_zm u_evol;
    matlib_create_zm( (input.x).len, 
                      (input.t).len, 
                      &(u_evol), 
                      MATLIB_COL_MAJOR, MATLIB_NO_TRANS);

    fem1d_ZILT2(input.N, data.IM, input.U_evol, u_evol);

    /* Write the final evolution data to a file 
     * */ 
    matlib_zmwrite_csv( "evol_data.dat", u_evol);
    /* Get the analytic evolution data
     * */ 
    pde1d_analytic_evol( params, pde1d_LSE_HermiteGaussian_WP_harmonic_potential,
                         input.x, input.t, u_evol);
    matlib_zmwrite_csv( "evol_data_analytic.dat", u_evol);


    matlib_xvwrite_csv("x_grid.dat", 1, &input.x);
    matlib_xvwrite_csv("t_grid.dat", 1, &input.t);

    CU_ASSERT_TRUE(true);
    pde1d_LSE_destroy_solverIVP(&input, &data);

    debug_exit("%s", "");
}

void test_pde1d_LSE_solve_IVP_evol3(void)
{
    debug_enter("%s", "");

    pde1d_LSE_data_t  input;
    pde1d_LSE_solver_t data;

    pde1d_LSE_set_defaultsIVP(&input);
    input.sol_mode = PDE1D_LSE_EVOLVE_ONLY;

    input.Nt = 2000;
    pde1d_LSE_init_solverIVP(&input, &data);
    pde1d_LSE_set_potential( &input, PDE1D_LSE_DYNAMIC, 
                             (void*)pde1d_LSE_timedependent_linear_potential);

    matlib_complex A_0 = 1.0 + I*0.0;
    matlib_complex a = 0.5 + I*0.5;
    matlib_real c = 0.5;
    matlib_complex phi_0 = 1.0 + I*0.0;
    matlib_real g_0 = 1.0;
    matlib_real mu  = 2.0*M_PI;

    void* params[6] = { (void*)&A_0, 
                        (void*)&a, 
                        (void*)&c, 
                        (void*)&phi_0,
                        (void*)&g_0,
                        (void*)&mu};

    input.params = params;
    /* Provide the initial condition 
     * */ 
    pde1d_LSE_Gaussian_WP_timedependent_linear_potential( params, 
                                              input.x, 
                                              (input.t.elem_p)[0],
                                              input.u_init);
    
    pde1d_LSE_solve_IVP2_evol(&input, &data);

    matlib_zm u_evol;
    matlib_create_zm( (input.x).len, 
                      (input.t).len, 
                      &(u_evol), 
                      MATLIB_COL_MAJOR, MATLIB_NO_TRANS);

    fem1d_ZILT2(input.N, data.IM, input.U_evol, u_evol);

    /* Write the final evolution data to a file 
     * */ 
    matlib_zmwrite_csv( "evol_data.dat", u_evol);
    /* Get the analytic evolution data
     * */ 
    pde1d_analytic_evol( params, pde1d_LSE_Gaussian_WP_timedependent_linear_potential,
                         input.x, input.t, u_evol);
    matlib_zmwrite_csv( "evol_data_analytic.dat", u_evol);


    matlib_xvwrite_csv("x_grid.dat", 1, &input.x);
    matlib_xvwrite_csv("t_grid.dat", 1, &input.t);

    CU_ASSERT_TRUE(true);
    pde1d_LSE_destroy_solverIVP(&input, &data);

    debug_exit("%s", "");
}
/*============================================================================*/
void test_pde1d_LSE_solve_IVP_error(void)
{
    debug_enter("%s", "");

    pde1d_LSE_data_t  input;
    pde1d_LSE_solver_t data;
    pde1d_LSE_set_defaultsIVP(&input);

    input.domain[0] = -30;
    input.domain[1] =  30;

    input.N  = 2000;
    input.Nt = 4000;

    input.sol_mode = PDE1D_LSE_ERROR_ONLY;
    pde1d_LSE_init_solverIVP(&input, &data);
    pde1d_LSE_set_potential( &input, PDE1D_LSE_STATIC, 
                             (void*)pde1d_LSE_constant_potential);
    
    matlib_complex A_0 = 1.0;
    matlib_complex a = 0.5 + I*0.5;
    matlib_real c = 0.5;
    matlib_complex phi_0 = 1.0;

    void* params[4] = { (void*)&A_0, 
                        (void*)&a, 
                        (void*)&c, 
                        (void*)&phi_0};

    input.params = params;
    input.u_analytic = pde1d_LSE_Gaussian_WP_constant_potential;

    pde1d_LSE_solve_IVP(&input, &data);


    matlib_xv tmp_v[3] = {(input.t), (input.e_rel), (input.e_abs)};
    matlib_xvwrite_csv("error_evol.dat", 3, tmp_v);
        
    matlib_real r = input.e_rel.elem_p[(input.Nt)];
    debug_body("Relative Error: %0.16g", r);

    CU_ASSERT_TRUE(r<4.0e-6);

    pde1d_LSE_destroy_solverIVP(&input, &data);
    debug_exit("%s", "");
}

void test_pde1d_LSE_solve_IVP_error3(void)
{
    debug_enter("%s", "");

    pde1d_LSE_data_t  input;
    pde1d_LSE_solver_t data;
    pde1d_LSE_set_defaultsIVP(&input);

    input.domain[0] = -30;
    input.domain[1] =  30;

    input.N  = 2000;
    input.Nt = 4000;

    input.sol_mode = PDE1D_LSE_ERROR_ONLY;
    pde1d_LSE_init_solverIVP(&input, &data);
    pde1d_LSE_set_potential( &input, PDE1D_LSE_DYNAMIC, 
                             (void*)pde1d_LSE_timedependent_linear_potential);
    
    matlib_complex A_0 = 1.0;
    matlib_complex a = 0.5 + I*0.5;
    matlib_real c = 0.5;
    matlib_complex phi_0 = 1.0;
    matlib_real g_0 = 1;
    matlib_real mu  = 2.0*M_PI;

    void* params[6] = { (void*)&A_0, 
                        (void*)&a, 
                        (void*)&c, 
                        (void*)&phi_0,
                        (void*)&g_0,
                        (void*)&mu};

    input.params = params;
    input.u_analytic = pde1d_LSE_Gaussian_WP_timedependent_linear_potential;

    pde1d_LSE_solve_IVP2_error(&input, &data);

    matlib_xv tmp_v[3] = {(input.t), (input.e_rel), (input.e_abs)};
    matlib_xvwrite_csv("error_evol.dat", 3, tmp_v);
        
    matlib_real r = input.e_rel.elem_p[(input.Nt)];
    debug_body("Relative Error: %0.16g", r);

    CU_ASSERT_TRUE(r<4.0e-6);

    pde1d_LSE_destroy_solverIVP(&input, &data);
    debug_exit("%s", "");
}
/*============================================================================+/
 | Test runner
 |
 |
 +============================================================================*/

int main(void)
{
    debug_enter("%s", "");
    mkl_domain_set_num_threads(4, MKL_BLAS);
    CU_pSuite pSuite = NULL;

    /* initialize the CUnit test registry */
    if (CUE_SUCCESS != CU_initialize_registry())
    {
        return CU_get_error();
    }

    /* Create a test array */
    CU_TestInfo test_array[] = 
    {
        //{ "Constant potential evolve", test_pde1d_LSE_solve_IVP_evol},
        //{ "Harmonic potential evolve", test_pde1d_LSE_solve_IVP_evol2},
        //{ "Linear time-dependent potential evolve", test_pde1d_LSE_solve_IVP_evol3},
        { "Constant potential error" , test_pde1d_LSE_solve_IVP_error},
        { "Linear time-dependent potential error" , test_pde1d_LSE_solve_IVP_error3},
        CU_TEST_INFO_NULL,
    };

    /* Create the test suite */ 
    CU_SuiteInfo suites[] = 
    {
        { "PDE-1D LSE", init_suite, clean_suite, test_array },
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

