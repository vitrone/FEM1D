/*============================================================================+/
 | File: Cunit_pthread_fem1d.c 
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
#include <unistd.h>

/* MKL */ 
#include "mkl.h"

#define NDEBUG

#include "legendre.h"
#include "fem1d.h"
#include "pfem1d.h"
#include "pthpool.h"
#include "assert.h"

static const double TOL = 1e-9;
#define MAX_NUM_CPU 7
/*============================================================================*/
void write_dm(char* file_name, matlib_dm M)
{
    matlib_index i, j, col_st, row_st;
    if(M.order == MATLIB_COL_MAJOR)
    {
        col_st = 1;
        row_st = M.lenc;
    }
    else if(M.order == MATLIB_ROW_MAJOR)
    {
        col_st = M.lenr;
        row_st = 1;
    }
    else
    {
        term_exec( "Storage order unknown (order: %d)", M.order);
    }

    FILE *fp = fopen(file_name, "w+");
    for (i=0; i<M.lenc; i++)
    {
        for (j=0; j<M.lenr-1; j++)
        {
            fprintf(fp, "% 0.16f\t", M.elem_p[i*col_st+j*row_st]);
        }
        fprintf(fp, "% 0.16f\n", M.elem_p[i*col_st+j*row_st]);
    }
    fclose(fp);
}

/*============================================================================*/
void Gaussian
( 
    matlib_dv x, 
    matlib_dv u
)
{

    debug_enter("%s", "");
    
    double *xptr;

    for (xptr = x.elem_p; xptr<(x.elem_p+x.len); xptr++)
    {
        *(u.elem_p) = exp(-*xptr**xptr/4.0);
        u.elem_p++;
    }
    debug_exit("%s", "");
}

void test_pfem1d_DFLT
(
    matlib_index p,
    matlib_index num_threads,
    matlib_dm    serial_time,
    matlib_dm    parallel_time
)
{
    debug_enter("polynomial degree: %d, nr. threads: %d", p, num_threads);
    matlib_index num_exp    = serial_time.lenc;
    matlib_index num_cycles = serial_time.lenr-1;

    struct timespec tb, te;
    double dt;

    clock_gettime(CLOCK_REALTIME, &tb);
    
    /* Create pthreads */
    pthpool_data_t mp[num_threads];
    
    pthpool_create_threads(num_threads, mp);

    clock_gettime(CLOCK_REALTIME, &te);
    dt = (double)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;
    debug_print( "Creation of threads[msec]: %0.6f", dt);

    matlib_index i, j;


    matlib_index N, N0 = 200;
    matlib_index P = 4*p;

    /* define the domain */ 
    double x_l = -5.0;
    double x_r =  5.0;

    matlib_dv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_dm FM, IM;

    matlib_create_dm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_dm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    legendre_LGLdataIM( xi, IM);

    matlib_dv x, u, U, V;
    double dim;
    double norm_actual, e_relative;
    double sum = 0;

    for(j=0; j<num_exp; j++)
    {
        /* generate the grid */ 
        N = (j+1)*N0;
        serial_time.elem_p[j]   = (double)N; 
        parallel_time.elem_p[j] = (double)N; 

        fem1d_ref2mesh (xi, N, x_l, x_r, &x);
        debug_body("length of x: %d", x.len);


        /* Provide the initial condition */
        matlib_create_dv( x.len, &u, MATLIB_COL_VECT);

        Gaussian(x, u);
        dim = N*(p+1);
        matlib_create_dv( dim, &U, MATLIB_COL_VECT);
        matlib_create_dv( dim, &V, MATLIB_COL_VECT);

        for(i=0; i<num_cycles; i++)
        {
            clock_gettime(CLOCK_REALTIME, &tb);
            fem1d_DFLT( N, FM, u, U);
            clock_gettime(CLOCK_REALTIME, &te);
            serial_time.elem_p[j+(i+1)*num_exp] = (double)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                  (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            clock_gettime(CLOCK_REALTIME, &tb);
            pfem1d_DFLT(N, FM, u, V, num_threads, mp);
            clock_gettime(CLOCK_REALTIME, &te);
            parallel_time.elem_p[j+(i+1)*num_exp] = (double)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                    (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            /* Analyze error */ 
            //BEGIN_DEBUG
            //    for(i=0; i<U.len; i++)
            //    {
            //        debug_print( "[%d] -> serial: %0.16f, parallel: %0.16f", 
            //                     i, U.elem_p[i], V.elem_p[i]);
            //    }
            //    norm_actual = matlib_dnrm2(U);
            //    matlib_daxpy(-1.0, U, V);
            //    e_relative = matlib_dnrm2(V)/norm_actual;
            //    debug_body("Relative error: % 0.16g", e_relative);
            //END_DEBUG

            
        }
        matlib_free((void*)x.elem_p);
        matlib_free((void*)u.elem_p);
        matlib_free((void*)U.elem_p);
        matlib_free((void*)V.elem_p);
    }
    
    //sleep(2);
    debug_body("%s", "signal threads to exit!");
    clock_gettime(CLOCK_REALTIME, &tb);
    pthpool_destroy_threads(num_threads, mp);
    clock_gettime(CLOCK_REALTIME, &te);
    dt = (double)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

    debug_print( "Exit threads[msec]: %0.6f", dt);
}
/*============================================================================*/

void test_pfem1d_DF2L
(
    matlib_index p,
    matlib_index num_threads,
    matlib_dm    serial_time,
    matlib_dm    parallel_time
)
{
    
    debug_enter("polynomial degree: %d, nr. threads: %d", p, num_threads);
    matlib_index num_exp    = serial_time.lenc;
    matlib_index num_cycles = serial_time.lenr-1;

    struct timespec tb, te;
    double dt;

    clock_gettime(CLOCK_REALTIME, &tb);
    
    /* Create pthreads */
    pthpool_data_t mp[num_threads];
    
    pthpool_create_threads(num_threads, mp);

    clock_gettime(CLOCK_REALTIME, &te);
    dt = (double)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;
    debug_print( "Creation of threads[msec]: %0.6f", dt);

    matlib_index i, j;


    matlib_index N, N0 = 200;
    matlib_index P = 4*p;

    /* define the domain */ 
    double x_l = -5.0;
    double x_r =  5.0;

    matlib_dv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_dm FM, IM;

    matlib_create_dm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_dm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    legendre_LGLdataIM( xi, IM);

    matlib_dv x, u, vb, U, V;
    double dim;
    double norm_actual, e_relative;

    for(j=0; j<num_exp; j++)
    {
        /* generate the grid */ 
        N = (j+1)*N0;
        serial_time.elem_p[j]   = (double)N; 
        parallel_time.elem_p[j] = (double)N; 

        fem1d_ref2mesh (xi, N, x_l, x_r, &x);
        debug_body("length of x: %d", x.len);


        /* Provide the initial condition */
        matlib_create_dv( x.len, &u, MATLIB_COL_VECT);

        Gaussian(x, u);
        dim = N*(p+1);
        matlib_create_dv( dim, &U, MATLIB_COL_VECT);
        matlib_create_dv( dim, &V, MATLIB_COL_VECT);
        matlib_create_dv( N*p+1, &vb, MATLIB_COL_VECT);
        fem1d_DFLT( N, FM, u, U);
        fem1d_DL2F( p, U, vb);

        for(i=0; i<num_cycles; i++)
        {
            clock_gettime(CLOCK_REALTIME, &tb);
            fem1d_DF2L( p, vb, U);
            clock_gettime(CLOCK_REALTIME, &te);
            serial_time.elem_p[j+(i+1)*num_exp] = (double)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                  (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            clock_gettime(CLOCK_REALTIME, &tb);
            pfem1d_DF2L(p, vb, V, num_threads, mp);
            clock_gettime(CLOCK_REALTIME, &te);
            parallel_time.elem_p[j+(i+1)*num_exp] = (double)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                    (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            //BEGIN_DEBUG
            //    for(i=0; i<U.len; i++)
            //    {
            //        debug_print( "[%d] -> serial: %0.16f, parallel: %0.16f", 
            //                     i, U.elem_p[i], V.elem_p[i]);
            //    }
            //    norm_actual = matlib_dnrm2(U);
            //    matlib_daxpy(-1.0, U, V);
            //    e_relative = matlib_dnrm2(V)/norm_actual;
            //    debug_body("Relative error: % 0.16g", e_relative);
            //END_DEBUG
            
        }
        matlib_free((void*)x.elem_p);
        matlib_free((void*)vb.elem_p);
        matlib_free((void*)u.elem_p);
        matlib_free((void*)U.elem_p);
        matlib_free((void*)V.elem_p);
    }
    
    debug_body("%s", "signal threads to exit!");
    clock_gettime(CLOCK_REALTIME, &tb);
    pthpool_destroy_threads(num_threads, mp);
    clock_gettime(CLOCK_REALTIME, &te);
    dt = (double)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

    debug_print( "Exit threads[msec]: %0.6f", dt);

    /* Analyze error */ 

}


void test_pfem1d_DPrjL2F
(
    matlib_index p,
    matlib_index num_threads,
    matlib_dm    serial_time,
    matlib_dm    parallel_time
)
{
    
    debug_enter("polynomial degree: %d, nr. threads: %d", p, num_threads);
    matlib_index num_exp    = serial_time.lenc;
    matlib_index num_cycles = serial_time.lenr-1;

    struct timespec tb, te;
    double dt;

    clock_gettime(CLOCK_REALTIME, &tb);
    
    /* Create pthreads */
    pthpool_data_t mp[num_threads];
    
    pthpool_create_threads(num_threads, mp);

    clock_gettime(CLOCK_REALTIME, &te);
    dt = (double)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;
    debug_print( "Creation of threads[msec]: %0.6f", dt);

    matlib_index i, j;


    matlib_index N, N0 = 200;
    matlib_index P = 4*p;

    /* define the domain */ 
    double x_l = -5.0;
    double x_r =  5.0;

    matlib_dv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_dm FM, IM;

    matlib_create_dm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_dm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    legendre_LGLdataIM( xi, IM);

    matlib_dv x, u, Pvb1, Pvb2, U;
    double dim;
    double norm_actual, e_relative;

    for(j=0; j<num_exp; j++)
    {
        /* generate the grid */ 
        N = (j+1)*N0;
        serial_time.elem_p[j]   = (double)N; 
        parallel_time.elem_p[j] = (double)N; 

        fem1d_ref2mesh (xi, N, x_l, x_r, &x);
        debug_body("length of x: %d", x.len);


        /* Provide the initial condition */
        matlib_create_dv( x.len, &u, MATLIB_COL_VECT);

        Gaussian(x, u);
        dim = N*(p+1);
        matlib_create_dv( dim, &U, MATLIB_COL_VECT);
        matlib_create_dv( N*p+1, &Pvb1, MATLIB_COL_VECT);
        matlib_create_dv( N*p+1, &Pvb2, MATLIB_COL_VECT);
        fem1d_DFLT( N, FM, u, U);

        for(i=0; i<num_cycles; i++)
        {
            clock_gettime(CLOCK_REALTIME, &tb);
            fem1d_DPrjL2F( p, U, Pvb1);
            clock_gettime(CLOCK_REALTIME, &te);
            serial_time.elem_p[j+(i+1)*num_exp] = (double)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                  (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            clock_gettime(CLOCK_REALTIME, &tb);
            pfem1d_DPrjL2F(p, U, Pvb2, num_threads, mp);
            clock_gettime(CLOCK_REALTIME, &te);
            parallel_time.elem_p[j+(i+1)*num_exp] = (double)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                    (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            /* Analyze error */ 
            //BEGIN_DEBUG
            //    for(i=0; i<Pvb1.len; i++)
            //    {
            //        debug_print( "[%d] -> serial: %0.16f, parallel: %0.16f", 
            //                     i, Pvb1.elem_p[i], Pvb2.elem_p[i]);
            //    }
            //    norm_actual = matlib_dnrm2(Pvb1);
            //    matlib_daxpy(-1.0, Pvb1, Pvb2);
            //    e_relative = matlib_dnrm2(Pvb2)/norm_actual;
            //    debug_body("Relative error: % 0.16g", e_relative);
            //END_DEBUG
            
        }
        matlib_free((void*)x.elem_p);
        matlib_free((void*)Pvb1.elem_p);
        matlib_free((void*)Pvb2.elem_p);
        matlib_free((void*)u.elem_p);
        matlib_free((void*)U.elem_p);
    }
    
    debug_body("%s", "signal threads to exit!");
    clock_gettime(CLOCK_REALTIME, &tb);
    pthpool_destroy_threads(num_threads, mp);
    clock_gettime(CLOCK_REALTIME, &te);
    dt = (double)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

    debug_print( "Exit threads[msec]: %0.6f", dt);


}

void test_pfem1d_DNorm2
(
    matlib_index p,
    matlib_index num_threads,
    matlib_dm    serial_time,
    matlib_dm    parallel_time
)
{
    
    debug_enter("polynomial degree: %d, nr. threads: %d", p, num_threads);
    matlib_index num_exp    = serial_time.lenc;
    matlib_index num_cycles = serial_time.lenr-1;

    struct timespec tb, te;
    double dt;

    clock_gettime(CLOCK_REALTIME, &tb);
    
    /* Create pthreads */
    pthpool_data_t mp[num_threads];
    
    pthpool_create_threads(num_threads, mp);

    clock_gettime(CLOCK_REALTIME, &te);
    dt = (double)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;
    debug_print( "Creation of threads[msec]: %0.6f", dt);

    matlib_index i, j;


    matlib_index N, N0 = 200;
    matlib_index P = 4*p;

    /* define the domain */ 
    double x_l = -5.0;
    double x_r =  5.0;

    matlib_dv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_dm FM, IM;

    matlib_create_dm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_dm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    legendre_LGLdataIM( xi, IM);

    matlib_dv x, u, U;
    double dim;
    double norm1, norm2, e_relative;

    for(j=0; j<num_exp; j++)
    {
        /* generate the grid */ 
        N = (j+1)*N0;
        serial_time.elem_p[j]   = (double)N; 
        parallel_time.elem_p[j] = (double)N; 

        fem1d_ref2mesh (xi, N, x_l, x_r, &x);
        debug_body("length of x: %d", x.len);


        /* Provide the initial condition */
        matlib_create_dv( x.len, &u, MATLIB_COL_VECT);

        Gaussian(x, u);
        dim = N*(p+1);
        matlib_create_dv( dim, &U, MATLIB_COL_VECT);
        fem1d_DFLT( N, FM, u, U);

        for(i=0; i<num_cycles; i++)
        {
            clock_gettime(CLOCK_REALTIME, &tb);
            norm1 = fem1d_DNorm2( p, N, U);
            clock_gettime(CLOCK_REALTIME, &te);
            serial_time.elem_p[j+(i+1)*num_exp] = (double)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                  (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            clock_gettime(CLOCK_REALTIME, &tb);
            norm2 = pfem1d_DNorm2(p, N, U, num_threads, mp);
            clock_gettime(CLOCK_REALTIME, &te);
            parallel_time.elem_p[j+(i+1)*num_exp] = (double)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                    (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            /* Analyze error */ 
            //BEGIN_DEBUG
            //    for(i=0; i<Pvb1.len; i++)
            //    {
            //        debug_print( "[%d] -> serial: %0.16f, parallel: %0.16f", 
            //                     i, Pvb1.elem_p[i], Pvb2.elem_p[i]);
            //    }
            //    norm_actual = matlib_dnrm2(Pvb1);
            //    matlib_daxpy(-1.0, Pvb1, Pvb2);
            //    e_relative = matlib_dnrm2(Pvb2)/norm_actual;
            //    debug_body("Relative error: % 0.16g", e_relative);
            //END_DEBUG
            
        }
        matlib_free((void*)x.elem_p);
        matlib_free((void*)u.elem_p);
        matlib_free((void*)U.elem_p);
    }
    
    debug_body("%s", "signal threads to exit!");
    clock_gettime(CLOCK_REALTIME, &tb);
    pthpool_destroy_threads(num_threads, mp);
    clock_gettime(CLOCK_REALTIME, &te);
    dt = (double)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

    debug_print( "Exit threads[msec]: %0.6f", dt);


}


/*============================================================================*/
void zGaussian
( 
    matlib_dv x, 
    matlib_zv u
)
{

    debug_enter("%s", "");
    
    double *xptr;

    for (xptr = x.elem_p; xptr<(x.elem_p+x.len); xptr++)
    {
        *(u.elem_p) = cexp(-(1.0+I*1.0)**xptr**xptr/4.0);
        u.elem_p++;
    }
    debug_exit("%s", "");
}


void test_pfem1d_ZFLT
(
    matlib_index p,
    matlib_index num_threads,
    matlib_dm    serial_time,
    matlib_dm    parallel_time
)
{
    
    debug_enter("polynomial degree: %d, nr. threads: %d", p, num_threads);
    matlib_index num_exp    = serial_time.lenc;
    matlib_index num_cycles = serial_time.lenr-1;

    struct timespec tb, te;
    double dt;

    clock_gettime(CLOCK_REALTIME, &tb);
    
    /* Create pthreads */
    pthpool_data_t mp[num_threads];
    
    pthpool_create_threads(num_threads, mp);

    clock_gettime(CLOCK_REALTIME, &te);
    dt = (double)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;
    debug_print( "Creation of threads[msec]: %0.6f", dt);

    matlib_index i, j;


    matlib_index N, N0 = 200;
    matlib_index P = 4*p;

    /* define the domain */ 
    double x_l = -5.0;
    double x_r =  5.0;

    matlib_dv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_dm FM, IM;

    matlib_create_dm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_dm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    legendre_LGLdataIM( xi, IM);

    matlib_dv x;
    matlib_zv u, U, V;
    double dim;
    double norm_actual, e_relative;

    for(j=0; j<num_exp; j++)
    {
        /* generate the grid */ 
        N = (j+1)*N0;
        serial_time.elem_p[j]   = (double)N; 
        parallel_time.elem_p[j] = (double)N; 

        fem1d_ref2mesh (xi, N, x_l, x_r, &x);
        debug_body("length of x: %d", x.len);


        /* Provide the initial condition */
        matlib_create_zv( x.len, &u, MATLIB_COL_VECT);

        zGaussian(x, u);
        dim = N*(p+1);
        matlib_create_zv( dim, &U, MATLIB_COL_VECT);
        matlib_create_zv( dim, &V, MATLIB_COL_VECT);

        for(i=0; i<num_cycles; i++)
        {
            clock_gettime(CLOCK_REALTIME, &tb);
            fem1d_ZFLT( N, FM, u, U);
            clock_gettime(CLOCK_REALTIME, &te);
            serial_time.elem_p[j+(i+1)*num_exp] = (double)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                  (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            clock_gettime(CLOCK_REALTIME, &tb);
            pfem1d_ZFLT(N, FM, u, V, num_threads, mp);
            clock_gettime(CLOCK_REALTIME, &te);
            parallel_time.elem_p[j+(i+1)*num_exp] = (double)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                    (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            /* Analyze error */ 
            //BEGIN_DEBUG
            //    for(i=0; i<U.len; i++)
            //    {
            //        debug_print( "[%d] -> serial: %0.16f, parallel: %0.16f", 
            //                     i, U.elem_p[i], V.elem_p[i]);
            //    }
            //    norm_actual = matlib_znrm2(U);
            //    matlib_zaxpy(-1.0, U, V);
            //    e_relative = matlib_znrm2(V)/norm_actual;
            //    debug_body("Relative error: % 0.16g", e_relative);
            //END_DEBUG
            
        }
        matlib_free((void*)x.elem_p);
        matlib_free((void*)u.elem_p);
        matlib_free((void*)U.elem_p);
        matlib_free((void*)V.elem_p);
    }
    
    debug_body("%s", "signal threads to exit!");
    clock_gettime(CLOCK_REALTIME, &tb);
    pthpool_destroy_threads(num_threads, mp);
    clock_gettime(CLOCK_REALTIME, &te);
    dt = (double)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

    debug_print( "Exit threads[msec]: %0.6f", dt);


}

void test_pfem1d_ZF2L
(
    matlib_index p,
    matlib_index num_threads,
    matlib_dm    serial_time,
    matlib_dm    parallel_time
)
{
    
    debug_enter("polynomial degree: %d, nr. threads: %d", p, num_threads);
    matlib_index num_exp    = serial_time.lenc;
    matlib_index num_cycles = serial_time.lenr-1;

    struct timespec tb, te;
    double dt;

    clock_gettime(CLOCK_REALTIME, &tb);
    
    /* Create pthreads */
    pthpool_data_t mp[num_threads];
    
    pthpool_create_threads(num_threads, mp);

    clock_gettime(CLOCK_REALTIME, &te);
    dt = (double)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;
    debug_print( "Creation of threads[msec]: %0.6f", dt);

    matlib_index i, j;

    matlib_index N, N0 = 200;
    matlib_index P = 4*p;

    /* define the domain */ 
    double x_l = -5.0;
    double x_r =  5.0;

    matlib_dv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_dm FM, IM;

    matlib_create_dm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_dm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    legendre_LGLdataIM( xi, IM);

    matlib_dv x;
    matlib_zv u, U, V, vb;
    double dim;
    double norm_actual, e_relative;

    for(j=0; j<num_exp; j++)
    {
        /* generate the grid */ 
        N = (j+1)*N0;
        serial_time.elem_p[j]   = (double)N; 
        parallel_time.elem_p[j] = (double)N; 

        fem1d_ref2mesh (xi, N, x_l, x_r, &x);
        debug_body("length of x: %d", x.len);


        /* Provide the initial condition */
        matlib_create_zv( x.len, &u, MATLIB_COL_VECT);

        zGaussian(x, u);
        //DEBUG_PRINT_ZV(u, "%s: ", "complex vector");
        dim = N*(p+1);
        matlib_create_zv( dim, &U, MATLIB_COL_VECT);
        matlib_create_zv( dim, &V, MATLIB_COL_VECT);
        matlib_create_zv( N*p+1, &vb, MATLIB_COL_VECT);
        fem1d_ZFLT( N, FM, u, U);
        fem1d_ZL2F( p, U, vb);

        for(i=0; i<num_cycles; i++)
        {
            clock_gettime(CLOCK_REALTIME, &tb);
            fem1d_ZF2L( p, vb, U);
            clock_gettime(CLOCK_REALTIME, &te);
            serial_time.elem_p[j+(i+1)*num_exp] = (double)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                  (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            clock_gettime(CLOCK_REALTIME, &tb);
            pfem1d_ZF2L(p, vb, V, num_threads, mp);
            clock_gettime(CLOCK_REALTIME, &te);
            parallel_time.elem_p[j+(i+1)*num_exp] = (double)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                    (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            //BEGIN_DEBUG
            //    for(i=0; i<U.len; i++)
            //    {
            //        debug_print( "[%d] -> serial: %0.16f% 0.16f, parallel: %0.16f% 0.16f", 
            //                     i, U.elem_p[i], V.elem_p[i]);
            //    }
            //    norm_actual = matlib_znrm2(U);
            //    matlib_zaxpy(-1.0, U, V);
            //    e_relative = matlib_znrm2(V)/norm_actual;
            //    debug_body("Relative error: % 0.16g", e_relative);
            //END_DEBUG
            
        }
        matlib_free((void*)x.elem_p);
        matlib_free((void*)u.elem_p);
        matlib_free((void*)vb.elem_p);
        matlib_free((void*)U.elem_p);
        matlib_free((void*)V.elem_p);
    }
    
    debug_body("%s", "signal threads to exit!");
    clock_gettime(CLOCK_REALTIME, &tb);
    pthpool_destroy_threads(num_threads, mp);
    clock_gettime(CLOCK_REALTIME, &te);
    dt = (double)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

    debug_print( "Exit threads[msec]: %0.6f", dt);

    /* Analyze error */ 

}

void test_pfem1d_ZPrjL2F
(
    matlib_index p,
    matlib_index num_threads,
    matlib_dm    serial_time,
    matlib_dm    parallel_time
)
{
    debug_enter("polynomial degree: %d, nr. threads: %d", p, num_threads);
    matlib_index num_exp    = serial_time.lenc;
    matlib_index num_cycles = serial_time.lenr-1;

    struct timespec tb, te;
    double dt;

    clock_gettime(CLOCK_REALTIME, &tb);
    
    /* Create pthreads */
    pthpool_data_t mp[num_threads];
    
    pthpool_create_threads(num_threads, mp);

    clock_gettime(CLOCK_REALTIME, &te);
    dt = (double)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;
    debug_print( "Creation of threads[msec]: %0.6f", dt);

    matlib_index i, j;


    matlib_index N, N0 = 200;
    matlib_index P = 4*p;

    /* define the domain */ 
    double x_l = -5.0;
    double x_r =  5.0;

    matlib_dv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_dm FM, IM;

    matlib_create_dm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_dm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    legendre_LGLdataIM( xi, IM);

    matlib_dv x;
    matlib_zv u, U, Pvb1, Pvb2;
    double dim;
    double norm_actual, e_relative;

    for(j=0; j<num_exp; j++)
    {
        /* generate the grid */ 
        N = (j+1)*N0;
        serial_time.elem_p[j]   = (double)N; 
        parallel_time.elem_p[j] = (double)N; 

        fem1d_ref2mesh (xi, N, x_l, x_r, &x);
        debug_body("length of x: %d", x.len);


        /* Provide the initial condition */
        matlib_create_zv( x.len, &u, MATLIB_COL_VECT);

        zGaussian(x, u);
        //DEBUG_PRINT_ZV(u, "%s: ", "complex vector");
        dim = N*(p+1);
        matlib_create_zv( dim, &U, MATLIB_COL_VECT);
        matlib_create_zv( N*p+1, &Pvb1, MATLIB_COL_VECT);
        matlib_create_zv( N*p+1, &Pvb2, MATLIB_COL_VECT);
        fem1d_ZFLT( N, FM, u, U);

        for(i=0; i<num_cycles; i++)
        {
            clock_gettime(CLOCK_REALTIME, &tb);
            fem1d_ZPrjL2F( p, U, Pvb1);
            clock_gettime(CLOCK_REALTIME, &te);
            serial_time.elem_p[j+(i+1)*num_exp] = (double)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                  (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            clock_gettime(CLOCK_REALTIME, &tb);
            pfem1d_ZPrjL2F(p, U, Pvb2, num_threads, mp);
            clock_gettime(CLOCK_REALTIME, &te);
            parallel_time.elem_p[j+(i+1)*num_exp] = (double)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                    (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            //BEGIN_DEBUG
            //    for(i=0; i<Pvb1.len; i++)
            //    {
            //        debug_print( "[%d] -> serial: %0.16f% 0.16f, parallel: %0.16f% 0.16f", 
            //                     i, Pvb1.elem_p[i], Pvb2.elem_p[i]);
            //    }
            //    norm_actual = matlib_znrm2(Pvb1);
            //    matlib_zaxpy(-1.0, Pvb1, Pvb2);
            //    e_relative = matlib_znrm2(Pvb2)/norm_actual;
            //    debug_body("Relative error: % 0.16g", e_relative);
            //END_DEBUG
            
        }
        matlib_free((void*)x.elem_p);
        matlib_free((void*)u.elem_p);
        matlib_free((void*)Pvb1.elem_p);
        matlib_free((void*)Pvb2.elem_p);
        matlib_free((void*)U.elem_p);
    }
    
    debug_body("%s", "signal threads to exit!");
    clock_gettime(CLOCK_REALTIME, &tb);
    pthpool_destroy_threads(num_threads, mp);
    clock_gettime(CLOCK_REALTIME, &te);
    dt = (double)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

    debug_print( "Exit threads[msec]: %0.6f", dt);

    /* Analyze error */ 

}

void test_pfem1d_ZNorm2
(
    matlib_index p,
    matlib_index num_threads,
    matlib_dm    serial_time,
    matlib_dm    parallel_time
)
{
    debug_enter("polynomial degree: %d, nr. threads: %d", p, num_threads);
    matlib_index num_exp    = serial_time.lenc;
    matlib_index num_cycles = serial_time.lenr-1;

    struct timespec tb, te;
    double dt;

    clock_gettime(CLOCK_REALTIME, &tb);
    
    /* Create pthreads */
    pthpool_data_t mp[num_threads];
    
    pthpool_create_threads(num_threads, mp);

    clock_gettime(CLOCK_REALTIME, &te);
    dt = (double)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;
    debug_print( "Creation of threads[msec]: %0.6f", dt);

    matlib_index i, j;


    matlib_index N, N0 = 200;
    matlib_index P = 4*p;

    /* define the domain */ 
    double x_l = -5.0;
    double x_r =  5.0;

    matlib_dv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_dm FM, IM;

    matlib_create_dm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_dm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    legendre_LGLdataIM( xi, IM);

    matlib_dv x;
    matlib_zv u, U;
    double dim;
    double norm1, norm2, e_relative;

    for(j=0; j<num_exp; j++)
    {
        /* generate the grid */ 
        N = (j+1)*N0;
        serial_time.elem_p[j]   = (double)N; 
        parallel_time.elem_p[j] = (double)N; 

        fem1d_ref2mesh (xi, N, x_l, x_r, &x);
        debug_body("length of x: %d", x.len);


        /* Provide the initial condition */
        matlib_create_zv( x.len, &u, MATLIB_COL_VECT);

        zGaussian(x, u);
        //DEBUG_PRINT_ZV(u, "%s: ", "complex vector");
        dim = N*(p+1);
        matlib_create_zv( dim, &U, MATLIB_COL_VECT);
        fem1d_ZFLT( N, FM, u, U);

        for(i=0; i<num_cycles; i++)
        {
            clock_gettime(CLOCK_REALTIME, &tb);
            norm1 = fem1d_ZNorm2( p, N, U);
            clock_gettime(CLOCK_REALTIME, &te);
            serial_time.elem_p[j+(i+1)*num_exp] = (double)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                  (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            clock_gettime(CLOCK_REALTIME, &tb);
            norm2 = pfem1d_ZNorm2( p, N, U, num_threads, mp);
            clock_gettime(CLOCK_REALTIME, &te);
            parallel_time.elem_p[j+(i+1)*num_exp] = (double)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                    (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            BEGIN_DEBUG
                e_relative = fabs(norm1-norm2)/norm1;
                debug_body("Relative error: % 0.16g", e_relative);
            END_DEBUG
            
        }
        matlib_free((void*)x.elem_p);
        matlib_free((void*)u.elem_p);
        matlib_free((void*)U.elem_p);
    }
    
    debug_body("%s", "signal threads to exit!");
    clock_gettime(CLOCK_REALTIME, &tb);
    pthpool_destroy_threads(num_threads, mp);
    clock_gettime(CLOCK_REALTIME, &te);
    dt = (double)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (double)(te.tv_nsec-tb.tv_nsec)/1.0e6;

    debug_print( "Exit threads[msec]: %0.6f", dt);

    /* Analyze error */ 

}

/*============================================================================*/

void test_performance
(
    matlib_index p,
    void         (*fp)(matlib_index, matlib_index, matlib_dm, matlib_dm),
    char*        path
) 
{
    /* define number of timing experiments and repeatitions */ 
    matlib_index num_exp    = 50;
    matlib_index num_cycles = 1000;
    
    matlib_index num_threads_max = 7;

    char file_name[80];

    matlib_dm serial_time, parallel_time;

    matlib_create_dm( num_exp, 
                      num_cycles+1, 
                      &serial_time, 
                      MATLIB_COL_MAJOR, 
                      MATLIB_NO_TRANS);
    matlib_create_dm( num_exp, 
                      num_cycles+1, 
                      &parallel_time, 
                      MATLIB_COL_MAJOR, 
                      MATLIB_NO_TRANS);

    for(matlib_index num_threads=2; num_threads < num_threads_max; num_threads++)
    {
        (*fp)(p, num_threads, serial_time, parallel_time);
        
        sprintf( file_name, "%s/serial_time_num_threads%d_p%02d.dat",
                 path, num_threads, p);
        debug_body("file name: %s", file_name);

        write_dm( file_name, serial_time);

        sprintf( file_name, "%s/parallel_time_num_threads%d_p%02d.dat",
                 path, num_threads, p);
        debug_body("file name: %s", file_name);

        write_dm( file_name, parallel_time);
        debug_body("%s","data written to files.");
    }
    
    matlib_free((void*)serial_time.elem_p);
    matlib_free((void*)parallel_time.elem_p);

    debug_exit("%s", "");

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

    char path[20];
    matlib_index p = 15;

    //strcpy( path, "pfem1d_DFLT");
    //test_performance(p, test_pfem1d_DFLT, path);

    //strcpy( path, "pfem1d_DF2L");
    //test_performance(p, test_pfem1d_DF2L, path);
    //
    //strcpy( path, "pfem1d_DPrjL2F");
    //test_performance(p, test_pfem1d_DPrjL2F, path);

    strcpy( path, "pfem1d_DNorm2");
    test_performance(p, test_pfem1d_DNorm2, path);
    //
    //strcpy( path, "pfem1d_ZFLT");
    //test_performance(p, test_pfem1d_ZFLT, path);

    //strcpy( path, "pfem1d_ZF2L");
    //test_performance(p, test_pfem1d_ZF2L, path);
    //
    //strcpy( path, "pfem1d_ZPrjL2F");
    //test_performance(p, test_pfem1d_ZPrjL2F, path);

    strcpy( path, "pfem1d_ZNorm2");
    test_performance(p, test_pfem1d_ZNorm2, path);
    return(0);
}

