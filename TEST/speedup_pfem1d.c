/*============================================================================+/
 | File: Cunit_pfem1d.c 
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

#include <sys/stat.h> /* for mkdir */ 

/* MKL */ 
#include "mkl.h"

#define NDEBUG
#define MATLIB_NTRACE_DATA

#include "legendre.h"
#include "fem1d.h"
#include "pfem1d.h"
#include "assert.h"

static const matlib_real TOL = 1e-9;

/*============================================================================*/
void Gaussian
( 
    matlib_xv x, 
    matlib_xv u
)
{

    debug_enter("%s", "");
    
    matlib_real *xptr;

    for (xptr = x.elem_p; xptr<(x.elem_p+x.len); xptr++)
    {
        *(u.elem_p) = exp(-*xptr**xptr/4.0);
        u.elem_p++;
    }
    debug_exit("%s", "");
}

void test_pfem1d_XFLT
(
    matlib_index p,
    matlib_index num_threads,
    matlib_xm    serial_time,
    matlib_xm    parallel_time
)
{
    debug_enter("polynomial degree: %d, nr. threads: %d", p, num_threads);
    matlib_index num_exp    = serial_time.lenc;
    matlib_index num_cycles = serial_time.lenr-1;

    struct timespec tb, te;
    matlib_real dt;

    clock_gettime(CLOCK_REALTIME, &tb);
    
    /* Create pthreads */
    pthpool_data_t mp[num_threads];
    
    pthpool_create_threads(num_threads, mp);

    clock_gettime(CLOCK_REALTIME, &te);
    dt = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;
    debug_print( "Creation of threads[msec]: %0.6f", dt);

    matlib_index i, j;


    matlib_index N, N0 = 200;
    matlib_index P = 4*p;

    /* define the domain */ 
    matlib_real x_l = -5.0;
    matlib_real x_r =  5.0;

    matlib_xv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_xm FM, IM;

    matlib_create_xm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_xm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    legendre_LGLdataIM( xi, IM);

    matlib_xv x, u, U, V;
    matlib_real dim;
    matlib_real norm_actual, e_relative;
    matlib_real sum = 0;

    for(j=0; j<num_exp; j++)
    {
        /* generate the grid */ 
        N = (j+1)*N0;
        serial_time.elem_p[j]   = (matlib_real)N; 
        parallel_time.elem_p[j] = (matlib_real)N; 

        fem1d_ref2mesh (xi, N, x_l, x_r, &x);
        debug_body("length of x: %d", x.len);


        /* Provide the initial condition */
        matlib_create_xv( x.len, &u, MATLIB_COL_VECT);

        Gaussian(x, u);
        dim = N*(p+1);
        matlib_create_xv( dim, &U, MATLIB_COL_VECT);
        matlib_create_xv( dim, &V, MATLIB_COL_VECT);

        for(i=0; i<num_cycles; i++)
        {
            clock_gettime(CLOCK_REALTIME, &tb);
            fem1d_XFLT( N, FM, u, U);
            clock_gettime(CLOCK_REALTIME, &te);
            serial_time.elem_p[j+(i+1)*num_exp] = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                  (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            clock_gettime(CLOCK_REALTIME, &tb);
            pfem1d_XFLT(N, FM, u, V, num_threads, mp);
            clock_gettime(CLOCK_REALTIME, &te);
            parallel_time.elem_p[j+(i+1)*num_exp] = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                    (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            BEGIN_DTRACE
                for(i=0; i<U.len; i++)
                {
                    debug_print( "[%d] -> serial: %0.16f, parallel: %0.16f", 
                                 i, U.elem_p[i], V.elem_p[i]);
                }
            END_DTRACE

            
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
    dt = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

    debug_print( "Exit threads[msec]: %0.6f", dt);
}
/*============================================================================*/

void test_pfem1d_XF2L
(
    matlib_index p,
    matlib_index num_threads,
    matlib_xm    serial_time,
    matlib_xm    parallel_time
)
{
    
    debug_enter("polynomial degree: %d, nr. threads: %d", p, num_threads);
    matlib_index num_exp    = serial_time.lenc;
    matlib_index num_cycles = serial_time.lenr-1;

    struct timespec tb, te;
    matlib_real dt;

    clock_gettime(CLOCK_REALTIME, &tb);
    
    /* Create pthreads */
    pthpool_data_t mp[num_threads];
    
    pthpool_create_threads(num_threads, mp);

    clock_gettime(CLOCK_REALTIME, &te);
    dt = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;
    debug_print( "Creation of threads[msec]: %0.6f", dt);

    matlib_index i, j;


    matlib_index N, N0 = 200;
    matlib_index P = 4*p;

    /* define the domain */ 
    matlib_real x_l = -5.0;
    matlib_real x_r =  5.0;

    matlib_xv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_xm FM, IM;

    matlib_create_xm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_xm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    legendre_LGLdataIM( xi, IM);

    matlib_xv x, u, vb, U, V;
    matlib_real dim;
    matlib_real norm_actual, e_relative;

    for(j=0; j<num_exp; j++)
    {
        /* generate the grid */ 
        N = (j+1)*N0;
        serial_time.elem_p[j]   = (matlib_real)N; 
        parallel_time.elem_p[j] = (matlib_real)N; 

        fem1d_ref2mesh (xi, N, x_l, x_r, &x);
        debug_body("length of x: %d", x.len);


        /* Provide the initial condition */
        matlib_create_xv( x.len, &u, MATLIB_COL_VECT);

        Gaussian(x, u);
        dim = N*(p+1);
        matlib_create_xv( dim, &U, MATLIB_COL_VECT);
        matlib_create_xv( dim, &V, MATLIB_COL_VECT);
        matlib_create_xv( N*p+1, &vb, MATLIB_COL_VECT);
        fem1d_XFLT( N, FM, u, U);
        fem1d_XL2F( p, U, vb);

        for(i=0; i<num_cycles; i++)
        {
            clock_gettime(CLOCK_REALTIME, &tb);
            fem1d_XF2L( p, vb, U);
            clock_gettime(CLOCK_REALTIME, &te);
            serial_time.elem_p[j+(i+1)*num_exp] = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                  (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            clock_gettime(CLOCK_REALTIME, &tb);
            pfem1d_XF2L(p, vb, V, num_threads, mp);
            clock_gettime(CLOCK_REALTIME, &te);
            parallel_time.elem_p[j+(i+1)*num_exp] = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                    (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            BEGIN_DTRACE
                for(i=0; i<U.len; i++)
                {
                    debug_print( "[%d] -> serial: %0.16f, parallel: %0.16f", 
                                 i, U.elem_p[i], V.elem_p[i]);
                }
            END_DTRACE
            
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
    dt = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

    debug_print( "Exit threads[msec]: %0.6f", dt);

    /* Analyze error */ 

}


void test_pfem1d_XPrjL2F
(
    matlib_index p,
    matlib_index num_threads,
    matlib_xm    serial_time,
    matlib_xm    parallel_time
)
{
    
    debug_enter("polynomial degree: %d, nr. threads: %d", p, num_threads);
    matlib_index num_exp    = serial_time.lenc;
    matlib_index num_cycles = serial_time.lenr-1;

    struct timespec tb, te;
    matlib_real dt;

    clock_gettime(CLOCK_REALTIME, &tb);
    
    /* Create pthreads */
    pthpool_data_t mp[num_threads];
    
    pthpool_create_threads(num_threads, mp);

    clock_gettime(CLOCK_REALTIME, &te);
    dt = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;
    debug_print( "Creation of threads[msec]: %0.6f", dt);

    matlib_index i, j;


    matlib_index N, N0 = 200;
    matlib_index P = 4*p;

    /* define the domain */ 
    matlib_real x_l = -5.0;
    matlib_real x_r =  5.0;

    matlib_xv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_xm FM, IM;

    matlib_create_xm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_xm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    legendre_LGLdataIM( xi, IM);

    matlib_xv x, u, Pvb1, Pvb2, U;
    matlib_real dim;
    matlib_real norm_actual, e_relative;

    for(j=0; j<num_exp; j++)
    {
        /* generate the grid */ 
        N = (j+1)*N0;
        serial_time.elem_p[j]   = (matlib_real)N; 
        parallel_time.elem_p[j] = (matlib_real)N; 

        fem1d_ref2mesh (xi, N, x_l, x_r, &x);
        debug_body("length of x: %d", x.len);


        /* Provide the initial condition */
        matlib_create_xv( x.len, &u, MATLIB_COL_VECT);

        Gaussian(x, u);
        dim = N*(p+1);
        matlib_create_xv( dim, &U, MATLIB_COL_VECT);
        matlib_create_xv( N*p+1, &Pvb1, MATLIB_COL_VECT);
        matlib_create_xv( N*p+1, &Pvb2, MATLIB_COL_VECT);
        fem1d_XFLT( N, FM, u, U);

        for(i=0; i<num_cycles; i++)
        {
            clock_gettime(CLOCK_REALTIME, &tb);
            fem1d_XPrjL2F( p, U, Pvb1);
            clock_gettime(CLOCK_REALTIME, &te);
            serial_time.elem_p[j+(i+1)*num_exp] = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                  (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            clock_gettime(CLOCK_REALTIME, &tb);
            pfem1d_XPrjL2F(p, U, Pvb2, num_threads, mp);
            clock_gettime(CLOCK_REALTIME, &te);
            parallel_time.elem_p[j+(i+1)*num_exp] = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                    (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            /* Analyze error */ 
            BEGIN_DTRACE
                for(i=0; i<Pvb1.len; i++)
                {
                    debug_print( "[%d] -> serial: %0.16f, parallel: %0.16f", 
                                 i, Pvb1.elem_p[i], Pvb2.elem_p[i]);
                }
            END_DTRACE
            
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
    dt = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

    debug_print( "Exit threads[msec]: %0.6f", dt);


}

void test_pfem1d_XNorm2
(
    matlib_index p,
    matlib_index num_threads,
    matlib_xm    serial_time,
    matlib_xm    parallel_time
)
{
    
    debug_enter("polynomial degree: %d, nr. threads: %d", p, num_threads);
    matlib_index num_exp    = serial_time.lenc;
    matlib_index num_cycles = serial_time.lenr-1;

    struct timespec tb, te;
    matlib_real dt;

    clock_gettime(CLOCK_REALTIME, &tb);
    
    /* Create pthreads */
    pthpool_data_t mp[num_threads];
    
    pthpool_create_threads(num_threads, mp);

    clock_gettime(CLOCK_REALTIME, &te);
    dt = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;
    debug_print( "Creation of threads[msec]: %0.6f", dt);

    matlib_index i, j;


    matlib_index N, N0 = 200;
    matlib_index P = 4*p;

    /* define the domain */ 
    matlib_real x_l = -5.0;
    matlib_real x_r =  5.0;

    matlib_xv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_xm FM, IM;

    matlib_create_xm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_xm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    legendre_LGLdataIM( xi, IM);

    matlib_xv x, u, U;
    matlib_real dim;
    matlib_real norm1, norm2, e_relative;

    for(j=0; j<num_exp; j++)
    {
        /* generate the grid */ 
        N = (j+1)*N0;
        serial_time.elem_p[j]   = (matlib_real)N; 
        parallel_time.elem_p[j] = (matlib_real)N; 

        fem1d_ref2mesh (xi, N, x_l, x_r, &x);
        debug_body("length of x: %d", x.len);


        /* Provide the initial condition */
        matlib_create_xv( x.len, &u, MATLIB_COL_VECT);

        Gaussian(x, u);
        dim = N*(p+1);
        matlib_create_xv( dim, &U, MATLIB_COL_VECT);
        fem1d_XFLT( N, FM, u, U);

        for(i=0; i<num_cycles; i++)
        {
            clock_gettime(CLOCK_REALTIME, &tb);
            norm1 = fem1d_XNorm2( p, N, U);
            clock_gettime(CLOCK_REALTIME, &te);
            serial_time.elem_p[j+(i+1)*num_exp] = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                  (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            clock_gettime(CLOCK_REALTIME, &tb);
            norm2 = pfem1d_XNorm2(p, N, U, num_threads, mp);
            clock_gettime(CLOCK_REALTIME, &te);
            parallel_time.elem_p[j+(i+1)*num_exp] = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                    (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            /* Analyze error */ 
            BEGIN_DTRACE
                e_relative = fabs(norm1-norm2)/norm1;
                debug_body("Relative error: % 0.16g", e_relative);
            END_DTRACE
            
        }
        matlib_free((void*)x.elem_p);
        matlib_free((void*)u.elem_p);
        matlib_free((void*)U.elem_p);
    }
    
    debug_body("%s", "signal threads to exit!");
    clock_gettime(CLOCK_REALTIME, &tb);
    pthpool_destroy_threads(num_threads, mp);
    clock_gettime(CLOCK_REALTIME, &te);
    dt = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

    debug_print( "Exit threads[msec]: %0.6f", dt);


}


/*============================================================================*/
void zGaussian
( 
    matlib_xv x, 
    matlib_zv u
)
{

    debug_enter("%s", "");
    
    matlib_real *xptr;

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
    matlib_xm    serial_time,
    matlib_xm    parallel_time
)
{
    
    debug_enter("polynomial degree: %d, nr. threads: %d", p, num_threads);
    matlib_index num_exp    = serial_time.lenc;
    matlib_index num_cycles = serial_time.lenr-1;

    struct timespec tb, te;
    matlib_real dt;

    clock_gettime(CLOCK_REALTIME, &tb);
    
    /* Create pthreads */
    pthpool_data_t mp[num_threads];
    
    pthpool_create_threads(num_threads, mp);

    clock_gettime(CLOCK_REALTIME, &te);
    dt = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;
    debug_print( "Creation of threads[msec]: %0.6f", dt);

    matlib_index i, j;


    matlib_index N, N0 = 200;
    matlib_index P = 4*p;

    /* define the domain */ 
    matlib_real x_l = -5.0;
    matlib_real x_r =  5.0;

    matlib_xv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_xm FM, IM;

    matlib_create_xm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_xm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    legendre_LGLdataIM( xi, IM);

    matlib_xv x;
    matlib_zv u, U, V;
    matlib_real dim;
    matlib_real norm_actual, e_relative;

    for(j=0; j<num_exp; j++)
    {
        /* generate the grid */ 
        N = (j+1)*N0;
        serial_time.elem_p[j]   = (matlib_real)N; 
        parallel_time.elem_p[j] = (matlib_real)N; 

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
            serial_time.elem_p[j+(i+1)*num_exp] = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                  (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            clock_gettime(CLOCK_REALTIME, &tb);
            pfem1d_ZFLT(N, FM, u, V, num_threads, mp);
            clock_gettime(CLOCK_REALTIME, &te);
            parallel_time.elem_p[j+(i+1)*num_exp] = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                    (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            /* Analyze error */ 
            BEGIN_DTRACE
                for(i=0; i<U.len; i++)
                {
                    debug_print( "[%d] -> serial: %0.16f, parallel: %0.16f", 
                                 i, U.elem_p[i], V.elem_p[i]);
                }
            END_DTRACE
            
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
    dt = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

    debug_print( "Exit threads[msec]: %0.6f", dt);


}

void test_pfem1d_ZF2L
(
    matlib_index p,
    matlib_index num_threads,
    matlib_xm    serial_time,
    matlib_xm    parallel_time
)
{
    
    debug_enter("polynomial degree: %d, nr. threads: %d", p, num_threads);
    matlib_index num_exp    = serial_time.lenc;
    matlib_index num_cycles = serial_time.lenr-1;

    struct timespec tb, te;
    matlib_real dt;

    clock_gettime(CLOCK_REALTIME, &tb);
    
    /* Create pthreads */
    pthpool_data_t mp[num_threads];
    
    pthpool_create_threads(num_threads, mp);

    clock_gettime(CLOCK_REALTIME, &te);
    dt = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;
    debug_print( "Creation of threads[msec]: %0.6f", dt);

    matlib_index i, j;

    matlib_index N, N0 = 200;
    matlib_index P = 4*p;

    /* define the domain */ 
    matlib_real x_l = -5.0;
    matlib_real x_r =  5.0;

    matlib_xv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_xm FM, IM;

    matlib_create_xm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_xm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    legendre_LGLdataIM( xi, IM);

    matlib_xv x;
    matlib_zv u, U, V, vb;
    matlib_real dim;
    matlib_real norm_actual, e_relative;

    for(j=0; j<num_exp; j++)
    {
        /* generate the grid */ 
        N = (j+1)*N0;
        serial_time.elem_p[j]   = (matlib_real)N; 
        parallel_time.elem_p[j] = (matlib_real)N; 

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
            serial_time.elem_p[j+(i+1)*num_exp] = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                  (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            clock_gettime(CLOCK_REALTIME, &tb);
            pfem1d_ZF2L(p, vb, V, num_threads, mp);
            clock_gettime(CLOCK_REALTIME, &te);
            parallel_time.elem_p[j+(i+1)*num_exp] = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                    (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            BEGIN_DTRACE
                for(i=0; i<U.len; i++)
                {
                    debug_print( "[%d] -> serial: %0.16f% 0.16f, parallel: %0.16f% 0.16f", 
                                 i, U.elem_p[i], V.elem_p[i]);
                }
            END_DTRACE
            
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
    dt = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

    debug_print( "Exit threads[msec]: %0.6f", dt);

    /* Analyze error */ 

}

void test_pfem1d_ZPrjL2F
(
    matlib_index p,
    matlib_index num_threads,
    matlib_xm    serial_time,
    matlib_xm    parallel_time
)
{
    debug_enter("polynomial degree: %d, nr. threads: %d", p, num_threads);
    matlib_index num_exp    = serial_time.lenc;
    matlib_index num_cycles = serial_time.lenr-1;

    struct timespec tb, te;
    matlib_real dt;

    clock_gettime(CLOCK_REALTIME, &tb);
    
    /* Create pthreads */
    pthpool_data_t mp[num_threads];
    
    pthpool_create_threads(num_threads, mp);

    clock_gettime(CLOCK_REALTIME, &te);
    dt = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;
    debug_print( "Creation of threads[msec]: %0.6f", dt);

    matlib_index i, j;


    matlib_index N, N0 = 200;
    matlib_index P = 4*p;

    /* define the domain */ 
    matlib_real x_l = -5.0;
    matlib_real x_r =  5.0;

    matlib_xv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_xm FM, IM;

    matlib_create_xm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_xm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    legendre_LGLdataIM( xi, IM);

    matlib_xv x;
    matlib_zv u, U, Pvb1, Pvb2;
    matlib_real dim;
    matlib_real norm_actual, e_relative;

    for(j=0; j<num_exp; j++)
    {
        /* generate the grid */ 
        N = (j+1)*N0;
        serial_time.elem_p[j]   = (matlib_real)N; 
        parallel_time.elem_p[j] = (matlib_real)N; 

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
            serial_time.elem_p[j+(i+1)*num_exp] = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                  (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            clock_gettime(CLOCK_REALTIME, &tb);
            pfem1d_ZPrjL2F(p, U, Pvb2, num_threads, mp);
            clock_gettime(CLOCK_REALTIME, &te);
            parallel_time.elem_p[j+(i+1)*num_exp] = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                    (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            BEGIN_DTRACE
                for(i=0; i<Pvb1.len; i++)
                {
                    debug_print( "[%d] -> serial: %0.16f% 0.16f, parallel: %0.16f% 0.16f", 
                                 i, Pvb1.elem_p[i], Pvb2.elem_p[i]);
                }
            END_DTRACE
            
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
    dt = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

    debug_print( "Exit threads[msec]: %0.6f", dt);
}

void test_pfem1d_ZNorm2
(
    matlib_index p,
    matlib_index num_threads,
    matlib_xm    serial_time,
    matlib_xm    parallel_time
)
{
    debug_enter("polynomial degree: %d, nr. threads: %d", p, num_threads);
    matlib_index num_exp    = serial_time.lenc;
    matlib_index num_cycles = serial_time.lenr-1;

    struct timespec tb, te;
    matlib_real dt;

    clock_gettime(CLOCK_REALTIME, &tb);
    
    /* Create pthreads */
    pthpool_data_t mp[num_threads];
    
    pthpool_create_threads(num_threads, mp);

    clock_gettime(CLOCK_REALTIME, &te);
    dt = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;
    debug_print( "Creation of threads[msec]: %0.6f", dt);

    matlib_index i, j;


    matlib_index N, N0 = 200;
    matlib_index P = 4*p;

    /* define the domain */ 
    matlib_real x_l = -5.0;
    matlib_real x_r =  5.0;

    matlib_xv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_xm FM, IM;

    matlib_create_xm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_xm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( xi, FM);
    legendre_LGLdataIM( xi, IM);

    matlib_xv x;
    matlib_zv u, U;
    matlib_real dim;
    matlib_real norm1, norm2, e_relative;

    for(j=0; j<num_exp; j++)
    {
        /* generate the grid */ 
        N = (j+1)*N0;
        serial_time.elem_p[j]   = (matlib_real)N; 
        parallel_time.elem_p[j] = (matlib_real)N; 

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
            serial_time.elem_p[j+(i+1)*num_exp] = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                  (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            clock_gettime(CLOCK_REALTIME, &tb);
            norm2 = pfem1d_ZNorm2( p, N, U, num_threads, mp);
            clock_gettime(CLOCK_REALTIME, &te);
            parallel_time.elem_p[j+(i+1)*num_exp] = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 +
                                                    (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

            BEGIN_DTRACE
                e_relative = fabs(norm1-norm2)/norm1;
                debug_body("Relative error: % 0.16g", e_relative);
            END_DTRACE
            
        }
        matlib_free((void*)x.elem_p);
        matlib_free((void*)u.elem_p);
        matlib_free((void*)U.elem_p);
    }
    
    debug_body("%s", "signal threads to exit!");
    clock_gettime(CLOCK_REALTIME, &tb);
    pthpool_destroy_threads(num_threads, mp);
    clock_gettime(CLOCK_REALTIME, &te);
    dt = (matlib_real)(te.tv_sec-tb.tv_sec)*1.0e3 + 
         (matlib_real)(te.tv_nsec-tb.tv_nsec)/1.0e6;

    debug_print( "Exit threads[msec]: %0.6f", dt);

    /* Analyze error */ 

}

/*============================================================================*/

void test_performance
(
    matlib_index p,
    void         (*fp)(matlib_index, matlib_index, matlib_xm, matlib_xm),
    char*        path
) 
{
    /* define number of timing experiments and repeatitions */ 
    matlib_index num_exp    = 50;
    matlib_index num_cycles = 1000;
    
    matlib_index num_threads_max = 7;

    char file_name[80];

    matlib_xm serial_time, parallel_time;

    matlib_create_xm( num_exp, 
                      num_cycles+1, 
                      &serial_time, 
                      MATLIB_COL_MAJOR, 
                      MATLIB_NO_TRANS);
    matlib_create_xm( num_exp, 
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

        matlib_xmwrite_csv( file_name, serial_time);

        sprintf( file_name, "%s/parallel_time_num_threads%d_p%02d.dat",
                 path, num_threads, p);
        debug_body("file name: %s", file_name);

        matlib_xmwrite_csv( file_name, parallel_time);
        debug_body("%s","data written to files.");
    }
    
    matlib_free((void*)serial_time.elem_p);
    matlib_free((void*)parallel_time.elem_p);

    debug_exit("%s", "");

}

void check_and_create_dir(char* path)
{
    struct stat sb = {0,}; /* stat buffer */ 
    if(stat(path, &sb)==-1)
    {
        mkdir(path, (S_IRWXU|S_IRWXG));
    }
    else
    {
        matlib_int uin = 1;
        printf("Directory exists.\n");
        printf("Enter [0] to overwrite (default: abort): ");
        scanf("%d", &uin);
        if(uin!=0)
        {
            printf("Aborting...\n");
            exit(EXIT_SUCCESS);
        }
    }
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

    /* Print options and get user input */ 
    char* options[] = { "pfem1d_XFLT",
                        "pfem1d_XF2L",
                        "pfem1d_XPrjL2F",
                        "pfem1d_XNorm2",
                        "pfem1d_ZFLT",
                        "pfem1d_ZF2L",
                        "pfem1d_ZPrjL2F",
                        "pfem1d_ZNorm2", 
                        NULL};

    void* fp[] = { test_pfem1d_XFLT,
                   test_pfem1d_XF2L,
                   test_pfem1d_XPrjL2F,
                   test_pfem1d_XNorm2,
                   test_pfem1d_ZFLT,
                   test_pfem1d_ZF2L,
                   test_pfem1d_ZPrjL2F,
                   test_pfem1d_ZNorm2, 
                   NULL};

    matlib_int i = 0;
    do 
    {
        printf("[%d] Test %s\n", i, options[i]);
        i++;
    } 
    while (options[i]!=NULL);

    matlib_int input, nr_try = 3;
    for(i=0; i< nr_try; i++ )
    {
        printf( "Enter option : ");
        scanf("%d", &input);
        if(options[input] != NULL)
        {
            char* path = options[input];
            printf("Running speedup test: %s\n", path);
            check_and_create_dir(path);

            matlib_index p = 4; /* default value */ 
            printf( "Enter polynomial degree (>1) : ");
            scanf("%d", &p);
            
            /* Carry out the test. */ 
            if((fp[input]!=NULL) && (p>1))
            {
                test_performance(p, fp[input], path);
            }
            break;
        }
        else
        {
            printf("Invalid input. Try again!\n");
        }
    }

    return(0);
}

