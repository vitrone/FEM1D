/*============================================================================+/
 | File: Cunit_fem1d.c 
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

//#define NDEBUG
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
void poly_func
( 
    const matlib_xv x, 
          matlib_xv y
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

void Gaussian_func_t
( 
    const matlib_xv x, 
    const matlib_xv t, 
          matlib_xm y
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

void Gaussian_zfunc_t
( 
    const matlib_xv x, 
    const matlib_xv t, 
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

/*============================================================================*/

matlib_real test_interpolation_general
(
    matlib_index  p,
    matlib_index  N,
    matlib_real domain[2],
    void   (*func_p)(matlib_xv, matlib_xv)
)
{
    debug_enter( "highest degree: %d, nr finite-elements: %d, domain: [% 0.3f, % 0.3f]",
                 p, N, domain[0], domain[1]);
    matlib_index i;

    /* Define computational domain */ 
    matlib_real x_l =  domain[0];
    matlib_real x_r =  domain[1];

    /* Define the number of finite-elements */ 

    /* Provide the degree of polynomial used for interpolation */ 
    matlib_index mesh_len = N*p+1;
    matlib_index nr_coeff = N*(p+1);

    matlib_xm FM;
    matlib_xv x, xi, quadW, u, U;

    legendre_LGLdataLT1( p, TOL, &xi, &quadW);

    matlib_create_xm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    legendre_LGLdataFM( xi, FM);

    matlib_create_xv( mesh_len, &u, MATLIB_COL_VECT);
    matlib_create_xv( nr_coeff, &U, MATLIB_COL_VECT);

    fem1d_ref2mesh (xi, N, x_l, x_r, &x);
    debug_body("length x: %d, length u: %d", x.len, u.len);
    
    (*func_p)(x, u);
    
    fem1d_XFLT( N, FM, u, U);
    DEBUG_PRINT_XV(U, "%s: ", "Transformed vector");
    
    /* Interpolation on arbitrary points: x1 */ 
    matlib_index P  = 2;
    mesh_len = N*P+1;
    matlib_index nr_samples = P+1;

    matlib_xm IMi;
    matlib_xv xii, u_interpol, u_actual, x1;

    matlib_create_xv( nr_samples, &xii,      MATLIB_COL_VECT);
    matlib_create_xv( mesh_len, &u_interpol, MATLIB_COL_VECT);
    matlib_create_xv( mesh_len, &u_actual,   MATLIB_COL_VECT);

    /* generate the grid on the reference interval */ 
    matlib_real dxii = 2.0/P;
    debug_body( "dxii=% 0.16f", dxii);

    *(xii.elem_p) = -1.0;
    for (i=1; i<P; i++)
    {
        *(xii.elem_p+i) = *(xii.elem_p+i-1) + dxii;
    }
    *(xii.elem_p+P) = 1.0;
    
    matlib_create_xm( xii.len, p+1, &IMi, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);
    legendre_LGLdataIM(xii, IMi);

    fem1d_XILT( N, IMi, U, u_interpol);


    /* find the actual values from the known function */ 
    fem1d_ref2mesh (xii, N, x_l, x_r, &x1);
    (*func_p)( x1, u_actual);

    BEGIN_DTRACE
        for (i=0; (i<x1.len); i++)
        {
            debug_print( "x=% 0.3f, interpolated "
                         "y:% 0.16f, actual y= % 0.16f",
                         *(x1.elem_p+i), 
                         *(u_interpol.elem_p+i), 
                         *(u_actual.elem_p+i));
        }
    END_DTRACE
    
    /* Compute the error */
    matlib_real norm_actual = matlib_xnrm2(u_actual);
    matlib_xaxpy(-1.0, u_actual, u_interpol);
    matlib_real e_relative = matlib_xnrm2(u_interpol)/norm_actual;
    
    /* Free allocated memory*/ 
    matlib_free( FM.elem_p);
    matlib_free( xi.elem_p);
    matlib_free( x.elem_p);
    matlib_free( u.elem_p);
    matlib_free( U.elem_p);
    matlib_free( quadW.elem_p);
    matlib_free( IMi.elem_p);
    matlib_free( xii.elem_p);
    matlib_free( x1.elem_p);
    matlib_free( u_actual.elem_p);
    matlib_free( u_interpol.elem_p);

    debug_exit("Relative interpolation error: % 0.16g", e_relative);
    return e_relative;

}
/*============================================================================*/

matlib_real test_interpolation_general2
(
    matlib_index     p,
    matlib_index     N,
    matlib_real    domain[2],
    matlib_xv t,
    void      (*func_p)(matlib_xv, matlib_xv, matlib_xm)
)
{
    debug_enter( "highest degree: %d, nr finite-elements: %d, domain: [% 0.3f, % 0.3f]",
                 p, N, domain[0], domain[1]);
    matlib_index i, j;

    /* Define computational domain */ 
    matlib_real x_l =  domain[0];
    matlib_real x_r =  domain[1];

    /* Define the number of finite-elements */ 

    /* Provide the degree of polynomial used for interpolation */ 
    matlib_index mesh_len = N*p+1;
    matlib_index nr_coeff = N*(p+1);

    matlib_xm FM, u, U;
    matlib_xv x, xi, quadW;

    legendre_LGLdataLT1( p, TOL, &xi, &quadW);

    matlib_create_xm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    legendre_LGLdataFM( xi, FM);

    matlib_create_xm( mesh_len, t.len, &u, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);
    matlib_create_xm( nr_coeff, t.len, &U, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);

    fem1d_ref2mesh (xi, N, x_l, x_r, &x);
    debug_body("length x: %d, size u: %d-by-%d", x.len, u.lenc, u.lenr);
    
    (*func_p)(x, t, u);
    
    fem1d_XFLT2( N, FM, u, U);
    DEBUG_PRINT_XM(U, "%s: ", "Transformed matrix");
    
    /* Interpolation on arbitrary points: x1 */ 
    matlib_index P  = 2;
    mesh_len = N*P+1;
    matlib_index nr_samples = P+1;

    matlib_xm IMi, u_interpol, u_actual;
    matlib_xv xii,  x1;

    matlib_create_xv( nr_samples, &xii,      MATLIB_COL_VECT);
    matlib_create_xm( mesh_len, t.len, &u_interpol, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);
    matlib_create_xm( mesh_len, t.len, &u_actual,   MATLIB_COL_MAJOR, MATLIB_NO_TRANS);

    /* generate the grid on the reference interval */ 
    matlib_real dxii = 2.0/P;
    debug_body( "dxii=% 0.16f", dxii);

    *(xii.elem_p) = -1.0;
    for (i=1; i<P; i++)
    {
        *(xii.elem_p+i) = *(xii.elem_p+i-1) + dxii;
    }
    *(xii.elem_p+P) = 1.0;
    
    matlib_create_xm( xii.len, p+1, &IMi, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);
    legendre_LGLdataIM(xii, IMi);

    fem1d_XILT2( N, IMi, U, u_interpol);


    /* find the actual values from the known function */ 
    fem1d_ref2mesh (xii, N, x_l, x_r, &x1);
    (*func_p)( x1, t, u_actual);
    
    /* Compute the error */
    matlib_real norm_actual, e_relative;

    matlib_index tmp_len = mesh_len * t.len; 

   /* 
    * matlib_xv u_tmp1 = MK_VM(u_interpol);
    * matlib_xv u_tmp2 = MK_VM(u_actual);
    * */

    matlib_xv u_tmp1 = { .len = tmp_len, 
                         .type = MATLIB_COL_VECT, 
                         .elem_p = u_interpol.elem_p};
    matlib_xv u_tmp2 = { .len = tmp_len, 
                         .type = MATLIB_COL_VECT, 
                         .elem_p = u_actual.elem_p};

    BEGIN_DTRACE
        for (j=0; (j<t.len); j++)
        {
            for (i=0; (i<x1.len); i++)
            {
                debug_print( "t=% 0.3f, x=% 0.3f, "
                             "interpolated y:% 0.16f, "
                             "actual y= % 0.16f",
                             *(t.elem_p+j),
                             *(x1.elem_p+i), 
                             *(u_interpol.elem_p+i+j*x1.len), 
                             *(u_actual.elem_p+i+j*x1.len));
            }
        }
    END_DTRACE

    norm_actual = matlib_xnrm2(u_tmp1);
    matlib_xaxpy(-1.0, u_tmp2, u_tmp1);
    e_relative = matlib_xnrm2(u_tmp1)/norm_actual;
    
    /* Free allocated memory*/ 
    matlib_free( FM.elem_p);
    matlib_free( xi.elem_p);
    matlib_free( x.elem_p);
    matlib_free( u.elem_p);
    matlib_free( U.elem_p);
    matlib_free( quadW.elem_p);
    matlib_free( IMi.elem_p);
    matlib_free( xii.elem_p);
    matlib_free( x1.elem_p);
    matlib_free( u_actual.elem_p);
    matlib_free( u_interpol.elem_p);

    debug_exit("Relative interpolation error: % 0.16g", e_relative);
    return e_relative;

}
/*============================================================================*/

void test_interpolation1(void)
{

    matlib_index p;
    matlib_index N = 500;
    matlib_real domain[2] = {-1.0, 1.0};
    matlib_real e_relative;

    for(p=2; p<12; p++)
    {
        e_relative = test_interpolation_general( p, N,
                                                 domain,
                                                 poly_func);

        CU_ASSERT_TRUE(e_relative<TOL);
    }

}
void test_interpolation2(void)
{

    matlib_index p;
    matlib_index N = 500;
    matlib_real domain[2] = {-5.0, 5.0};
    matlib_real e_relative;

    for(p=2; p<12; p++)
    {
        e_relative = test_interpolation_general( p, N,
                                                 domain,
                                                 Gaussian_func);

        CU_ASSERT_TRUE(e_relative<TOL);
    }
}

void test_interpolation3(void)
{
    debug_enter("%s", "");

    matlib_index p;
    matlib_index N = 400;
    matlib_real domain[2] = {-5.0, 5.0};
    matlib_real ta[4] = {0, 0.5, 1.0, 1.5};
    matlib_xv t = {.len = 4, .type = MATLIB_COL_VECT, .elem_p = ta};

    matlib_real e_relative;
    for(p=2; p<12; p++)
    {
        e_relative = test_interpolation_general2( p, N,
                                                 domain,
                                                 t,
                                                 Gaussian_func_t);

        CU_ASSERT_TRUE(e_relative<TOL);
    }

}
/*============================================================================*/

matlib_real test_interpolation_general3
(
    matlib_index p,
    matlib_index N
)
{
    matlib_index i;
    
    /* Define computational domain */ 
    matlib_real x_l = -5.0;
    matlib_real x_r = 5.0;

    /* Define the number of finite-elements */ 

    /* Provide the degree of polynomial used for interpolation */ 
    matlib_index mesh_len = N*p+1;
    matlib_index nr_coeff = N*(p+1);

    matlib_xm FM;
    matlib_xv x, xi, quadW;
    matlib_zv u, U;

    legendre_LGLdataLT1( p, TOL, &xi, &quadW);

    matlib_create_xm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    legendre_LGLdataFM( xi, FM);

    matlib_create_zv( mesh_len, &u, MATLIB_COL_VECT);
    matlib_create_zv( nr_coeff, &U, MATLIB_COL_VECT);

    fem1d_ref2mesh (xi, N, x_l, x_r, &x);
    debug_body("length x: %d, length u: %d", x.len, u.len);
    
    Gaussian_zfunc(x, u);
    
    fem1d_ZFLT( N, FM, u, U);
    DEBUG_PRINT_ZV(U, "%s: ", "Transformed complex-vector");

    /* Interpolation on arbitrary points: x1 */ 
    matlib_index P  = 2;
    mesh_len = N*P+1;
    matlib_index nr_samples = P+1;

    matlib_xm IMi;
    matlib_xv xii, x1;
    matlib_zv u_interpol, u_actual;

    matlib_create_xv( nr_samples, &xii,      MATLIB_COL_VECT);
    matlib_create_zv( mesh_len, &u_interpol, MATLIB_COL_VECT);
    matlib_create_zv( mesh_len, &u_actual,   MATLIB_COL_VECT);

    /* generate the grid on the reference interval */ 
    matlib_real dxii = 2.0/P;
    debug_body( "dxii=% 0.16f", dxii);

    *(xii.elem_p) = -1.0;
    for (i=1; i<P; i++)
    {
        *(xii.elem_p+i) = *(xii.elem_p+i-1) + dxii;
    }
    *(xii.elem_p+P) = 1.0;
    
    matlib_create_xm( xii.len, p+1, &IMi, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);
    legendre_LGLdataIM(xii, IMi);

    fem1d_ZILT( N, IMi, U, u_interpol);

    /* find the actual values from the known function */ 
    fem1d_ref2mesh (xii, N, x_l, x_r, &x1);
    Gaussian_zfunc( x1, u_actual);

    BEGIN_DTRACE
        for (i=0; (i<x1.len); i++)
        {
            debug_print( "x=% 0.3f, "
                         "interpolated y: % 0.16f %+0.16fi, "
                         "actual y = % 0.16f %+0.16fi",
                         *(x1.elem_p+i), 
                         *(u_interpol.elem_p+i), 
                         *(u_actual.elem_p+i));
        }
    END_DTRACE
    
    /* Compute the error */
    matlib_real norm_actual = matlib_znrm2(u_actual);
    matlib_zaxpy(-1.0, u_actual, u_interpol);
    matlib_real e_relative = matlib_znrm2(u_interpol)/norm_actual;
    
    /* Free allocated memory*/ 
    matlib_free( FM.elem_p);
    matlib_free( xi.elem_p);
    matlib_free( x.elem_p);
    matlib_free( u.elem_p);
    matlib_free( U.elem_p);
    matlib_free( quadW.elem_p);
    matlib_free( IMi.elem_p);
    matlib_free( xii.elem_p);
    matlib_free( x1.elem_p);
    matlib_free( u_actual.elem_p);
    matlib_free( u_interpol.elem_p);

    debug_exit("Relative interpolation error: % 0.16g", e_relative);
    return(e_relative);
}
/*============================================================================*/

matlib_real test_interpolation_general4
(
    matlib_index p,
    matlib_index N
)
{
    matlib_index i, j;
    
    /* Define computational domain */ 
    matlib_real x_l = -5.0;
    matlib_real x_r = 5.0;

    matlib_real ta[4] = {0, 0.5, 1.0, 1.5};
    matlib_xv t = { .len = 4, 
                    .type = MATLIB_COL_VECT, 
                    .elem_p = ta};

    /* Define the number of finite-elements */ 

    /* Provide the degree of polynomial used for interpolation */ 
    matlib_index mesh_len = N*p+1;
    matlib_index nr_coeff = N*(p+1);

    matlib_xm FM;
    matlib_xv x, xi, quadW;
    matlib_zm u, U;

    legendre_LGLdataLT1( p, TOL, &xi, &quadW);

    matlib_create_xm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    legendre_LGLdataFM( xi, FM);

    matlib_create_zm( mesh_len, t.len, &u, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);
    matlib_create_zm( nr_coeff, t.len, &U, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);

    fem1d_ref2mesh (xi, N, x_l, x_r, &x);
    debug_body("length x: %d, size u: %d-by-%d", x.len, u.lenc, u.lenr);
    
    Gaussian_zfunc_t(x, t, u);
    
    fem1d_ZFLT2( N, FM, u, U);
    DEBUG_PRINT_ZM(U, "%s: ", "Transformed complex-matrix");

    /* Interpolation on arbitrary points: x1 */ 
    matlib_index P  = 2;
    mesh_len = N*P+1;
    matlib_index nr_samples = P+1;

    matlib_xm IMi;
    matlib_xv xii, x1;
    matlib_zm u_interpol, u_actual;

    matlib_create_xv( nr_samples, &xii, MATLIB_COL_VECT);
    matlib_create_zm( mesh_len, t.len, &u_interpol, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);
    matlib_create_zm( mesh_len, t.len, &u_actual,   MATLIB_COL_MAJOR, MATLIB_NO_TRANS);

    /* generate the grid on the reference interval */ 
    matlib_real dxii = 2.0/P;
    debug_body( "dxii=% 0.16f", dxii);

    *(xii.elem_p) = -1.0;
    for (i=1; i<P; i++)
    {
        *(xii.elem_p+i) = *(xii.elem_p+i-1) + dxii;
    }
    *(xii.elem_p+P) = 1.0;
    
    matlib_create_xm( xii.len, p+1, &IMi, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);
    legendre_LGLdataIM(xii, IMi);

    fem1d_ZILT2( N, IMi, U, u_interpol);

    /* find the actual values from the known function */ 
    fem1d_ref2mesh (xii, N, x_l, x_r, &x1);
    Gaussian_zfunc_t( x1, t, u_actual);

    BEGIN_DTRACE
        for (j=0; (j<t.len); j++)
        {
            for (i=0; (i<x1.len); i++)
            {
                debug_print( "t=% 0.3f, x=% 0.3f, "
                             "interpolated y: % 0.16f %+0.16fi, "
                             "actual y = % 0.16f %+0.16fi",
                             *(t.elem_p+j),
                             *(x1.elem_p+i), 
                             *(u_interpol.elem_p+i+j*x1.len), 
                             *(u_actual.elem_p+i+j*x1.len));
            }
        }
    END_DTRACE
    
    /* Compute the error */
    matlib_zv u_tmp1 = MK_VM(u_interpol);
    matlib_zv u_tmp2 = MK_VM(u_actual);

    matlib_real norm_actual = matlib_znrm2(u_tmp1);
    matlib_zaxpy(-1.0, u_tmp2, u_tmp1);
    matlib_real e_relative = matlib_znrm2(u_tmp1)/norm_actual;
    
    /* Free allocated memory*/ 
    matlib_free(FM.elem_p);
    matlib_free(xi.elem_p);
    matlib_free(x.elem_p);
    matlib_free(u.elem_p);
    matlib_free(U.elem_p);
    matlib_free(quadW.elem_p);
    matlib_free(IMi.elem_p);
    matlib_free(xii.elem_p);
    matlib_free(x1.elem_p);
    matlib_free(u_actual.elem_p);
    matlib_free(u_interpol.elem_p);

    debug_exit("Relative interpolation error: % 0.16g", e_relative);
    return(e_relative);
}
/*============================================================================*/

void test_interpolation4(void)
{

    matlib_index p;
    matlib_index N = 800;
    matlib_real e_relative;

    for(p=2; p<12; p++)
    {
        e_relative = test_interpolation_general3( p, N);

        CU_ASSERT_TRUE(e_relative<TOL);
    }
}

void test_interpolation5(void)
{
    debug_enter("%s", "");

    matlib_index p;
    matlib_index N = 800;

    matlib_real e_relative;
    for(p=2; p<12; p++)
    {
        e_relative = test_interpolation_general4( p, N);

        CU_ASSERT_TRUE(e_relative<TOL);
    }

}


/*============================================================================*/

matlib_real test_L2_dnorm_general
(
    matlib_index p,
    matlib_index N
)
{
    debug_enter("highest degree: %d, nr. fem-elements: %d", p, N);
    matlib_real x_l = -5.0;
    matlib_real x_r =  5.0;

    matlib_real J = (x_r-x_l)/(2.0*N);
    
    matlib_index mesh_len = N*p+1;
    matlib_index nr_coeff = N*(p+1);

    matlib_xv x, xi, quadW;
    matlib_xv u, U;
    matlib_xm FM;

    legendre_LGLdataLT1( p, TOL, &xi, &quadW);

    matlib_create_xm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    legendre_LGLdataFM( xi, FM);

    matlib_create_xv( mesh_len, &u, MATLIB_COL_VECT);
    matlib_create_xv( nr_coeff, &U, MATLIB_COL_VECT);

    fem1d_ref2mesh (xi, N, x_l, x_r, &x);
    debug_body("length x: %d, length u: %d", x.len, u.len);
    
    Gaussian_func(x, u);
    
    fem1d_XFLT( N, FM, u, U);
    DEBUG_PRINT_XV(U, "%s: ", "Transformed vector");

    matlib_real expected_norm = pow(M_PI, 0.25);
    matlib_real lp_norm = sqrt(J)*fem1d_XNorm2(p, N, U);
    matlib_real e_relative = fabs(expected_norm-lp_norm)/expected_norm;

    /* Free allocated memory*/ 
    matlib_free(FM.elem_p);
    matlib_free(xi.elem_p);
    matlib_free(x.elem_p);
    matlib_free(u.elem_p);
    matlib_free(U.elem_p);
    matlib_free(quadW.elem_p);

    debug_exit("Relative error in norm2: % 0.16g", e_relative);
    return e_relative;

}

matlib_real test_L2_znorm_general
(
    matlib_index p,
    matlib_index N
)
{
    debug_enter("highest degree: %d, nr. fem-elements: %d", p, N);
    matlib_real x_l = -5.0;
    matlib_real x_r =  5.0;

    matlib_real J = (x_r-x_l)/(2.0*N);
    
    matlib_index mesh_len = N*p+1;
    matlib_index nr_coeff = N*(p+1);

    matlib_xv x, xi, quadW;
    matlib_zv u, U;
    matlib_xm FM;

    legendre_LGLdataLT1( p, TOL, &xi, &quadW);

    matlib_create_xm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    legendre_LGLdataFM( xi, FM);

    matlib_create_zv( mesh_len, &u, MATLIB_COL_VECT);
    matlib_create_zv( nr_coeff, &U, MATLIB_COL_VECT);

    fem1d_ref2mesh (xi, N, x_l, x_r, &x);
    debug_body("length x: %d, length u: %d", x.len, u.len);
    
    Gaussian_zfunc(x, u);
    
    fem1d_ZFLT( N, FM, u, U);
    DEBUG_PRINT_ZV(U, "%s: ", "Transformed vector");

    matlib_real expected_norm = pow(M_PI, 0.25);
    matlib_real lp_norm = sqrt(J)*fem1d_ZNorm2(p, N, U);
    matlib_real e_relative = fabs(expected_norm-lp_norm)/expected_norm;

    /* Free allocated memory*/ 
    matlib_free(FM.elem_p);
    matlib_free(xi.elem_p);
    matlib_free(x.elem_p);
    matlib_free(u.elem_p);
    matlib_free(U.elem_p);
    matlib_free(quadW.elem_p);

    debug_exit("Relative error in norm2: % 0.16g", e_relative);
    return e_relative;

}
/*============================================================================*/

void test_L2_dnorm(void)
{
    matlib_index N = 300;
    matlib_index p_max = 15;
    matlib_real e_relative;
    for(matlib_index p=2; p<p_max; p++)
    {
        e_relative = test_L2_dnorm_general(p, N);
        CU_ASSERT_TRUE(e_relative<TOL);
    }
}
void test_L2_znorm(void)
{
    matlib_index N = 600;
    matlib_index p_max = 15;
    matlib_real e_relative;
    for(matlib_index p=2; p<p_max; p++)
    {
        e_relative = test_L2_znorm_general(p, N);
        CU_ASSERT_TRUE(e_relative<TOL);
    }
}
/*============================================================================*/
matlib_real test_fem1d_XL2F_general
(
    matlib_index p,
    matlib_index N
)
{
    matlib_real x_l =  -5.0;
    matlib_real x_r =   5.0;

    debug_enter( "polynomial degree: %d, nr. of LGL points: %d", p, p+1 );

    matlib_xv x, xi, quadW;
    legendre_LGLdataLT1( p, TOL, &xi, &quadW);
    fem1d_ref2mesh (xi, N, x_l, x_r, &x);
    
    matlib_xm FM;
    matlib_create_xm( xi.len, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    legendre_LGLdataFM( xi, FM);
    
    
    matlib_xv u, U, U1;
    matlib_create_xv( x.len,    &u, MATLIB_COL_VECT);
    matlib_create_xv( N*(p+1),  &U, MATLIB_COL_VECT);
    matlib_create_xv( N*(p+1), &U1, MATLIB_COL_VECT);

    Gaussian_func(x, u);
    DEBUG_PRINT_XV(u, "%s: ", "value of Gaussian");
    fem1d_XFLT( N, FM, u, U); /* Legendre transform */ 

    matlib_xv vb;
    matlib_create_xv( N*p+1, &vb, MATLIB_COL_VECT);
    
    /* Transformation from Lgendre basis to FEM-shape func basis */ 
    fem1d_XL2F(p, U, vb);  
    DEBUG_PRINT_XV(vb, "%s: ", "coeff. fem basis func.");

    fem1d_XF2L(p, vb, U1);  
    BEGIN_DTRACE
        matlib_index i;
        for(i=0;i<U.len; i++)
        {
            debug_print( "[%d]->U: % 0.16f, U1: % 0.16f", 
                         i, *(U.elem_p+i), *(U1.elem_p+i));
        }
    END_DTRACE
    matlib_real norm_actual = matlib_xnrm2(U);
    matlib_xaxpy(-1.0, U, U1);
    matlib_real e_relative = matlib_xnrm2(U1)/norm_actual;

    return(e_relative);
    debug_exit("Relative interpolation error: % 0.16g", e_relative);
}

matlib_real test_fem1d_ZL2F_general
(
    matlib_index p,
    matlib_index N
)
{
    matlib_real x_l =  -5.0;
    matlib_real x_r =   5.0;

    debug_enter( "polynomial degree: %d, nr. of LGL points: %d", p, p+1 );

    matlib_xv x, xi, quadW;
    legendre_LGLdataLT1( p, TOL, &xi, &quadW);
    fem1d_ref2mesh (xi, N, x_l, x_r, &x);
    
    matlib_xm FM;
    matlib_create_xm( xi.len, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    legendre_LGLdataFM( xi, FM);
    
    
    matlib_zv u, U, U1;
    matlib_create_zv( x.len,    &u, MATLIB_COL_VECT);
    matlib_create_zv( N*(p+1),  &U, MATLIB_COL_VECT);
    matlib_create_zv( N*(p+1), &U1, MATLIB_COL_VECT);

    Gaussian_zfunc(x, u);
    DEBUG_PRINT_ZV(u, "%s: ", "value of Gaussian");
    fem1d_ZFLT( N, FM, u, U); /* Legendre transform */ 

    matlib_zv vb;
    matlib_create_zv( N*p+1, &vb, MATLIB_COL_VECT);
    
    /* Transformation from Lgendre basis to FEM-shape func basis */ 
    fem1d_ZL2F(p, U, vb);  
    DEBUG_PRINT_ZV(vb, "%s: ", "coeff. fem basis func.");

    fem1d_ZF2L(p, vb, U1);  
    BEGIN_DTRACE
        matlib_index i;
        for(i=0;i<U.len; i++)
        {
            debug_print( "[%d]->U: % 0.16f %+0.16f, U1: % 0.16f %+0.16f", 
                         i, *(U.elem_p+i), *(U1.elem_p+i));
        }
    END_DTRACE
    matlib_real norm_actual = matlib_znrm2(U);
    matlib_zaxpy(-1.0, U, U1);
    matlib_real e_relative = matlib_znrm2(U1)/norm_actual;

    return(e_relative);
    debug_exit("Relative interpolation error: % 0.16g", e_relative);
}
/*============================================================================*/

void test_fem1d_XL2F1(void)
{
    matlib_index N = 40;
    matlib_index p_max = 15;
    matlib_real e_relative;
    for(matlib_index p=2; p<p_max; p++)
    {
        e_relative = test_fem1d_XL2F_general(p, N);
        CU_ASSERT_TRUE(e_relative<TOL);
    }

}
void test_fem1d_ZL2F1(void)
{
    matlib_index N = 40;
    matlib_index p_max = 15;
    matlib_real e_relative;
    for(matlib_index p=2; p<p_max; p++)
    {
        e_relative = test_fem1d_ZL2F_general(p, N);
        CU_ASSERT_TRUE(e_relative<TOL);
    }

}

/*============================================================================*/

void test_fem1d_quadM1(void)
{
    matlib_index p = 11;
    matlib_index P = 2*p;

    debug_enter( "polynomial degree: %d, nr. of LGL points: %d", p, P+1 );

    matlib_xv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);

    matlib_xm IMi, Q, IMii, Q1;

    matlib_create_xm( xi.len, p+1,  &IMi, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_xm( xi.len, p+1, &IMii, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataIM(xi, IMi);
    legendre_LGLdataIM(xi, IMii);

    DEBUG_PRINT_XM(IMi,  "%s: ", "inverse transform matrix (IM: col major)");
    DEBUG_PRINT_XM(IMii, "%s: ", "inverse transform matrix (IM: row major)");
    
    fem1d_quadM(quadW,  IMi,  &Q);
    fem1d_quadM(quadW, IMii, &Q1);
    
    DEBUG_PRINT_XM(Q,  "%s: ", "Quadrature matrix (IM: col major)");
    DEBUG_PRINT_XM(Q1, "%s: ", "Quadrature matrix (IM: row major)");

    matlib_xv u_tmp1 = MK_VM(Q);
    matlib_xv u_tmp2 = MK_VM(Q1);

    matlib_real norm_actual = matlib_xnrm2(u_tmp1);
    
    matlib_xaxpy(-1.0, u_tmp2, u_tmp1);
    DEBUG_PRINT_XV(u_tmp1, "%s: ", "Q-Q1");
    
    matlib_real e_relative = matlib_xnrm2(u_tmp1)/norm_actual;
    debug_exit("Relative error: % 0.16g", e_relative);


    CU_ASSERT_TRUE(e_relative<TOL);
    debug_exit("%s","");
}

/*============================================================================*/

matlib_real test_fem1d_MEMI_general( matlib_index p)
{
    debug_enter("polynomial degree: %d", p);
    matlib_index i;
    matlib_index P = p+1;

    matlib_xv q;
    fem1d_MEMI(p, &q);
    DEBUG_PRINT_XV(q, "%s :", "Master element mass integrals");

    matlib_xv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_xm FM, IM, Q;
    matlib_create_xm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataIM( xi, IM);
    fem1d_quadM( quadW, IM, &Q);
    DEBUG_PRINT_XM(Q,  "%s: ", "Quadrature matrix (IM: col major)");

    matlib_xv q1, ones;
    matlib_create_xv( Q.lenc,   &q1, MATLIB_COL_VECT);
    matlib_create_xv( xi.len, &ones, MATLIB_COL_VECT);
    
    for(i=0; i<xi.len; i++)
    {
        ones.elem_p[i] = 1.0;
    }

    matlib_xgemv( 1.0, Q, ones, 0, q1);

    BEGIN_DTRACE
        for(i=0;i<q.len; i++)
        {
            debug_print( "[%d]-> q: % 0.16f, q1: % 0.16f", 
                         i, *(q.elem_p+i), *(q1.elem_p+i));
        }
    END_DTRACE

    matlib_real norm_actual = matlib_xnrm2(q);
    
    matlib_xaxpy(-1.0, q, q1);
    matlib_real e_relative = matlib_xnrm2(q1)/norm_actual;
    debug_exit("Relative error: % 0.16g", e_relative);

    return(e_relative);
}
void test_fem1d_MEMI(void)
{
    matlib_index p[11] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    matlib_real e_relative;
    for(matlib_index i=0; i<11; i++)
    {
        e_relative = test_fem1d_MEMI_general(p[i]);
        CU_ASSERT_TRUE(e_relative<TOL);
    }
}
/*============================================================================*/
void full_sym_mat
(
    matlib_xm_sparse sM,
    matlib_xm*      M
)
{
    matlib_index i, j;
    
    matlib_create_xm( sM.lenc, sM.lenr, M, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    

    for(i=0; i<sM.lenc; i++)
    {
        j = sM.rowIn[i];
        M->elem_p[i*sM.lenr+sM.colIn[j]] = sM.elem_p[j];
        for(j=sM.rowIn[i]+1; j<sM.rowIn[i+1]; j++)
        {
            M->elem_p[i*sM.lenr+sM.colIn[j]] = sM.elem_p[j];
            M->elem_p[sM.colIn[j]*sM.lenr+i] = sM.elem_p[j];
        }
    }

}
matlib_real test_fem1d_XGMM_general
(
    matlib_index p,
    matlib_index nr_LGL,
    matlib_index N,
    matlib_real domain[2],
    void  (*func_p)(matlib_xv, matlib_xv), 
    void  (*potential_p)(matlib_xv, matlib_xv)
)
/* 
 * potential function: phi(x)
 * Compute mass matrix for phi(x)u(x)
 * */ 
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

    matlib_index dim = N*p+1, nnz = N*p*(p+3)/2+1;
    matlib_index row[dim+1], col[nnz];
    matlib_real ugpmm[nnz];

    /* Transform it to FEM-basis: vertex and bubble basis functions */  
    matlib_xv Uvb;
    matlib_create_xv( dim, &Uvb, MATLIB_COL_VECT);

    fem1d_XL2F(p, U, Uvb);


    DEBUG_PRINT_XV(U,  "%s: ", "Legendre Transform");
    DEBUG_PRINT_XV(Uvb, "%s: ", "coeff. fem basis func.");

    /* Aseemble the global mass matrix */ 
    matlib_xv Y, y;
    matlib_create_xv( x.len,   &y, MATLIB_COL_VECT);
    matlib_create_xv( N*(p+1), &Y, MATLIB_COL_VECT);
   
    (*potential_p)(x,y);

    debug_body("nnz: %d", nnz);

    matlib_xv q;
    matlib_create_xv( Q.lenc*N, &q, MATLIB_COL_VECT);
    fem1d_XFLT( N, Q, y, q);
    DEBUG_PRINT_XV(q, "%s: ", "inner product");

    fem1d_XCSRGMM(p, N, q, row, col, ugpmm);

    debug_body("nnz: %d, nnz actual: %d", row[dim], nnz);

    matlib_xm_sparse M = { .lenc = dim, 
                          .lenr = dim, 
                          .rowIn = row, 
                          .colIn = col, 
                          .elem_p = ugpmm};

    BEGIN_DTRACE
        debug_print("dimension of the sparse square matrix: %d", dim);
        for(i=0; i<dim; i++)
        {
            for(j=row[i]; j<row[i+1]; j++)
            {
                debug_print("UGPMM(%d,%d): % 0.16f",i, col[j], ugpmm[j]);
            }
        }
    END_DTRACE
    
    /* Compute projection onto FEM-basis using the GMM */ 
    matlib_xv Pvb;
    matlib_create_xv( dim, &Pvb, MATLIB_COL_VECT);

    matlib_xcsrsymv(MATLIB_UPPER, M, Uvb, Pvb);
       
    /* compute projection from the LP representation of phi(x)*u(x)*/
    matlib_xv Pvb1;
    matlib_create_xv( dim, &Pvb1, MATLIB_COL_VECT);

    for(i=0; i<x.len; i++)
    {
        y.elem_p[i] = u.elem_p[i] * y.elem_p[i];
    
    }
    fem1d_XFLT( N, FM, y, Y);
    fem1d_XPrjL2F(p, Y, Pvb1);

    BEGIN_DTRACE
        for(i=0;i<dim; i++)
        {
            debug_print( "[%d]-> Pvb: % 0.16f, Pvb1: % 0.16f", 
                         i, *(Pvb.elem_p+i), *(Pvb1.elem_p+i));
        }
    END_DTRACE

    matlib_real norm_actual = matlib_xnrm2(Pvb1);
    
    matlib_xaxpy(-1.0, Pvb1, Pvb);

    matlib_real e_relative = matlib_xnrm2(Pvb)/norm_actual;

    debug_exit("Relative error: % 0.16g", e_relative);

    return(e_relative);
}
void constant_dpotential( const matlib_xv x, matlib_xv y)
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
void harmonic_dpotential( const matlib_xv x, matlib_xv y)
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
void sinusoidal_dpotential( const matlib_xv x, matlib_xv y)
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
void test_fem1d_XGMM1(void)
{
    matlib_index p, nr_LGL;
    matlib_index N = 600;
    matlib_real domain[2] = {-5.0, 5.0};
    matlib_real e_relative;
    for(matlib_index p=2; p<13; p++)
    {
        nr_LGL = 2*p+1;
        e_relative = test_fem1d_XGMM_general( p, nr_LGL, N, domain,
                                              Gaussian_func, 
                                              constant_dpotential);
        CU_ASSERT_TRUE(e_relative<TOL);

        e_relative = test_fem1d_XGMM_general( p, nr_LGL, N, domain,
                                              Gaussian_func, 
                                              harmonic_dpotential);
        CU_ASSERT_TRUE(e_relative<TOL);

        e_relative = test_fem1d_XGMM_general( p, nr_LGL, N, domain,
                                              Gaussian_func, 
                                              sinusoidal_dpotential);
        CU_ASSERT_TRUE(e_relative<TOL);
    } 
}
/*============================================================================*/
matlib_real test_fem1d_ZGMM_general
(
    matlib_index p,
    matlib_index nr_LGL,
    matlib_index N,
    matlib_real domain[2],
    void  (*func_p)(matlib_xv, matlib_zv), 
    void  (*potential_p)(matlib_xv, matlib_zv)
)
/* 
 * potential function: phi(x)
 * Compute mass matrix for phi(x)u(x)
 * */ 
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

    matlib_zm_sparse M = { .lenc = dim, 
                          .lenr = dim, 
                          .rowIn = row, 
                          .colIn = col, 
                          .elem_p = ugpmm};

    BEGIN_DTRACE
        debug_print("dimension of the sparse square matrix: %d", dim);
        for(i=0; i<dim; i++)
        {
            for(j=row[i]; j<row[i+1]; j++)
            {
                debug_print("UGPMM(%d,%d): % 0.16f %+0.16fi", i, col[j], ugpmm[j]);
            }
        }
    END_DTRACE
    
    /* Compute projection onto FEM-basis using the GMM */ 
    matlib_zv Pvb;
    matlib_create_zv( dim, &Pvb, MATLIB_COL_VECT);

    matlib_zcsrsymv(MATLIB_UPPER, M, Uvb, Pvb);
       
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

    BEGIN_DTRACE
        for(i=0;i<dim; i++)
        {
            debug_print( "[%d]-> Pvb-Pvb1: % 0.16f%+0.16fi", 
                         i, *(Pvb.elem_p+i)-*(Pvb1.elem_p+i));
        }
    END_DTRACE

    matlib_real norm_actual = matlib_znrm2(Pvb1);
    
    matlib_zaxpy(-1.0, Pvb1, Pvb);

    matlib_real e_relative = matlib_znrm2(Pvb)/norm_actual;

    debug_exit("Relative error: % 0.16g", e_relative);

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

void test_fem1d_ZGMM1(void)
{
    matlib_index p, nr_LGL;
    matlib_index N = 1600;
    matlib_real domain[2] = {-5.0, 5.0};
    matlib_real e_relative;

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
        nr_LGL = 2*p+1;
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
/*============================================================================*/

void test_fem1d_zm_nsparse_GMM_general
(
    matlib_index p,
    matlib_index nr_LGL,
    matlib_index N,
    matlib_real  domain[2],
    void  (*potential_p)(matlib_xv, matlib_real, matlib_real, matlib_zm)
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

    matlib_real t = 0, dt = 1.0e-3;
    matlib_zm phi, q;
    matlib_zm_nsparse M;
    matlib_zm_sparse M1;
    matlib_index nsparse = 10;
    fem1d_zm_nsparse_GMM(p, N, nsparse, Q, &phi, &q, &M, FEM1D_GMM_INIT);
    fem1d_zm_nsparse_GMM(p, N, nsparse, Q, NULL, NULL, &M, FEM1D_GET_SPARSITY_ONLY);

    matlib_zv v1 = {.len = M.rowIn[M.lenc]};
    matlib_zv v2 = {.len = M.rowIn[M.lenc]};
    matlib_zv phi1 = {.len = x.len};

    matlib_real norm_actual, e_relative;
    for(i=0; i<5; i++)
    {
        (*potential_p)(x, t, dt, phi);
        fem1d_zm_nsparse_GMM(p, N, nsparse, Q, &phi, &q, &M, FEM1D_GET_NZE_ONLY);
        for(j=0; j<nsparse; j++)
        {
            phi1.elem_p = phi.elem_p + j*x.len;
            fem1d_zm_sparse_GMM(p, Q, phi1, &M1);

            v1.elem_p = *(M.elem_p + j);
            v2.elem_p = M1.elem_p;
            BEGIN_DTRACE
                for(i=0; i<v1.len; i++)
                {
                    debug_print( "[%d]-> v1: % 0.16f%+0.16f, v2: % 0.16f%+0.16f",
                                 i, v1.elem_p[i], v2.elem_p[i]);
                }
            END_DTRACE
            
            norm_actual = matlib_znrm2(v2);
            matlib_zaxpy(-1.0, v2, v1);
            e_relative = matlib_znrm2(v1)/norm_actual;

            debug_exit("Relative error(%d): % 0.16g", j, e_relative);
            CU_ASSERT_TRUE(e_relative<TOL);

            matlib_free(M1.rowIn);
            matlib_free(M1.colIn);
            matlib_free(M1.elem_p);
            t += dt;
        }
    }
    fem1d_zm_nsparse_GMM(p, N, nsparse, Q, &phi, &q, &M, FEM1D_GMM_FREE);

}

void linear_timedependent_zpotential
( 
    matlib_xv      x, 
    matlib_real    t0, 
    matlib_real    dt, 
    matlib_zm      y
)
{
    debug_enter("%s", "");
    matlib_index i, j;
    for(j=0; j<y.lenr; j++)
    {
        for (i=0; i<x.len; i++)
        {
            *(y.elem_p) =  *(x.elem_p+i)*cos(2*M_PI*(t0));
            y.elem_p++;
        }
        t0 += dt;
    }
    debug_exit("%s", "");
}

void test_fem1d_zm_nsparse_GMM(void)
{
    matlib_index p, nr_LGL;
    matlib_index N = 50;
    matlib_real domain[2] = {-5.0, 5.0};
    matlib_real e_relative;

    p = 4;
    nr_LGL = 2*p+1;
    test_fem1d_zm_nsparse_GMM_general( p, nr_LGL, N, domain,
                                       linear_timedependent_zpotential);
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
        //{ "Interpolation for polynomial"           , test_interpolation1 },
        //{ "Interpolation for Gaussian"             , test_interpolation2 },
        //{ "Interpolation for Gaussian(t)"          , test_interpolation3 },
        //{ "Interpolation for complex Gaussian"     , test_interpolation4 },
        //{ "Interpolation for complex Gaussian(t)"  , test_interpolation5 },
        //{ "L2-Norm for Real"                       , test_L2_dnorm       },
        //{ "L2-Norm for Complex"                    , test_L2_znorm       },
        //{ "Transformation L2F, F2L for Real"       , test_fem1d_XL2F1    },
        //{ "Transformation L2F, F2L for Complex"    , test_fem1d_ZL2F1    },
        //{ "Quadrature Matrix"                      , test_fem1d_quadM1   },
        //{ "MEMI"                                   , test_fem1d_MEMI     },
        //{ "Global mass matrix for Gaussian real"   , test_fem1d_XGMM1    },
        //{ "Global mass matrix for Gaussian complex", test_fem1d_ZGMM1    },
        { "N-Sparse"                          , test_fem1d_zm_nsparse_GMM},
        CU_TEST_INFO_NULL,
    };

    /* Create the test suite */ 
    CU_SuiteInfo suites[] = 
    {
        { "Fem1d Lobatto basis", init_suite, clean_suite, test_array },
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

