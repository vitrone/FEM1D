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

/*============================================================================*/

double test_interpolation_general
(
    matlib_index  p,
    matlib_index  N,
    double domain[2],
    void   (*func_p)(matlib_dv, matlib_dv)
)
{
    debug_enter( "highest degree: %d, nr finite-elements: %d, domain: [% 0.3f, % 0.3f]",
                 p, N, domain[0], domain[1]);
    matlib_index i;

    /* Define computational domain */ 
    double x_l =  domain[0];
    double x_r =  domain[1];

    /* Define the number of finite-elements */ 

    /* Provide the degree of polynomial used for interpolation */ 
    matlib_index mesh_len = N*p+1;
    matlib_index nr_coeff = N*(p+1);

    matlib_dm FM;
    matlib_dv x, xi, quadW, u, U;

    legendre_LGLdataLT1( p, TOL, &xi, &quadW);

    matlib_create_dm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    legendre_LGLdataFM( xi, FM);

    matlib_create_dv( mesh_len, &u, MATLIB_COL_VECT);
    matlib_create_dv( nr_coeff, &U, MATLIB_COL_VECT);

    fem1d_ref2mesh (xi, N, x_l, x_r, &x);
    debug_body("length x: %d, length u: %d", x.len, u.len);
    
    (*func_p)(x, u);
    
    fem1d_DFLT( N, FM, u, U);
    DEBUG_PRINT_DV(U, "%s: ", "Transformed vector");
    
    /* Interpolation on arbitrary points: x1 */ 
    matlib_index P  = 2;
    mesh_len = N*P+1;
    matlib_index nr_samples = P+1;

    matlib_dm IMi;
    matlib_dv xii, u_interpol, u_actual, x1;

    matlib_create_dv( nr_samples, &xii,      MATLIB_COL_VECT);
    matlib_create_dv( mesh_len, &u_interpol, MATLIB_COL_VECT);
    matlib_create_dv( mesh_len, &u_actual,   MATLIB_COL_VECT);

    /* generate the grid on the reference interval */ 
    double dxii = 2.0/P;
    debug_body( "dxii=% 0.16f", dxii);

    *(xii.elem_p) = -1.0;
    for (i=1; i<P; i++)
    {
        *(xii.elem_p+i) = *(xii.elem_p+i-1) + dxii;
    }
    *(xii.elem_p+P) = 1.0;
    
    matlib_create_dm( xii.len, p+1, &IMi, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);
    legendre_LGLdataIM(xii, IMi);

    fem1d_DILT( N, IMi, U, u_interpol);


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
    double norm_actual = matlib_dnrm2(u_actual);
    matlib_daxpy(-1.0, u_actual, u_interpol);
    double e_relative = matlib_dnrm2(u_interpol)/norm_actual;
    
    /* Free allocated memory*/ 
    matlib_free( (void*) FM.elem_p);
    matlib_free( (void*) xi.elem_p);
    matlib_free( (void*) x.elem_p);
    matlib_free( (void*) u.elem_p);
    matlib_free( (void*) U.elem_p);
    matlib_free( (void*) quadW.elem_p);
    matlib_free( (void*) IMi.elem_p);
    matlib_free( (void*) xii.elem_p);
    matlib_free( (void*) x1.elem_p);
    matlib_free( (void*) u_actual.elem_p);
    matlib_free( (void*) u_interpol.elem_p);

    debug_exit("Relative interpolation error: % 0.16g", e_relative);
    return e_relative;

}
/*============================================================================*/

double test_interpolation_general2
(
    matlib_index     p,
    matlib_index     N,
    double    domain[2],
    matlib_dv t,
    void      (*func_p)(matlib_dv, matlib_dv, matlib_dm)
)
{
    debug_enter( "highest degree: %d, nr finite-elements: %d, domain: [% 0.3f, % 0.3f]",
                 p, N, domain[0], domain[1]);
    matlib_index i, j;

    /* Define computational domain */ 
    double x_l =  domain[0];
    double x_r =  domain[1];

    /* Define the number of finite-elements */ 

    /* Provide the degree of polynomial used for interpolation */ 
    matlib_index mesh_len = N*p+1;
    matlib_index nr_coeff = N*(p+1);

    matlib_dm FM, u, U;
    matlib_dv x, xi, quadW;

    legendre_LGLdataLT1( p, TOL, &xi, &quadW);

    matlib_create_dm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    legendre_LGLdataFM( xi, FM);

    matlib_create_dm( mesh_len, t.len, &u, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);
    matlib_create_dm( nr_coeff, t.len, &U, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);

    fem1d_ref2mesh (xi, N, x_l, x_r, &x);
    debug_body("length x: %d, size u: %d-by-%d", x.len, u.lenc, u.lenr);
    
    (*func_p)(x, t, u);
    
    fem1d_DFLT2( N, FM, u, U);
    DEBUG_PRINT_DM(U, "%s: ", "Transformed matrix");
    
    /* Interpolation on arbitrary points: x1 */ 
    matlib_index P  = 2;
    mesh_len = N*P+1;
    matlib_index nr_samples = P+1;

    matlib_dm IMi, u_interpol, u_actual;
    matlib_dv xii,  x1;

    matlib_create_dv( nr_samples, &xii,      MATLIB_COL_VECT);
    matlib_create_dm( mesh_len, t.len, &u_interpol, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);
    matlib_create_dm( mesh_len, t.len, &u_actual,   MATLIB_COL_MAJOR, MATLIB_NO_TRANS);

    /* generate the grid on the reference interval */ 
    double dxii = 2.0/P;
    debug_body( "dxii=% 0.16f", dxii);

    *(xii.elem_p) = -1.0;
    for (i=1; i<P; i++)
    {
        *(xii.elem_p+i) = *(xii.elem_p+i-1) + dxii;
    }
    *(xii.elem_p+P) = 1.0;
    
    matlib_create_dm( xii.len, p+1, &IMi, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);
    legendre_LGLdataIM(xii, IMi);

    fem1d_DILT2( N, IMi, U, u_interpol);


    /* find the actual values from the known function */ 
    fem1d_ref2mesh (xii, N, x_l, x_r, &x1);
    (*func_p)( x1, t, u_actual);
    
    /* Compute the error */
    double norm_actual, e_relative;

    matlib_index tmp_len = mesh_len * t.len; 

   /* 
    * matlib_dv u_tmp1 = MK_VM(u_interpol);
    * matlib_dv u_tmp2 = MK_VM(u_actual);
    * */

    matlib_dv u_tmp1 = { .len = tmp_len, 
                         .type = MATLIB_COL_VECT, 
                         .elem_p = u_interpol.elem_p};
    matlib_dv u_tmp2 = { .len = tmp_len, 
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

    norm_actual = matlib_dnrm2(u_tmp1);
    matlib_daxpy(-1.0, u_tmp2, u_tmp1);
    e_relative = matlib_dnrm2(u_tmp1)/norm_actual;
    
    /* Free allocated memory*/ 
    matlib_free( (void*) FM.elem_p);
    matlib_free( (void*) xi.elem_p);
    matlib_free( (void*) x.elem_p);
    matlib_free( (void*) u.elem_p);
    matlib_free( (void*) U.elem_p);
    matlib_free( (void*) quadW.elem_p);
    matlib_free( (void*) IMi.elem_p);
    matlib_free( (void*) xii.elem_p);
    matlib_free( (void*) x1.elem_p);
    matlib_free( (void*) u_actual.elem_p);
    matlib_free( (void*) u_interpol.elem_p);

    debug_exit("Relative interpolation error: % 0.16g", e_relative);
    return e_relative;

}
/*============================================================================*/

void test_interpolation1(void)
{

    matlib_index p;
    matlib_index N = 500;
    double domain[2] = {-1.0, 1.0};
    double e_relative;

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
    double domain[2] = {-5.0, 5.0};
    double e_relative;

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
    double domain[2] = {-5.0, 5.0};
    double ta[4] = {0, 0.5, 1.0, 1.5};
    matlib_dv t = {.len = 4, .type = MATLIB_COL_VECT, .elem_p = ta};

    double e_relative;
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

double test_interpolation_general3
(
    matlib_index p,
    matlib_index N
)
{
    matlib_index i;
    
    /* Define computational domain */ 
    double x_l = -5.0;
    double x_r = 5.0;

    /* Define the number of finite-elements */ 

    /* Provide the degree of polynomial used for interpolation */ 
    matlib_index mesh_len = N*p+1;
    matlib_index nr_coeff = N*(p+1);

    matlib_dm FM;
    matlib_dv x, xi, quadW;
    matlib_zv u, U;

    legendre_LGLdataLT1( p, TOL, &xi, &quadW);

    matlib_create_dm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
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

    matlib_dm IMi;
    matlib_dv xii, x1;
    matlib_zv u_interpol, u_actual;

    matlib_create_dv( nr_samples, &xii,      MATLIB_COL_VECT);
    matlib_create_zv( mesh_len, &u_interpol, MATLIB_COL_VECT);
    matlib_create_zv( mesh_len, &u_actual,   MATLIB_COL_VECT);

    /* generate the grid on the reference interval */ 
    double dxii = 2.0/P;
    debug_body( "dxii=% 0.16f", dxii);

    *(xii.elem_p) = -1.0;
    for (i=1; i<P; i++)
    {
        *(xii.elem_p+i) = *(xii.elem_p+i-1) + dxii;
    }
    *(xii.elem_p+P) = 1.0;
    
    matlib_create_dm( xii.len, p+1, &IMi, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);
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
    double norm_actual = matlib_znrm2(u_actual);
    matlib_zaxpy(-1.0, u_actual, u_interpol);
    double e_relative = matlib_znrm2(u_interpol)/norm_actual;
    
    /* Free allocated memory*/ 
    matlib_free( (void*) FM.elem_p);
    matlib_free( (void*) xi.elem_p);
    matlib_free( (void*) x.elem_p);
    matlib_free( (void*) u.elem_p);
    matlib_free( (void*) U.elem_p);
    matlib_free( (void*) quadW.elem_p);
    matlib_free( (void*) IMi.elem_p);
    matlib_free( (void*) xii.elem_p);
    matlib_free( (void*) x1.elem_p);
    matlib_free( (void*) u_actual.elem_p);
    matlib_free( (void*) u_interpol.elem_p);

    debug_exit("Relative interpolation error: % 0.16g", e_relative);
    return(e_relative);
}
/*============================================================================*/

double test_interpolation_general4
(
    matlib_index p,
    matlib_index N
)
{
    matlib_index i, j;
    
    /* Define computational domain */ 
    double x_l = -5.0;
    double x_r = 5.0;

    double ta[4] = {0, 0.5, 1.0, 1.5};
    matlib_dv t = { .len = 4, 
                    .type = MATLIB_COL_VECT, 
                    .elem_p = ta};

    /* Define the number of finite-elements */ 

    /* Provide the degree of polynomial used for interpolation */ 
    matlib_index mesh_len = N*p+1;
    matlib_index nr_coeff = N*(p+1);

    matlib_dm FM;
    matlib_dv x, xi, quadW;
    matlib_zm u, U;

    legendre_LGLdataLT1( p, TOL, &xi, &quadW);

    matlib_create_dm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
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

    matlib_dm IMi;
    matlib_dv xii, x1;
    matlib_zm u_interpol, u_actual;

    matlib_create_dv( nr_samples, &xii, MATLIB_COL_VECT);
    matlib_create_zm( mesh_len, t.len, &u_interpol, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);
    matlib_create_zm( mesh_len, t.len, &u_actual,   MATLIB_COL_MAJOR, MATLIB_NO_TRANS);

    /* generate the grid on the reference interval */ 
    double dxii = 2.0/P;
    debug_body( "dxii=% 0.16f", dxii);

    *(xii.elem_p) = -1.0;
    for (i=1; i<P; i++)
    {
        *(xii.elem_p+i) = *(xii.elem_p+i-1) + dxii;
    }
    *(xii.elem_p+P) = 1.0;
    
    matlib_create_dm( xii.len, p+1, &IMi, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);
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

    double norm_actual = matlib_znrm2(u_tmp1);
    matlib_zaxpy(-1.0, u_tmp2, u_tmp1);
    double e_relative = matlib_znrm2(u_tmp1)/norm_actual;
    
    /* Free allocated memory*/ 
    matlib_free( (void*) FM.elem_p);
    matlib_free( (void*) xi.elem_p);
    matlib_free( (void*) x.elem_p);
    matlib_free( (void*) u.elem_p);
    matlib_free( (void*) U.elem_p);
    matlib_free( (void*) quadW.elem_p);
    matlib_free( (void*) IMi.elem_p);
    matlib_free( (void*) xii.elem_p);
    matlib_free( (void*) x1.elem_p);
    matlib_free( (void*) u_actual.elem_p);
    matlib_free( (void*) u_interpol.elem_p);

    debug_exit("Relative interpolation error: % 0.16g", e_relative);
    return(e_relative);
}
/*============================================================================*/

void test_interpolation4(void)
{

    matlib_index p;
    matlib_index N = 800;
    double e_relative;

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

    double e_relative;
    for(p=2; p<12; p++)
    {
        e_relative = test_interpolation_general4( p, N);

        CU_ASSERT_TRUE(e_relative<TOL);
    }

}


/*============================================================================*/

double test_L2_dnorm_general
(
    matlib_index p,
    matlib_index N
)
{
    debug_enter("highest degree: %d, nr. fem-elements: %d", p, N);
    double x_l = -5.0;
    double x_r =  5.0;

    double J = (x_r-x_l)/(2.0*N);
    
    matlib_index mesh_len = N*p+1;
    matlib_index nr_coeff = N*(p+1);

    matlib_dv x, xi, quadW;
    matlib_dv u, U;
    matlib_dm FM;

    legendre_LGLdataLT1( p, TOL, &xi, &quadW);

    matlib_create_dm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    legendre_LGLdataFM( xi, FM);

    matlib_create_dv( mesh_len, &u, MATLIB_COL_VECT);
    matlib_create_dv( nr_coeff, &U, MATLIB_COL_VECT);

    fem1d_ref2mesh (xi, N, x_l, x_r, &x);
    debug_body("length x: %d, length u: %d", x.len, u.len);
    
    Gaussian_func(x, u);
    
    fem1d_DFLT( N, FM, u, U);
    DEBUG_PRINT_DV(U, "%s: ", "Transformed vector");

    double expected_norm = pow(M_PI, 0.25);
    double lp_norm = sqrt(J)*fem1d_DNorm2(p, N, U);
    double e_relative = fabs(expected_norm-lp_norm)/expected_norm;

    /* Free allocated memory*/ 
    matlib_free( (void*) FM.elem_p);
    matlib_free( (void*) xi.elem_p);
    matlib_free( (void*) x.elem_p);
    matlib_free( (void*) u.elem_p);
    matlib_free( (void*) U.elem_p);
    matlib_free( (void*) quadW.elem_p);

    debug_exit("Relative error in norm2: % 0.16g", e_relative);
    return e_relative;

}

double test_L2_znorm_general
(
    matlib_index p,
    matlib_index N
)
{
    debug_enter("highest degree: %d, nr. fem-elements: %d", p, N);
    double x_l = -5.0;
    double x_r =  5.0;

    double J = (x_r-x_l)/(2.0*N);
    
    matlib_index mesh_len = N*p+1;
    matlib_index nr_coeff = N*(p+1);

    matlib_dv x, xi, quadW;
    matlib_zv u, U;
    matlib_dm FM;

    legendre_LGLdataLT1( p, TOL, &xi, &quadW);

    matlib_create_dm( p+1, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    legendre_LGLdataFM( xi, FM);

    matlib_create_zv( mesh_len, &u, MATLIB_COL_VECT);
    matlib_create_zv( nr_coeff, &U, MATLIB_COL_VECT);

    fem1d_ref2mesh (xi, N, x_l, x_r, &x);
    debug_body("length x: %d, length u: %d", x.len, u.len);
    
    Gaussian_zfunc(x, u);
    
    fem1d_ZFLT( N, FM, u, U);
    DEBUG_PRINT_ZV(U, "%s: ", "Transformed vector");

    double expected_norm = pow(M_PI, 0.25);
    double lp_norm = sqrt(J)*fem1d_ZNorm2(p, N, U);
    double e_relative = fabs(expected_norm-lp_norm)/expected_norm;

    /* Free allocated memory*/ 
    matlib_free( (void*) FM.elem_p);
    matlib_free( (void*) xi.elem_p);
    matlib_free( (void*) x.elem_p);
    matlib_free( (void*) u.elem_p);
    matlib_free( (void*) U.elem_p);
    matlib_free( (void*) quadW.elem_p);

    debug_exit("Relative error in norm2: % 0.16g", e_relative);
    return e_relative;

}
/*============================================================================*/

void test_L2_dnorm(void)
{
    matlib_index N = 300;
    matlib_index p_max = 15;
    double e_relative;
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
    double e_relative;
    for(matlib_index p=2; p<p_max; p++)
    {
        e_relative = test_L2_znorm_general(p, N);
        CU_ASSERT_TRUE(e_relative<TOL);
    }
}
/*============================================================================*/
double test_fem1d_DL2F_general
(
    matlib_index p,
    matlib_index N
)
{
    double x_l =  -5.0;
    double x_r =   5.0;

    debug_enter( "polynomial degree: %d, nr. of LGL points: %d", p, p+1 );

    matlib_dv x, xi, quadW;
    legendre_LGLdataLT1( p, TOL, &xi, &quadW);
    fem1d_ref2mesh (xi, N, x_l, x_r, &x);
    
    matlib_dm FM;
    matlib_create_dm( xi.len, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
    legendre_LGLdataFM( xi, FM);
    
    
    matlib_dv u, U, U1;
    matlib_create_dv( x.len,    &u, MATLIB_COL_VECT);
    matlib_create_dv( N*(p+1),  &U, MATLIB_COL_VECT);
    matlib_create_dv( N*(p+1), &U1, MATLIB_COL_VECT);

    Gaussian_func(x, u);
    DEBUG_PRINT_DV(u, "%s: ", "value of Gaussian");
    fem1d_DFLT( N, FM, u, U); /* Legendre transform */ 

    matlib_dv vb;
    matlib_create_dv( N*p+1, &vb, MATLIB_COL_VECT);
    
    /* Transformation from Lgendre basis to FEM-shape func basis */ 
    fem1d_DL2F(p, U, vb);  
    DEBUG_PRINT_DV(vb, "%s: ", "coeff. fem basis func.");

    fem1d_DF2L(p, vb, U1);  
    BEGIN_DTRACE
        matlib_index i;
        for(i=0;i<U.len; i++)
        {
            debug_print( "[%d]->U: % 0.16f, U1: % 0.16f", 
                         i, *(U.elem_p+i), *(U1.elem_p+i));
        }
    END_DTRACE
    double norm_actual = matlib_dnrm2(U);
    matlib_daxpy(-1.0, U, U1);
    double e_relative = matlib_dnrm2(U1)/norm_actual;

    return(e_relative);
    debug_exit("Relative interpolation error: % 0.16g", e_relative);
}

double test_fem1d_ZL2F_general
(
    matlib_index p,
    matlib_index N
)
{
    double x_l =  -5.0;
    double x_r =   5.0;

    debug_enter( "polynomial degree: %d, nr. of LGL points: %d", p, p+1 );

    matlib_dv x, xi, quadW;
    legendre_LGLdataLT1( p, TOL, &xi, &quadW);
    fem1d_ref2mesh (xi, N, x_l, x_r, &x);
    
    matlib_dm FM;
    matlib_create_dm( xi.len, xi.len, &FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    
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
    double norm_actual = matlib_znrm2(U);
    matlib_zaxpy(-1.0, U, U1);
    double e_relative = matlib_znrm2(U1)/norm_actual;

    return(e_relative);
    debug_exit("Relative interpolation error: % 0.16g", e_relative);
}
/*============================================================================*/

void test_fem1d_DL2F1(void)
{
    matlib_index N = 40;
    matlib_index p_max = 15;
    double e_relative;
    for(matlib_index p=2; p<p_max; p++)
    {
        e_relative = test_fem1d_DL2F_general(p, N);
        CU_ASSERT_TRUE(e_relative<TOL);
    }

}
void test_fem1d_ZL2F1(void)
{
    matlib_index N = 40;
    matlib_index p_max = 15;
    double e_relative;
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

    matlib_dv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);

    matlib_dm IMi, Q, IMii, Q1;

    matlib_create_dm( xi.len, p+1,  &IMi, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    
    matlib_create_dm( xi.len, p+1, &IMii, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataIM(xi, IMi);
    legendre_LGLdataIM(xi, IMii);

    DEBUG_PRINT_DM(IMi,  "%s: ", "inverse transform matrix (IM: col major)");
    DEBUG_PRINT_DM(IMii, "%s: ", "inverse transform matrix (IM: row major)");
    
    fem1d_quadM(quadW,  IMi,  &Q);
    fem1d_quadM(quadW, IMii, &Q1);
    
    DEBUG_PRINT_DM(Q,  "%s: ", "Quadrature matrix (IM: col major)");
    DEBUG_PRINT_DM(Q1, "%s: ", "Quadrature matrix (IM: row major)");

    matlib_dv u_tmp1 = MK_VM(Q);
    matlib_dv u_tmp2 = MK_VM(Q1);

    double norm_actual = matlib_dnrm2(u_tmp1);
    
    matlib_daxpy(-1.0, u_tmp2, u_tmp1);
    DEBUG_PRINT_DV(u_tmp1, "%s: ", "Q-Q1");
    
    double e_relative = matlib_dnrm2(u_tmp1)/norm_actual;
    debug_exit("Relative error: % 0.16g", e_relative);


    CU_ASSERT_TRUE(e_relative<TOL);
    debug_exit("%s","");
}

/*============================================================================*/

double test_fem1d_MEMI_general( matlib_index p)
{
    debug_enter("polynomial degree: %d", p);
    matlib_index i;
    matlib_index P = p+1;

    matlib_dv q;
    fem1d_MEMI(p, &q);
    DEBUG_PRINT_DV(q, "%s :", "Master element mass integrals");

    matlib_dv xi, quadW;
    legendre_LGLdataLT1( P, TOL, &xi, &quadW);
    
    matlib_dm FM, IM, Q;
    matlib_create_dm( xi.len, p+1, &IM, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataIM( xi, IM);
    fem1d_quadM( quadW, IM, &Q);
    DEBUG_PRINT_DM(Q,  "%s: ", "Quadrature matrix (IM: col major)");

    matlib_dv q1, ones;
    matlib_create_dv( Q.lenc,   &q1, MATLIB_COL_VECT);
    matlib_create_dv( xi.len, &ones, MATLIB_COL_VECT);
    
    for(i=0; i<xi.len; i++)
    {
        ones.elem_p[i] = 1.0;
    }

    matlib_dgemv( 1.0, Q, ones, 0, q1);

    BEGIN_DTRACE
        for(i=0;i<q.len; i++)
        {
            debug_print( "[%d]-> q: % 0.16f, q1: % 0.16f", 
                         i, *(q.elem_p+i), *(q1.elem_p+i));
        }
    END_DTRACE

    double norm_actual = matlib_dnrm2(q);
    
    matlib_daxpy(-1.0, q, q1);
    double e_relative = matlib_dnrm2(q1)/norm_actual;
    debug_exit("Relative error: % 0.16g", e_relative);

    return(e_relative);
}
void test_fem1d_MEMI(void)
{
    matlib_index p[11] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    double e_relative;
    for(matlib_index i=0; i<11; i++)
    {
        e_relative = test_fem1d_MEMI_general(p[i]);
        CU_ASSERT_TRUE(e_relative<TOL);
    }
}
/*============================================================================*/
void full_sym_mat
(
    matlib_dsparsem sM,
    matlib_dm*      M
)
{
    matlib_index i, j;
    
    matlib_create_dm( sM.lenc, sM.lenr, M, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    

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
    matlib_dv Pvb;
    matlib_create_dv( dim, &Pvb, MATLIB_COL_VECT);

    matlib_dcsrsymv(MATLIB_UPPER, M, Uvb, Pvb);
       
    /* compute projection from the LP representation of phi(x)*u(x)*/
    matlib_dv Pvb1;
    matlib_create_dv( dim, &Pvb1, MATLIB_COL_VECT);

    for(i=0; i<x.len; i++)
    {
        y.elem_p[i] = u.elem_p[i] * y.elem_p[i];
    
    }
    fem1d_DFLT( N, FM, y, Y);
    fem1d_DPrjL2F(p, Y, Pvb1);

    BEGIN_DTRACE
        for(i=0;i<dim; i++)
        {
            debug_print( "[%d]-> Pvb: % 0.16f, Pvb1: % 0.16f", 
                         i, *(Pvb.elem_p+i), *(Pvb1.elem_p+i));
        }
    END_DTRACE

    double norm_actual = matlib_dnrm2(Pvb1);
    
    matlib_daxpy(-1.0, Pvb1, Pvb);

    double e_relative = matlib_dnrm2(Pvb)/norm_actual;

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
    matlib_index N = 600;
    double domain[2] = {-5.0, 5.0};
    double e_relative;
    for(matlib_index p=2; p<13; p++)
    {
        nr_LGL = 2*p+1;
        e_relative = test_fem1d_DGMM_general( p, nr_LGL, N, domain,
                                              Gaussian_func, 
                                              constant_dpotential);
        CU_ASSERT_TRUE(e_relative<TOL);

        e_relative = test_fem1d_DGMM_general( p, nr_LGL, N, domain,
                                              Gaussian_func, 
                                              harmonic_dpotential);
        CU_ASSERT_TRUE(e_relative<TOL);

        e_relative = test_fem1d_DGMM_general( p, nr_LGL, N, domain,
                                              Gaussian_func, 
                                              sinusoidal_dpotential);
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

    double norm_actual = matlib_znrm2(Pvb1);
    
    matlib_zaxpy(-1.0, Pvb1, Pvb);

    double e_relative = matlib_znrm2(Pvb)/norm_actual;

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
    matlib_index N = 1600;
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
        { "Interpolation for polynomial"           , test_interpolation1 },
        { "Interpolation for Gaussian"             , test_interpolation2 },
        { "Interpolation for Gaussian(t)"          , test_interpolation3 },
        { "Interpolation for complex Gaussian"     , test_interpolation4 },
        { "Interpolation for complex Gaussian(t)"  , test_interpolation5 },
        { "L2-Norm for Real"                       , test_L2_dnorm       },
        { "L2-Norm for Complex"                    , test_L2_znorm       },
        { "Transformation L2F, F2L for Real"       , test_fem1d_DL2F1    },
        { "Transformation L2F, F2L for Complex"    , test_fem1d_ZL2F1    },
        { "Quadrature Matrix"                      , test_fem1d_quadM1   },
        { "MEMI"                                   , test_fem1d_MEMI     },
        { "Global mass matrix for Gaussian real"   , test_fem1d_DGMM1    },
        { "Global mass matrix for Gaussian complex", test_fem1d_ZGMM1    },
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
