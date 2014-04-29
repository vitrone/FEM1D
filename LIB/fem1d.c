/*============================================================================+/
 | 
 |  Name   : fem1d.c                                                      
 |  Author : Vishal Vaibhav                                               
 |                                                                        
 |  Description : This module defines all the functions for implementing  
 |  finite element method. See the documentation for mathematical details.
 |                                                                        
 |  History : Created 21 July 2013.                                       
 |
/+============================================================================*/

#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mkl.h"

#define NDEBUG
#define MATLIB_NTRACE_DATA

#include "fem1d.h"
#include "assert.h"

/*============================================================================*/
void fem1d_ref2mesh
(
    matlib_xv    xi,
    matlib_index N,
    matlib_real  x_l,
    matlib_real  x_r,
    matlib_xv*   x
)
/*============================================================================+/
 | Description: This function computes the transformation of points from      | 
 | reference domain to the computational domain.                              |
 |                                                                            |
 | Input(s) :                                                                 |
 |  xi  - array of points in the reference domian including 1 and -1          |
 |  N   - number of finite elements in the computational domain               |
 |  x_l - left boundary                                                       |
 |  x_r - right boundary                                                      |
 |                                                                            |
 | Output(s):                                                                 |
 |  x - transformed coordinate points of the from the reference domain        |
 |                                                                            |
/+============================================================================*/

{
    debug_enter( "nr. of points of reference domain: %d,"
                 "nr. of finite-elements: %d", xi.len, N );
    matlib_index p = xi.len-1;
    matlib_index mesh_len = N*p+1;

    matlib_create_xv( mesh_len, x, MATLIB_COL_VECT);
    
    matlib_real dx = (x_r-x_l)/N;
    matlib_real J  = dx/2.0;
    matlib_real tmp;
    matlib_index ptmp, i, j;
    
    /* xm_i=x_l+J */
    tmp  = x_l+J;
    ptmp = 0;
    for (i=0; i<N; i++)
    {
        for(j=0; j<p; j++)
        {
            *(x->elem_p+ptmp+j) = J**(xi.elem_p+j)+tmp;
        }
        tmp  = tmp+dx;
        ptmp = ptmp+p;
    }
    *(x->elem_p+ptmp) = x_r;
    
    /* length of x: (n_xi-1)*N + 1 */ 
    debug_exit("length of x: %d", mesh_len);
}

/*============================================================================+/
 | Legendre transformation routines
/+============================================================================*/
void fem1d_XFLT
(
    const matlib_index N,
    const matlib_xm    FM,
          matlib_xv    u,
          matlib_xv    U
)
/* 
 * N : can be deduced from the length of the matrices provided but
 *     this leads to untraceble bugs. Hence this is made as an argument.
 * */ 
{
    debug_enter( "nr. finite-elements: %d, "
                 "matrix FM: %d-by-%d, "
                 "vectors u: %d, U:%d", 
                 N, FM.lenc, FM.lenr, u.len, U.len );

    /* highest degree of Legendre polynomials retained 
     * */
    matlib_index p = FM.lenc-1; 
    
    /* nr. LGL-points -1 used for sampling in FEM-element 
     * */ 
    matlib_index P = FM.lenr-1; 
    
    debug_body( "degree of polynomial: %d,"
                "nr. sampling point: %d", 
                p, FM.lenr );

    assert(((FM.elem_p!=NULL) && (u.elem_p !=NULL)) && (U.elem_p!=NULL));

    if(P>1)
    {
        debug_body("nr. computed finite-elements: %d", (u.len-1)/P);
        assert(N == (u.len-1)/P);
        if (U.len == N*(p+1))
        {
            matlib_index strideFM;
            _CBLAS_ORDER order_enum;
            GET_CBLASORDER_AND_STRIDE(order_enum, strideFM, FM);

            matlib_real one = 1.0, zero = 0;
            matlib_index incu = 1;
            matlib_index i;
            for(i=0; i<N; i++)
            {
                cblas_dgemv( order_enum, 
                             CblasNoTrans, 
                             FM.lenc, 
                             FM.lenr, 
                             one, 
                             FM.elem_p, 
                             strideFM, 
                             u.elem_p,
                             incu,
                             zero, 
                             U.elem_p, 
                             incu);
                (u.elem_p) += P;
                (U.elem_p) += (FM.lenc);
            }
        }
        else
        {
            term_execb( "size of vectors/matrices incorrect: matrix "
                        "FM: %d-by-%d, vectors u: %d, U:%d", 
                        FM.lenc, FM.lenr, u.len, U.len );
        }
    }
    else
    {
        term_exec( "incorrect nr. of sampling points in ref. domain: %d",
                    FM.lenr );
    }
    debug_exit("%s","");
}

void fem1d_ZFLT
(
    const matlib_index N,
    const matlib_xm    FM,
          matlib_zv    u,
          matlib_zv    U
)
/* 
 * N : can be deduced from the length of the matrices provided but
 *     this leads to untraceble bugs. Hence this is made as an argument.
 * */ 
{
    debug_enter( "nr. finite-elements: %d, "
                 "matrix FM: %d-by-%d, "
                 "vectors u: %d, U:%d", 
                 N, FM.lenc, FM.lenr, u.len, U.len );

    /* highest degree of Legendre polynomials retained 
     * */
    matlib_index p = FM.lenc-1; 
    
    /* nr. LGL-points -1 used for sampling in FEM-element 
     * */ 
    matlib_index P = FM.lenr-1; 
    
    debug_body( "degree of polynomial: %d,"
                "nr. sampling point: %d", 
                p, FM.lenr );

    assert(((FM.elem_p!=NULL) && (u.elem_p !=NULL)) && (U.elem_p!=NULL));

    if(P>1)
    {
        debug_body("nr. computed finite-elements: %d", (u.len-1)/P);
        assert(N == (u.len-1)/P);
        if (U.len == N*(p+1))
        {
            matlib_index strideFM;
            _CBLAS_ORDER order_enum;
            GET_CBLASORDER_AND_STRIDE(order_enum, strideFM, FM);

            matlib_real one = 1.0, zero = 0;
            matlib_index incu = 2;
            matlib_index i;
            for(i=0; i<N; i++)
            {
                cblas_dgemv( order_enum, 
                             CblasNoTrans, 
                             FM.lenc, 
                             FM.lenr, 
                             one, 
                             FM.elem_p, 
                             strideFM, 
                             (matlib_real*)u.elem_p,
                             incu,
                             zero, 
                             (matlib_real*)U.elem_p, 
                             incu);
                cblas_dgemv( order_enum, 
                             CblasNoTrans, 
                             FM.lenc, 
                             FM.lenr, 
                             one, 
                             FM.elem_p, 
                             strideFM, 
                             ((matlib_real*)u.elem_p+1),
                             incu,
                             zero, 
                             ((matlib_real*)U.elem_p+1), 
                             incu);
                (u.elem_p) += P;
                (U.elem_p) += (FM.lenc);
            }
        }
        else
        {
            term_execb( "size of vectors/matrices incorrect: matrix "
                        "FM: %d-by-%d, vectors u: %d, U:%d", 
                        FM.lenc, FM.lenr, u.len, U.len );
        }
    }
    else
    {
        term_exec( "incorrect nr. of sampling points in ref. domain: %d",
                    FM.lenr );
    }
    debug_exit("%s","");
}

/*============================================================================*/
void fem1d_XILT
(
    const matlib_index N,
    const matlib_xm    IM,
          matlib_xv    U,
          matlib_xv    u
)
/* 
 * N : can be deduced from the length of the matrices provided but
 *     this leads to untraceble bugs. Hence this is made as an argument.
 * */ 
{
    debug_enter( "nr. finite-elements: %d, "
                 "matrix IM: %d-by-%d, "
                 "vectors u: %d, U:%d", 
                 N, IM.lenc, IM.lenr, U.len, u.len );

    /* nr. LGL-points-1 used for sampling in FEM-element 
     * */
    matlib_index P = IM.lenc-1; 
    
    /* Highest degree of Legendre polynomials
     * */
    matlib_index p = IM.lenr-1; 

    debug_body( "degree of polynomial: %d, "
                "nr. sampling point: %d", 
                p, IM.lenc );

    assert(((IM.elem_p!=NULL) && (u.elem_p !=NULL)) && (U.elem_p!=NULL));

    if(p>1)
    {
        debug_body("nr. computed finite-elements: %d", U.len/(p+1));
        assert(N == U.len/(p+1));
        if (u.len == N*P+1)
        {
            matlib_index strideIM;
            _CBLAS_ORDER order_enum;
            
            GET_CBLASORDER_AND_STRIDE(order_enum, strideIM, IM);
            
            matlib_index incU = 1;
            matlib_index incu = 1;
            matlib_index i;
            for(i=0; i<N; i++)
            {
                cblas_dgemv( order_enum, 
                             CblasNoTrans, 
                             IM.lenc, 
                             IM.lenr,
                             1.0, 
                             IM.elem_p, 
                             strideIM, 
                             U.elem_p, 
                             incU, 
                             0,
                             u.elem_p,
                             incu);
                (U.elem_p) += (IM.lenr);
                (u.elem_p) += P;
            }
        }
        else
        {
            term_execb( "size of vectors/matrices incorrect, "
                        "IM: %d-by-%d, U: %d, u: %d",
                        IM.lenc, IM.lenr, U.len, u.len );
        }
    }
    else
    {
        term_exec( "highest degree of polynomials incorrect: %d",
                    IM.lenr );
    }
    debug_exit("%s","");
}

void fem1d_ZILT
(
    const matlib_index N,
    const matlib_xm    IM,
          matlib_zv    U,
          matlib_zv    u
)
/* 
 * N : can be deduced from the length of the matrices provided but
 *     this leads to untraceble bugs. Hence this is made as an argument.
 * */ 
{
    debug_enter( "nr. finite-elements: %d, "
                 "matrix IM: %d-by-%d, "
                 "vectors u: %d, U:%d", 
                 N, IM.lenc, IM.lenr, U.len, u.len );

    /* nr. LGL-points-1 used for sampling in FEM-element 
     * */
    matlib_index P = IM.lenc-1; 
    
    /* Highest degree of Legendre polynomials
     * */
    matlib_index p = IM.lenr-1; 

    debug_body( "degree of polynomial: %d, "
                "nr. sampling point: %d", 
                p, IM.lenc );

    assert(((IM.elem_p!=NULL) && (u.elem_p !=NULL)) && (U.elem_p!=NULL));

    if(p>1)
    {
        debug_body("nr. computed finite-elements: %d", U.len/(p+1));
        assert(N == U.len/(p+1));
        if (u.len == N*P+1)
        {
            matlib_index strideIM;
            _CBLAS_ORDER order_enum;
            
            GET_CBLASORDER_AND_STRIDE(order_enum, strideIM, IM);
            
            matlib_index incu = 2;
            matlib_index i;
            for(i=0; i<N; i++)
            {
                cblas_dgemv( order_enum, 
                             CblasNoTrans, 
                             IM.lenc, 
                             IM.lenr,
                             1.0, 
                             IM.elem_p, 
                             strideIM, 
                             (matlib_real*)U.elem_p, 
                             incu, 
                             0,
                             (matlib_real*)u.elem_p,
                             incu);
                cblas_dgemv( order_enum, 
                             CblasNoTrans, 
                             IM.lenc, 
                             IM.lenr,
                             1.0, 
                             IM.elem_p, 
                             strideIM, 
                             ((matlib_real*)U.elem_p+1), 
                             incu, 
                             0,
                             ((matlib_real*)u.elem_p+1),
                             incu);
                (U.elem_p) += (IM.lenr);
                (u.elem_p) += P;
            }
        }
        else
        {
            term_execb( "size of vectors/matrices incorrect: ",
                        IM.lenc, IM.lenr, U.len, u.len );
        }
    }
    else
    {
        term_exec( "highest degree of polynomials incorrect: %d",
                    IM.lenr );
    }
    debug_exit("%s","");
}
/*============================================================================*/

void fem1d_XFLT2
(
    const matlib_index N,
    const matlib_xm    FM,
          matlib_xm    u,
          matlib_xm    U
)
/* 
 * N : can be deduced from the length of the matrices provided but
 *     this leads to untraceble bugs. Hence this is made as an argument.
 * */ 
{
    debug_enter( "nr. finite-elements: %d, matrices FM: %d-by-%d, "
                 "u: %d-by-%d, U: %d-by-%d", N, 
                 FM.lenc, FM.lenr, 
                 u.lenc , u.lenr,
                 U.lenc , U.lenr );

   /* highest degree of Legendre polynomials retained 
    * */ 
    matlib_index p = FM.lenc-1; 
    
    /* nr. LGL-points -1 used for sampling in FEM-element 
     * */
    matlib_index P = FM.lenr-1;  


    bool order_OK = (u.order == MATLIB_COL_MAJOR) && 
                    (U.order == MATLIB_COL_MAJOR);

    if(!order_OK)
    {
        term_exec( "Storage order of matrices incorrect: %s", "u and U");
    }

    assert(((FM.elem_p!=NULL) && (u.elem_p !=NULL)) && (U.elem_p!=NULL));
    assert(P>1);

    if(P>1)
    {
        debug_body("nr. computed finite-elements: %d", (u.lenc-1)/P);
        assert(N == (u.lenc-1)/P);

        if ((U.lenc == N*(p+1)) && (u.lenr == U.lenr))
        {
            debug_body( "highest degree: %d, "
                         "nr. samples on ref. domain: %d, "
                         "nr. of transform vectors: %d",
                         p, FM.lenr, u.lenr);

            matlib_index strideFM;
            _CBLAS_ORDER order_enum;

            GET_CBLASORDER_AND_STRIDE(order_enum, strideFM, FM);

            matlib_index incu = 1;
            matlib_index incU = 1;
            matlib_index i, j;
            
            for (j=0; j<u.lenr; j++)
            {
                for(i=0; i<N; i++)
                {
                    cblas_dgemv( order_enum, 
                                 CblasNoTrans, 
                                 FM.lenc, 
                                 FM.lenr, 
                                 1.0, 
                                 FM.elem_p, 
                                 strideFM,
                                 u.elem_p, 
                                 incu, 
                                 0, 
                                 U.elem_p, 
                                 incU);

                    (u.elem_p) += P;
                    (U.elem_p) += FM.lenc;
                }
                (u.elem_p)++;
            }
        }
        else
        {
            term_exec( "size of matrices incorrect: "
                       "FM: %d-by-%d, u: %d-by-%d, U: %d-by-%d", N, 
                       FM.lenc, FM.lenr, 
                       u.lenc , u.lenr,
                       U.lenc , U.lenr );
        }
    }
    else
    {
        term_exec( "incorrect nr. of sampling points in ref. domain: %d",
                    FM.lenr );
    }

    debug_exit("%s","");
}

void fem1d_ZFLT2
(
    const matlib_index N,
    const matlib_xm    FM,
          matlib_zm    u,
          matlib_zm    U
)
/* 
 * N : can be deduced from the length of the matrices provided but
 *     this leads to untraceble bugs. Hence this is made as an argument.
 * */ 
{
    debug_enter( "nr. finite-elements: %d, matrices FM: %d-by-%d, "
                 "u: %d-by-%d, U: %d-by-%d", N, 
                 FM.lenc, FM.lenr, 
                 u.lenc , u.lenr,
                 U.lenc , U.lenr );

   /* highest degree of Legendre polynomials retained 
    * */ 
    matlib_index p = FM.lenc-1; 
    
    /* nr. LGL-points -1 used for sampling in FEM-element 
     * */
    matlib_index P = FM.lenr-1;  


    bool order_OK = (u.order == MATLIB_COL_MAJOR) && 
                    (U.order == MATLIB_COL_MAJOR);

    if(!order_OK)
    {
        term_exec( "Storage order of matrices incorrect: %s", "u and U");
    }

    assert(((FM.elem_p!=NULL) && (u.elem_p !=NULL)) && (U.elem_p!=NULL));
    assert(P>1);

    if(P>1)
    {
        debug_body("nr. computed finite-elements: %d", (u.lenc-1)/P);
        assert(N == (u.lenc-1)/P);

        if ((U.lenc == N*(p+1)) && (u.lenr == U.lenr))
        {
            debug_body( "highest degree: %d, "
                         "nr. samples on ref. domain: %d, "
                         "nr. of transform vectors: %d",
                         p, FM.lenr, u.lenr);

            matlib_index strideFM;
            _CBLAS_ORDER order_enum;

            GET_CBLASORDER_AND_STRIDE(order_enum, strideFM, FM);

            matlib_index incu = 2;
            matlib_index incU = 2;
            matlib_index i, j;
            
            for (j=0; j<u.lenr; j++)
            {
                for(i=0; i<N; i++)
                {
                    cblas_dgemv( order_enum, 
                                 CblasNoTrans, 
                                 FM.lenc, 
                                 FM.lenr, 
                                 1.0, 
                                 FM.elem_p, 
                                 strideFM,
                                 (matlib_real*)u.elem_p, 
                                 incu, 
                                 0, 
                                 (matlib_real*)U.elem_p, 
                                 incU);
                    cblas_dgemv( order_enum, 
                                 CblasNoTrans, 
                                 FM.lenc, 
                                 FM.lenr, 
                                 1.0, 
                                 FM.elem_p, 
                                 strideFM,
                                 (matlib_real*)u.elem_p+1, 
                                 incu, 
                                 0, 
                                 (matlib_real*)U.elem_p+1, 
                                 incU);

                    (u.elem_p) += P;
                    (U.elem_p) += FM.lenc;
                }
                (u.elem_p)++;
            }
        }
        else
        {
            term_exec( "size of matrices incorrect: "
                       "FM: %d-by-%d, u: %d-by-%d, U: %d-by-%d", N, 
                       FM.lenc, FM.lenr, 
                       u.lenc , u.lenr,
                       U.lenc , U.lenr );
        }
    }
    else
    {
        term_exec( "incorrect nr. of sampling points in ref. domain: %d",
                    FM.lenr );
    }

    debug_exit("%s","");
}
/*======================================================================*/
void fem1d_XILT2
(
    const matlib_index N,
    const matlib_xm    IM,
          matlib_xm    U,
          matlib_xm    u
)
/* 
 * N : can be deduced from the length of the matrices provided but
 *     this leads to untraceble bugs. Hence this is made as an argument.
 * */ 
{
    debug_enter( "nr. finite-elements: %d, matrices FM: %d-by-%d, "
                 "u: %d-by-%d, U: %d-by-%d", N, 
                 IM.lenc, IM.lenr, 
                 u.lenc , u.lenr,
                 U.lenc , U.lenr );


    /* nr. LGL-points-1 used for sampling in FEM-element 
     * */
    matlib_index P = IM.lenc-1; 
    
    /* Highest degree of Legendre polynomials 
     * */
    matlib_index p = IM.lenr-1; 

    bool order_OK = (U.order == MATLIB_COL_MAJOR) && 
                    (u.order == MATLIB_COL_MAJOR);
    if(!order_OK)
    {
        term_exec( "Storage order of matrices incorrect: %s", "u and U");
    }

    assert(((IM.elem_p != NULL) && (u.elem_p != NULL)) && (U.elem_p != NULL));

    if(p>1)
    {
        debug_body( "nr. computed finite-elements: %d", U.lenc/(p+1));
        assert(N == U.lenc/(p+1));
        if ((u.lenc == N*P+1) && (U.lenr == u.lenr))
        {
            debug_enter( "highest degree: %d, "
                         "nr. samples on ref. domain: %d, "
                         "nr. of transform vectors: %d",
                         p, P+1, u.lenr);

            matlib_index strideIM;
            _CBLAS_ORDER order_enum;

            GET_CBLASORDER_AND_STRIDE(order_enum, strideIM, IM);

            matlib_index incU = 1;
            matlib_index incu = 1;
            matlib_index i, j;

            for(i=0; i<N; i++)
            {
                cblas_dgemv( order_enum, 
                             CblasNoTrans, 
                             IM.lenc, 
                             IM.lenr, 
                             1.0, 
                             IM.elem_p, 
                             strideIM, 
                             U.elem_p, 
                             incU, 
                             0, 
                             u.elem_p, 
                             incu);
                (U.elem_p) += IM.lenr;
                (u.elem_p) += P;
            }
            
            for (j=1; j<u.lenr; j++)
            {
                (u.elem_p)++;
                for(i=0; i<N; i++)
                {
                    cblas_dgemv( order_enum, 
                                 CblasNoTrans, 
                                 IM.lenc, 
                                 IM.lenr, 
                                 1.0, 
                                 IM.elem_p, 
                                 strideIM, 
                                 U.elem_p, 
                                 incU, 
                                 0, 
                                 u.elem_p, 
                                 incu);

                    (U.elem_p) += IM.lenr;
                    (u.elem_p) += P;
                }
            }
        }
        else
        {
            term_execb( "length of matrices incorrect: matrix "
                        "IM: %d-by-%d, U: %d-by-%d, u: %d-by-%d",
                        IM.lenc, IM.lenr, 
                        U.lenc, U.lenr, 
                        u.lenc, u.lenr );
        }
    }
    else
    {
        term_exec( "highest degree of polynomials incorrect: %d",
                    IM.lenr );
    }
    debug_exit("%s", "");

}

void fem1d_ZILT2
(
    const matlib_index N,
    const matlib_xm    IM,
          matlib_zm    U,
          matlib_zm    u
)
/* 
 * N : can be deduced from the length of the matrices provided but
 *     this leads to untraceble bugs. Hence this is made as an argument.
 * */ 
{
    debug_enter( "nr. finite-elements: %d, matrices IM: %d-by-%d, "
                 "u: %d-by-%d, U: %d-by-%d", N, 
                 IM.lenc, IM.lenr, 
                 u.lenc , u.lenr,
                 U.lenc , U.lenr );


    /* nr. LGL-points-1 used for sampling in FEM-element 
     * */
    matlib_index P = IM.lenc-1; 
    
    /* Highest degree of Legendre polynomials 
     * */
    matlib_index p = IM.lenr-1; 

    bool order_OK = (U.order == MATLIB_COL_MAJOR) && 
                    (u.order == MATLIB_COL_MAJOR);
    if(!order_OK)
    {
        term_exec( "Storage order of matrices incorrect: %s", "u and U");
    }

    assert(((IM.elem_p != NULL) && (u.elem_p != NULL)) && (U.elem_p != NULL));

    if(p>1)
    {
        debug_body( "nr. computed finite-elements: %d", U.lenc/(p+1));
        assert(N == U.lenc/(p+1));
        if ((u.lenc == N*P+1) && (U.lenr == u.lenr))
        {
            debug_enter( "highest degree: %d, "
                         "nr. samples on ref. domain: %d, "
                         "nr. of transform vectors: %d",
                         p, P+1, u.lenr);

            matlib_index strideIM;
            _CBLAS_ORDER order_enum;

            GET_CBLASORDER_AND_STRIDE(order_enum, strideIM, IM);

            matlib_index incU = 2;
            matlib_index incu = 2;
            matlib_index i, j;

            for(i=0; i<N; i++)
            {
                cblas_dgemv( order_enum, 
                             CblasNoTrans, 
                             IM.lenc, 
                             IM.lenr, 
                             1.0, 
                             IM.elem_p, 
                             strideIM, 
                             (matlib_real*)U.elem_p, 
                             incU, 
                             0, 
                             (matlib_real*)u.elem_p, 
                             incu);
                cblas_dgemv( order_enum, 
                             CblasNoTrans, 
                             IM.lenc, 
                             IM.lenr, 
                             1.0, 
                             IM.elem_p, 
                             strideIM, 
                             (matlib_real*)U.elem_p+1, 
                             incU, 
                             0, 
                             (matlib_real*)u.elem_p+1, 
                             incu);
                (U.elem_p) += IM.lenr;
                (u.elem_p) += P;
            }
            
            for (j=1; j<u.lenr; j++)
            {
                (u.elem_p)++;
                for(i=0; i<N; i++)
                {
                    cblas_dgemv( order_enum, 
                                 CblasNoTrans, 
                                 IM.lenc, 
                                 IM.lenr, 
                                 1.0, 
                                 IM.elem_p, 
                                 strideIM, 
                                 (matlib_real*)U.elem_p, 
                                 incU, 
                                 0, 
                                 (matlib_real*)u.elem_p, 
                                 incu);
                    cblas_dgemv( order_enum, 
                                 CblasNoTrans, 
                                 IM.lenc, 
                                 IM.lenr, 
                                 1.0, 
                                 IM.elem_p, 
                                 strideIM, 
                                 (matlib_real*)U.elem_p+1, 
                                 incU, 
                                 0, 
                                 (matlib_real*)u.elem_p+1, 
                                 incu);

                    (U.elem_p) += IM.lenr;
                    (u.elem_p) += P;
                }
            }
        }
        else
        {
            term_execb( "length of matrices incorrect: matrix "
                        "IM: %d-by-%d, U: %d-by-%d, u: %d-by-%d",
                        IM.lenc, IM.lenr, 
                        U.lenc, U.lenr, 
                        u.lenc, u.lenr );
        }
    }
    else
    {
        term_exec( "highest degree of polynomials incorrect: %d",
                    IM.lenr );
    }
    debug_exit("%s", "");

}
/*============================================================================+/
 | Transformation from Legendre basis to FEM-basis and vice-versa
/+============================================================================*/
 
void fem1d_XF2L
(
    const matlib_index p, 
    const matlib_xv    vb,
          matlib_xv    u
)
/* 
 * p : highest degree of polynomials involved, deducebale from length of vectors
 *     but provide here to make debugging easy
 *
 * */ 
{
    debug_enter( "highest polynomial degree: %d "
                 "length of vectors vb: %d, u: %d", 
                 p, vb.len, u.len );

    /* vb.len = N*p+1, u.len = N*(p+1)
     * vb.len-1+N = u.len
     * */ 
    matlib_index N = (vb.len-1)/p;
    debug_body( "nr finite elements: %d ", N);

    assert((vb.elem_p != NULL) && (u.elem_p != NULL));
    
    void (*fp[9])( matlib_index, 
                   matlib_real*, 
                   matlib_real*, 
                   matlib_real*) = { fem1d_xshapefunc2lp_2,
                                     fem1d_xshapefunc2lp_3, 
                                     fem1d_xshapefunc2lp_4, 
                                     fem1d_xshapefunc2lp_5, 
                                     fem1d_xshapefunc2lp_6, 
                                     fem1d_xshapefunc2lp_7, 
                                     fem1d_xshapefunc2lp_8, 
                                     fem1d_xshapefunc2lp_9, 
                                     fem1d_xshapefunc2lp_10};


    if(u.len == vb.len+(N-1))
    {
        if(p>10)
        {
            fem1d_xshapefunc2lp( p, N, vb.elem_p, vb.elem_p+N+1, u.elem_p);
        }
        else
        {
            matlib_index func_index = p-2;
            (*fp[func_index])(N, vb.elem_p, vb.elem_p+N+1, u.elem_p);
        }
    }
    else
    {
        term_execb( "length of vectors incorrect: "
                    "vb: %d, u: %d",
                    vb.len, u.len);
    }

    debug_exit("%s", "");
}

void fem1d_ZF2L
(
    const matlib_index p, 
    const matlib_zv    vb,
          matlib_zv    u
)
/* 
 * p : highest degree of polynomials involved, deducebale from length of vectors
 *     but provide here to make debugging easy
 *
 * */ 
{
    debug_enter( "highest polynomial degree: %d "
                 "length of vectors vb: %d, u: %d", 
                 p, vb.len, u.len );

    matlib_index N = (vb.len-1)/p;
    debug_body( "nr finite elements: %d ", N);

    assert((vb.elem_p != NULL) && (u.elem_p != NULL));
    
    void (*fp[9])( matlib_index, 
                   matlib_complex*, 
                   matlib_complex*, 
                   matlib_complex*) = { fem1d_zshapefunc2lp_2,
                                        fem1d_zshapefunc2lp_3, 
                                        fem1d_zshapefunc2lp_4, 
                                        fem1d_zshapefunc2lp_5, 
                                        fem1d_zshapefunc2lp_6, 
                                        fem1d_zshapefunc2lp_7, 
                                        fem1d_zshapefunc2lp_8, 
                                        fem1d_zshapefunc2lp_9, 
                                        fem1d_zshapefunc2lp_10};


    if(u.len == vb.len+N-1)
    {
        if(p>10)
        {
            fem1d_zshapefunc2lp( p, N, vb.elem_p, vb.elem_p+N+1, u.elem_p);
        }
        else
        {
            matlib_index func_index = p-2;
            (*fp[func_index])(N, vb.elem_p, vb.elem_p+N+1, u.elem_p);
        }
    }
    else
    {
        term_execb( "length of vectors incorrect: "
                    "vb: %d, u: %d",
                    vb.len, u.len);
    }

    debug_exit("%s", "");
}

/*============================================================================*/

void fem1d_XL2F
(
    const matlib_index p, 
    const matlib_xv    u,
          matlib_xv    vb
)
/* 
 * p : highest degree of polynomials involved, deducebale from length of vectors
 *     but provide here to make debugging easy
 *
 * */ 
{
    debug_enter( "highest polynomial degree: %d "
                 "length of vectors vb: %d, u: %d", 
                 p, vb.len, u.len );

    matlib_index N = u.len/(p+1);
    debug_body( "nr finite elements: %d ", N);

    assert((vb.elem_p != NULL) && (u.elem_p != NULL));
    
    void (*fp[9])( matlib_index, 
                   matlib_real*, 
                   matlib_real*, 
                   matlib_real* ) = { lp2fem1d_xshapefunc_2,
                                      lp2fem1d_xshapefunc_3, 
                                      lp2fem1d_xshapefunc_4, 
                                      lp2fem1d_xshapefunc_5, 
                                      lp2fem1d_xshapefunc_6, 
                                      lp2fem1d_xshapefunc_7, 
                                      lp2fem1d_xshapefunc_8, 
                                      lp2fem1d_xshapefunc_9, 
                                      lp2fem1d_xshapefunc_10};


    if(u.len == vb.len+(N-1))
    {
        if(p>10)
        {
            lp2fem1d_xshapefunc( p, N, u.elem_p, vb.elem_p, vb.elem_p+N+1);
        }
        else
        {
            matlib_index func_index = p-2;
            (*fp[func_index])(N, u.elem_p, vb.elem_p, vb.elem_p+N+1);
        }
    }
    else
    {
        term_execb( "length of vectors incorrect: "
                    "vb: %d, u: %d",
                    vb.len, u.len);
    }

    debug_exit("%s", "");
}

void fem1d_ZL2F
(
    const matlib_index p, 
    const matlib_zv    u,
          matlib_zv    vb
)
/* 
 * p : highest degree of polynomials involved, deducebale from length of vectors
 *     but provide here to make debugging easy
 *
 * */ 
{
    debug_enter( "highest polynomial degree: %d "
                 "length of vectors vb: %d, u: %d", 
                 p, vb.len, u.len );

    matlib_index N = u.len/(p+1);
    debug_body( "nr finite elements: %d ", N);

    assert((vb.elem_p != NULL) && (u.elem_p != NULL));
    
    void (*fp[9])( matlib_index, 
                   matlib_complex*, 
                   matlib_complex*, 
                   matlib_complex*) = { lp2fem1d_zshapefunc_2,
                                        lp2fem1d_zshapefunc_3, 
                                        lp2fem1d_zshapefunc_4, 
                                        lp2fem1d_zshapefunc_5, 
                                        lp2fem1d_zshapefunc_6, 
                                        lp2fem1d_zshapefunc_7, 
                                        lp2fem1d_zshapefunc_8, 
                                        lp2fem1d_zshapefunc_9, 
                                        lp2fem1d_zshapefunc_10};


    if(u.len == vb.len+(N-1))
    {
        if(p>10)
        {
            lp2fem1d_zshapefunc( p, N, u.elem_p, vb.elem_p, vb.elem_p+N+1);
        }
        else
        {
            matlib_index func_index = p-2;
            (*fp[func_index])(N, u.elem_p, vb.elem_p, vb.elem_p+N+1);
        }
    }
    else
    {
        term_execb( "length of vectors incorrect: "
                    "vb: %d, u: %d",
                    vb.len, u.len);
    }

    debug_exit("%s", "");
}

/*============================================================================*/


void fem1d_xshapefunc2lp
(
    matlib_index p, 
    matlib_index N, 
    matlib_real *v,               
    matlib_real *b,               
    matlib_real *u                
)
/* (N+1)-by-1 vector, vertex function basis */
/* (p-1)-by-1 vector, bubble function basis      */
/* (p+1)*N-by-1 vector in FEM-basis         */

{
    /* p must be greater than 3 */ 
    matlib_index i,j;

    matlib_real A[4];
    matlib_real B[p-3], C[p-3], *pB, *pC;

    A[0] = -1.0/sqrt(6);
    A[1] = -1.0/sqrt(10);
    A[2] =  1.0/sqrt(2*(2*p-3));
    A[3] =  1.0/sqrt(2*(2*p-1));

    for 
    (
        pB = &B[0], pC = &C[0], j=0;
        j<(p-3);
        pB++, pC++, j++
    )
    {
        *pB =  1.0/sqrt(2*(2*j+3)); 
        *pC = -1.0/sqrt(2*(2*j+7)); 
    }

    for (i=0; i<N; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v)+A[0]**b    ; u++;
        *u = 0.5*(-*(v-1) + *v)+A[1]**(b+1); u++;

        for 
        (
            pB = &B[0], pC = &C[0];
            (pB < &B[p-3]);
            pB++, pC++, b++, u++
        )
        {
            *u = *pB**b + *pC**(b+2);
        }
        *u = A[2]**b, u++, b++; 
        *u = A[3]**b, u++, b++; 
    }
}
void fem1d_zshapefunc2lp
(
    matlib_index    p, 
    matlib_index    N, 
    matlib_complex* v,               
    matlib_complex* b,               
    matlib_complex* u                
)
/* (N+1)-by-1 vector, vertex function basis */
/* (p-1)-by-1 vector, bubble function basis      */
/* (p+1)*N-by-1 vector in FEM-basis         */

{
    /* p must be greater than 3 */ 
    matlib_index i,j;

    matlib_real A[4];
    matlib_real B[p-3], C[p-3], *pB, *pC;

    A[0] = -1.0/sqrt(6);
    A[1] = -1.0/sqrt(10);
    A[2] =  1.0/sqrt(2*(2*p-3));
    A[3] =  1.0/sqrt(2*(2*p-1));

    for 
    (
        pB = &B[0], pC = &C[0], j=0;
        j<(p-3);
        pB++, pC++, j++
    )
    {
        *pB =  1.0/sqrt(2*(2*j+3)); 
        *pC = -1.0/sqrt(2*(2*j+7)); 
    }

    for (i=0; i<N; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v)+A[0]**b    ; u++;
        *u = 0.5*(-*(v-1) + *v)+A[1]**(b+1); u++;

        for 
        (
            pB = &B[0], pC = &C[0];
            (pB < &B[p-3]);
            pB++, pC++, b++, u++
        )
        {
            *u = *pB**b + *pC**(b+2);
        }
        *u = A[2]**b, u++, b++; 
        *u = A[3]**b, u++, b++; 
    }
}
/*============================================================================*/
void lp2fem1d_xshapefunc
(
    matlib_index p,
    matlib_index N,
    matlib_real *u,               /* (p+1)*N-by-1 vector in FEM-basis         */
    matlib_real *v,               /* (N+1)-by-1 vector, vertex function basis */
    matlib_real *b                /* (p-1)*N-by-1 vector, bubble func basis   */
)
{
    /* p must be greater than 3 */ 
    matlib_index i, j;
    matlib_real A[5], B[p-3], C[p-3];
    matlib_real *pA = &A[0], *pB = &B[0], *pC = &C[0];

    *pA = sqrt(2*(2*p-1));   pA++; /* A[0] */
    *pA = sqrt(2*(2*p-3));   pA++; /* A[1] */
    *pA = 1.0/sqrt(6);       pA++; /* A[2] */
    *pA = 1.0/sqrt(10);      pA++; /* A[3] */
    *pA = 2.0**(pA-1);             /* A[4] */

    pA = &A[0];
    /* Start filling the vector from the last fem-element */
    b  = b + (p-1)*N -1;
    u  = u + (p+1)*N -1;
    v  = v + N;

    *b = *pA**u;      /* A[0], b-->(p-2)-th elem, u-->p-th elem     */ 
    pA++; u--; b--;
    *b = *pA**u;      /* A[1], b-->(p-3)-th elem, u-->(p-1)-th elem */
    pA++; u--;

    /* 
     * fill the coefficients to be used later
     * loop p-3 times 
     *
     * */ 
    for 
    (
        j=2, pB = &B[0], pC = &C[0]; 
        j<p-1; 
        j++, pC++, pB++, u--
    )
    {   
        *pB = sqrt(2*(2*(p-j)-1));
        *pC = *pB/sqrt(2*(2*(p-j)+3)); 
        b--;
        *b  = *pB**u + *pC**(b+2);
    }
    /* 
     * u  --> 1st elem, 
     * b  --> 0th elem, 
     * v  --> N
     * pA --> A[2] 
     *
     * */
    *v = *pA**b + *(pA+1)**(b+1) + *u + *(u-1), v--; 
    *v = *pA**b - *(pA+1)**(b+1) - *u + *(u-1);

    for (i=1; i<N; i++)
    {
        u -= 2;
        b--; 
        
        pA = &A[0];
        *b = *pA**u;      /* A[0], b-->(p-2)-th elem, u-->p-th elem     */ 
        pA++; u--; b--; 
        *b = *pA**u;      /* A[1], b-->(p-3)-th elem, u-->(p-1)-th elem */
        pA += 3; u--;

        /* loop p-3 times */ 
        for 
        (
            j=2, pB = &B[0], pC = &C[0]; 
            pB < &B[p-3];
            pC++, pB++, u--
        )
        {   
            b--;
            *b = *pB**u + *pC**(b+2);
        }
        /* u --> 1st elem, b --> 0th elem, */
        v--; 
        *v = *(v+1) - *pA**(b+1) - 2.0**u;
    }     
}
void lp2fem1d_zshapefunc
(
    matlib_index p,
    matlib_index N,
    matlib_complex *u,               
    matlib_complex *v,               
    matlib_complex *b                
)
{
    /* p must be greater than 3 */ 
    matlib_index i, j;
    matlib_real A[5], B[p-3], C[p-3];
    matlib_real *pA = &A[0], *pB = &B[0], *pC = &C[0];

    *pA = sqrt(2*(2*p-1)); pA++; /* A[0] */
    *pA = sqrt(2*(2*p-3)); pA++; /* A[1] */
    *pA = 1.0/sqrt(6);       pA++; /* A[2] */
    *pA = 1.0/sqrt(10);      pA++; /* A[3] */
    *pA = 2.0**(pA-1);             /* A[4] */

    pA = &A[0];
    /* Start filling the vector from the last fem-element */
    b  = b + (p-1)*N -1;
    u  = u + (p+1)*N -1;
    v  = v + N;

    *b = *pA**u;      /* A[0], b-->(p-2)-th elem, u-->p-th elem     */ 
    pA++; u--; b--;
    *b = *pA**u;      /* A[1], b-->(p-3)-th elem, u-->(p-1)-th elem */
    pA++; u--;

    /* 
     * fill the coefficients to be used later
     * loop p-3 times 
     *
     * */ 
    for 
    (
        j=2, pB = &B[0], pC = &C[0]; 
        j<p-1; 
        j++, pC++, pB++, u--
    )
    {   
        *pB = sqrt(2*(2*(p-j)-1));
        *pC = *pB/sqrt(2*(2*(p-j)+3)); 
        b--;
        *b  = *pB**u + *pC**(b+2);
    }
    /* 
     * u  --> 1st elem, 
     * b  --> 0th elem, 
     * v  --> N-1  
     * pA --> A[2] 
     *
     * */
    *v = *pA**b + *(pA+1)**(b+1) + *u + *(u-1), v--; 
    *v = *pA**b - *(pA+1)**(b+1) - *u + *(u-1);

    for (i=1; i<N; i++)
    {
        u -= 2;
        b--; 
        
        pA = &A[0];
        *b = *pA**u;      /* A[0], b-->(p-2)-th elem, u-->p-th elem     */ 
        pA++; u--; b--; 
        *b = *pA**u;      /* A[1], b-->(p-3)-th elem, u-->(p-1)-th elem */
        pA += 3; u--;

        /* loop p-3 times */ 
        for 
        (
            j=2, pB = &B[0], pC = &C[0]; 
            pB < &B[p-3];
            pC++, pB++, u--
        )
        {   
            b--;
            *b = *pB**u + *pC**(b+2);
        }
        /* u --> 1st elem, b --> 0th elem, */
        v--; 
        *v = *(v+1) - *pA**(b+1) - 2**u;
    }     
}

/* Unrolled version of the conversion routines for given degree         */ 
/*======================================================================*/

 

void fem1d_xshapefunc2lp_2
(
    matlib_index N, 
    matlib_real *v,
    matlib_real *b,
    matlib_real *u
)
{

    matlib_index i;
    for (i=0; i<N; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b; u++;
        *u = 0.5*(-*(v-1) + *v)          ; u++;

        *u = _A03**b, u++, b++;
    }
}
 

 

void fem1d_xshapefunc2lp_3
(
    matlib_index N, 
    matlib_real *v,
    matlib_real *b,
    matlib_real *u
)
{

    matlib_index i;
    for (i=0; i<N; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;

        *u = _A03**b, u++, b++;
        *u = _A04**b, u++, b++;
    }
}
 

 

void fem1d_xshapefunc2lp_4
(
    matlib_index N, 
    matlib_real *v,
    matlib_real *b,
    matlib_real *u
)
{

    matlib_index i;
    for (i=0; i<N; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;

        *u = _A04**b, u++, b++;
        *u = _A05**b, u++, b++;
    }
}
 

 

void fem1d_xshapefunc2lp_5
(
    matlib_index N, 
    matlib_real *v,
    matlib_real *b,
    matlib_real *u
)
{

    matlib_index i;
    for (i=0; i<N; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;

        *u = _A05**b, u++, b++;
        *u = _A06**b, u++, b++;
    }
}
 

 

void fem1d_xshapefunc2lp_6
(
    matlib_index N, 
    matlib_real *v,
    matlib_real *b,
    matlib_real *u
)
{

    matlib_index i;
    for (i=0; i<N; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;
        *u = _A05**b + _B05**(b+2), u++, b++;

        *u = _A06**b, u++, b++;
        *u = _A07**b, u++, b++;
    }
}
 

 

void fem1d_xshapefunc2lp_7
(
    matlib_index N, 
    matlib_real *v,
    matlib_real *b,
    matlib_real *u
)
{

    matlib_index i;
    for (i=0; i<N; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;
        *u = _A05**b + _B05**(b+2), u++, b++;
        *u = _A06**b + _B06**(b+2), u++, b++;

        *u = _A07**b, u++, b++;
        *u = _A08**b, u++, b++;
    }
}
 

 

void fem1d_xshapefunc2lp_8
(
    matlib_index N, 
    matlib_real *v,
    matlib_real *b,
    matlib_real *u
)
{

    matlib_index i;
    for (i=0; i<N; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;
        *u = _A05**b + _B05**(b+2), u++, b++;
        *u = _A06**b + _B06**(b+2), u++, b++;
        *u = _A07**b + _B07**(b+2), u++, b++;

        *u = _A08**b, u++, b++;
        *u = _A09**b, u++, b++;
    }
}
 

 

void fem1d_xshapefunc2lp_9
(
    matlib_index N, 
    matlib_real *v,
    matlib_real *b,
    matlib_real *u
)
{

    matlib_index i;
    for (i=0; i<N; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;
        *u = _A05**b + _B05**(b+2), u++, b++;
        *u = _A06**b + _B06**(b+2), u++, b++;
        *u = _A07**b + _B07**(b+2), u++, b++;
        *u = _A08**b + _B08**(b+2), u++, b++;

        *u = _A09**b, u++, b++;
        *u = _A10**b, u++, b++;
    }
}
 

 

void fem1d_xshapefunc2lp_10
(
    matlib_index N, 
    matlib_real *v,
    matlib_real *b,
    matlib_real *u
)
{

    matlib_index i;
    for (i=0; i<N; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;
        *u = _A05**b + _B05**(b+2), u++, b++;
        *u = _A06**b + _B06**(b+2), u++, b++;
        *u = _A07**b + _B07**(b+2), u++, b++;
        *u = _A08**b + _B08**(b+2), u++, b++;
        *u = _A09**b + _B09**(b+2), u++, b++;

        *u = _A10**b, u++, b++;
        *u = _A11**b, u++, b++;
    }
}
 

/* COMPLEX VERSION */ 

 

void fem1d_zshapefunc2lp_2
(
    matlib_index N, 
    matlib_complex *v,
    matlib_complex *b,
    matlib_complex *u
)
{

    matlib_index i;
    for (i=0; i<N; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b; u++;
        *u = 0.5*(-*(v-1) + *v)          ; u++;

        *u = _A03**b, u++, b++;
    }
}
 

 

void fem1d_zshapefunc2lp_3
(
    matlib_index N, 
    matlib_complex *v,
    matlib_complex *b,
    matlib_complex *u
)
{

    matlib_index i;
    for (i=0; i<N; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;

        *u = _A03**b, u++, b++;
        *u = _A04**b, u++, b++;
    }
}
 

 

void fem1d_zshapefunc2lp_4
(
    matlib_index N, 
    matlib_complex *v,
    matlib_complex *b,
    matlib_complex *u
)
{

    matlib_index i;
    for (i=0; i<N; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;

        *u = _A04**b, u++, b++;
        *u = _A05**b, u++, b++;
    }
}
 

 

void fem1d_zshapefunc2lp_5
(
    matlib_index N, 
    matlib_complex *v,
    matlib_complex *b,
    matlib_complex *u
)
{

    matlib_index i;
    for (i=0; i<N; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;

        *u = _A05**b, u++, b++;
        *u = _A06**b, u++, b++;
    }
}
 

 

void fem1d_zshapefunc2lp_6
(
    matlib_index N, 
    matlib_complex *v,
    matlib_complex *b,
    matlib_complex *u
)
{

    matlib_index i;
    for (i=0; i<N; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;
        *u = _A05**b + _B05**(b+2), u++, b++;

        *u = _A06**b, u++, b++;
        *u = _A07**b, u++, b++;
    }
}
 

 

void fem1d_zshapefunc2lp_7
(
    matlib_index N, 
    matlib_complex *v,
    matlib_complex *b,
    matlib_complex *u
)
{

    matlib_index i;
    for (i=0; i<N; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;
        *u = _A05**b + _B05**(b+2), u++, b++;
        *u = _A06**b + _B06**(b+2), u++, b++;

        *u = _A07**b, u++, b++;
        *u = _A08**b, u++, b++;
    }
}
 

 

void fem1d_zshapefunc2lp_8
(
    matlib_index N, 
    matlib_complex *v,
    matlib_complex *b,
    matlib_complex *u
)
{

    matlib_index i;
    for (i=0; i<N; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;
        *u = _A05**b + _B05**(b+2), u++, b++;
        *u = _A06**b + _B06**(b+2), u++, b++;
        *u = _A07**b + _B07**(b+2), u++, b++;

        *u = _A08**b, u++, b++;
        *u = _A09**b, u++, b++;
    }
}
 

 

void fem1d_zshapefunc2lp_9
(
    matlib_index N, 
    matlib_complex *v,
    matlib_complex *b,
    matlib_complex *u
)
{

    matlib_index i;
    for (i=0; i<N; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;
        *u = _A05**b + _B05**(b+2), u++, b++;
        *u = _A06**b + _B06**(b+2), u++, b++;
        *u = _A07**b + _B07**(b+2), u++, b++;
        *u = _A08**b + _B08**(b+2), u++, b++;

        *u = _A09**b, u++, b++;
        *u = _A10**b, u++, b++;
    }
}
 

 

void fem1d_zshapefunc2lp_10
(
    matlib_index N, 
    matlib_complex *v,
    matlib_complex *b,
    matlib_complex *u
)
{

    matlib_index i;
    for (i=0; i<N; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;
        *u = _A05**b + _B05**(b+2), u++, b++;
        *u = _A06**b + _B06**(b+2), u++, b++;
        *u = _A07**b + _B07**(b+2), u++, b++;
        *u = _A08**b + _B08**(b+2), u++, b++;
        *u = _A09**b + _B09**(b+2), u++, b++;

        *u = _A10**b, u++, b++;
        *u = _A11**b, u++, b++;
    }
}
 

/*======================================================================*/
 

void lp2fem1d_xshapefunc_2
(
    matlib_index N,
    matlib_real *u,   
    matlib_real *v,   
    matlib_real *b    
)
{
    matlib_index i;
    b  = b + 1*N -1;
    u  = u + 3*N -1;
    v  = v + N;

    *b = _C02**u; u--;
    *v = _D00**b + *u + *(u-1), v--;
    *v = _D00**b - *u + *(u-1); v--;

    for (i=1; i<N; i++)
    {
        u -=2; b--;
        *b = _C02**u; u--;
        *v = -*(v+1) + 2**(u-1) + 2*_D00**b, v--;
    }
}
 
 

void lp2fem1d_xshapefunc_3
(
    matlib_index N,
    matlib_real *u,   
    matlib_real *v,   
    matlib_real *b    
)
{
    matlib_index i;
    b  = b + 2*N -1;
    u  = u + 4*N -1;
    v  = v + N;

    *b = _C03**u; u--; b--;
    *b = _C02**u; u--; b--;
    *v = _D00**(b+1) + _D01**(b+2) + *u + *(u-1), v--;
    *v = _D00**(b+1) - _D01**(b+2) - *u + *(u-1); v--;

    for (i=1; i<N; i++)
    {

        u -= 2;

        *b = _C03**u; u--; b--;
        *b = _C02**u; u--; b--;
        *v = *(v+1) - 2*_D01**(b+2) - 2**u, v--;
    }
}
 

 

void lp2fem1d_xshapefunc_4
(
    matlib_index N,
    matlib_real *u,   
    matlib_real *v,   
    matlib_real *b    
)
{
    matlib_index i;
    b  = b + 3*N -1;
    u  = u + 5*N -1;
    v  = v + N;

    *b = _C04**u; u--; b--;
    *b = _C03**u; u--; b--;
    *b = _C02**u + _D02**(b+2); u--; b--;
    *v = _D00**(b+1) + _D01**(b+2) + *u + *(u-1), v--;
    *v = _D00**(b+1) - _D01**(b+2) - *u + *(u-1); v--;

    for (i=1; i<N; i++)
    {

        u -= 2;

        *b = _C04**u; u--; b--;
        *b = _C03**u; u--; b--;
        *b = _C02**u + _D02**(b+2); u--; b--;
        *v = *(v+1) - 2*_D01**(b+2) - 2**u, v--;
    }
}
 

 

void lp2fem1d_xshapefunc_5
(
    matlib_index N,
    matlib_real *u,   
    matlib_real *v,   
    matlib_real *b    
)
{
    matlib_index i;
    b  = b + 4*N -1;
    u  = u + 6*N -1;
    v  = v + N;

    *b = _C05**u; u--; b--;
    *b = _C04**u; u--; b--;
    *b = _C03**u + _D03**(b+2); u--; b--;
    *b = _C02**u + _D02**(b+2); u--; b--;
    *v = _D00**(b+1) + _D01**(b+2) + *u + *(u-1), v--;
    *v = _D00**(b+1) - _D01**(b+2) - *u + *(u-1); v--;

    for (i=1; i<N; i++)
    {

        u -= 2;

        *b = _C05**u; u--; b--;
        *b = _C04**u; u--; b--;
        *b = _C03**u + _D03**(b+2); u--; b--;
        *b = _C02**u + _D02**(b+2); u--; b--;
        *v = *(v+1) - 2*_D01**(b+2) - 2**u, v--;
    }
}
 

 

void lp2fem1d_xshapefunc_6
(
    matlib_index N,
    matlib_real *u,   
    matlib_real *v,   
    matlib_real *b    
)
{
    matlib_index i;
    b  = b + 5*N -1;
    u  = u + 7*N -1;
    v  = v + N;

    *b = _C06**u; u--; b--;
    *b = _C05**u; u--; b--;
    *b = _C04**u + _D04**(b+2); u--; b--;
    *b = _C03**u + _D03**(b+2); u--; b--;
    *b = _C02**u + _D02**(b+2); u--; b--;
    *v = _D00**(b+1) + _D01**(b+2) + *u + *(u-1), v--;
    *v = _D00**(b+1) - _D01**(b+2) - *u + *(u-1); v--;

    for (i=1; i<N; i++)
    {

        u -= 2;

        *b = _C06**u; u--; b--;
        *b = _C05**u; u--; b--;
        *b = _C04**u + _D04**(b+2); u--; b--;
        *b = _C03**u + _D03**(b+2); u--; b--;
        *b = _C02**u + _D02**(b+2); u--; b--;
        *v = *(v+1) - 2*_D01**(b+2) - 2**u, v--;
    }
}
 

 

void lp2fem1d_xshapefunc_7
(
    matlib_index N,
    matlib_real *u,   
    matlib_real *v,   
    matlib_real *b    
)
{
    matlib_index i;
    b  = b + 6*N -1;
    u  = u + 8*N -1;
    v  = v + N;

    *b = _C07**u; u--; b--;
    *b = _C06**u; u--; b--;
    *b = _C05**u + _D05**(b+2); u--; b--;
    *b = _C04**u + _D04**(b+2); u--; b--;
    *b = _C03**u + _D03**(b+2); u--; b--;
    *b = _C02**u + _D02**(b+2); u--; b--;
    *v = _D00**(b+1) + _D01**(b+2) + *u + *(u-1), v--;
    *v = _D00**(b+1) - _D01**(b+2) - *u + *(u-1); v--;

    for (i=1; i<N; i++)
    {

        u -= 2;

        *b = _C07**u; u--; b--;
        *b = _C06**u; u--; b--;
        *b = _C05**u + _D05**(b+2); u--; b--;
        *b = _C04**u + _D04**(b+2); u--; b--;
        *b = _C03**u + _D03**(b+2); u--; b--;
        *b = _C02**u + _D02**(b+2); u--; b--;
        *v = *(v+1) - 2*_D01**(b+2) - 2**u, v--;
    }
}
 

 

void lp2fem1d_xshapefunc_8
(
    matlib_index N,
    matlib_real *u,   
    matlib_real *v,   
    matlib_real *b    
)
{
    matlib_index i;
    b  = b + 7*N -1;
    u  = u + 9*N -1;
    v  = v + N;

    *b = _C08**u; u--; b--;
    *b = _C07**u; u--; b--;
    *b = _C06**u + _D06**(b+2); u--; b--;
    *b = _C05**u + _D05**(b+2); u--; b--;
    *b = _C04**u + _D04**(b+2); u--; b--;
    *b = _C03**u + _D03**(b+2); u--; b--;
    *b = _C02**u + _D02**(b+2); u--; b--;
    *v = _D00**(b+1) + _D01**(b+2) + *u + *(u-1), v--;
    *v = _D00**(b+1) - _D01**(b+2) - *u + *(u-1); v--;

    for (i=1; i<N; i++)
    {

        u -= 2;

        *b = _C08**u; u--; b--;
        *b = _C07**u; u--; b--;
        *b = _C06**u + _D06**(b+2); u--; b--;
        *b = _C05**u + _D05**(b+2); u--; b--;
        *b = _C04**u + _D04**(b+2); u--; b--;
        *b = _C03**u + _D03**(b+2); u--; b--;
        *b = _C02**u + _D02**(b+2); u--; b--;
        *v = *(v+1) - 2*_D01**(b+2) - 2**u, v--;
    }
}
 

 

void lp2fem1d_xshapefunc_9
(
    matlib_index N,
    matlib_real *u,   
    matlib_real *v,   
    matlib_real *b    
)
{
    matlib_index i;
    b  = b + 8*N -1;
    u  = u + 10*N -1;
    v  = v + N;

    *b = _C09**u; u--; b--;
    *b = _C08**u; u--; b--;
    *b = _C07**u + _D07**(b+2); u--; b--;
    *b = _C06**u + _D06**(b+2); u--; b--;
    *b = _C05**u + _D05**(b+2); u--; b--;
    *b = _C04**u + _D04**(b+2); u--; b--;
    *b = _C03**u + _D03**(b+2); u--; b--;
    *b = _C02**u + _D02**(b+2); u--; b--;
    *v = _D00**(b+1) + _D01**(b+2) + *u + *(u-1), v--;
    *v = _D00**(b+1) - _D01**(b+2) - *u + *(u-1); v--;

    for (i=1; i<N; i++)
    {

        u -= 2;

        *b = _C09**u; u--; b--;
        *b = _C08**u; u--; b--;
        *b = _C07**u + _D07**(b+2); u--; b--;
        *b = _C06**u + _D06**(b+2); u--; b--;
        *b = _C05**u + _D05**(b+2); u--; b--;
        *b = _C04**u + _D04**(b+2); u--; b--;
        *b = _C03**u + _D03**(b+2); u--; b--;
        *b = _C02**u + _D02**(b+2); u--; b--;
        *v = *(v+1) - 2*_D01**(b+2) - 2**u, v--;
    }
}
 

 

void lp2fem1d_xshapefunc_10
(
    matlib_index N,
    matlib_real *u,   
    matlib_real *v,   
    matlib_real *b    
)
{
    matlib_index i;
    b  = b + 9*N -1;
    u  = u + 11*N -1;
    v  = v + N;

    *b = _C10**u; u--; b--;
    *b = _C09**u; u--; b--;
    *b = _C08**u + _D08**(b+2); u--; b--;
    *b = _C07**u + _D07**(b+2); u--; b--;
    *b = _C06**u + _D06**(b+2); u--; b--;
    *b = _C05**u + _D05**(b+2); u--; b--;
    *b = _C04**u + _D04**(b+2); u--; b--;
    *b = _C03**u + _D03**(b+2); u--; b--;
    *b = _C02**u + _D02**(b+2); u--; b--;
    *v = _D00**(b+1) + _D01**(b+2) + *u + *(u-1), v--;
    *v = _D00**(b+1) - _D01**(b+2) - *u + *(u-1); v--;

    for (i=1; i<N; i++)
    {

        u -= 2;

        *b = _C10**u; u--; b--;
        *b = _C09**u; u--; b--;
        *b = _C08**u + _D08**(b+2); u--; b--;
        *b = _C07**u + _D07**(b+2); u--; b--;
        *b = _C06**u + _D06**(b+2); u--; b--;
        *b = _C05**u + _D05**(b+2); u--; b--;
        *b = _C04**u + _D04**(b+2); u--; b--;
        *b = _C03**u + _D03**(b+2); u--; b--;
        *b = _C02**u + _D02**(b+2); u--; b--;
        *v = *(v+1) - 2*_D01**(b+2) - 2**u, v--;
    }
}
 

/* COMPLEX VERSION */ 


 

void lp2fem1d_zshapefunc_2
(
    matlib_index   N,
    matlib_complex *u,   
    matlib_complex *v,   
    matlib_complex *b    
)
{
    matlib_index i;
    b  = b + 1*N -1;
    u  = u + 3*N -1;
    v  = v + N;

    *b = _C02**u; u--;
    *v = _D00**b + *u + *(u-1), v--;
    *v = _D00**b - *u + *(u-1); v--;

    for (i=1; i<N; i++)
    {
        u -=2; b--;
        *b = _C02**u; u--;
        *v = -*(v+1) + 2**(u-1) + 2*_D00**b, v--;
    }
}
 
 

void lp2fem1d_zshapefunc_3
(
    matlib_index   N,
    matlib_complex *u,   
    matlib_complex *v,   
    matlib_complex *b    
)
{
    matlib_index i;
    b  = b + 2*N -1;
    u  = u + 4*N -1;
    v  = v + N;

    *b = _C03**u; u--; b--;
    *b = _C02**u; u--; b--;
    *v = _D00**(b+1) + _D01**(b+2) + *u + *(u-1), v--;
    *v = _D00**(b+1) - _D01**(b+2) - *u + *(u-1); v--;

    for (i=1; i<N; i++)
    {

        u -= 2;

        *b = _C03**u; u--; b--;
        *b = _C02**u; u--; b--;
        *v = *(v+1) - 2*_D01**(b+2) - 2**u, v--;
    }
}
 

 

void lp2fem1d_zshapefunc_4
(
    matlib_index   N,
    matlib_complex *u,   
    matlib_complex *v,   
    matlib_complex *b    
)
{
    matlib_index i;
    b  = b + 3*N -1;
    u  = u + 5*N -1;
    v  = v + N;

    *b = _C04**u; u--; b--;
    *b = _C03**u; u--; b--;
    *b = _C02**u + _D02**(b+2); u--; b--;
    *v = _D00**(b+1) + _D01**(b+2) + *u + *(u-1), v--;
    *v = _D00**(b+1) - _D01**(b+2) - *u + *(u-1); v--;

    for (i=1; i<N; i++)
    {

        u -= 2;

        *b = _C04**u; u--; b--;
        *b = _C03**u; u--; b--;
        *b = _C02**u + _D02**(b+2); u--; b--;
        *v = *(v+1) - 2*_D01**(b+2) - 2**u, v--;
    }
}
 

 

void lp2fem1d_zshapefunc_5
(
    matlib_index   N,
    matlib_complex *u,   
    matlib_complex *v,   
    matlib_complex *b    
)
{
    matlib_index i;
    b  = b + 4*N -1;
    u  = u + 6*N -1;
    v  = v + N;

    *b = _C05**u; u--; b--;
    *b = _C04**u; u--; b--;
    *b = _C03**u + _D03**(b+2); u--; b--;
    *b = _C02**u + _D02**(b+2); u--; b--;
    *v = _D00**(b+1) + _D01**(b+2) + *u + *(u-1), v--;
    *v = _D00**(b+1) - _D01**(b+2) - *u + *(u-1); v--;

    for (i=1; i<N; i++)
    {

        u -= 2;

        *b = _C05**u; u--; b--;
        *b = _C04**u; u--; b--;
        *b = _C03**u + _D03**(b+2); u--; b--;
        *b = _C02**u + _D02**(b+2); u--; b--;
        *v = *(v+1) - 2*_D01**(b+2) - 2**u, v--;
    }
}
 

 

void lp2fem1d_zshapefunc_6
(
    matlib_index   N,
    matlib_complex *u,   
    matlib_complex *v,   
    matlib_complex *b    
)
{
    matlib_index i;
    b  = b + 5*N -1;
    u  = u + 7*N -1;
    v  = v + N;

    *b = _C06**u; u--; b--;
    *b = _C05**u; u--; b--;
    *b = _C04**u + _D04**(b+2); u--; b--;
    *b = _C03**u + _D03**(b+2); u--; b--;
    *b = _C02**u + _D02**(b+2); u--; b--;
    *v = _D00**(b+1) + _D01**(b+2) + *u + *(u-1), v--;
    *v = _D00**(b+1) - _D01**(b+2) - *u + *(u-1); v--;

    for (i=1; i<N; i++)
    {

        u -= 2;

        *b = _C06**u; u--; b--;
        *b = _C05**u; u--; b--;
        *b = _C04**u + _D04**(b+2); u--; b--;
        *b = _C03**u + _D03**(b+2); u--; b--;
        *b = _C02**u + _D02**(b+2); u--; b--;
        *v = *(v+1) - 2*_D01**(b+2) - 2**u, v--;
    }
}
 

 

void lp2fem1d_zshapefunc_7
(
    matlib_index    N,
    matlib_complex *u,   
    matlib_complex *v,   
    matlib_complex *b    
)
{
    matlib_index i;
    b  = b + 6*N -1;
    u  = u + 8*N -1;
    v  = v + N;

    *b = _C07**u; u--; b--;
    *b = _C06**u; u--; b--;
    *b = _C05**u + _D05**(b+2); u--; b--;
    *b = _C04**u + _D04**(b+2); u--; b--;
    *b = _C03**u + _D03**(b+2); u--; b--;
    *b = _C02**u + _D02**(b+2); u--; b--;
    *v = _D00**(b+1) + _D01**(b+2) + *u + *(u-1), v--;
    *v = _D00**(b+1) - _D01**(b+2) - *u + *(u-1); v--;

    for (i=1; i<N; i++)
    {

        u -= 2;

        *b = _C07**u; u--; b--;
        *b = _C06**u; u--; b--;
        *b = _C05**u + _D05**(b+2); u--; b--;
        *b = _C04**u + _D04**(b+2); u--; b--;
        *b = _C03**u + _D03**(b+2); u--; b--;
        *b = _C02**u + _D02**(b+2); u--; b--;
        *v = *(v+1) - 2*_D01**(b+2) - 2**u, v--;
    }
}
 

 

void lp2fem1d_zshapefunc_8
(
    matlib_index   N,
    matlib_complex *u,   
    matlib_complex *v,   
    matlib_complex *b    
)
{
    matlib_index i;
    b  = b + 7*N -1;
    u  = u + 9*N -1;
    v  = v + N;

    *b = _C08**u; u--; b--;
    *b = _C07**u; u--; b--;
    *b = _C06**u + _D06**(b+2); u--; b--;
    *b = _C05**u + _D05**(b+2); u--; b--;
    *b = _C04**u + _D04**(b+2); u--; b--;
    *b = _C03**u + _D03**(b+2); u--; b--;
    *b = _C02**u + _D02**(b+2); u--; b--;
    *v = _D00**(b+1) + _D01**(b+2) + *u + *(u-1), v--;
    *v = _D00**(b+1) - _D01**(b+2) - *u + *(u-1); v--;

    for (i=1; i<N; i++)
    {

        u -= 2;

        *b = _C08**u; u--; b--;
        *b = _C07**u; u--; b--;
        *b = _C06**u + _D06**(b+2); u--; b--;
        *b = _C05**u + _D05**(b+2); u--; b--;
        *b = _C04**u + _D04**(b+2); u--; b--;
        *b = _C03**u + _D03**(b+2); u--; b--;
        *b = _C02**u + _D02**(b+2); u--; b--;
        *v = *(v+1) - 2*_D01**(b+2) - 2**u, v--;
    }
}
 

 

void lp2fem1d_zshapefunc_9
(
    matlib_index   N,
    matlib_complex *u,   
    matlib_complex *v,   
    matlib_complex *b    
)
{
    matlib_index i;
    b  = b + 8*N -1;
    u  = u + 10*N -1;
    v  = v + N;

    *b = _C09**u; u--; b--;
    *b = _C08**u; u--; b--;
    *b = _C07**u + _D07**(b+2); u--; b--;
    *b = _C06**u + _D06**(b+2); u--; b--;
    *b = _C05**u + _D05**(b+2); u--; b--;
    *b = _C04**u + _D04**(b+2); u--; b--;
    *b = _C03**u + _D03**(b+2); u--; b--;
    *b = _C02**u + _D02**(b+2); u--; b--;
    *v = _D00**(b+1) + _D01**(b+2) + *u + *(u-1), v--;
    *v = _D00**(b+1) - _D01**(b+2) - *u + *(u-1); v--;

    for (i=1; i<N; i++)
    {

        u -= 2;

        *b = _C09**u; u--; b--;
        *b = _C08**u; u--; b--;
        *b = _C07**u + _D07**(b+2); u--; b--;
        *b = _C06**u + _D06**(b+2); u--; b--;
        *b = _C05**u + _D05**(b+2); u--; b--;
        *b = _C04**u + _D04**(b+2); u--; b--;
        *b = _C03**u + _D03**(b+2); u--; b--;
        *b = _C02**u + _D02**(b+2); u--; b--;
        *v = *(v+1) - 2*_D01**(b+2) - 2**u, v--;
    }
}
 

 

void lp2fem1d_zshapefunc_10
(
    matlib_index   N,
    matlib_complex *u,   
    matlib_complex *v,   
    matlib_complex *b    
)
{
    matlib_index i;
    b  = b + 9*N -1;
    u  = u + 11*N -1;
    v  = v + N;

    *b = _C10**u; u--; b--;
    *b = _C09**u; u--; b--;
    *b = _C08**u + _D08**(b+2); u--; b--;
    *b = _C07**u + _D07**(b+2); u--; b--;
    *b = _C06**u + _D06**(b+2); u--; b--;
    *b = _C05**u + _D05**(b+2); u--; b--;
    *b = _C04**u + _D04**(b+2); u--; b--;
    *b = _C03**u + _D03**(b+2); u--; b--;
    *b = _C02**u + _D02**(b+2); u--; b--;
    *v = _D00**(b+1) + _D01**(b+2) + *u + *(u-1), v--;
    *v = _D00**(b+1) - _D01**(b+2) - *u + *(u-1); v--;

    for (i=1; i<N; i++)
    {

        u -= 2;

        *b = _C10**u; u--; b--;
        *b = _C09**u; u--; b--;
        *b = _C08**u + _D08**(b+2); u--; b--;
        *b = _C07**u + _D07**(b+2); u--; b--;
        *b = _C06**u + _D06**(b+2); u--; b--;
        *b = _C05**u + _D05**(b+2); u--; b--;
        *b = _C04**u + _D04**(b+2); u--; b--;
        *b = _C03**u + _D03**(b+2); u--; b--;
        *b = _C02**u + _D02**(b+2); u--; b--;
        *v = *(v+1) - 2*_D01**(b+2) - 2**u, v--;
    }
}
 

/*============================================================================+/
 | Projection from LP basis representation to FEM-basis
/+============================================================================*/

void fem1d_XPrjL2F
(
    matlib_index p,
    matlib_xv    u,
    matlib_xv    Pvb
)
{
    debug_enter( "highest polynomial degree: %d "
                 "length of vectors u: %d, Pbv: %d", 
                 p, u.len, Pvb.len );

    matlib_index N = u.len/(p+1);
    debug_body( "nr finite elements: %d ", N);

    assert((u.elem_p != NULL) && (Pvb.elem_p != NULL));

    void (*fp[9])( matlib_index, 
                   matlib_real*, 
                   matlib_real*, 
                   matlib_real*) = { fem1d_xprjLP2FEM_ShapeFunc_2,
                                fem1d_xprjLP2FEM_ShapeFunc_3, 
                                fem1d_xprjLP2FEM_ShapeFunc_4, 
                                fem1d_xprjLP2FEM_ShapeFunc_5, 
                                fem1d_xprjLP2FEM_ShapeFunc_6, 
                                fem1d_xprjLP2FEM_ShapeFunc_7, 
                                fem1d_xprjLP2FEM_ShapeFunc_8, 
                                fem1d_xprjLP2FEM_ShapeFunc_9, 
                                fem1d_xprjLP2FEM_ShapeFunc_10};

    if(u.len == Pvb.len+N-1)
    {
        if(p>10)
        {
            fem1d_xprjLP2FEM_ShapeFunc( p, N, u.elem_p, Pvb.elem_p, Pvb.elem_p+N+1);
        }
        else
        {
            matlib_index func_index = p-2;
            (*fp[func_index])(N, u.elem_p, Pvb.elem_p, Pvb.elem_p+N+1);
        }
    }
    else
    {
        term_execb( "length of vectors incorrect: "
                    "u: %d, Pvb: %d",
                    u.len, Pvb.len);
    }

    debug_exit("%s", "");

}

void fem1d_ZPrjL2F
(
    matlib_index     p,
    matlib_zv u,
    matlib_zv Pvb
)
{
    debug_enter( "highest polynomial degree: %d "
                 "length of vectors u: %d, Pbv: %d", 
                 p, u.len, Pvb.len );

    matlib_index N = u.len/(p+1);
    debug_body( "nr finite elements: %d ", N);

    assert((u.elem_p != NULL) && (Pvb.elem_p != NULL));

    void (*fp[9])( matlib_index, 
                   matlib_complex*, 
                   matlib_complex*, 
                   matlib_complex*) = { fem1d_zprjLP2FEM_ShapeFunc_2,
                                        fem1d_zprjLP2FEM_ShapeFunc_3, 
                                        fem1d_zprjLP2FEM_ShapeFunc_4, 
                                        fem1d_zprjLP2FEM_ShapeFunc_5, 
                                        fem1d_zprjLP2FEM_ShapeFunc_6, 
                                        fem1d_zprjLP2FEM_ShapeFunc_7, 
                                        fem1d_zprjLP2FEM_ShapeFunc_8, 
                                        fem1d_zprjLP2FEM_ShapeFunc_9, 
                                        fem1d_zprjLP2FEM_ShapeFunc_10};

    if(u.len == Pvb.len+N-1)
    {
        if(p>10)
        {
            fem1d_zprjLP2FEM_ShapeFunc( p, N, u.elem_p, Pvb.elem_p, Pvb.elem_p+N+1);
        }
        else
        {
            matlib_index func_index = p-2;
            (*fp[func_index])(N, u.elem_p, Pvb.elem_p, Pvb.elem_p+N+1);
        }
    }
    else
    {
        term_execb( "length of vectors incorrect: "
                    "u: %d, Pvb: %d",
                    u.len, Pvb.len);
    }

    debug_exit("%s", "");

}
/*============================================================================*/

void fem1d_xprjLP2FEM_ShapeFunc
(
    matlib_index p, 
    matlib_index N, 
    matlib_real *u, 
    matlib_real *Pv,                     
    matlib_real *Pb                      
)
{
    matlib_index i;
    matlib_real B[p-1], C[p-1], tmp, *pB, *pC;

    for 
    ( 
        i=0, pB = &B[0], pC = &C[0]; 
        i<p-1; 
        i++, pB++, pC++ 
    )
    {
        tmp =  1.0/sqrt(4*i+6);
        *pB =  tmp/(i+2.5);
        *pC = -tmp/(i+0.5);
    }

    tmp = 0;
    for (i=0; i<N; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        for 
        (
            pB = &B[0], pC = &C[0];
            pB < &B[p-1];
            pB++, pC++, Pb++, u++
        )
        {
            *Pb = *pB**(u+2) + *pC**u;
        }
        u += 2;
    }
    *Pv = tmp;

}

void fem1d_zprjLP2FEM_ShapeFunc
(
    matlib_index    p, 
    matlib_index    N, 
    matlib_complex* u, 
    matlib_complex* Pv,
    matlib_complex* Pb 
)
/* Pv : vector of size (N+1)
 * Pb : vector of size N*(p-1)
 * 
 * */
{
    matlib_index i;
    matlib_real B[p-1], C[p-1], tmp, *pB, *pC;
    matlib_complex tmpc;

    for 
    ( 
        i=0, pB = &B[0], pC = &C[0]; 
        i<p-1; 
        i++, pB++, pC++ 
    )
    {
        tmp =  1.0/sqrt(4*i+6);
        *pB =  tmp/(i+2.5);
        *pC = -tmp/(i+0.5);
    }

    tmpc = 0;
    for (i=0; i<N; i++, Pv++)
    {
        *Pv  = (*u - *(u+1)/3) + tmpc;
        tmpc = (*u + *(u+1)/3);
        for 
        (
            pB = &B[0], pC = &C[0];
            pB < &B[p-1];
            pB++, pC++, Pb++, u++
        )
        {
            *Pb = *pB**(u+2) + *pC**u;
        }
        u += 2;
    }
    *Pv = tmpc;

}
/*============================================================================*/

void fem1d_xprjLP2FEM_ShapeFunc_2
(
    matlib_index N,
    matlib_real *u, 
    matlib_real *Pv,   
    matlib_real *Pb    
)
{
    matlib_index i;
    matlib_real tmp = 0;
    for (i=0; i<N; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        u += 2;

    }
    *Pv = tmp;
}
 

 

void fem1d_xprjLP2FEM_ShapeFunc_3
(
    matlib_index N,
    matlib_real *u, 
    matlib_real *Pv,   
    matlib_real *Pb    
)
{
    matlib_index i;
    matlib_real tmp = 0;
    for (i=0; i<N; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        u += 2;

    }
    *Pv = tmp;
}
 

 

void fem1d_xprjLP2FEM_ShapeFunc_4
(
    matlib_index N,
    matlib_real *u, 
    matlib_real *Pv,   
    matlib_real *Pb    
)
{
    matlib_index i;
    matlib_real tmp = 0;
    for (i=0; i<N; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        u += 2;

    }
    *Pv = tmp;
}
 

 

void fem1d_xprjLP2FEM_ShapeFunc_5
(
    matlib_index N,
    matlib_real *u, 
    matlib_real *Pv,   
    matlib_real *Pb    
)
{
    matlib_index i;
    matlib_real tmp = 0;
    for (i=0; i<N; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        u += 2;

    }
    *Pv = tmp;
}
 

 

void fem1d_xprjLP2FEM_ShapeFunc_6
(
    matlib_index N,
    matlib_real *u, 
    matlib_real *Pv,   
    matlib_real *Pb    
)
{
    matlib_index i;
    matlib_real tmp = 0;
    for (i=0; i<N; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        *Pb = _E04**(u+2) + _F04**u; Pb++; u++;
        u += 2;

    }
    *Pv = tmp;
}
 

 

void fem1d_xprjLP2FEM_ShapeFunc_7
(
    matlib_index N,
    matlib_real *u, 
    matlib_real *Pv,   
    matlib_real *Pb    
)
{
    matlib_index i;
    matlib_real tmp = 0;
    for (i=0; i<N; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        *Pb = _E04**(u+2) + _F04**u; Pb++; u++;
        *Pb = _E05**(u+2) + _F05**u; Pb++; u++;
        u += 2;

    }
    *Pv = tmp;
}
 

 

void fem1d_xprjLP2FEM_ShapeFunc_8
(
    matlib_index N,
    matlib_real *u, 
    matlib_real *Pv,   
    matlib_real *Pb    
)
{
    matlib_index i;
    matlib_real tmp = 0;
    for (i=0; i<N; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        *Pb = _E04**(u+2) + _F04**u; Pb++; u++;
        *Pb = _E05**(u+2) + _F05**u; Pb++; u++;
        *Pb = _E06**(u+2) + _F06**u; Pb++; u++;
        u += 2;

    }
    *Pv = tmp;
}
 

 

void fem1d_xprjLP2FEM_ShapeFunc_9
(
    matlib_index N,
    matlib_real *u, 
    matlib_real *Pv,   
    matlib_real *Pb    
)
{
    matlib_index i;
    matlib_real tmp = 0;
    for (i=0; i<N; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        *Pb = _E04**(u+2) + _F04**u; Pb++; u++;
        *Pb = _E05**(u+2) + _F05**u; Pb++; u++;
        *Pb = _E06**(u+2) + _F06**u; Pb++; u++;
        *Pb = _E07**(u+2) + _F07**u; Pb++; u++;
        u += 2;

    }
    *Pv = tmp;
}
 

 

void fem1d_xprjLP2FEM_ShapeFunc_10
(
    matlib_index N,
    matlib_real *u, 
    matlib_real *Pv,   
    matlib_real *Pb    
)
{
    matlib_index i;
    matlib_real tmp = 0;
    for (i=0; i<N; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        *Pb = _E04**(u+2) + _F04**u; Pb++; u++;
        *Pb = _E05**(u+2) + _F05**u; Pb++; u++;
        *Pb = _E06**(u+2) + _F06**u; Pb++; u++;
        *Pb = _E07**(u+2) + _F07**u; Pb++; u++;
        *Pb = _E08**(u+2) + _F08**u; Pb++; u++;
        u += 2;

    }
    *Pv = tmp;
}
 

/* COMPLEX VERSION */ 

 

void fem1d_zprjLP2FEM_ShapeFunc_2
(
    matlib_index           N,
    matlib_complex* u, 
    matlib_complex* Pv,   
    matlib_complex* Pb    
)
{
    matlib_index i;
    matlib_complex tmp = 0;
    for (i=0; i<N; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        u += 2;

    }
    *Pv = tmp;
}
 

 

void fem1d_zprjLP2FEM_ShapeFunc_3
(
    matlib_index           N,
    matlib_complex* u, 
    matlib_complex* Pv,   
    matlib_complex* Pb    
)
{
    matlib_index i;
    matlib_complex tmp = 0;
    for (i=0; i<N; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        u += 2;

    }
    *Pv = tmp;
}
 

 

void fem1d_zprjLP2FEM_ShapeFunc_4
(
    matlib_index           N,
    matlib_complex* u, 
    matlib_complex* Pv,   
    matlib_complex* Pb    
)
{
    matlib_index i;
    matlib_complex tmp = 0;
    for (i=0; i<N; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        u += 2;

    }
    *Pv = tmp;
}
 

 

void fem1d_zprjLP2FEM_ShapeFunc_5
(
    matlib_index           N,
    matlib_complex* u, 
    matlib_complex* Pv,   
    matlib_complex* Pb    
)
{
    matlib_index i;
    matlib_complex tmp = 0;
    for (i=0; i<N; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        u += 2;

    }
    *Pv = tmp;
}
 

 

void fem1d_zprjLP2FEM_ShapeFunc_6
(
    matlib_index           N,
    matlib_complex* u, 
    matlib_complex* Pv,   
    matlib_complex* Pb    
)
{
    matlib_index i;
    matlib_complex tmp = 0;
    for (i=0; i<N; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        *Pb = _E04**(u+2) + _F04**u; Pb++; u++;
        u += 2;

    }
    *Pv = tmp;
}
 

 

void fem1d_zprjLP2FEM_ShapeFunc_7
(
    matlib_index           N,
    matlib_complex* u, 
    matlib_complex* Pv,   
    matlib_complex* Pb    
)
{
    matlib_index i;
    matlib_complex tmp = 0;
    for (i=0; i<N; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        *Pb = _E04**(u+2) + _F04**u; Pb++; u++;
        *Pb = _E05**(u+2) + _F05**u; Pb++; u++;
        u += 2;

    }
    *Pv = tmp;
}
 

 

void fem1d_zprjLP2FEM_ShapeFunc_8
(
    matlib_index           N,
    matlib_complex* u, 
    matlib_complex* Pv,   
    matlib_complex* Pb    
)
{
    matlib_index i;
    matlib_complex tmp = 0;
    for (i=0; i<N; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        *Pb = _E04**(u+2) + _F04**u; Pb++; u++;
        *Pb = _E05**(u+2) + _F05**u; Pb++; u++;
        *Pb = _E06**(u+2) + _F06**u; Pb++; u++;
        u += 2;

    }
    *Pv = tmp;
}
 

 

void fem1d_zprjLP2FEM_ShapeFunc_9
(
    matlib_index           N,
    matlib_complex* u, 
    matlib_complex* Pv,   
    matlib_complex* Pb    
)
{
    matlib_index i;
    matlib_complex tmp = 0;
    for (i=0; i<N; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        *Pb = _E04**(u+2) + _F04**u; Pb++; u++;
        *Pb = _E05**(u+2) + _F05**u; Pb++; u++;
        *Pb = _E06**(u+2) + _F06**u; Pb++; u++;
        *Pb = _E07**(u+2) + _F07**u; Pb++; u++;
        u += 2;

    }
    *Pv = tmp;
}
 

 

void fem1d_zprjLP2FEM_ShapeFunc_10
(
    matlib_index           N,
    matlib_complex* u, 
    matlib_complex* Pv,   
    matlib_complex* Pb    
)
{
    matlib_index i;
    matlib_complex tmp = 0;
    for (i=0; i<N; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        *Pb = _E04**(u+2) + _F04**u; Pb++; u++;
        *Pb = _E05**(u+2) + _F05**u; Pb++; u++;
        *Pb = _E06**(u+2) + _F06**u; Pb++; u++;
        *Pb = _E07**(u+2) + _F07**u; Pb++; u++;
        *Pb = _E08**(u+2) + _F08**u; Pb++; u++;
        u += 2;

    }
    *Pv = tmp;
}
 



/*============================================================================+/
 | L2 Norm in Legendre Basis
/+============================================================================*/

matlib_real fem1d_XNorm2
(
    matlib_index p,
    matlib_index N,
    matlib_xv    u
)
{
    debug_enter( "Highest degree of polynomial: %d "
                 "nr. fnite-elements : %d "
                 "length of u: %d", p, N, u.len);

    matlib_real (*fp[9])( matlib_index, 
                          matlib_real* ) = { fem1d_xlp_snorm2_d_2, 
                                             fem1d_xlp_snorm2_d_3, 
                                             fem1d_xlp_snorm2_d_4, 
                                             fem1d_xlp_snorm2_d_5, 
                                             fem1d_xlp_snorm2_d_6, 
                                             fem1d_xlp_snorm2_d_7, 
                                             fem1d_xlp_snorm2_d_8, 
                                             fem1d_xlp_snorm2_d_9, 
                                             fem1d_xlp_snorm2_d_10};
 
    matlib_real snorm2;
    assert(u.elem_p!=NULL);
    if(u.len == (p+1)*N)
    {
        if(p>10)
        {
            snorm2 = fem1d_xlp_snorm2_d( p, N, u.elem_p);
        }
        else
        {
            matlib_index func_index = p-2;
            snorm2 = (*fp[func_index])( N, u.elem_p);
        }
    
    }
    else
    {
        term_execb( "length of vector incorrect: "
                    "u: %d (degree : %d, nr. fem-elements: %d)",
                    u.len, p, N);
    }

    debug_exit("L2 norm: %0.16f", sqrt(snorm2));
    return sqrt(snorm2);

}
matlib_real fem1d_ZNorm2
(
    matlib_index     p,
    matlib_index     N,
    matlib_zv u
)
{
    debug_enter( "Highest degree of polynomial: %d "
                 "nr. fnite-elements : %d "
                 "length of u: %d", p, N, u.len);

    matlib_real (*fp[9])(matlib_index, matlib_complex*) = { fem1d_zlp_snorm2_d_2, 
                                                 fem1d_zlp_snorm2_d_3, 
                                                 fem1d_zlp_snorm2_d_4, 
                                                 fem1d_zlp_snorm2_d_5, 
                                                 fem1d_zlp_snorm2_d_6, 
                                                 fem1d_zlp_snorm2_d_7, 
                                                 fem1d_zlp_snorm2_d_8, 
                                                 fem1d_zlp_snorm2_d_9, 
                                                 fem1d_zlp_snorm2_d_10};
 
    matlib_real snorm2;
    assert(u.elem_p!=NULL);
    if(u.len == (p+1)*N)
    {
        if(p>10)
        {
            snorm2 = fem1d_zlp_snorm2_d( p, N, u.elem_p);
        }
        else
        {
            matlib_index func_index = p-2;
            snorm2 = (*fp[func_index])( N, u.elem_p);
        }
    
    }
    else
    {
        term_execb( "length of vector incorrect: "
                    "u: %d (degree : %d, nr. fem-elements: %d)",
                    u.len, p, N);
    }

    debug_exit("L2 norm: %0.16f", sqrt(snorm2));
    return sqrt(snorm2);

}


matlib_real fem1d_xlp_snorm2_d(matlib_index p, matlib_index N, matlib_real *u)
{
    matlib_real snorm = 0;
    matlib_index i;
    matlib_real coeff[p+1], *pcoeff;

    for 
    ( 
        pcoeff = &coeff[0], i=0; 
        i<p+1; 
        i++, pcoeff++
    )
    {
        *pcoeff = 1.0/(i+0.5); 
    }

    for (i=0; i<N; i++)
    { 
        for 
        ( 
            pcoeff = &coeff[0]; 
            pcoeff < &coeff[0]+p+1; 
            u++, pcoeff++
        )
        {
            snorm += *u**u**pcoeff;
        }
    }
    return snorm;
}

matlib_real fem1d_zlp_snorm2_d(matlib_index p, matlib_index N, matlib_complex *u)
{
    matlib_real snorm = 0;
    matlib_index i;
    matlib_real coeff[p+1], *pcoeff;

    for 
    ( 
        pcoeff = &coeff[0], i=0; 
        i<p+1; 
        i++, pcoeff++
    )
    {
        *pcoeff = 1.0/(i+0.5); 
    }

    for (i=0; i<N; i++)
    {
        for 
        ( 
            pcoeff = &coeff[0]; 
            pcoeff < &coeff[0]+p+1; 
            u++, pcoeff++
        )
        {
            snorm += *u*conj(*u)**pcoeff;
        }
    }
    return snorm;
}
/*============================================================================*/

 
matlib_real fem1d_xlp_snorm2_d_2(matlib_index N, matlib_real *u)
{
    matlib_real snorm = 0, tmp;
    matlib_index i;
    for (i=0; i<N; i++)
    {
        tmp = 	  _N0 * *(u+0) * *(u+0)
        	+ _N1 * *(u+1) * *(u+1)
        	+ _N2 * *(u+2) * *(u+2);
        u += 3;
        snorm += tmp;
    }
    return snorm;
}
 
 
matlib_real fem1d_xlp_snorm2_d_3(matlib_index N, matlib_real *u)
{
    matlib_real snorm = 0, tmp;
    matlib_index i;
    for (i=0; i<N; i++)
    {
        tmp = 	  _N0 * *(u+0) * *(u+0)
        	+ _N1 * *(u+1) * *(u+1)
        	+ _N2 * *(u+2) * *(u+2)
        	+ _N3 * *(u+3) * *(u+3);
        u += 4;
        snorm += tmp;
    }
    return snorm;
}
 
 
matlib_real fem1d_xlp_snorm2_d_4(matlib_index N, matlib_real *u)
{
    matlib_real snorm = 0, tmp;
    matlib_index i;
    for (i=0; i<N; i++)
    {
        tmp = 	  _N0 * *(u+0) * *(u+0)
        	+ _N1 * *(u+1) * *(u+1)
        	+ _N2 * *(u+2) * *(u+2)
        	+ _N3 * *(u+3) * *(u+3)
        	+ _N4 * *(u+4) * *(u+4);
        u += 5;
        snorm += tmp;
    }
    return snorm;
}
 
 
matlib_real fem1d_xlp_snorm2_d_5(matlib_index N, matlib_real *u)
{
    matlib_real snorm = 0, tmp;
    matlib_index i;
    for (i=0; i<N; i++)
    {
        tmp = 	  _N0 * *(u+0) * *(u+0)
        	+ _N1 * *(u+1) * *(u+1)
        	+ _N2 * *(u+2) * *(u+2)
        	+ _N3 * *(u+3) * *(u+3)
        	+ _N4 * *(u+4) * *(u+4)
        	+ _N5 * *(u+5) * *(u+5);
        u += 6;
        snorm += tmp;
    }
    return snorm;
}
 
 
matlib_real fem1d_xlp_snorm2_d_6(matlib_index N, matlib_real *u)
{
    matlib_real snorm = 0, tmp;
    matlib_index i;
    for (i=0; i<N; i++)
    {
        tmp = 	  _N0 * *(u+0) * *(u+0)
        	+ _N1 * *(u+1) * *(u+1)
        	+ _N2 * *(u+2) * *(u+2)
        	+ _N3 * *(u+3) * *(u+3)
        	+ _N4 * *(u+4) * *(u+4)
        	+ _N5 * *(u+5) * *(u+5)
        	+ _N6 * *(u+6) * *(u+6);
        u += 7;
        snorm += tmp;
    }
    return snorm;
}
 
 
matlib_real fem1d_xlp_snorm2_d_7(matlib_index N, matlib_real *u)
{
    matlib_real snorm = 0, tmp;
    matlib_index i;
    for (i=0; i<N; i++)
    {
        tmp = 	  _N0 * *(u+0) * *(u+0)
        	+ _N1 * *(u+1) * *(u+1)
        	+ _N2 * *(u+2) * *(u+2)
        	+ _N3 * *(u+3) * *(u+3)
        	+ _N4 * *(u+4) * *(u+4)
        	+ _N5 * *(u+5) * *(u+5)
        	+ _N6 * *(u+6) * *(u+6)
        	+ _N7 * *(u+7) * *(u+7);
        u += 8;
        snorm += tmp;
    }
    return snorm;
}
 
 
matlib_real fem1d_xlp_snorm2_d_8(matlib_index N, matlib_real *u)
{
    matlib_real snorm = 0, tmp;
    matlib_index i;
    for (i=0; i<N; i++)
    {
        tmp = 	  _N0 * *(u+0) * *(u+0)
        	+ _N1 * *(u+1) * *(u+1)
        	+ _N2 * *(u+2) * *(u+2)
        	+ _N3 * *(u+3) * *(u+3)
        	+ _N4 * *(u+4) * *(u+4)
        	+ _N5 * *(u+5) * *(u+5)
        	+ _N6 * *(u+6) * *(u+6)
        	+ _N7 * *(u+7) * *(u+7)
        	+ _N8 * *(u+8) * *(u+8);
        u += 9;
        snorm += tmp;
    }
    return snorm;
}
 
 
matlib_real fem1d_xlp_snorm2_d_9(matlib_index N, matlib_real *u)
{
    matlib_real snorm = 0, tmp;
    matlib_index i;
    for (i=0; i<N; i++)
    {
        tmp = 	  _N0 * *(u+0) * *(u+0)
        	+ _N1 * *(u+1) * *(u+1)
        	+ _N2 * *(u+2) * *(u+2)
        	+ _N3 * *(u+3) * *(u+3)
        	+ _N4 * *(u+4) * *(u+4)
        	+ _N5 * *(u+5) * *(u+5)
        	+ _N6 * *(u+6) * *(u+6)
        	+ _N7 * *(u+7) * *(u+7)
        	+ _N8 * *(u+8) * *(u+8)
        	+ _N9 * *(u+9) * *(u+9);
        u += 10;
        snorm += tmp;
    }
    return snorm;
}
 
 
matlib_real fem1d_xlp_snorm2_d_10(matlib_index N, matlib_real *u)
{
    matlib_real snorm = 0, tmp;
    matlib_index i;
    for (i=0; i<N; i++)
    {
        tmp = 	  _N0 * *(u+0) * *(u+0)
        	+ _N1 * *(u+1) * *(u+1)
        	+ _N2 * *(u+2) * *(u+2)
        	+ _N3 * *(u+3) * *(u+3)
        	+ _N4 * *(u+4) * *(u+4)
        	+ _N5 * *(u+5) * *(u+5)
        	+ _N6 * *(u+6) * *(u+6)
        	+ _N7 * *(u+7) * *(u+7)
        	+ _N8 * *(u+8) * *(u+8)
        	+ _N9 * *(u+9) * *(u+9)
        	+ _N10 * *(u+10) * *(u+10);
        u += 11;
        snorm += tmp;
    }
    return snorm;
}
 

/* COMPLEX VERSION */ 
 
matlib_real fem1d_zlp_snorm2_d_2(matlib_index N, matlib_complex *u)
{
    matlib_real snorm = 0, tmp;
    matlib_index i;
    for (i=0; i<N; i++)
    {
        tmp = 	  _N0 * *(u+0) * conj(*(u+0))
        	+ _N1 * *(u+1) * conj(*(u+1))
        	+ _N2 * *(u+2) * conj(*(u+2));
        u += 3;
        snorm += tmp;
    }
    return snorm;
}
 
 
matlib_real fem1d_zlp_snorm2_d_3(matlib_index N, matlib_complex *u)
{
    matlib_real snorm = 0, tmp;
    matlib_index i;
    for (i=0; i<N; i++)
    {
        tmp = 	  _N0 * *(u+0) * conj(*(u+0))
        	+ _N1 * *(u+1) * conj(*(u+1))
        	+ _N2 * *(u+2) * conj(*(u+2))
        	+ _N3 * *(u+3) * conj(*(u+3));
        u += 4;
        snorm += tmp;
    }
    return snorm;
}
 
 
matlib_real fem1d_zlp_snorm2_d_4(matlib_index N, matlib_complex *u)
{
    matlib_real snorm = 0, tmp;
    matlib_index i;
    for (i=0; i<N; i++)
    {
        tmp = 	  _N0 * *(u+0) * conj(*(u+0))
        	+ _N1 * *(u+1) * conj(*(u+1))
        	+ _N2 * *(u+2) * conj(*(u+2))
        	+ _N3 * *(u+3) * conj(*(u+3))
        	+ _N4 * *(u+4) * conj(*(u+4));
        u += 5;
        snorm += tmp;
    }
    return snorm;
}
 
 
matlib_real fem1d_zlp_snorm2_d_5(matlib_index N, matlib_complex *u)
{
    matlib_real snorm = 0, tmp;
    matlib_index i;
    for (i=0; i<N; i++)
    {
        tmp = 	  _N0 * *(u+0) * conj(*(u+0))
        	+ _N1 * *(u+1) * conj(*(u+1))
        	+ _N2 * *(u+2) * conj(*(u+2))
        	+ _N3 * *(u+3) * conj(*(u+3))
        	+ _N4 * *(u+4) * conj(*(u+4))
        	+ _N5 * *(u+5) * conj(*(u+5));
        u += 6;
        snorm += tmp;
    }
    return snorm;
}
 
 
matlib_real fem1d_zlp_snorm2_d_6(matlib_index N, matlib_complex *u)
{
    matlib_real snorm = 0, tmp;
    matlib_index i;
    for (i=0; i<N; i++)
    {
        tmp = 	  _N0 * *(u+0) * conj(*(u+0))
        	+ _N1 * *(u+1) * conj(*(u+1))
        	+ _N2 * *(u+2) * conj(*(u+2))
        	+ _N3 * *(u+3) * conj(*(u+3))
        	+ _N4 * *(u+4) * conj(*(u+4))
        	+ _N5 * *(u+5) * conj(*(u+5))
        	+ _N6 * *(u+6) * conj(*(u+6));
        u += 7;
        snorm += tmp;
    }
    return snorm;
}
 
 
matlib_real fem1d_zlp_snorm2_d_7(matlib_index N, matlib_complex *u)
{
    matlib_real snorm = 0, tmp;
    matlib_index i;
    for (i=0; i<N; i++)
    {
        tmp = 	  _N0 * *(u+0) * conj(*(u+0))
        	+ _N1 * *(u+1) * conj(*(u+1))
        	+ _N2 * *(u+2) * conj(*(u+2))
        	+ _N3 * *(u+3) * conj(*(u+3))
        	+ _N4 * *(u+4) * conj(*(u+4))
        	+ _N5 * *(u+5) * conj(*(u+5))
        	+ _N6 * *(u+6) * conj(*(u+6))
        	+ _N7 * *(u+7) * conj(*(u+7));
        u += 8;
        snorm += tmp;
    }
    return snorm;
}
 
 
matlib_real fem1d_zlp_snorm2_d_8(matlib_index N, matlib_complex *u)
{
    matlib_real snorm = 0, tmp;
    matlib_index i;
    for (i=0; i<N; i++)
    {
        tmp = 	  _N0 * *(u+0) * conj(*(u+0))
        	+ _N1 * *(u+1) * conj(*(u+1))
        	+ _N2 * *(u+2) * conj(*(u+2))
        	+ _N3 * *(u+3) * conj(*(u+3))
        	+ _N4 * *(u+4) * conj(*(u+4))
        	+ _N5 * *(u+5) * conj(*(u+5))
        	+ _N6 * *(u+6) * conj(*(u+6))
        	+ _N7 * *(u+7) * conj(*(u+7))
        	+ _N8 * *(u+8) * conj(*(u+8));
        u += 9;
        snorm += tmp;
    }
    return snorm;
}
 
 
matlib_real fem1d_zlp_snorm2_d_9(matlib_index N, matlib_complex *u)
{
    matlib_real snorm = 0, tmp;
    matlib_index i;
    for (i=0; i<N; i++)
    {
        tmp = 	  _N0 * *(u+0) * conj(*(u+0))
        	+ _N1 * *(u+1) * conj(*(u+1))
        	+ _N2 * *(u+2) * conj(*(u+2))
        	+ _N3 * *(u+3) * conj(*(u+3))
        	+ _N4 * *(u+4) * conj(*(u+4))
        	+ _N5 * *(u+5) * conj(*(u+5))
        	+ _N6 * *(u+6) * conj(*(u+6))
        	+ _N7 * *(u+7) * conj(*(u+7))
        	+ _N8 * *(u+8) * conj(*(u+8))
        	+ _N9 * *(u+9) * conj(*(u+9));
        u += 10;
        snorm += tmp;
    }
    return snorm;
}
 
 
matlib_real fem1d_zlp_snorm2_d_10(matlib_index N, matlib_complex *u)
{
    matlib_real snorm = 0, tmp;
    matlib_index i;
    for (i=0; i<N; i++)
    {
        tmp = 	  _N0 * *(u+0) * conj(*(u+0))
        	+ _N1 * *(u+1) * conj(*(u+1))
        	+ _N2 * *(u+2) * conj(*(u+2))
        	+ _N3 * *(u+3) * conj(*(u+3))
        	+ _N4 * *(u+4) * conj(*(u+4))
        	+ _N5 * *(u+5) * conj(*(u+5))
        	+ _N6 * *(u+6) * conj(*(u+6))
        	+ _N7 * *(u+7) * conj(*(u+7))
        	+ _N8 * *(u+8) * conj(*(u+8))
        	+ _N9 * *(u+9) * conj(*(u+9))
        	+ _N10 * *(u+10) * conj(*(u+10));
        u += 11;
        snorm += tmp;
    }
    return snorm;
}
 


/*============================================================================+/
 | Dot product of two vectors
/+============================================================================*/


matlib_real fem1d_xdot
( 
    matlib_index l,
    matlib_real* x,
    matlib_real* y
)
{
    matlib_index i;
    matlib_real sum = *x**y;

    for( i=1; i<l; i++)
    {
        x++;
        y++;
        sum += *x**y;
    
    }
    return sum;
}

matlib_complex fem1d_zdot
( 
    matlib_index           l,
    matlib_complex* x,
    matlib_complex* y
)
{
    matlib_index i;
    matlib_real sum = *x**y;

    for( i=1; i<l; i++)
    {
        x++;
        y++;
        sum += *x**y;
    
    }
    return sum;
}
/*============================================================================*/

matlib_real fem1d_xdot2(matlib_index l, matlib_real *u, matlib_real *v)
{
    matlib_real tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1];
    return tmp;
}

matlib_real fem1d_xdot3(matlib_index l, matlib_real *u, matlib_real *v)
{
    matlib_real tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2];
    return tmp;
}

matlib_real fem1d_xdot4(matlib_index l, matlib_real *u, matlib_real *v)
{
    matlib_real tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3];
    return tmp;
}

matlib_real fem1d_xdot5(matlib_index l, matlib_real *u, matlib_real *v)
{
    matlib_real tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4];
    return tmp;
}
 
 
matlib_real fem1d_xdot6(matlib_index l, matlib_real *u, matlib_real *v)
{
    matlib_real tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4]
        	+ u[5] * v[5];
    return tmp;
}

matlib_real fem1d_xdot7(matlib_index l, matlib_real *u, matlib_real *v)
{
    matlib_real tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4]
        	+ u[5] * v[5]
        	+ u[6] * v[6];
    return tmp;
}

matlib_real fem1d_xdot8(matlib_index l, matlib_real *u, matlib_real *v)
{
    matlib_real tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4]
        	+ u[5] * v[5]
        	+ u[6] * v[6]
        	+ u[7] * v[7];
    return tmp;
}

matlib_real fem1d_xdot9(matlib_index l, matlib_real *u, matlib_real *v)
{
    matlib_real tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4]
        	+ u[5] * v[5]
        	+ u[6] * v[6]
        	+ u[7] * v[7]
        	+ u[8] * v[8];
    return tmp;
}
 
 
matlib_real fem1d_xdot10(matlib_index l, matlib_real *u, matlib_real *v)
{
    matlib_real tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4]
        	+ u[5] * v[5]
        	+ u[6] * v[6]
        	+ u[7] * v[7]
        	+ u[8] * v[8]
        	+ u[9] * v[9];
    return tmp;
}
 
 
matlib_real fem1d_xdot11(matlib_index l, matlib_real *u, matlib_real *v)
{
    matlib_real tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4]
        	+ u[5] * v[5]
        	+ u[6] * v[6]
        	+ u[7] * v[7]
        	+ u[8] * v[8]
        	+ u[9] * v[9]
        	+ u[10] * v[10];
    return tmp;
}
 
 
matlib_real fem1d_xdot12(matlib_index l, matlib_real *u, matlib_real *v)
{
    matlib_real tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4]
        	+ u[5] * v[5]
        	+ u[6] * v[6]
        	+ u[7] * v[7]
        	+ u[8] * v[8]
        	+ u[9] * v[9]
        	+ u[10] * v[10]
        	+ u[11] * v[11];
    return tmp;
}
 
 
matlib_real fem1d_xdot13(matlib_index l, matlib_real *u, matlib_real *v)
{
    matlib_real tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4]
        	+ u[5] * v[5]
        	+ u[6] * v[6]
        	+ u[7] * v[7]
        	+ u[8] * v[8]
        	+ u[9] * v[9]
        	+ u[10] * v[10]
        	+ u[11] * v[11]
        	+ u[12] * v[12];
    return tmp;
}
 
 
matlib_real fem1d_xdot14(matlib_index l, matlib_real *u, matlib_real *v)
{
    matlib_real tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4]
        	+ u[5] * v[5]
        	+ u[6] * v[6]
        	+ u[7] * v[7]
        	+ u[8] * v[8]
        	+ u[9] * v[9]
        	+ u[10] * v[10]
        	+ u[11] * v[11]
        	+ u[12] * v[12]
        	+ u[13] * v[13];
    return tmp;
}
 
 
matlib_real fem1d_xdot15(matlib_index l, matlib_real *u, matlib_real *v)
{
    matlib_real tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4]
        	+ u[5] * v[5]
        	+ u[6] * v[6]
        	+ u[7] * v[7]
        	+ u[8] * v[8]
        	+ u[9] * v[9]
        	+ u[10] * v[10]
        	+ u[11] * v[11]
        	+ u[12] * v[12]
        	+ u[13] * v[13]
        	+ u[14] * v[14];
    return tmp;
}
 
 
matlib_real fem1d_xdot16(matlib_index l, matlib_real *u, matlib_real *v)
{
    matlib_real tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4]
        	+ u[5] * v[5]
        	+ u[6] * v[6]
        	+ u[7] * v[7]
        	+ u[8] * v[8]
        	+ u[9] * v[9]
        	+ u[10] * v[10]
        	+ u[11] * v[11]
        	+ u[12] * v[12]
        	+ u[13] * v[13]
        	+ u[14] * v[14]
        	+ u[15] * v[15];
    return tmp;
}
 
/* COMPLEX VERSION */
 
matlib_complex fem_zdot2(matlib_index l, matlib_complex *u, matlib_complex *v)
{
    matlib_complex tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1];
    return tmp;
}
 
 
matlib_complex fem_zdot3(matlib_index l, matlib_complex *u, matlib_complex *v)
{
    matlib_complex tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2];
    return tmp;
}
 
 
matlib_complex fem_zdot4(matlib_index l, matlib_complex *u, matlib_complex *v)
{
    matlib_complex tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3];
    return tmp;
}
 
 
matlib_complex fem_zdot5(matlib_index l, matlib_complex *u, matlib_complex *v)
{
    matlib_complex tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4];
    return tmp;
}
 
 
matlib_complex fem_zdot6(matlib_index l, matlib_complex *u, matlib_complex *v)
{
    matlib_complex tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4]
        	+ u[5] * v[5];
    return tmp;
}
 
 
matlib_complex fem_zdot7(matlib_index l, matlib_complex *u, matlib_complex *v)
{
    matlib_complex tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4]
        	+ u[5] * v[5]
        	+ u[6] * v[6];
    return tmp;
}
 
 
matlib_complex fem_zdot8(matlib_index l, matlib_complex *u, matlib_complex *v)
{
    matlib_complex tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4]
        	+ u[5] * v[5]
        	+ u[6] * v[6]
        	+ u[7] * v[7];
    return tmp;
}
 
 
matlib_complex fem_zdot9(matlib_index l, matlib_complex *u, matlib_complex *v)
{
    matlib_complex tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4]
        	+ u[5] * v[5]
        	+ u[6] * v[6]
        	+ u[7] * v[7]
        	+ u[8] * v[8];
    return tmp;
}
 
 
matlib_complex fem_zdot10(matlib_index l, matlib_complex *u, matlib_complex *v)
{
    matlib_complex tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4]
        	+ u[5] * v[5]
        	+ u[6] * v[6]
        	+ u[7] * v[7]
        	+ u[8] * v[8]
        	+ u[9] * v[9];
    return tmp;
}
 
 
matlib_complex fem_zdot11(matlib_index l, matlib_complex *u, matlib_complex *v)
{
    matlib_complex tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4]
        	+ u[5] * v[5]
        	+ u[6] * v[6]
        	+ u[7] * v[7]
        	+ u[8] * v[8]
        	+ u[9] * v[9]
        	+ u[10] * v[10];
    return tmp;
}
 
 
matlib_complex fem_zdot12(matlib_index l, matlib_complex *u, matlib_complex *v)
{
    matlib_complex tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4]
        	+ u[5] * v[5]
        	+ u[6] * v[6]
        	+ u[7] * v[7]
        	+ u[8] * v[8]
        	+ u[9] * v[9]
        	+ u[10] * v[10]
        	+ u[11] * v[11];
    return tmp;
}
 
 
matlib_complex fem_zdot13(matlib_index l, matlib_complex *u, matlib_complex *v)
{
    matlib_complex tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4]
        	+ u[5] * v[5]
        	+ u[6] * v[6]
        	+ u[7] * v[7]
        	+ u[8] * v[8]
        	+ u[9] * v[9]
        	+ u[10] * v[10]
        	+ u[11] * v[11]
        	+ u[12] * v[12];
    return tmp;
}
 
 
matlib_complex fem_zdot14(matlib_index l, matlib_complex *u, matlib_complex *v)
{
    matlib_complex tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4]
        	+ u[5] * v[5]
        	+ u[6] * v[6]
        	+ u[7] * v[7]
        	+ u[8] * v[8]
        	+ u[9] * v[9]
        	+ u[10] * v[10]
        	+ u[11] * v[11]
        	+ u[12] * v[12]
        	+ u[13] * v[13];
    return tmp;
}
 
 
matlib_complex fem_zdot15(matlib_index l, matlib_complex *u, matlib_complex *v)
{
    matlib_complex tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4]
        	+ u[5] * v[5]
        	+ u[6] * v[6]
        	+ u[7] * v[7]
        	+ u[8] * v[8]
        	+ u[9] * v[9]
        	+ u[10] * v[10]
        	+ u[11] * v[11]
        	+ u[12] * v[12]
        	+ u[13] * v[13]
        	+ u[14] * v[14];
    return tmp;
}

matlib_complex fem_zdot16(matlib_index l, matlib_complex *u, matlib_complex *v)
{
    matlib_complex tmp;
    tmp =	  u[0] * v[0]
        	+ u[1] * v[1]
        	+ u[2] * v[2]
        	+ u[3] * v[3]
        	+ u[4] * v[4]
        	+ u[5] * v[5]
        	+ u[6] * v[6]
        	+ u[7] * v[7]
        	+ u[8] * v[8]
        	+ u[9] * v[9]
        	+ u[10] * v[10]
        	+ u[11] * v[11]
        	+ u[12] * v[12]
        	+ u[13] * v[13]
        	+ u[14] * v[14]
        	+ u[15] * v[15];
    return tmp;
}

/*============================================================================+/
 | Building Global Mass Matrix 
/+============================================================================*/

void fem1d_quadM
(
    matlib_xv  quadW,
    matlib_xm  IM,
    matlib_xm* Q
)
/* 
 * IM: Inverse transform matrix, nr. LP basis elements: IM.lenr,
 * sampling points: IM.lenc, 
 * 
 * IM is needed in a colum major format
 * quadW: quadrature weight
 * quadW.len = IM.lenc
 * 
 * The quadrature matrix is defined as 
 *              Q = Transpose(Diag(quadW)*(IM*B))
 * This function outputs Q in row major format so that samples of the basis
 * functions are adjacent in memory.
 *
 * */
{
    debug_enter( "size of quadW: %d, "
                 "size of IM: %d-by-%d ",
                 quadW.len, 
                 IM.lenc, IM.lenr);

    

    warn_if(!(IM.lenc>IM.lenr), "%s", "Numerical quadrature becomes inaccurate" );

    if(quadW.len==IM.lenc)
    {
        matlib_index i, j, k;
        /* Highest polynomial degree: p = IM.lenr -1
         * nr_combi = 3+(p-1)*(p+4)/2 
         * */ 
        matlib_index nr_combi = 3+(IM.lenr-2)*(IM.lenr+3)/2;

        debug_body( "nr. of basis elements: %d, "
                    "nr. of combination of basis-elements: %d", 
                    IM.lenr, nr_combi);
        
        matlib_create_xm( nr_combi, IM.lenc, Q, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);

        matlib_xm B, M;
        matlib_create_xm( IM.lenr, IM.lenr, &B, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);
        matlib_create_xm( IM.lenc, IM.lenr, &M, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);

        /* zeroth column */ 
        *B.elem_p           =  1.0/2;
        *(B.elem_p+1)       = -1.0/2; 

        /* first column */ 
        *(B.elem_p+0+B.lenc) =  1.0/2;
        *(B.elem_p+1+B.lenc) =  1.0/2;

        for (i=2; i<B.lenr; i++)
        {
            *(B.elem_p+i   + i*B.lenc)   =  1.0/sqrt(2.0*(2.0*i-1.0));
            *(B.elem_p+i-2 + i*B.lenc)   = -1.0/sqrt(2.0*(2.0*i-1.0));
        }
        DEBUG_PRINT_XM(B, "%s: ", "FEM-basis in terms of LP");

        if(IM.order==MATLIB_COL_MAJOR)
        {
            /* 
             * M <--one*IMi*B + zero*M 
             *
             * jth basis function values (vertex/bubble): {M_ij}_{0<=i<=lenc}
             *
             * */ 
            matlib_xgemm(1.0, IM, B, 0, M);
            DEBUG_PRINT_XM(M, "%s (IM: col major): ", "sampling FEM-basis functions");
        }
        else if (IM.order==MATLIB_ROW_MAJOR)
        {
            /* Define the transpose and put the 'op' info as 'MATLIB_TRANS' so
             * that it can serve the purpose of IM in COL_MAJOR format.
             * */ 
            matlib_xm IM_T = { .lenc  = IM.lenr,
                               .lenr  = IM.lenc,
                               .order = MATLIB_COL_MAJOR,
                               .op    = MATLIB_TRANS,     
                               .elem_p = IM.elem_p};

            matlib_xgemm(1.0, IM_T, B, 0, M);
            DEBUG_PRINT_XM(M, "%s (IM: row major): ", "sampling FEM-basis functions");
        }
        else
        {
            term_execb("Storage order of the matrix unknown: (order: %d)", M.order);        
        }
        DEBUG_PRINT_XM(M, "%s: ", "sampling FEM-basis functions");

        for (i=0; i<M.lenc; i++)
        {
            Q->elem_p[0*M.lenc + i] = *(M.elem_p+i       )**(M.elem_p+i       ) * quadW.elem_p[i];
            Q->elem_p[1*M.lenc + i] = *(M.elem_p+i       )**(M.elem_p+i+M.lenc) * quadW.elem_p[i];
            Q->elem_p[2*M.lenc + i] = *(M.elem_p+i+M.lenc)**(M.elem_p+i+M.lenc) * quadW.elem_p[i];
        }
        
        for (j=2; j<M.lenr; j++)
        {
            for (i=0; i<M.lenc; i++)
            {
                Q->elem_p[(j+1)*M.lenc + i] = 
                    *(M.elem_p+i)**(M.elem_p+i+j*M.lenc) * quadW.elem_p[i];
                
                Q->elem_p[(j+M.lenr-1)*M.lenc + i] = 
                    *(M.elem_p+i+M.lenc)**(M.elem_p+i+j*M.lenc) * quadW.elem_p[i];
            }
        }
        
        matlib_index shiftIn = 2*(IM.lenr-2)+3;

        for (k=2; k<M.lenr; k++) /* runs over higher order LGL-functions, total nr.: (p-1)*/ 
        {
            for (j=k; j<M.lenr; j++)
            {
                for (i=0; i<M.lenc; i++)
                {
                    Q->elem_p[shiftIn*M.lenc+i] = 
                        *(M.elem_p+i+k*M.lenc)**(M.elem_p+i+j*M.lenc)*quadW.elem_p[i];
                }
                shiftIn++;
            }
        }
        BEGIN_DTRACE
            matlib_real ONES[quadW.len];
            for (i=0; i<quadW.len; i++){ONES[i]=1.0; }

            matlib_real value;
            
            value = fem1d_xdot(Q->lenr, Q->elem_p, ONES);
            debug_body("inner product (%d,%d): sum_iQ(0,i) = % 0.16f", 0, 0, value);
            
            value = fem1d_xdot(Q->lenr, Q->elem_p+Q->lenr, ONES);
            debug_body("inner product (%d,%d): sum_iQ(1,i) = % 0.16f", 0, 1, value);
            
            value = fem1d_xdot(Q->lenr, Q->elem_p+2*Q->lenr, ONES);
            debug_body("inner product (%d,%d): sum_iQ(2,i) = % 0.16f", 1, 1, value);

            for (j=2; j<M.lenr; j++)
            {
                value = fem1d_xdot(Q->lenr, Q->elem_p+(j+1)*Q->lenr, ONES);
                debug_body( "inner product (%d,%d): sum_iQ(%d,i) = % 0.16f",
                            0, j, j+1, value);

            }
            for (j=2; j<M.lenr; j++)
            {
                value = fem1d_xdot(Q->lenr, Q->elem_p+(j+M.lenr-1)*Q->lenr, ONES);
                debug_body( "inner product (%d,%d): sum_iQ(%d,i) = % 0.16f", 
                            1, j, j+M.lenr-1, value);
            }
        
            shiftIn = 2*(IM.lenr-2)+3;

            for (k=2; k<M.lenr; k++)
            {
                for (j=k; j<M.lenr; j++)
                {
                    value = fem1d_xdot(Q->lenr, Q->elem_p+shiftIn*Q->lenr, ONES);
                    debug_body( "inner product (%d,%d): sum_iQ(%d,i) = % 0.16f",
                                k, j, shiftIn, value);
                    shiftIn++;
                }
            }
        END_DTRACE
    }
    else
    {
        term_exec( "dimension incorrect: length of quadW: %d, "
                   "size of IM: %d-by-%d", quadW.len, IM.lenc, IM.lenr );
    }
    debug_exit("%s", "");

}/* fem1d_quadM */ 

/*============================================================================*/

void fem1d_MEMI
/* Master Element Mass Integral */ 
(
    matlib_index p,
    matlib_xv*   q
)
{
    debug_enter(" degree of polynomial: %d", p);
    matlib_index nr_combi = 3+(p-1)*(p+4)/2; 

    debug_body(" nr. of combinations of basis elements: %d", nr_combi);

    matlib_create_xv( nr_combi, q, MATLIB_COL_VECT);
    
    matlib_index i, j, shiftq;
    /* l_i --> Lobatto-heirarchic fem-shape functions 
     * */ 
    q->elem_p[0] = 2.0/3; /* (l0, l0) */ 
    q->elem_p[1] = 1.0/3; /* (l0, l1) */
    q->elem_p[2] = 2.0/3; /* (l1, l1) */

    /* higher order basis elements
     * */ 
    q->elem_p[3] = -sqrt(6.0)/6; /* (l0, l2) */ 
    q->elem_p[4] =  sqrt(10)/30; /* (l0, l3) */
    for(i=2; i<p-1; i++)
    {
        q->elem_p[i+3] = 0; /* (l0, l_i+3)*/ 
    }
    shiftq = 3 +(p-1);
    q->elem_p[shiftq]   = -sqrt(6.0)/6; /* (l1, l2) */ 
    q->elem_p[shiftq+1] = -sqrt(10)/30; /* (l1, l3) */
    for(i=2; i<p-1; i++)
    {
        q->elem_p[i+shiftq] = 0; /* (l1, l_i+3)*/ 
    }
    shiftq = 3 + 2*(p-1);
    for (i=2; i<p+1; i++)
    {
        q->elem_p[shiftq] = 2.0/((2.0*i+1)*(2.0*i-3));
        shiftq ++;

        for(j=i+1; j<p+1; j++)
        {
            if(j==i+2)
            {
                q->elem_p[shiftq] = -1.0/((2.0*i+1)*sqrt((2.0*i-1)*(2.0*i+3)));
            }
            else
            {
                q->elem_p[shiftq] = 0;
            }
            shiftq ++;
        }
    }
    debug_exit("%s", "");
}
/*============================================================================*/

void fem1d_GMMSparsity
/* Determine the sparsity structure in CSR format */ 
(
    matlib_index            p,
    matlib_index            N,                        
    matlib_index*           row,                     
    matlib_index*           col                     
)
{

    debug_enter( "degree of polynomial: %d, "
                 "Number of finite-elements: %d",
                 p, N) ;

    /* number of nonzero elements nnz = N*p*(p+3)/2+1;
     * */
    matlib_index nr_combi = (p-1)*(p+4)/2+3;
    debug_body("nr. of combinations of basis functions: %d", nr_combi);

    matlib_index s, i, l, l0, m, st;

    i = 0;
    s = 0;
    
    /* In this routine, we traverse the mass matrix row-wise
     * and record the row entries and column index
     * the mass matrix is symmetric therefore it suffices to just store the
     * upper-triangular elements
     * i-th row: 
     * matrix elements: ugpmm[row[i]]... ugpmm[row[i+1]-1]
     * column matlib_index   : col[row[i]] < col[row[i]+1] <...< col[row[i+1]-1]
     * For symmetric matrices, col[row[i]]=i
     *
     * */ 

    /* For quadrature it is important to note that 
     * K_1 = (x_0, x_1), K_2 = (x_1, x_2),...
     * basis element v_1 is peaked at x_1, v_2 at x_2,...
     * */ 

    /* zero-th row 
     *
     * (v_0,v_0)_{K_1}
     * (v_0,v_1)_{K_1}
     * */ 
    row[i] = s;
    col[s] = i; s++;
    
    col[s] = i+1; s++;
    
    for (l=N+1; l<(N+p); l++)
    {
        col[s] = l;
        //debug_body("s=%d", s);
        s++;
    }

    l0 = N+1;

    for( i=1; i<N; i++)
    {
        /* i-th row: start from the diagonal 
         * (v_i,v_i)_{K_i} + (v_i,v_i)_{K_{i+1}}
         * (v_i,v_{i+1})_{K_{i+1}}
         * */ 
        row[i] = s;
        col[s] = i;
        //debug_body("s=%d", s);
        s++;
        col[s] = i + 1; 
        //debug_body("s=%d", s);
        s++;
        /* v_i interacts with [b_i?], [b_{i+1,?}]
         * (p-1) bubble functions in each element
         */
        for (l=l0; l<(l0+p-1); l++)
        {
            col[s] = l;
            //debug_body("s=%d", s);
            s++;
        }
        l0 = l;

        for (l=l0; l<(l0+p-1); l++)
        {
            col[s] = l;
            //debug_body("s=%d", s);
            s++;
        
        }
    }
    /* last row: only the diagonal 
     * (v_N,v_N)_{K_N}
     * */ 
    row[i] = s;
    col[s] = i;
    //debug_body("s=%d", s);
    s++;
    for (l=l0; l<(l0+p-1); l++)
    {
        col[s] = l;
        //debug_body("s=%d", s);
        s++;
    }
    i++;

    /* vertex-vertex and vertex-bubble combinations are done
     * do the bubble-bubble combination
     * */ 

    for( st=0; st<N*nr_combi; st +=nr_combi)
    {
        for( l=0; l<p-1; l++)
        {
            row[i] = s;
            for (m=i; m<(i+p-1-l); m++)
            {
                col[s] = m;
                //debug_body("s=%d", s);
                s++;
            }
            i++;
        }
    }
    row[i] = s;
    debug_exit("%s", "");
}
/*============================================================================*/

void fem1d_XCSRGMM
/* Double CSR - Assemble Global Mass Matrix*/ 
(
    matlib_index  p,
    matlib_index  N,                        
    matlib_xv     q,
    matlib_index* row,                     
    matlib_index* col,                     
    matlib_real*  ugpmm                   
                                    
)
{

    debug_enter( "degree of polynomial: %d, "
                 "Number of finite-elements: %d, "
                 "length of q: %d", 
                 p, N, q.len) ;

   matlib_index nr_combi = q.len/N; 
   debug_body("nr. of combinations of basis functions: %d", nr_combi);
   

    /* number of nonzero elements nnz = N*p*(p+3)/2+1;
     * */
    matlib_index s, i, l, l0, m, st;

    i = 0;
    s = 0;
    
    /* In this routine, we traverse the mass matrix row-wise
     * and record the row entries and column index
     * the mass matrix is symmetric therefore it suffices to just store the
     * upper-triangular elements
     * i-th row: 
     * matrix elements: ugpmm[row[i]]... ugpmm[row[i+1]-1]
     * column matlib_index   : col[row[i]] < col[row[i]+1] <...< col[row[i+1]-1]
     * For symmetric matrices, col[row[i]]=i
     *
     * */ 

    /* For quadrature it is important to note that 
     * K_1 = (x_0, x_1), K_2 = (x_1, x_2),...
     * basis element v_1 is peaked at x_1, v_2 at x_2,...
     * */ 

    /* zero-th row 
     *
     * (v_0,v_0)_{K_1}
     * (v_0,v_1)_{K_1}
     * */ 
    row[i] = s;
    col[s] = i;   ugpmm[s] = *(q.elem_p);
    debug_body( "row[%d]=%d, col[%d] = %d, ugpmm[%d] = % 0.16f",
                i, row[i], s, col[s], s, ugpmm[s]);
    s++;

    col[s] = i+1; ugpmm[s] = *(q.elem_p+1);
    debug_body( "row[%d]=%d, col[%d] = %d, ugpmm[%d] = % 0.16f",
                i, row[i], s, col[s], s, ugpmm[s]);
    s++;

    matlib_index shiftq = 3;
    
    for (l=N+1; l<(N+p); l++)
    {
        col[s] = l;
        ugpmm[s] = (*(q.elem_p + shiftq));
        debug_body( "row[%d]=%d, col[%d] = %d, ugpmm[%d] = % 0.16f",
                    i, row[i], s, col[s], s, ugpmm[s]);
        shiftq ++;
        s++;
    }

    l0 = N+1;
    st = 0;
    for( i=1; i<N; i++)
    {
        /* i-th row: start from the diagonal 
         * (v_i,v_i)_{K_i} + (v_i,v_i)_{K_{i+1}}
         * (v_i,v_{i+1})_{K_{i+1}}
         * */ 
        row[i] = s;
        col[s] = i;
        ugpmm[s] = *(q.elem_p+st+2) + *(q.elem_p+st+nr_combi);
        debug_body( "row[%d]=%d, col[%d] = %d, ugpmm[%d] = % 0.16f",
                    i, row[i], s, col[s], s, ugpmm[s]);
        s++;
        col[s] = i + 1; 
        ugpmm[s] = *(q.elem_p+st+nr_combi+1); 
        debug_body( "row[%d]=%d, col[%d] = %d, ugpmm[%d] = % 0.16f",
                    i, row[i], s, col[s], s, ugpmm[s]);
        s++;
        /* v_i interacts with [b_i?], [b_{i+1,?}]
         * (p-1) bubble functions in each element
         */
        shiftq = st+3+p-1;
        for (l=l0; l<(l0+p-1); l++)
        {
            col[s] = l;
            ugpmm[s] = *(q.elem_p + shiftq);
            debug_body( "row[%d]=%d, col[%d] = %d, ugpmm[%d] = % 0.16f",
                         i, row[i], s, col[s], s, ugpmm[s]);
            shiftq ++;
            s++;
        }
        l0 = l;

        shiftq = st + nr_combi+3;
        for (l=l0; l<(l0+p-1); l++)
        {
            col[s] = l;
            ugpmm[s] = *(q.elem_p+shiftq);
            debug_body( "row[%d]=%d, col[%d] = %d, ugpmm[%d] = % 0.16f",
                         i, row[i], s, col[s], s, ugpmm[s]);
            shiftq++;
            s++;
        
        }
        st = st+nr_combi;
    }
    debug_body("s=%d", s);
    /* last row: only the diagonal 
     * (v_N,v_N)_{K_N}
     * */ 
    row[i] = s;
    col[s] = i;
    ugpmm[s] = *(q.elem_p+st+2);
    debug_body( "row[%d]=%d, col[%d] = %d, ugpmm[%d] = % 0.16f",
                 i, row[i], s, col[s], s, ugpmm[s]);
    s++;
    shiftq = st+3+(p-1);
    for (l=l0; l<(l0+p-1); l++)
    {
        col[s] = l;
        ugpmm[s] = *(q.elem_p+shiftq);
        debug_body( "row[%d]=%d, col[%d] = %d, ugpmm[%d] = % 0.16f",
                     i, row[i], s, col[s], s, ugpmm[s]);
        shiftq++;
        s++;
    }
    i++;

    /* vertex-vertex and vertex-bubble combinations are done
     * do the bubble-bubble combination
     * */ 
    matlib_index shiftq0 = 3+2*(p-1);

    for( st=0; st<N*nr_combi; st +=nr_combi)
    {
        shiftq = shiftq0;
        for( l=0; l<p-1; l++)
        {
            row[i] = s;
            for (m=i; m<(i+p-1-l); m++)
            {
                col[s] = m;
                ugpmm[s] = *(q.elem_p+st+shiftq);
                debug_body( "row[%d]=%d, col[%d] = %d, ugpmm[%d] = % 0.16f",
                            i, row[i], s, col[s], s, ugpmm[s]);
                shiftq++;
                s++;
            }
            i++;
        }
    }
    row[i] = s;
    debug_exit("%s", "");
}
/*============================================================================*/
void fem1d_XCSRGMM2
/* Double CSR - Assemble Global Mass Matrix*/ 
(
    matlib_index     p,
    matlib_index     N,                        
    matlib_xm q,
    matlib_real**  ugpmm_p                   
                                    
)
{

    debug_enter( "degree of polynomial: %d, "
                 "Number of finite-elements: %d, "
                 "length of q: %d-by-%d", 
                 p, N, q.lenc, q.lenr);

   matlib_index nr_combi = q.lenc/N; 
   debug_body("nr. of combinations of basis functions: %d", nr_combi);
   

    /* number of nonzero elements nnz = N*p*(p+3)/2+1;
     * */
    matlib_index s, i, l, l0, m, st;

    i = 0;
    s = 0;
    
    /* In this routine, we traverse the mass matrix row-wise
     * and record the row entries and column index
     * the mass matrix is symmetric therefore it suffices to just store the
     * upper-triangular elements
     * i-th row: 
     * matrix elements: ugpmm[row[i]]... ugpmm[row[i+1]-1]
     * column matlib_index   : col[row[i]] < col[row[i]+1] <...< col[row[i+1]-1]
     * For symmetric matrices, col[row[i]]=i
     *
     * */ 

    /* For quadrature it is important to note that 
     * K_1 = (x_0, x_1), K_2 = (x_1, x_2),...
     * basis element v_1 is peaked at x_1, v_2 at x_2,...
     * */ 

    /* zero-th row 
     *
     * (v_0,v_0)_{K_1}
     * (v_0,v_1)_{K_1}
     * */ 
    matlib_real **ugpmm;

    for(ugpmm=ugpmm_p; ugpmm<ugpmm_p+q.lenr; ugpmm++)
    {
        (*ugpmm)[s] = *(q.elem_p);
        s++;
    
        (*ugpmm)[s] = *(q.elem_p+1);
        s++;

        matlib_index shiftq = 3;
        
        for (l=N+1; l<(N+p); l++)
        {
            (*ugpmm)[s] = (*(q.elem_p + shiftq));
            shiftq ++;
            s++;
        }

        l0 = N+1;
        st = 0;
        for( i=1; i<N; i++)
        {
            /* i-th row: start from the diagonal 
             * (v_i,v_i)_{K_i} + (v_i,v_i)_{K_{i+1}}
             * (v_i,v_{i+1})_{K_{i+1}}
             * */ 
            (*ugpmm)[s] = *(q.elem_p+st+2) + *(q.elem_p+st+nr_combi);
            s++;
            (*ugpmm)[s] = *(q.elem_p+st+nr_combi+1); 
            s++;
            /* v_i interacts with [b_i?], [b_{i+1,?}]
             * (p-1) bubble functions in each element
             */
            shiftq = st+3+p-1;
            for (l=l0; l<(l0+p-1); l++)
            {
                (*ugpmm)[s] = *(q.elem_p + shiftq);
                shiftq ++;
                s++;
            }
            l0 = l;

            shiftq = st + nr_combi+3;
            for (l=l0; l<(l0+p-1); l++)
            {
                (*ugpmm)[s] = *(q.elem_p+shiftq);
                shiftq++;
                s++;
            
            }
            st = st+nr_combi;
        }
        debug_body("s=%d", s);
        /* last row: only the diagonal 
         * (v_N,v_N)_{K_N}
         * */ 
        (*ugpmm)[s] = *(q.elem_p+st+2);
        s++;
        shiftq = st+3+(p-1);
        for (l=l0; l<(l0+p-1); l++)
        {
            (*ugpmm)[s] = *(q.elem_p+shiftq);
            shiftq++;
            s++;
        }       
        i++;

        /* vertex-vertex and vertex-bubble combinations are done
         * do the bubble-bubble combination
         * */ 
        matlib_index shiftq0 = 3+2*(p-1);

        for( st=0; st<N*nr_combi; st +=nr_combi)
        {
            shiftq = shiftq0;
            for( l=0; l<p-1; l++)
            {
                for (m=i; m<(i+p-1-l); m++)
                {
                    (*ugpmm)[s] = *(q.elem_p+st+shiftq);
                    shiftq++;
                    s++;
                }
                i++;
            }
        }
        q.elem_p += q.lenc;
    }

    debug_exit("%s", "");
}
/*============================================================================*/

void fem1d_ZCSRGMM
/* Complex CSR - Assemble Global Mass Matrix*/ 
(
    matlib_index     p,
    matlib_index     N,                        
    matlib_zv        q,
    matlib_index*    row,                     
    matlib_index*    col,                     
    matlib_complex*  ugpmm                   
                                    
)
{

    debug_enter( "degree of polynomial: %d, "
                 "Number of finite-elements: %d, "
                 "length of q: %d", 
                 p, N, q.len) ;

   matlib_index nr_combi = q.len/N; 
   debug_body("nr. of combinations of basis functions: %d", nr_combi);
   

    /* number of nonzero elements nnz = N*p*(p+3)/2+1;
     * */
    matlib_index s, i, l, l0, m, st;

    i = 0;
    s = 0;
    
    /* In this routine, we traverse the mass matrix row-wise
     * and record the row entries and column index
     * the mass matrix is symmetric therefore it suffices to just store the
     * upper-triangular elements
     * i-th row: 
     * matrix elements: ugpmm[row[i]]... ugpmm[row[i+1]-1]
     * column matlib_index   : col[row[i]] < col[row[i]+1] <...< col[row[i+1]-1]
     * For symmetric matrices, col[row[i]]=i
     *
     * */ 

    /* For quadrature it is important to note that 
     * K_1 = (x_0, x_1), K_2 = (x_1, x_2),...
     * basis element v_1 is peaked at x_1, v_2 at x_2,...
     * */ 

    /* zero-th row 
     *
     * (v_0,v_0)_{K_1}
     * (v_0,v_1)_{K_1}
     * */ 
    row[i] = s;
    col[s] = i;   ugpmm[s] = *(q.elem_p);      s++;
    
    col[s] = i+1; ugpmm[s] = *(q.elem_p+1);    s++;
    
    matlib_index shiftq = 3;
    
    for (l=N+1; l<(N+p); l++)
    {
        col[s] = l;
        ugpmm[s] = (*(q.elem_p + shiftq));
        debug_body( "row[%d]=%d, col[%d] = %d, ugpmm[%d] = % 0.16f% 0.16fi",
                    i, row[i], s, col[s], s, ugpmm[s]);
        shiftq ++;
        s++;
    }

    l0 = N+1;
    st = 0;
    for( i=1; i<N; i++)
    {
        /* i-th row: start from the diagonal 
         * (v_i,v_i)_{K_i} + (v_i,v_i)_{K_{i+1}}
         * (v_i,v_{i+1})_{K_{i+1}}
         * */ 
        row[i] = s;
        col[s] = i;
        ugpmm[s] = *(q.elem_p+st+2) + *(q.elem_p+st+nr_combi);
        debug_body( "row[%d]=%d, col[%d] = %d, ugpmm[%d] = % 0.16f%+0.16fi",
                    i, row[i], s, col[s], s, ugpmm[s]);
        s++;
        col[s] = i + 1; 
        ugpmm[s] = *(q.elem_p+st+nr_combi+1); 
        debug_body( "row[%d]=%d, col[%d] = %d, ugpmm[%d] = % 0.16f%+0.16fi",
                    i, row[i], s, col[s], s, ugpmm[s]);
        s++;
        /* v_i interacts with [b_i?], [b_{i+1,?}]
         * (p-1) bubble functions in each element
         */
        shiftq = st+3+p-1;
        for (l=l0; l<(l0+p-1); l++)
        {
            col[s] = l;
            ugpmm[s] = *(q.elem_p + shiftq);
            debug_body( "row[%d]=%d, col[%d] = %d, ugpmm[%d] = % 0.16f%+0.16fi",
                         i, row[i], s, col[s], s, ugpmm[s]);
            shiftq ++;
            s++;
        }
        l0 = l;

        shiftq = st + nr_combi+3;
        for (l=l0; l<(l0+p-1); l++)
        {
            col[s] = l;
            ugpmm[s] = *(q.elem_p+shiftq);
            debug_body( "row[%d]=%d, col[%d] = %d, ugpmm[%d] = % 0.16f%+0.16fi",
                         i, row[i], s, col[s], s, ugpmm[s]);
            shiftq++;
            s++;
        
        }
        st = st+nr_combi;
    }
    debug_body("s=%d", s);
    /* last row: only the diagonal 
     * (v_N,v_N)_{K_N}
     * */ 
    row[i] = s;
    col[s] = i;
    ugpmm[s] = *(q.elem_p+st+2);
    debug_body( "row[%d]=%d, col[%d] = %d, ugpmm[%d] = % 0.16f%+0.16fi",
                 i, row[i], s, col[s], s, ugpmm[s]);
    s++;
    shiftq = st+3+(p-1);
    for (l=l0; l<(l0+p-1); l++)
    {
        col[s] = l;
        ugpmm[s] = *(q.elem_p+shiftq);
        debug_body( "row[%d]=%d, col[%d] = %d, ugpmm[%d] = % 0.16f%+0.16fi",
                     i, row[i], s, col[s], s, ugpmm[s]);
        shiftq++;
        s++;
    }
    i++;

    /* vertex-vertex and vertex-bubble combinations are done
     * do the bubble-bubble combination
     * */ 
    matlib_index shiftq0 = 3+2*(p-1);

    for( st=0; st<N*nr_combi; st +=nr_combi)
    {
        shiftq = shiftq0;
        for( l=0; l<p-1; l++)
        {
            row[i] = s;
            for (m=i; m<(i+p-1-l); m++)
            {
                col[s] = m;
                ugpmm[s] = *(q.elem_p+st+shiftq);
                debug_body( "row[%d]=%d, col[%d] = %d, ugpmm[%d] = % 0.16f%+0.16fi",
                            i, row[i], s, col[s], s, ugpmm[s]);
                shiftq++;
                s++;
            }
            i++;
        }
    }
    row[i] = s;
    debug_exit("%s", "");
}

void fem1d_ZCSRGMM2
/* Complex CSR - Assemble Global Mass Matrix*/ 
(
    matlib_index           p,
    matlib_index           N,                        
    matlib_zm       q,
    matlib_complex** ugpmm_p                  
)
{

    debug_enter( "degree of polynomial: %d, "
                 "Number of finite-elements: %d, "
                 "length of q: %d-by-%d", 
                 p, N, q.lenc, q.lenr);

   matlib_index nr_combi = q.lenc/N; 
   debug_body("nr. of combinations of basis functions: %d", nr_combi);
   

    /* number of nonzero elements nnz = N*p*(p+3)/2+1;
     * */
    matlib_index s, i, l, l0, m, st;

    
    /* In this routine, we traverse the mass matrix row-wise
     * and record the row entries and column index
     * the mass matrix is symmetric therefore it suffices to just store the
     * upper-triangular elements
     * i-th row: 
     * matrix elements: ugpmm[row[i]]... ugpmm[row[i+1]-1]
     * column matlib_index   : col[row[i]] < col[row[i]+1] <...< col[row[i+1]-1]
     * For symmetric matrices, col[row[i]]=i
     *
     * */ 

    /* For quadrature it is important to note that 
     * K_1 = (x_0, x_1), K_2 = (x_1, x_2),...
     * basis element v_1 is peaked at x_1, v_2 at x_2,...
     * */ 

    /* zero-th row 
     *
     * (v_0,v_0)_{K_1}
     * (v_0,v_1)_{K_1}
     * */ 
    matlib_complex **ugpmm;

    for(ugpmm=ugpmm_p; ugpmm<ugpmm_p+q.lenr; ugpmm++)
    {
        i = 0;
        s = 0;
        (*ugpmm)[s] = *(q.elem_p);
        s++;
    
        (*ugpmm)[s] = *(q.elem_p+1);
        s++;

        matlib_index shiftq = 3;
        
        for (l=N+1; l<(N+p); l++)
        {
            (*ugpmm)[s] = (*(q.elem_p + shiftq));
            shiftq ++;
            s++;
        }

        l0 = N+1;
        st = 0;
        for( i=1; i<N; i++)
        {
            /* i-th row: start from the diagonal 
             * (v_i,v_i)_{K_i} + (v_i,v_i)_{K_{i+1}}
             * (v_i,v_{i+1})_{K_{i+1}}
             * */ 
            (*ugpmm)[s] = *(q.elem_p+st+2) + *(q.elem_p+st+nr_combi);
            s++;
            (*ugpmm)[s] = *(q.elem_p+st+nr_combi+1); 
            s++;
            /* v_i interacts with [b_i?], [b_{i+1,?}]
             * (p-1) bubble functions in each element
             */
            shiftq = st+3+p-1;
            for (l=l0; l<(l0+p-1); l++)
            {
                (*ugpmm)[s] = *(q.elem_p + shiftq);
                shiftq ++;
                s++;
            }
            l0 = l;

            shiftq = st + nr_combi+3;
            for (l=l0; l<(l0+p-1); l++)
            {
                (*ugpmm)[s] = *(q.elem_p+shiftq);
                shiftq++;
                s++;
            
            }
            st = st+nr_combi;
        }
        debug_body("s=%d", s);
        /* last row: only the diagonal 
         * (v_N,v_N)_{K_N}
         * */ 
        (*ugpmm)[s] = *(q.elem_p+st+2);
        s++;
        shiftq = st+3+(p-1);
        for (l=l0; l<(l0+p-1); l++)
        {
            (*ugpmm)[s] = *(q.elem_p+shiftq);
            shiftq++;
            s++;
        }
        i++;

        /* vertex-vertex and vertex-bubble combinations are done
         * do the bubble-bubble combination
         * */ 
        matlib_index shiftq0 = 3+2*(p-1);

        for( st=0; st<N*nr_combi; st +=nr_combi)
        {
            shiftq = shiftq0;
            for( l=0; l<p-1; l++)
            {
                for (m=i; m<(i+p-1-l); m++)
                {
                    (*ugpmm)[s] = *(q.elem_p+st+shiftq);
                    shiftq++;
                    s++;
                }
                i++;
            }
        }
        q.elem_p += q.lenc;
    }

    debug_exit("%s", "");
}
/*============================================================================*/
void fem1d_xm_sparse_GMM
/* Complex - Assemble Global Mass Matrix*/ 
(
    matlib_index      p,
    matlib_xm         Q,
    matlib_xv         phi,
    matlib_xm_sparse* M
)
{
    debug_enter( "poynomial degree: %d, "
                 "size of Q: %d-by-%d), "
                 "length of phi: %d",
                 p, Q.lenc, Q.lenr, phi.len);

    matlib_index N = (phi.len-1)/(Q.lenr-1);
    matlib_index nnz = N*(Q.lenc-1)+1;

    debug_body( "nr. finite-elements: %d, "
                "nr. of non-zero elements: %d", N, nnz);
    assert(Q.lenc==(p-1)*(p+4)/2+3);

    matlib_index dim = N*p+1;
    M->lenc   = dim;
    M->lenr   = dim;
    M->rowIn  = calloc( dim+1, sizeof(matlib_index));
    M->colIn  = calloc(   nnz, sizeof(matlib_index));
    M->elem_p = calloc(   nnz, sizeof(double));

    matlib_xv q;
    matlib_create_xv( Q.lenc*N, &q, MATLIB_COL_VECT);

    fem1d_XFLT( N, Q, phi, q);
    DEBUG_PRINT_XV(q, "%s: ", "inner product");

    fem1d_XCSRGMM(p, N, q, M->rowIn, M->colIn, M->elem_p);
    
    matlib_free((void*) q.elem_p);
    debug_exit("%s", "");

}

void fem1d_zm_sparse_GMM
/* Complex - Assemble Global Mass Matrix*/ 
(
    matlib_index      p,
    matlib_xm         Q,
    matlib_zv         phi,
    matlib_zm_sparse* M
)
{
    debug_enter( "poynomial degree: %d, "
                 "size of Q: %d-by-%d), "
                 "length of phi: %d",
                 p, Q.lenc, Q.lenr, phi.len);
    matlib_index N = (phi.len-1)/(Q.lenr-1);
    matlib_index nnz = N*(Q.lenc-1)+1;

    debug_body( "nr. finite-elements: %d, "
                "nr. of non-zero elements: %d", N, nnz);
    assert(Q.lenc==(p-1)*(p+4)/2+3);

    matlib_index dim = N*p+1;
    M->lenc   = dim;
    M->lenr   = dim;
    M->rowIn  = calloc( dim+1, sizeof(matlib_index));
    M->colIn  = calloc(   nnz, sizeof(matlib_index));
    M->elem_p = calloc(   nnz, sizeof(matlib_complex));

    matlib_zv q;
    matlib_create_zv( Q.lenc*N, &q, MATLIB_COL_VECT);

    fem1d_ZFLT( N, Q, phi, q);
    DEBUG_PRINT_ZV(q, "%s: ", "inner product");

    fem1d_ZCSRGMM(p, N, q, M->rowIn, M->colIn, M->elem_p);
    
    matlib_free((void*) q.elem_p);
    debug_exit("%s", "");

}
/*============================================================================*/

void fem1d_xm_nsparse_GMM
/* Double - Assemble Global Mass Matrix*/ 
(
    matlib_index       p,
    matlib_index       N,
    matlib_index       nsparse,
    matlib_xm          Q,
    matlib_xm*         phi,
    matlib_xm*         q,
    matlib_xm_nsparse* M,
    FEM1D_OP_GMM       op_enum /* option */ 
)
{
    debug_enter( "poynomial degree: %d, "
                 "nr. of fem-elements: %d, "
                 "nr. of sparse matrices: %d, "
                 "size of Q: %d-by-%d), ",
                 p, N, nsparse, Q.lenc, Q.lenr);

    assert(Q.lenc==(p-1)*(p+4)/2+3);
    if(op_enum == FEM1D_GMM_INIT)
    {
        matlib_create_xm( Q.lenr, nsparse, phi, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);
        matlib_create_xm( Q.lenc*N, nsparse, q, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);
        
        matlib_index dim = N*p+1;
        matlib_index nnz = N*(Q.lenc-1)+1;
        debug_body( "nr. of non-zero elements of M: %d", nnz);
        M->lenc   = dim;
        M->lenr   = dim;
        M->rowIn  = calloc( dim+1, sizeof(matlib_index));
        M->colIn  = calloc(   nnz, sizeof(matlib_index));
        M->elem_p = (matlib_real**)calloc( nsparse, sizeof(matlib_real *));

        for(matlib_index i = 0; i<nsparse; i++)
        {
            M->elem_p[i] = calloc( nnz, sizeof(double));
        }
    }
    else if((op_enum == FEM1D_GET_NZE_ONLY) || 
            (op_enum == FEM1D_GET_SPARSITY_NZE))
    {
        fem1d_XFLT2( N, Q, *phi, *q);
        fem1d_XCSRGMM2(p, N, *q, M->elem_p);
    }
    else if((op_enum == FEM1D_GET_SPARSITY_ONLY) ||
            (op_enum == FEM1D_GET_SPARSITY_NZE))
    {
        fem1d_GMMSparsity(p, N, M->rowIn, M->colIn);
    }
    else if(op_enum == FEM1D_GMM_FREE)
    {
        matlib_free(phi->elem_p);
        matlib_free(q->elem_p);
        matlib_free(M->colIn);
        matlib_free(M->rowIn);

        for(matlib_index i = 0; i<nsparse; i++)
        {
            matlib_free(M->elem_p[i]);
        }
    }

    debug_exit("%s", "");
}

void fem1d_zm_nsparse_GMM
/* ZNSparseM - Global Mass Matrix */ 
(
    matlib_index       p,
    matlib_index       N,
    matlib_index       nsparse,
    matlib_xm          Q,
    matlib_zm*         phi,
    matlib_zm*         q,
    matlib_zm_nsparse* M,
    FEM1D_OP_GMM       op_enum /* option */ 
)
/* 
 * This function assembles 'n' sparse matrices which have the same sparsity
 * structure.
 *
 *
 * */ 
{
    debug_enter( "poynomial degree: %d, "
                 "nr. of fem-elements: %d, "
                 "nr. of sparse matrices: %d, "
                 "size of Q: %d-by-%d)",
                 p, N, nsparse, Q.lenc, Q.lenr);

    assert(Q.lenc==(p-1)*(p+4)/2+3);
    if(op_enum == FEM1D_GMM_INIT)
    {
        debug_body("%s", "Initializing data structures");
        matlib_index dim = N*p+1;
        matlib_index nnz = N*(Q.lenc-1)+1;
        matlib_create_zm( (Q.lenr-1)*N+1, nsparse, phi, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);
        matlib_create_zm( Q.lenc*N,       nsparse,   q, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);
        
        debug_body( "nr. of non-zero elements of M: %d", nnz);
        M->lenc    = dim;
        M->lenr    = dim;
        M->nsparse = nsparse;
        M->rowIn   = calloc( dim+1, sizeof(matlib_index));
        M->colIn   = calloc(   nnz, sizeof(matlib_index));
        M->elem_p  = (matlib_complex**)calloc( nsparse, sizeof(matlib_complex *));

        for(matlib_index i = 0; i<nsparse; i++)
        {
            M->elem_p[i] = calloc( nnz, sizeof(matlib_complex));
        }
    }
    if((op_enum == FEM1D_GET_SPARSITY_ONLY) ||
            (op_enum == FEM1D_GET_SPARSITY_NZE))
    {
        fem1d_GMMSparsity(p, N, M->rowIn, M->colIn);
        debug_body( "nr. of non-zero elements: %d", 
                    M->rowIn[M->lenc]);
    }
    if((op_enum == FEM1D_GET_NZE_ONLY) || 
            (op_enum == FEM1D_GET_SPARSITY_NZE))
    {
        debug_body( "size of phi: %d-by-%d, "
                    "size of q: %d-by-%d, "
                    "size of M: %d-by-%d ",
                    phi->lenc, phi->lenr,
                    q->lenc, q->lenr, 
                    M->lenc, M->lenr);
        fem1d_ZFLT2( N, Q, *phi, *q);
        fem1d_ZCSRGMM2(p, N, *q, M->elem_p);
    }
    if(op_enum == FEM1D_GMM_FREE)
    {
        matlib_free(phi->elem_p);
        matlib_free(q->elem_p);
        matlib_free(M->colIn);
        matlib_free(M->rowIn);
        for(matlib_index i = 0; i<nsparse; i++)
        {
            matlib_free(M->elem_p[i]);
        }
    }

    debug_exit("%s", "");
}

/*============================================================================*/
