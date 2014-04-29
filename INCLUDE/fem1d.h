#ifndef FEM1D_H
#define FEM1D_H

/*============================================================================+/
 | Include all the dependencies
/+============================================================================*/
#include "basic.h"
#include "matlib.h"
#include "debug.h"
#include "ehandler.h"

/*============================================================================+/
 |DATA STRUCTURES AND ENUMS
/+============================================================================*/
/* NZE : Non-Zero Elements */ 
typedef enum
{
    FEM1D_GMM_INIT,
    FEM1D_GET_SPARSITY_ONLY,
    FEM1D_GET_NZE_ONLY,
    FEM1D_GET_SPARSITY_NZE,
    FEM1D_GMM_FREE

} FEM1D_OP_GMM; /* OPTIONS GMM */ 


/*============================================================================*/

void fem1d_ref2mesh
(
    matlib_xv    xi,
    matlib_index N,
    matlib_real  x_l,
    matlib_real  x_r,
    matlib_xv*   x
);

/*============================================================================+/
 | Legendre transformation routines
 | Naming convention: fem1d_<data-type><function descriptor>
 | X: Real data
 | Z: Complex data
/+============================================================================*/


void fem1d_XFLT
(
    const matlib_index N,
    const matlib_xm    FM,
          matlib_xv    u,
          matlib_xv    U
);

void fem1d_ZFLT
(
    const matlib_index N,
    const matlib_xm    FM,
          matlib_zv    u,
          matlib_zv    U
);

void fem1d_XILT
(
    const matlib_index N,
    const matlib_xm    IM,
          matlib_xv    U,
          matlib_xv    u
);

void fem1d_ZILT
(
    const matlib_index N,
    const matlib_xm    IM,
          matlib_zv    U,
          matlib_zv    u
);

void fem1d_XFLT2
(
    const matlib_index N,
    const matlib_xm    FM,
    const matlib_xm    u,
          matlib_xm    U
);
void fem1d_ZFLT2
(
    const matlib_index N,
    const matlib_xm FM,
          matlib_zm u,
          matlib_zm U
);
void fem1d_XILT2
(
    const matlib_index N,
    const matlib_xm IM,
    const matlib_xm U,
          matlib_xm u
);
void fem1d_ZILT2
(
    const matlib_index N,
    const matlib_xm IM,
          matlib_zm U,
          matlib_zm u
);

/*======================================================================*/

/* 
 * Transformation from Legendre basis to FEM-basis
 *
 * */
void fem1d_XF2L
(
    const matlib_index p, 
    const matlib_xv    vb,
          matlib_xv    u
);
void fem1d_ZF2L
(
    const matlib_index p, 
    const matlib_zv vb,
          matlib_zv u
);

void fem1d_XL2F
(
    const matlib_index p, 
    const matlib_xv u,
          matlib_xv vb
);
void fem1d_ZL2F
(
    const matlib_index p, 
    const matlib_zv u,
          matlib_zv vb
);

void fem1d_xshapefunc2lp
(
    matlib_index p, 
    matlib_index N, 
    matlib_real *v, 
    matlib_real *b,
    matlib_real *u
);
void fem1d_zshapefunc2lp
(
    matlib_index    p, 
    matlib_index    N, 
    matlib_complex* v,               
    matlib_complex* b,               
    matlib_complex* u                
);

void lp2fem1d_xshapefunc
(
    matlib_index p, 
    matlib_index N, 
    matlib_real *u, 
    matlib_real *v, 
    matlib_real *b
);
void lp2fem1d_zshapefunc
(
    matlib_index p, 
    matlib_index N, 
    matlib_complex *u, 
    matlib_complex *v, 
    matlib_complex *b
);

/* Printed with precision 0.20f.*/
#define _A01 -0.40824829046386307274
#define _A02 -0.31622776601683794118
#define _A03 0.40824829046386307274
#define _B03 -0.26726124191242439654
#define _A04 0.31622776601683794118
#define _B04 -0.23570226039551586683
#define _A05 0.26726124191242439654
#define _B05 -0.21320071635561041457
#define _A06 0.23570226039551586683
#define _B06 -0.19611613513818404453
#define _A07 0.21320071635561041457
#define _B07 -0.18257418583505535814
#define _A08 0.19611613513818404453
#define _B08 -0.17149858514250881925
#define _A09 0.18257418583505535814
#define _B09 -0.16222142113076254422
#define _A10 0.17149858514250881925
#define _A11 0.16222142113076254422

void fem1d_xshapefunc2lp_2 ( matlib_index N, matlib_real *v, matlib_real *b, matlib_real *u );
void fem1d_xshapefunc2lp_3 ( matlib_index N, matlib_real *v, matlib_real *b, matlib_real *u );
void fem1d_xshapefunc2lp_4 ( matlib_index N, matlib_real *v, matlib_real *b, matlib_real *u );
void fem1d_xshapefunc2lp_5 ( matlib_index N, matlib_real *v, matlib_real *b, matlib_real *u );
void fem1d_xshapefunc2lp_6 ( matlib_index N, matlib_real *v, matlib_real *b, matlib_real *u );
void fem1d_xshapefunc2lp_7 ( matlib_index N, matlib_real *v, matlib_real *b, matlib_real *u );
void fem1d_xshapefunc2lp_8 ( matlib_index N, matlib_real *v, matlib_real *b, matlib_real *u );
void fem1d_xshapefunc2lp_9 ( matlib_index N, matlib_real *v, matlib_real *b, matlib_real *u );
void fem1d_xshapefunc2lp_10( matlib_index N, matlib_real *v, matlib_real *b, matlib_real *u );

void fem1d_zshapefunc2lp_2 ( matlib_index N, matlib_complex *v, matlib_complex *b, matlib_complex *u );
void fem1d_zshapefunc2lp_3 ( matlib_index N, matlib_complex *v, matlib_complex *b, matlib_complex *u );
void fem1d_zshapefunc2lp_4 ( matlib_index N, matlib_complex *v, matlib_complex *b, matlib_complex *u );
void fem1d_zshapefunc2lp_5 ( matlib_index N, matlib_complex *v, matlib_complex *b, matlib_complex *u );
void fem1d_zshapefunc2lp_6 ( matlib_index N, matlib_complex *v, matlib_complex *b, matlib_complex *u );
void fem1d_zshapefunc2lp_7 ( matlib_index N, matlib_complex *v, matlib_complex *b, matlib_complex *u );
void fem1d_zshapefunc2lp_8 ( matlib_index N, matlib_complex *v, matlib_complex *b, matlib_complex *u );
void fem1d_zshapefunc2lp_9 ( matlib_index N, matlib_complex *v, matlib_complex *b, matlib_complex *u );
void fem1d_zshapefunc2lp_10( matlib_index N, matlib_complex *v, matlib_complex *b, matlib_complex *u );

/* Printed with precision 0.20f.*/
#define _C10 6.16441400296897601407
#define _C09 5.83095189484530074253
#define _C08 5.47722557505166118830
#define _C07 5.09901951359278449161
#define _C06 4.69041575982342973106
#define _C05 4.24264068711928477029
#define _C04 3.74165738677394132949
#define _C03 3.16227766016837952279
#define _C02 2.44948974278317788134
#define _D08 0.88852331663863859390
#define _D07 0.87447463219520615851
#define _D06 0.85634883857767529758
#define _D05 0.83205029433784361004
#define _D04 0.79772403521746559907
#define _D03 0.74535599249993000903
#define _D02 0.65465367070797708671
#define _D01 0.31622776601683794118
#define _D00 0.40824829046386307274
void lp2fem1d_xshapefunc_2 ( matlib_index N, matlib_real *u, matlib_real *v, matlib_real *b );
void lp2fem1d_xshapefunc_3 ( matlib_index N, matlib_real *u, matlib_real *v, matlib_real *b );
void lp2fem1d_xshapefunc_4 ( matlib_index N, matlib_real *u, matlib_real *v, matlib_real *b );
void lp2fem1d_xshapefunc_5 ( matlib_index N, matlib_real *u, matlib_real *v, matlib_real *b );
void lp2fem1d_xshapefunc_6 ( matlib_index N, matlib_real *u, matlib_real *v, matlib_real *b );
void lp2fem1d_xshapefunc_7 ( matlib_index N, matlib_real *u, matlib_real *v, matlib_real *b );
void lp2fem1d_xshapefunc_8 ( matlib_index N, matlib_real *u, matlib_real *v, matlib_real *b );
void lp2fem1d_xshapefunc_9 ( matlib_index N, matlib_real *u, matlib_real *v, matlib_real *b );
void lp2fem1d_xshapefunc_10( matlib_index N, matlib_real *u, matlib_real *v, matlib_real *b );

void lp2fem1d_zshapefunc_2 ( matlib_index N, matlib_complex *u, matlib_complex *v, matlib_complex *b );
void lp2fem1d_zshapefunc_3 ( matlib_index N, matlib_complex *u, matlib_complex *v, matlib_complex *b );
void lp2fem1d_zshapefunc_4 ( matlib_index N, matlib_complex *u, matlib_complex *v, matlib_complex *b );
void lp2fem1d_zshapefunc_5 ( matlib_index N, matlib_complex *u, matlib_complex *v, matlib_complex *b );
void lp2fem1d_zshapefunc_6 ( matlib_index N, matlib_complex *u, matlib_complex *v, matlib_complex *b );
void lp2fem1d_zshapefunc_7 ( matlib_index N, matlib_complex *u, matlib_complex *v, matlib_complex *b );
void lp2fem1d_zshapefunc_8 ( matlib_index N, matlib_complex *u, matlib_complex *v, matlib_complex *b );
void lp2fem1d_zshapefunc_9 ( matlib_index N, matlib_complex *u, matlib_complex *v, matlib_complex *b );
void lp2fem1d_zshapefunc_10( matlib_index N, matlib_complex *u, matlib_complex *v, matlib_complex *b );
/* Projection subroutines                                               */
/*======================================================================*/
/* Printed with precision 0.20f.*/
#define _E00 0.16329931618554521799
#define _F00 -0.81649658092772614548
#define _E01 0.09035079029052513200
#define _F01 -0.21081851067789195153
#define _E02 0.05939138709164986513
#define _F02 -0.10690449676496975584
#define _E03 0.04285495643554834005
#define _F03 -0.06734350297014739251
#define _E04 0.03280011020855545106
#define _F04 -0.04737793696791342546
#define _E05 0.02614881801842454043
#define _F05 -0.03565747911603345949
#define _E06 0.02147931598059474659
#define _F06 -0.02808833628231621055
#define _E07 0.01805248264657987461
#define _F07 -0.02286647801900117474
#define _E08 0.01544965915531071841
#define _F08 -0.01908487307420735773

void fem1d_XPrjL2F
(
    matlib_index p,
    matlib_xv    u,
    matlib_xv    Pvb
);

void fem1d_xprjLP2FEM_ShapeFunc
(
    matlib_index p, 
    matlib_index N, 
    matlib_real* u, 
    matlib_real* Pv, 
    matlib_real* Pb
);

void fem1d_ZPrjL2F
(
    matlib_index p,
    matlib_zv    u,
    matlib_zv    Pvb
);
void fem1d_zprjLP2FEM_ShapeFunc
(
    matlib_index    p, 
    matlib_index    N, 
    matlib_complex* u, 
    matlib_complex* Pv,              
    matlib_complex* Pb               
);


void fem1d_xprjLP2FEM_ShapeFunc_2 ( matlib_index N, matlib_real *u, matlib_real *Pv, matlib_real *Pb);
void fem1d_xprjLP2FEM_ShapeFunc_3 ( matlib_index N, matlib_real *u, matlib_real *Pv, matlib_real *Pb);
void fem1d_xprjLP2FEM_ShapeFunc_4 ( matlib_index N, matlib_real *u, matlib_real *Pv, matlib_real *Pb);
void fem1d_xprjLP2FEM_ShapeFunc_5 ( matlib_index N, matlib_real *u, matlib_real *Pv, matlib_real *Pb);
void fem1d_xprjLP2FEM_ShapeFunc_6 ( matlib_index N, matlib_real *u, matlib_real *Pv, matlib_real *Pb);
void fem1d_xprjLP2FEM_ShapeFunc_7 ( matlib_index N, matlib_real *u, matlib_real *Pv, matlib_real *Pb);
void fem1d_xprjLP2FEM_ShapeFunc_8 ( matlib_index N, matlib_real *u, matlib_real *Pv, matlib_real *Pb);
void fem1d_xprjLP2FEM_ShapeFunc_9 ( matlib_index N, matlib_real *u, matlib_real *Pv, matlib_real *Pb);
void fem1d_xprjLP2FEM_ShapeFunc_10( matlib_index N, matlib_real *u, matlib_real *Pv, matlib_real *Pb);

void fem1d_zprjLP2FEM_ShapeFunc_2 ( matlib_index N, matlib_complex *u, matlib_complex *Pv, matlib_complex *Pb);
void fem1d_zprjLP2FEM_ShapeFunc_3 ( matlib_index N, matlib_complex *u, matlib_complex *Pv, matlib_complex *Pb);
void fem1d_zprjLP2FEM_ShapeFunc_4 ( matlib_index N, matlib_complex *u, matlib_complex *Pv, matlib_complex *Pb);
void fem1d_zprjLP2FEM_ShapeFunc_5 ( matlib_index N, matlib_complex *u, matlib_complex *Pv, matlib_complex *Pb);
void fem1d_zprjLP2FEM_ShapeFunc_6 ( matlib_index N, matlib_complex *u, matlib_complex *Pv, matlib_complex *Pb);
void fem1d_zprjLP2FEM_ShapeFunc_7 ( matlib_index N, matlib_complex *u, matlib_complex *Pv, matlib_complex *Pb);
void fem1d_zprjLP2FEM_ShapeFunc_8 ( matlib_index N, matlib_complex *u, matlib_complex *Pv, matlib_complex *Pb);
void fem1d_zprjLP2FEM_ShapeFunc_9 ( matlib_index N, matlib_complex *u, matlib_complex *Pv, matlib_complex *Pb);
void fem1d_zprjLP2FEM_ShapeFunc_10( matlib_index N, matlib_complex *u, matlib_complex *Pv, matlib_complex *Pb);
/* Squared L2 norm                                                      */ 
/*======================================================================*/
/* Printed with precision 0.20f.*/
#define _N0 2.00000000000000000000
#define _N1 0.66666666666666662966
#define _N2 0.40000000000000002220
#define _N3 0.28571428571428569843
#define _N4 0.22222222222222220989
#define _N5 0.18181818181818182323
#define _N6 0.15384615384615385469
#define _N7 0.13333333333333333148
#define _N8 0.11764705882352941013
#define _N9 0.10526315789473683626
#define _N10 0.09523809523809523281

matlib_real fem1d_XNorm2
(
    matlib_index p,
    matlib_index N,
    matlib_xv    u
);
matlib_real fem1d_ZNorm2
(
    matlib_index p,
    matlib_index N,
    matlib_zv u
);

matlib_real fem1d_xlp_snorm2_d(matlib_index p, matlib_index N, matlib_real *u);
matlib_real fem1d_zlp_snorm2_d(matlib_index p, matlib_index N, matlib_complex *u);

matlib_real fem1d_xlp_snorm2_d_2 ( matlib_index N, matlib_real *u);
matlib_real fem1d_xlp_snorm2_d_3 ( matlib_index N, matlib_real *u);
matlib_real fem1d_xlp_snorm2_d_4 ( matlib_index N, matlib_real *u);
matlib_real fem1d_xlp_snorm2_d_5 ( matlib_index N, matlib_real *u);
matlib_real fem1d_xlp_snorm2_d_6 ( matlib_index N, matlib_real *u);
matlib_real fem1d_xlp_snorm2_d_7 ( matlib_index N, matlib_real *u);
matlib_real fem1d_xlp_snorm2_d_8 ( matlib_index N, matlib_real *u);
matlib_real fem1d_xlp_snorm2_d_9 ( matlib_index N, matlib_real *u);
matlib_real fem1d_xlp_snorm2_d_10( matlib_index N, matlib_real *u);

matlib_real fem1d_zlp_snorm2_d_2 ( matlib_index N, matlib_complex *u);
matlib_real fem1d_zlp_snorm2_d_3 ( matlib_index N, matlib_complex *u);
matlib_real fem1d_zlp_snorm2_d_4 ( matlib_index N, matlib_complex *u);
matlib_real fem1d_zlp_snorm2_d_5 ( matlib_index N, matlib_complex *u);
matlib_real fem1d_zlp_snorm2_d_6 ( matlib_index N, matlib_complex *u);
matlib_real fem1d_zlp_snorm2_d_7 ( matlib_index N, matlib_complex *u);
matlib_real fem1d_zlp_snorm2_d_8 ( matlib_index N, matlib_complex *u);
matlib_real fem1d_zlp_snorm2_d_9 ( matlib_index N, matlib_complex *u);
matlib_real fem1d_zlp_snorm2_d_10( matlib_index N, matlib_complex *u);
/*======================================================================*/

matlib_real fem1d_xdot
( 
    matlib_index l,                        /* Length of the vectors            */
    matlib_real *x,                      /* vector of length l               */
    matlib_real *y
);
matlib_complex fem1d_zdot
( 
    matlib_index    l,
    matlib_complex* x,
    matlib_complex* y
);
matlib_real fem1d_xdot2 ( matlib_index l, matlib_real *u, matlib_real *v);
matlib_real fem1d_xdot3 ( matlib_index l, matlib_real *u, matlib_real *v);
matlib_real fem1d_xdot4 ( matlib_index l, matlib_real *u, matlib_real *v);
matlib_real fem1d_xdot5 ( matlib_index l, matlib_real *u, matlib_real *v);
matlib_real fem1d_xdot6 ( matlib_index l, matlib_real *u, matlib_real *v);
matlib_real fem1d_xdot7 ( matlib_index l, matlib_real *u, matlib_real *v);
matlib_real fem1d_xdot8 ( matlib_index l, matlib_real *u, matlib_real *v);
matlib_real fem1d_xdot9 ( matlib_index l, matlib_real *u, matlib_real *v);
matlib_real fem1d_xdot10( matlib_index l, matlib_real *u, matlib_real *v);
matlib_real fem1d_xdot11( matlib_index l, matlib_real *u, matlib_real *v);
matlib_real fem1d_xdot12( matlib_index l, matlib_real *u, matlib_real *v);
matlib_real fem1d_xdot13( matlib_index l, matlib_real *u, matlib_real *v);
matlib_real fem1d_xdot14( matlib_index l, matlib_real *u, matlib_real *v);
matlib_real fem1d_xdot15( matlib_index l, matlib_real *u, matlib_real *v);
matlib_real fem1d_xdot16( matlib_index l, matlib_real *u, matlib_real *v);

matlib_complex fem1d_zdot2( matlib_index l, matlib_complex *u, matlib_complex *v);
matlib_complex fem1d_zdot3( matlib_index l, matlib_complex *u, matlib_complex *v);
matlib_complex fem1d_zdot4( matlib_index l, matlib_complex *u, matlib_complex *v);
matlib_complex fem1d_zdot5( matlib_index l, matlib_complex *u, matlib_complex *v);
matlib_complex fem1d_zdot6( matlib_index l, matlib_complex *u, matlib_complex *v);
matlib_complex fem1d_zdot7( matlib_index l, matlib_complex *u, matlib_complex *v);
matlib_complex fem1d_zdot8( matlib_index l, matlib_complex *u, matlib_complex *v);
matlib_complex fem1d_zdot9( matlib_index l, matlib_complex *u, matlib_complex *v);
matlib_complex fem1d_zdot10( matlib_index l, matlib_complex *u, matlib_complex *v);
matlib_complex fem1d_zdot11( matlib_index l, matlib_complex *u, matlib_complex *v);
matlib_complex fem1d_zdot12( matlib_index l, matlib_complex *u, matlib_complex *v);
matlib_complex fem1d_zdot13( matlib_index l, matlib_complex *u, matlib_complex *v);
matlib_complex fem1d_zdot14( matlib_index l, matlib_complex *u, matlib_complex *v);
matlib_complex fem1d_zdot15( matlib_index l, matlib_complex *u, matlib_complex *v);
matlib_complex fem1d_zdot16( matlib_index l, matlib_complex *u, matlib_complex *v);

/*============================================================================+/
 | Building Global Mass Matrix 
/+============================================================================*/
void fem1d_quadM
(
    matlib_xv  quadW,
    matlib_xm  IM,
    matlib_xm* Q
);
void fem1d_MEMI
(
    matlib_index p,
    matlib_xv*   q
);
void fem1d_GMMSparsity
/* Determine the sparsity structure in CSR format */ 
(
    matlib_index  p,
    matlib_index  N,                        
    matlib_index* row,                     
    matlib_index* col                     
);

void fem1d_XCSRGMM
/* Double CSR - Assemble Global Mass Matrix*/ 
(
    matlib_index  p,
    matlib_index  N,                        
    matlib_xv     q,
    matlib_index* row,                     
    matlib_index* col,                     
    matlib_real*  ugpmm                   
                                    
);

void fem1d_XCSRGMM2
/* Double CSR - Assemble Global Mass Matrix*/ 
(
    matlib_index  p,
    matlib_index  N,                        
    matlib_xm     q,
    matlib_real** ugpmm_p              
                                    
);

void fem1d_ZCSRGMM
/* Complex CSR - Assemble Global Mass Matrix*/ 
(
    matlib_index    p,
    matlib_index    N,                        
    matlib_zv       q,
    matlib_index*   row,                     
    matlib_index*   col,                     
    matlib_complex* ugpmm                   
                                    
);

void fem1d_ZCSRGMM2
/* Double CSR - Assemble Global Mass Matrix*/ 
(
    matlib_index     p,
    matlib_index     N,                        
    matlib_zm        q,
    matlib_complex** ugpmm_p                
);


void fem1d_xm_sparse_GMM
/* Complex - Assemble Global Mass Matrix*/ 
(
    matlib_index      p,
    matlib_xm         Q,
    matlib_xv         phi,
    matlib_xm_sparse* M
);

void fem1d_zm_sparse_GMM
/* Complex - Assemble Global Mass Matrix*/ 
(
    matlib_index      p,
    matlib_xm         Q,
    matlib_zv         phi,
    matlib_zm_sparse* M
);

void fem1d_xm_nsparse_GMM
/* Real - Assemble Global Mass Matrix*/ 
(
    matlib_index       p,
    matlib_index       N,
    matlib_index       nsparse,
    matlib_xm          Q,
    matlib_xm*         phi,
    matlib_xm*         q,
    matlib_xm_nsparse* M,
    FEM1D_OP_GMM       op_enum /* option */ 
);

void fem1d_zm_nsparse_GMM
/* Complex - Global Mass Matrix */ 
(
    matlib_index       p,
    matlib_index       N,
    matlib_index       nsparse,
    matlib_xm          Q,
    matlib_zm*         phi,
    matlib_zm*         q,
    matlib_zm_nsparse* M,
    FEM1D_OP_GMM       op_enum /* option */ 
);

#endif
