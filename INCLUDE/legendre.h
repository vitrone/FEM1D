#ifndef LEGENDRE_H
#define LEGENDRE_H

/*============================================================================+/
 | Include all the dependencies
/+============================================================================*/
#include "basic.h"
#include "matlib.h"
#include "debug.h"
#include "ehandler.h"

/*============================================================================*/

double Legendre_LGL( double x, int n);

void CalcLP
( 
    int n, 
    double *LP, 
    double x, 
    double C[n-1],
    double D[n-1]
);

void find_Gauss_points
( 
    matlib_index n, 
    double tol, 
    double* zeros
);

void find_LGL_points
( 
    matlib_index n, 
    double tol, 
    double* zeros, 
    double* quadW,
    double *gzeros
);

void backward_transform_matrix
(
    const matlib_index   P,            
    const matlib_index   p,            
    const double* x,          
          double* pILTM
);
void forward_transform_matrix
(                                                                 
    const matlib_index   p,            
    const double* zeros,      
          double* pFLTM
);

void forward_transform_matrix2
(
    const matlib_index   p,                        
    const matlib_index   P,                        
    const double* zeros,                  
          double* pFLTM                   
);

/*======================================================================*/
/************************************************************************
 *                                                                      *
 * Column major version of the above functions for MATLAB               *
 *                                                                      *
 ************************************************************************/
void backward_transform_matrix_colmajor
( 
    matlib_index P, 
    matlib_index p, 
    double *x, 
    double *pILTM
);
void forward_transform_matrix_colmajor
( 
    matlib_index p, 
    double *zeros, 
    double *pFLTM
);
void forward_transform_matrix2_colmajor
( 
    matlib_index p, 
    matlib_index P, 
    double *zeros, 
    double *pFLTM                   /* (p+1)-by-(P+1) matrix, col. major*/
);
/*============================================================================*/

void legendre_LGLdataLT1
( 
    const matlib_index      p, 
    const double     tol,
          matlib_dv* zeros,
          matlib_dv* quadW
);
void legendre_LGLdataLT2
( 
    const matlib_index      p, 
    const double     tol,
          matlib_dv* zeros,
          matlib_dv* quadW,
          matlib_dm* FM,
          matlib_dm* IM
);
void legendre_LGLdataFM
( 
    const matlib_dv xi,
          matlib_dm FM
);
void legendre_LGLdataIM
( 
    const matlib_dv xi,
          matlib_dm IM
);
#endif
