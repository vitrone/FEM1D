/*============================================================================+/
 | mxLGLdataIM.c
/+============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mex.h"
#include "matrix.h"
#include "legendre.h"


void mexFunction
(
          int      nlhs, 
          mxArray* plhs[],
          int      nrhs, 
    const mxArray* prhs[]
)
{

    if(nrhs!=2) 
    {
        mexErrMsgTxt("Two inputs required.");
    }
    if(nlhs!=1) 
    {
        mexErrMsgTxt("One output required.");
    } 
    if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) 
    {
        mexErrMsgTxt("Input vector must be type double.");
    }
    
    matlib_index p;
    if(mxGetScalar(prhs[0])>1)
    {
        p = (matlib_index)floor(mxGetScalar(prhs[0]));
    }
    else
    {
        mexErrMsgTxt("Polynomial degree must be>1.");
    } 

    if (mxGetNumberOfDimensions(prhs[1]) != 2)
    {
        mexErrMsgTxt("Vector input required for LGL-points.");
    }
    double* xi = mxGetPr(prhs[1]);

    matlib_index P;
    if(mxGetNumberOfElements(prhs[1])>0)
    {
        P = (matlib_index)mxGetNumberOfElements(prhs[1])-1;
    }
    else
    {
        mexErrMsgTxt("Input vector is empty.");
    } 
    
    plhs[0] = mxCreateDoubleMatrix( P+1, p+1, mxREAL);
    double *IM = mxGetPr(plhs[0]);
    backward_transform_matrix_colmajor(P, p, xi, IM);
}


