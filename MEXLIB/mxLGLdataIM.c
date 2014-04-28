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
        mexErrMsgIdAndTxt("LGLdata:IM:nrhs","Two inputs required.");
    }
    if(nlhs!=1) 
    {
        mexErrMsgIdAndTxt("LGLdata:IM:nlhs","One output required.");
    } 
    if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) 
    {
        mexErrMsgIdAndTxt("LGLdata:IM","Input vector must be type double.");
    }
    
    matlib_index p = (matlib_index)mxGetScalar(prhs[0]);

    double* xi = mxGetPr(prhs[1]);
    matlib_index P = (matlib_index)mxGetNumberOfElements(prhs[1])-1;
    
    
    plhs[0] = mxCreateDoubleMatrix( P+1, p+1, mxREAL);
    double *IM = mxGetPr(plhs[0]);
    backward_transform_matrix_colmajor(P, p, xi, IM);
}


