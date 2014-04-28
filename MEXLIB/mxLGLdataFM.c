/*============================================================================+/
 | mxLGLdataFM.c
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
        mexErrMsgIdAndTxt("LGLdata:FM:nrhs","Two inputs required.");
    }
    if(nlhs!=1) 
    {
        mexErrMsgIdAndTxt("LGLdata:FM:nlhs","One output required.");
    } 
    if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) 
    {
        mexErrMsgIdAndTxt("LGLdata:FM","Input vector must be type double.");
    }
    
    matlib_index p = (matlib_index)mxGetScalar(prhs[0]);

    double* xi = mxGetPr(prhs[1]);
    matlib_index P = (matlib_index)mxGetNumberOfElements(prhs[1])-1;
    
    
    plhs[0] = mxCreateDoubleMatrix( p+1, P+1, mxREAL);
    double *FM = mxGetPr(plhs[0]);
    forward_transform_matrix2_colmajor( p, P, xi, FM);
}


