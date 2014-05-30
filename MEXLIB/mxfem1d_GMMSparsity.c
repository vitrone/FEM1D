/*============================================================================+/
 | mxfem1d_GMMSparsity.c
/+============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mex.h"
#include "matrix.h"
#include "fem1d.h"


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
        mexErrMsgIdAndTxt("FEM1D:GMMSparsity:nrhs","Two inputs required.");
    }
    if(nlhs!=2) 
    {
        mexErrMsgIdAndTxt("FEM1D:GMMSparsity:nlhs","Two outputs required.");
    } 
    
    
    if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) 
    {
        mexErrMsgIdAndTxt("FEM1D:quadM:quadW","Input vector must be type double.");
    }

    matlib_index p = (matlib_index)mxGetScalar(prhs[0]);
    matlib_index N = (matlib_index)mxGetScalar(prhs[1]);

    plhs[2] = mxCreateDoubleMatrix( 1, N*p+1, mxREAL);

}

                 
                 
                 
                 

