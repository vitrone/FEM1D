/*============================================================================+/
 | mxLGLdataLT1.c
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
        mexErrMsgIdAndTxt("LGLdataLT1:nrhs","Two inputs required.");
    }
    if(nlhs!=2) 
    {
        mexErrMsgIdAndTxt("LGLdataLT1:nlhs","Two outputs required.");
    } 
    
    matlib_index p = (matlib_index)mxGetScalar(prhs[0]);
    double tol     = (double)mxGetScalar(prhs[1]);

    plhs[0] = mxCreateDoubleMatrix( p+1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix( p+1, 1, mxREAL);

    double *zeros = mxGetPr(plhs[0]);
    double *quadW = mxGetPr(plhs[1]);
    double gzeros[p];

    find_Gauss_points( p, tol, gzeros);
    find_LGL_points( p, tol, zeros, quadW, gzeros);

}
