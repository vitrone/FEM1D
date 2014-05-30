/*============================================================================+/
 | mxfem1d_Norm2.c
 | Handles complex as well as real inut data.
 | Complex mxArray is not compatible with C complex types. Hence, a complex 
 | number is treated as a 2-tuple.
 | 
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

    if(nrhs!=3) 
    {
        mexErrMsgTxt("Three inputs required.");
    }
    if(nlhs!=1) 
    {
        mexErrMsgTxt("One output required.");
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
    matlib_index N;
    if(mxGetScalar(prhs[1])>0) 
    {
        N = (matlib_index)floor(mxGetScalar(prhs[1]));
    }
    else
    {
        mexErrMsgTxt("Number of finite-elements must be a positive integer.");
    } 
    
    matlib_index dim = N*(p+1);
    if(mxGetNumberOfElements(prhs[2])!= dim)
    {
        mexErrMsgTxt("Size of the input vector is inconsistent.");

    }
    matlib_xv ur = { .len = dim, 
                     .elem_p = mxGetPr(prhs[2]),
                     .type   = MATLIB_COL_VECT };

    double rnorm2, inorm2 = 0;
    rnorm2 = fem1d_XNorm2(p, N, ur);

    if(mxIsComplex(prhs[2]))
    {

        matlib_xv ui = { .len = dim, 
                         .elem_p = mxGetPi(prhs[2]),
                         .type   = MATLIB_COL_VECT };
        inorm2 = fem1d_XNorm2(p, N, ui);
    }
    plhs[0] = mxCreateDoubleScalar(sqrt(rnorm2*rnorm2 + inorm2*inorm2));
}
