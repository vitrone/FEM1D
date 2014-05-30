/*============================================================================+/
 | mxfem1d_F2L.c
 | Handles complex as well as real inut data.
 | Complex mxArray is not compatible with C complex types. Hence, a complex 
 | number is treated as a 2-tuple.
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
        mexErrMsgTxt("Two inputs required.");
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

    matlib_int len = mxGetNumberOfElements(prhs[1]);
    if((len<3) || (((len-1) % p)!= 0))
    {
        mexErrMsgTxt("Length of input vector is incorrect.");
    }
    matlib_xv vbr = { .len    = len, 
                      .elem_p = mxGetPr(prhs[1]),
                      .type   = MATLIB_COL_VECT};

    matlib_index dim = (vbr.len-1)*(p+1)/p;
    

    if(mxIsComplex(prhs[1]))
    {
        matlib_xv vbi = { .len    = vbr.len, 
                          .elem_p = mxGetPr(prhs[1]),
                          .type   = MATLIB_COL_VECT};
        plhs[0] = mxCreateDoubleMatrix( dim, 1, mxCOMPLEX);
        matlib_xv ur = { .len    = dim, 
                         .elem_p = mxGetPr(plhs[0]),
                         .type   = MATLIB_COL_VECT};
        matlib_xv ui = { .len    = dim, 
                         .elem_p = mxGetPr(plhs[0]),
                         .type   = MATLIB_COL_VECT};

        fem1d_XF2L(p, vbr, ur);
        fem1d_XF2L(p, vbi, ur);
    }
    else
    {
        plhs[0] = mxCreateDoubleMatrix( dim, 1, mxREAL);
        matlib_xv ur = { .len    = dim,
                         .elem_p = mxGetPr(plhs[0]),
                         .type   = MATLIB_COL_VECT};
        fem1d_XF2L(p, vbr, ur);
    
    }
}
