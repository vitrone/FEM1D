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
        mexErrMsgIdAndTxt("FEM1D:F2L:nrhs","Two inputs required.");
    }
    if(nlhs!=1) 
    {
        mexErrMsgIdAndTxt("FEM1D:F2L:nlhs","One output required.");
    } 
    
    matlib_index p = (matlib_index)mxGetScalar(prhs[0]);

    matlib_dv vbr = { .len    = mxGetN(prhs[1]), 
                      .elem_p = mxGetPr(prhs[1]),
                      .type   = MATLIB_COL_VECT};

    matlib_index dim = (vbr.len-1)*(p+1)/p;
    

    if(mxIsComplex(prhs[2]))
    {
        matlib_dv vbi = { .len    = vbr.len, 
                          .elem_p = mxGetPr(prhs[1]),
                          .type   = MATLIB_COL_VECT};
        plhs[0] = mxCreateDoubleMatrix( dim, 1, mxCOMPLEX);
        matlib_dv ur = { .len    = dim, 
                         .elem_p = mxGetPr(plhs[0]),
                         .type   = MATLIB_COL_VECT};
        matlib_dv ui = { .len    = dim, 
                         .elem_p = mxGetPr(plhs[0]),
                         .type   = MATLIB_COL_VECT};

        fem1d_DF2L(p, vbr, ur);
        fem1d_DF2L(p, vbi, ur);
    }
    else
    {
        plhs[0] = mxCreateDoubleMatrix( dim, 1, mxREAL);
        matlib_dv ur = { .len    = dim,
                         .elem_p = mxGetPr(plhs[0]),
                         .type   = MATLIB_COL_VECT};
        fem1d_DF2L(p, vbr, ur);
    
    }
}
