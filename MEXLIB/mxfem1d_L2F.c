/*============================================================================+/
 | mxfem1d_L2F.c
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
        mexErrMsgIdAndTxt("FEM1D:L2F:nrhs","Two inputs required.");
    }
    if(nlhs!=1) 
    {
        mexErrMsgIdAndTxt("FEM1D:L2F:nlhs","One output required.");
    } 
    
    matlib_index p = (matlib_index)mxGetScalar(prhs[0]);

    matlib_xv ur = { .len    = mxGetN(prhs[1]), 
                     .elem_p = mxGetPr(prhs[1]),
                     .type   = MATLIB_COL_VECT};

    matlib_index dim = (ur.len)*p/(p+1)+1;
    

    if(mxIsComplex(prhs[1]))
    {
        matlib_xv ui = {  .len    = ur.len,
                          .elem_p = mxGetPr(prhs[1]),
                          .type   = MATLIB_COL_VECT};
        plhs[0] = mxCreateDoubleMatrix( dim, 1, mxCOMPLEX);
        matlib_xv vbr = { .len    = dim, 
                          .elem_p = mxGetPr(plhs[0]),
                          .type   = MATLIB_COL_VECT};
        matlib_xv vbi = { .len    = dim, 
                          .elem_p = mxGetPr(plhs[0]),
                          .type   = MATLIB_COL_VECT};

        fem1d_XL2F(p, ur, vbr);
        fem1d_XL2F(p, ui, vbi);
    }
    else
    {
        plhs[0] = mxCreateDoubleMatrix( dim, 1, mxREAL);
        matlib_xv vbr = { .len    = dim,
                          .elem_p = mxGetPr(plhs[0]),
                          .type   = MATLIB_COL_VECT};
        fem1d_XL2F(p, ur, vbr);
    
    }
}

