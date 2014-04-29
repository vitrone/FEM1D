/*============================================================================+/
 | mxfem1d_PrjL2F.c
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
        mexErrMsgIdAndTxt("FEM1D:PrjL2F:nrhs","Two inputs required.");
    }
    if(nlhs!=1) 
    {
        mexErrMsgIdAndTxt("FEM1D:PrjL2F:nlhs","One output required.");
    } 
    
    matlib_index p = (matlib_index)mxGetScalar(prhs[0]);

    matlib_dv ur = { .len    = mxGetN(prhs[1]), 
                     .elem_p = mxGetPr(prhs[1]),
                     .type   = MATLIB_COL_VECT};

    matlib_index dim = (ur.len)*p/(p+1)+1;

    if(mxIsComplex(prhs[1]))
    {
        matlib_dv ui = {  .len    = ur.len,
                          .elem_p = mxGetPr(prhs[1]),
                          .type   = MATLIB_COL_VECT};
        plhs[0] = mxCreateDoubleMatrix( dim, 1, mxCOMPLEX);
        matlib_dv Pvbr = { .len    = dim, 
                           .elem_p = mxGetPr(plhs[0]),
                           .type   = MATLIB_COL_VECT};
        matlib_dv Pvbi = { .len    = dim, 
                           .elem_p = mxGetPr(plhs[0]),
                           .type   = MATLIB_COL_VECT};

        fem1d_DPrjL2F(p, ur, Pvbr);
        fem1d_DPrjL2F(p, ui, Pvbi);
    }
    else
    {
        plhs[0] = mxCreateDoubleMatrix( dim, 1, mxREAL);
        matlib_dv Pvbr = { .len    = dim,
                           .elem_p = mxGetPr(plhs[0]),
                           .type   = MATLIB_COL_VECT};
        fem1d_DPrjL2F(p, ur, Pvbr);
    
    }
    

}
