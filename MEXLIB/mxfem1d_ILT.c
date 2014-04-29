/*============================================================================+/
 | mxfem1d_ILT.c
 | Handles complex as well as real inut data.
 | Complex mxArray is not compatible with C complex types. Hence, a complex 
 | number is treated as a 2-tuple.
 | Usage: 
 |   u = mxfem1d_ILT(N, IM, U)
 |   N : number of finite elements
 |   U is the Legendre transform of u.
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
        mexErrMsgIdAndTxt("FEM1D:ILT:nrhs","Three inputs required.");
    }
    if(nlhs!=1) 
    {
        mexErrMsgIdAndTxt("FEM1D:ILT:nlhs","One output required.");
    } 
    
    matlib_index N = (matlib_index)mxGetScalar(prhs[0]);
    
    if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) 
    {
        mexErrMsgIdAndTxt("FEM1D:ILT:IM","Input matrix must be type double.");
    }

    matlib_xm IM = { .lenc = mxGetM(prhs[1]), 
                     .lenr = mxGetN(prhs[1]), 
                     .elem_p = mxGetPr(prhs[1]),
                     .op     = MATLIB_NO_TRANS,
                     .order  = MATLIB_COL_MAJOR };

    matlib_xm Ur = { .lenc = mxGetM(prhs[2]), 
                     .lenr = mxGetN(prhs[2]), 
                     .elem_p = mxGetPr(prhs[2]), 
                     .op     = MATLIB_NO_TRANS,
                     .order  = MATLIB_COL_MAJOR};
    
    matlib_index dim = N*(IM.lenc-1)+1;

    if(mxIsComplex(prhs[2]))
    {
        matlib_xm Ui = { .lenc = Ur.lenc, 
                         .lenr = Ur.lenr, 
                         .elem_p = mxGetPi(prhs[2]),
                         .op     = MATLIB_NO_TRANS,
                         .order  = MATLIB_COL_MAJOR};

        plhs[0] = mxCreateDoubleMatrix( dim, Ur.lenr, mxCOMPLEX);
        matlib_xm ur = { .lenc = dim, 
                         .lenr = Ur.lenr, 
                         .elem_p = mxGetPr(plhs[0]), 
                         .op     = MATLIB_NO_TRANS,
                         .order  = MATLIB_COL_MAJOR};

        matlib_xm ui = { .lenc = dim, 
                         .lenr = Ur.lenr, 
                         .elem_p = mxGetPi(plhs[0]), 
                         .op     = MATLIB_NO_TRANS,
                         .order  = MATLIB_COL_MAJOR};
        fem1d_XILT2(N, IM, Ur, ur);
        fem1d_XILT2(N, IM, Ui, ui);
    }
    else
    {
        plhs[0] = mxCreateDoubleMatrix( dim, Ur.lenr, mxREAL);
        matlib_xm ur = { .lenc = dim, 
                         .lenr = Ur.lenr, 
                         .elem_p = mxGetPr(plhs[0]), 
                         .op     = MATLIB_NO_TRANS,
                         .order  = MATLIB_COL_MAJOR};

        fem1d_XILT2(N, IM, Ur, ur);
    }
}
