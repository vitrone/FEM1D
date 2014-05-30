/*============================================================================+/
 | mxfem1d_FLT.c
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

    if(nrhs!=3) 
    {
        mexErrMsgTxt("Three inputs required.");
    }
    if(nlhs!=1) 
    {
        mexErrMsgTxt("One output required.");
    } 
    matlib_index N;
    if(mxGetScalar(prhs[0])>0) 
    {
        N = (matlib_index)floor(mxGetScalar(prhs[0]));
    }
    else
    {
        mexErrMsgTxt("Number of finite-elements must be a positive integer.");
    } 
    
    if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) 
    {
        mexErrMsgTxt("Input matrix must be type double.");
    }
    if (mxGetNumberOfDimensions(prhs[1]) != 2)
    {
        mexErrMsgTxt("Transformation matrix must be two dimensional.");
    }
    
    if (mxGetNumberOfDimensions(prhs[2]) != 2)
    {
        mexErrMsgTxt("Third argument must be a vector or a matrix.");
    }

    matlib_xm FM = { .lenc = mxGetM(prhs[1]), 
                     .lenr = mxGetN(prhs[1]), 
                     .elem_p = mxGetPr(prhs[1]),
                     .op     = MATLIB_NO_TRANS,
                     .order  = MATLIB_COL_MAJOR };

    matlib_xm ur = { .lenc = mxGetM(prhs[2]), 
                     .lenr = mxGetN(prhs[2]), 
                     .elem_p = mxGetPr(prhs[2]), 
                     .op     = MATLIB_NO_TRANS,
                     .order  = MATLIB_COL_MAJOR};
    
    matlib_index dim = N*FM.lenc;
    matlib_index P = FM.lenr - 1;
    if(N*P != (ur.lenc-1)) 
    {
        mexErrMsgTxt("Dimension mismatch.");
    } 

    if(mxIsComplex(prhs[2]))
    {
        matlib_xm ui = { .lenc = ur.lenc, 
                         .lenr = ur.lenr, 
                         .elem_p = mxGetPi(prhs[2]),
                         .op     = MATLIB_NO_TRANS,
                         .order  = MATLIB_COL_MAJOR};

        plhs[0] = mxCreateDoubleMatrix( dim, ur.lenr, mxCOMPLEX);
        matlib_xm Ur = { .lenc = dim, 
                         .lenr = ur.lenr, 
                         .elem_p = mxGetPr(plhs[0]), 
                         .op     = MATLIB_NO_TRANS,
                         .order  = MATLIB_COL_MAJOR};

        matlib_xm Ui = { .lenc = dim, 
                         .lenr = ur.lenr, 
                         .elem_p = mxGetPi(plhs[0]), 
                         .op     = MATLIB_NO_TRANS,
                         .order  = MATLIB_COL_MAJOR};
        fem1d_XFLT2(N, FM, ur, Ur);
        fem1d_XFLT2(N, FM, ui, Ui);
    }
    else
    {
        plhs[0] = mxCreateDoubleMatrix( dim, ur.lenr, mxREAL);
        matlib_xm Ur = { .lenc = dim, 
                         .lenr = ur.lenr, 
                         .elem_p = mxGetPr(plhs[0]), 
                         .op     = MATLIB_NO_TRANS,
                         .order  = MATLIB_COL_MAJOR};

        fem1d_XFLT2(N, FM, ur, Ur);
    }
}
