#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mex.h"
#include "matrix.h"
#include "legendre.h"

void mexFunction(int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
    double *x;
    Index p, P;
    
    P = (Index) mxGetScalar(prhs[0]);
    p = (Index) mxGetScalar(prhs[1]);
    x = mxGetPr(prhs[2]);
    
    /* IM is (P+1)-by-(p+1) matrix */

    plhs[0] = mxCreateDoubleMatrix(P+1, p+1, mxREAL);
    double *pIM = mxGetPr(plhs[0]);
    backward_transform_matrix_colmajor(P, p, x, pIM);
}


