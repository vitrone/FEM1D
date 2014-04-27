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
    
    p = (Index) mxGetScalar(prhs[0]);
    P = (Index) mxGetScalar(prhs[1]);
    x = mxGetPr(prhs[2]);
    
    plhs[0] = mxCreateDoubleMatrix( p+1, P+1, mxREAL);
    double *pFM = mxGetPr(plhs[0]);
    forward_transform_matrix2_colmajor( p, P, x, pFM);
}


