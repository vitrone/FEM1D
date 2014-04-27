#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mex.h"
#include "matrix.h"
#include "jacobi.h"

void mexFunction(int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
    double tol, frac_pow, *NUM, *DENOM;

    int order;
    order = (Index)mxGetScalar(prhs[0]);
    frac_pow = mxGetScalar(prhs[1]);
    tol = mxGetScalar(prhs[2]);
    
    plhs[0] = mxCreateDoubleMatrix(1, order+1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, order, mxREAL);
    NUM = mxGetPr(plhs[0]);
    DENOM = mxGetPr(plhs[1]);
    double poles[order];
    
    find_zeros( order, -frac_pow, frac_pow, poles, tol);
    diag_pade_pf( order, frac_pow, NUM, DENOM, poles);
}