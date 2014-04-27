#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mex.h"
#include "matrix.h"
#include "legendre.h"

void mexFunction(int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
    double tol, *zeros, *quadW;
    Index p;
    
    p = (Index)mxGetScalar(prhs[0]);
    tol = mxGetScalar(prhs[1]);
    
    plhs[0] = mxCreateDoubleMatrix( p+1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix( p+1, 1, mxREAL);
    
    zeros = mxGetPr(plhs[0]);
    quadW = mxGetPr(plhs[1]);
    double gzeros[p];
        
    find_Gauss_points( p, tol, gzeros);
    find_LGL_points( p, tol, zeros, quadW, gzeros);
}