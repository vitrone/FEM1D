#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mex.h"
#include "matrix.h"
#include "legendre.h"
#include "fem1d.h"
#include "debug.h"


void mexFunction(int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
    double *sM, *IMi, *quadW;

    Index p, P;    
    
    p = (Index)mxGetScalar(prhs[0]);
    P = (Index)mxGetScalar(prhs[1]);    
    IMi = mxGetPr(prhs[2]);
    quadW = mxGetPr(prhs[3]);
    
    Index lenc_sM = P+1; 
    Index lenr_sM = 3+(p-1)*(p+4)/2;

    plhs[0] = mxCreateDoubleMatrix( lenc_sM, lenr_sM, mxREAL);
    sM = mxGetPr(plhs[0]);

    fem1d_superM_colMajor(p, P, IMi, quadW, sM);
}
