#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mex.h"
#include "matrix.h"
#include "legendre.h"
#include "fem1d.h"


void mexFunction(int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
    double *FM, *ur, *vr;

    Index p, N, P;    
    
    p = (Index)mxGetScalar(prhs[0]);
    N = (Index)mxGetScalar(prhs[1]);
    P = (Index)mxGetScalar(prhs[2]);
    FM = mxGetPr(prhs[3]);
    ur = mxGetPr(prhs[4]);    
    
    Index tmp = (p+1)*N;
    bool s = mxIsComplex(prhs[4]);
    if (s)
    {
        double* ui  =   mxGetPi(prhs[4]);
        plhs[0] = mxCreateDoubleMatrix( tmp, 1, mxCOMPLEX);
        vr = mxGetPr(plhs[0]);
        double* vi = mxGetPi(plhs[0]);
        fem_dflt_ColMajor(p, N, P, FM, ur, vr);
        fem_dflt_ColMajor(p, N, P, FM, ui, vi);
    } 
    else 
    {
        plhs[0] = mxCreateDoubleMatrix( tmp, 1, mxREAL);
        vr = mxGetPr(plhs[0]);
        fem_dflt_ColMajor(p, N, P, FM, ur, vr);
    }    
    
}
