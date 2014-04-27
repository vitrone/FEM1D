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

    Index p, N, P, nu, i;    
    
    p = (Index)mxGetScalar(prhs[0]);
    N = (Index)mxGetScalar(prhs[1]);
    P = (Index)mxGetScalar(prhs[2]);
    FM = mxGetPr(prhs[3]);
    ur = mxGetPr(prhs[4]);    
    nu = (Index)mxGetScalar(prhs[5]);    
    
    Index tmp = (p+1)*N;
    bool s = mxIsComplex(prhs[4]);
    if (s)
    {
        double* ui  =   mxGetPi(prhs[4]);
        plhs[0] = mxCreateDoubleMatrix( tmp, nu, mxCOMPLEX);
        vr = mxGetPr(plhs[0]);
        double* vi = mxGetPi(plhs[0]);
        fem_dflt_ColMajor2(p, N, P, FM, ur, vr, nu);
        fem_dflt_ColMajor2(p, N, P, FM, ui, vi, nu);
        
    } 
    else 
    {
        plhs[0] = mxCreateDoubleMatrix( tmp, nu, mxREAL);
        vr = mxGetPr(plhs[0]);
        fem_dflt_ColMajor2(p, N, P, FM, ur, vr, nu);
    }    
    
}
