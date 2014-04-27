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
    double *ur, *vr, *norm, norm2r;

    Index p, N;    
    
    p = (Index)mxGetScalar(prhs[0]);
    N = (Index)mxGetScalar(prhs[1]);
    ur = mxGetPr(prhs[2]);
    double (*fp[10])(Index, double*) = { fem1d_dlp_snorm2_d_2, 
                                         fem1d_dlp_snorm2_d_3, 
                                         fem1d_dlp_snorm2_d_4, 
                                         fem1d_dlp_snorm2_d_5, 
                                         fem1d_dlp_snorm2_d_6, 
                                         fem1d_dlp_snorm2_d_7, 
                                         fem1d_dlp_snorm2_d_8, 
                                         fem1d_dlp_snorm2_d_9, 
                                         fem1d_dlp_snorm2_d_10};
    
    bool s = (mxIsComplex(prhs[2]));
    plhs[0] = mxCreateDoubleMatrix( 1, 1, mxREAL);
    norm = mxGetPr(plhs[0]);
    if (p>10)
    {
        if (s)
        {
            double* ui = mxGetPi(prhs[2]);
            norm2r = fem1d_dlp_snorm2_d( p, N, ur);
            double norm2i = fem1d_dlp_snorm2_d( p, N, ui);
            *norm = sqrt(norm2r+norm2i);
        } 
        else 
        {
            norm2r = fem1d_dlp_snorm2_d( p, N, ur);
            *norm = sqrt(norm2r);
        }
    } 
    else 
    {
        Index tmp = p-2;
        if (s)
        {
            double* ui = mxGetPi(prhs[2]);
            norm2r = (*fp[tmp])( N, ur);
            double norm2i = (*fp[tmp])( N, ui);
            *norm = sqrt(norm2r+norm2i);
        } 
        else 
        {
            norm2r = (*fp[tmp])( N, ur);
            *norm = sqrt(norm2r);
        }
    }    
}
