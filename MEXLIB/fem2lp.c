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
    double *ur, *vr, *br;

    matlib_index p, N;    
    
    p = (matlib_index)mxGetScalar(prhs[0]);
    N = (matlib_index)mxGetScalar(prhs[1]);
    vr = mxGetPr(prhs[2]);
    br = mxGetPr(prhs[3]);
    void (*fp[9])(matlib_index, double*, double*, double*) = { fem1d_dshapefunc2lp_2,
                                                        fem1d_dshapefunc2lp_3, 
                                                        fem1d_dshapefunc2lp_4, 
                                                        fem1d_dshapefunc2lp_5, 
                                                        fem1d_dshapefunc2lp_6, 
                                                        fem1d_dshapefunc2lp_7, 
                                                        fem1d_dshapefunc2lp_8, 
                                                        fem1d_dshapefunc2lp_9, 
                                                        fem1d_dshapefunc2lp_10};
    
    matlib_index tmp1 = N*(p+1);
    bool s = (mxIsComplex(prhs[2]));
    bool t = (mxIsComplex(prhs[3]));
    if (p>10)
    {
        if (s||t)
        {
            double *bi, *vi;
            if (s) 
            {
                vi = mxGetPi(prhs[2]);
            }
            else 
            {
                double tmpvi[N+1]; 
                vi = &tmpvi[0];
            }
            if (t) 
            {
                bi = mxGetPi(prhs[3]);
            }
            else 
            {
                double tmpbi[(p-1)*N]; 
                bi = &tmpbi[0];
            }
            plhs[0] = mxCreateDoubleMatrix( tmp1, 1, mxCOMPLEX);
            ur = mxGetPr(plhs[0]);
            double* ui = mxGetPi(plhs[0]);
            fem1d_dshapefunc2lp( p, N, vr, br, ur);
            fem1d_dshapefunc2lp( p, N, vi, bi, ui);
        } 
        else 
        {
            plhs[0] = mxCreateDoubleMatrix( tmp1, 1, mxREAL);
            ur = mxGetPr(plhs[0]);
            fem1d_dshapefunc2lp( p, N, vr, br, ur);
        }
    } 
    else 
    {
        matlib_index tmp2 = p-2;
        if (s||t)
        {
            double *bi, *vi;
            if (s) 
            {
                vi = mxGetPi(prhs[2]);
            }
            else 
            {
                double tmpvi[N+1]; 
                vi = &tmpvi[0];
            }
            if (t) 
            {
                bi = mxGetPi(prhs[3]);
            }
            else 
            {
                double tmpbi[(p-1)*N]; 
                bi = &tmpbi[0];
            }
            plhs[0] = mxCreateDoubleMatrix( tmp1, 1, mxCOMPLEX);
            ur = mxGetPr(plhs[0]);
            double* ui = mxGetPi(plhs[0]);
            (*fp[tmp2])( N, vr, br, ur);
            (*fp[tmp2])( N, vi, bi, ui);
        } 
        else 
        {
            plhs[0] = mxCreateDoubleMatrix( tmp1, 1, mxREAL);
            ur = mxGetPr(plhs[0]);
            (*fp[tmp2])( N, vr, br, ur);
        }
    }    
}
