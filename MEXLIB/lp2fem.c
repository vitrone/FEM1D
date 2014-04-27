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

    Index p, N;    
    
    p = (Index)mxGetScalar(prhs[0]);
    N = (Index)mxGetScalar(prhs[1]);
    ur = mxGetPr(prhs[2]);
    void (*fp[9])(Index, double*, double*, double*) = { lp2fem1d_dshapefunc_2,
                                                        lp2fem1d_dshapefunc_3, 
                                                        lp2fem1d_dshapefunc_4, 
                                                        lp2fem1d_dshapefunc_5, 
                                                        lp2fem1d_dshapefunc_6, 
                                                        lp2fem1d_dshapefunc_7, 
                                                        lp2fem1d_dshapefunc_8, 
                                                        lp2fem1d_dshapefunc_9, 
                                                        lp2fem1d_dshapefunc_10 };
    
    Index tmp1 = N+1, tmp2 = (p-1)*N;
    bool s = mxIsComplex(prhs[2]);
    if (p>10)
    {
        if (s)
        {
            double* ui = mxGetPi(prhs[2]);
            plhs[0] = mxCreateDoubleMatrix( tmp1, 1, mxCOMPLEX);
            plhs[1] = mxCreateDoubleMatrix( tmp2, 1, mxCOMPLEX);
            vr = mxGetPr(plhs[0]);
            br = mxGetPr(plhs[1]);
            double* vi = mxGetPi(plhs[0]);
            double* bi = mxGetPi(plhs[1]);
            lp2fem1d_dshapefunc( p, N, ur, vr, br);
            lp2fem1d_dshapefunc( p, N, ui, vi, bi);
        } 
        else 
        {
            plhs[0] = mxCreateDoubleMatrix( tmp1, 1, mxREAL);
            plhs[1] = mxCreateDoubleMatrix( tmp2, 1, mxREAL);
            vr = mxGetPr(plhs[0]);
            br = mxGetPr(plhs[1]);
            lp2fem1d_dshapefunc( p, N, ur, vr, br);
        }
    } 
    else 
    {
        Index tmp3 = p-2;
        if (s)
        {
            double* ui = mxGetPi(prhs[2]);
            plhs[0] = mxCreateDoubleMatrix( tmp1, 1, mxCOMPLEX);
            plhs[1] = mxCreateDoubleMatrix( tmp2, 1, mxCOMPLEX);
            vr = mxGetPr(plhs[0]);
            br = mxGetPr(plhs[1]);
            double* vi = mxGetPi(plhs[0]);
            double* bi = mxGetPi(plhs[1]);
            (*fp[tmp3])( N, ur, vr, br);
            (*fp[tmp3])( N, ui, vi, bi);
        } 
        else 
        {
            plhs[0] = mxCreateDoubleMatrix( tmp1, 1, mxREAL);
            plhs[1] = mxCreateDoubleMatrix( tmp2, 1, mxREAL);
            vr = mxGetPr(plhs[0]);
            br = mxGetPr(plhs[1]);
            (*fp[tmp3])( N, ur, vr, br);
        }
    }    
}
