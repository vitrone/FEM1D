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
    void (*fp[9])(Index, double*, double*, double*) = { fem1d_dprjLP2FEM_ShapeFunc_2,
                                                        fem1d_dprjLP2FEM_ShapeFunc_3, 
                                                        fem1d_dprjLP2FEM_ShapeFunc_4, 
                                                        fem1d_dprjLP2FEM_ShapeFunc_5, 
                                                        fem1d_dprjLP2FEM_ShapeFunc_6, 
                                                        fem1d_dprjLP2FEM_ShapeFunc_7, 
                                                        fem1d_dprjLP2FEM_ShapeFunc_8, 
                                                        fem1d_dprjLP2FEM_ShapeFunc_9, 
                                                        fem1d_dprjLP2FEM_ShapeFunc_10};
    
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
            fem1d_dprjLP2FEM_ShapeFunc( p, N, ur, vr, br);
            fem1d_dprjLP2FEM_ShapeFunc( p, N, ui, vi, bi);
        } 
        else 
        {
            plhs[0] = mxCreateDoubleMatrix( tmp1, 1, mxREAL);
            plhs[1] = mxCreateDoubleMatrix( tmp2, 1, mxREAL);
            vr = mxGetPr(plhs[0]);
            br = mxGetPr(plhs[1]);
            fem1d_dprjLP2FEM_ShapeFunc( p, N, ur, vr, br);
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
