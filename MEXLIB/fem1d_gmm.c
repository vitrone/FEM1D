#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mex.h"
#include "matrix.h"
#include "legendre.h"
#include "fem1d.h"


void mexFunction
(
          int      nlhs, 
          mxArray* plhs[],
          int      nrhs, 
    const mxArray *prhs[]
)
{
    double *M, *pot;
    double  *ugpmm, *lgpmm;

    Index p, N, P, nnz;    
    
    p = (Index)mxGetScalar(prhs[0]);
    N = (Index)mxGetScalar(prhs[1]);
    P = (Index)mxGetScalar(prhs[2]);
    pot = mxGetPr(prhs[3]);
    M = mxGetPr(prhs[4]);    
    nnz = (Index)mxGetScalar(prhs[5]);
    
    /* Number of non-zero elements nnz = N*p*(p+3)/2+1;*/

    plhs[0] = mxCreateDoubleMatrix( nnz, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix( nnz, 1, mxREAL);
    
    ugpmm = mxGetPr(plhs[0]);
    lgpmm = mxGetPr(plhs[1]);
    fem1d_gmm( p, N, P, M, pot, ugpmm, lgpmm);
    
}
