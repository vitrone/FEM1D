#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mex.h"
#include "matrix.h"
#include "legendre.h"
#include "fem1d.h"
#include "debug.h"


void mexFunction
(
          int      nlhs, 
          mxArray* plhs[],
          int      nrhs, 
    const mxArray *prhs[]
)
{
    double *M, *pot;

    Index p, N, P, nnz1, nnz2;    
    
    p = (Index)mxGetScalar(prhs[0]);
    N = (Index)mxGetScalar(prhs[1]);
    P = (Index)mxGetScalar(prhs[2]);
    pot = mxGetPr(prhs[3]);
    M = mxGetPr(prhs[4]); 

    debug_enter("degree of polynomial: %d, nr. LGL points: %d", p, P+1);
    
    nnz1 = N*p*(p+3)/2+1;
    nnz2 = N*p*(p+2)+1;
    
    /* Number of non-zero elements nnz = N*p*(p+3)/2+1;*/
    Index dim = N*p+1;
    Index row[dim], col[nnz1];
    double ugpmm[nnz1];

    plhs[0] = mxCreateNumericMatrix( nnz2, 1, mxUINT64_CLASS, mxREAL);
    plhs[1] = mxCreateNumericMatrix( nnz2, 1, mxUINT64_CLASS, mxREAL);
    plhs[2] = mxCreateDoubleMatrix( nnz2, 1, mxREAL);
    
    Index* row1 = (Index*)mxGetPr(plhs[0]);
    Index* col1 = (Index*)mxGetPr(plhs[1]);
    double* gpmm = mxGetPr(plhs[2]);

    debug_body("ugpmm: Length of columns = %d, nr of non-zero elements = %d", dim, nnz1 );
    debug_body("gpmm: Length of columns = %d, nr of non-zero elements = %d", dim, nnz2 );
    
    fem1d_ugpmm_CSR( p, N, P, M, pot, row, col, ugpmm);
    
    if ((row[dim-1]+1) != nnz1)
    {
        debug_body( "incorrect number of non-zero elements of ugpmm: %i != %i", row[dim-1], nnz1);
    
    }
    else
    {
        fem1d_gpmm_COO( dim, row, col, ugpmm, row1, col1, gpmm);
    }
    debug_exit("%s", "");
}
