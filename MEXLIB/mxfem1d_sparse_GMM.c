/*============================================================================+/
 | mxfem1d_sparse_GMM.c
 | Symmetric sparse matrix
 | CSC3 format with Lower triangluar part
/+============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mex.h"
#include "matrix.h"
#include "fem1d.h"


void mexFunction
(
          int      nlhs, 
          mxArray* plhs[],
          int      nrhs, 
    const mxArray* prhs[]
)
{
    if(nrhs!=3) 
    {
        mexErrMsgTxt("Two inputs required.");
    }
    if(nlhs!=1) 
    {
        mexErrMsgTxt("One output required.");
    } 

    matlib_index p;
    if(mxGetScalar(prhs[0])>1)
    {
        p = (matlib_index)floor(mxGetScalar(prhs[0]));
    }
    else
    {
        mexErrMsgTxt("Polynomial degree must be>1.");
    } 

    if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) 
    {
        mexErrMsgTxt("Input matrix must be type double.");
    }
    if (mxGetNumberOfDimensions(prhs[1]) != 2)
    {
        mexErrMsgTxt("Second argument must be a vector.");
    }
    if (mxGetNumberOfDimensions(prhs[2]) != 2)
    {
        mexErrMsgTxt("Third argument must be a matrix.");
    }
    matlib_index nr_combi = 3+(p-1)*(p+4)/2;
    matlib_xm Q = { .lenc = mxGetM(prhs[1]), 
                    .lenr = mxGetN(prhs[1]), 
                    .elem_p = mxGetPr(prhs[1]),
                    .op     = MATLIB_NO_TRANS,
                    .order  = MATLIB_COL_MAJOR };
    if((Q.lenc!=nr_combi) ||(Q.lenr<(p+1)))
    {
        mexErrMsgTxt("Size of quadrature matrix incorrect.");
    }
    matlib_xv phi_r = { .len = mxGetNumberOfElements(prhs[2]),
                        .elem_p = mxGetPr(prhs[2])};

    matlib_index P = Q.lenr - 1;
    if((phi_r.len<P) || (((phi_r.len-1)%P)!=0)) 
    {
        mexErrMsgTxt("Size of the potential vector is inconsistent.");
    } 
    matlib_index N = (phi_r.len-1)/P;
    matlib_index nnz = N*(Q.lenc-1)+1;

    matlib_index dim = N*p+1;
    matlib_index* rowIn = calloc( dim+1, sizeof(matlib_index));
    matlib_index* colIn = calloc(   nnz, sizeof(matlib_index));
    matlib_real* elem_p = calloc(   nnz, sizeof(matlib_real));

    matlib_xv q;
    matlib_create_xv( Q.lenc*N, &q, MATLIB_COL_VECT);

    matlib_index i, j;
 
    fem1d_XFLT( N, Q, phi_r, q);
    fem1d_XCSRGMM(p, N, q, rowIn, colIn, elem_p);

    if(mxIsComplex(prhs[2]))
    {
        matlib_xv phi_i = { .len = phi_r.len,
                            .elem_p = mxGetPi(prhs[2])};

        plhs[0] = mxCreateSparse( dim, dim, nnz, mxCOMPLEX); 

        mwSize* sir = mxGetIr(plhs[0]);
        mwSize* sjc = mxGetJc(plhs[0]);

        matlib_real* sr = mxGetPr(plhs[0]);
        matlib_real* si = mxGetPi(plhs[0]);

        for(i=0; i<dim; i++)
        {
            sjc[i] = rowIn[i];
            for(j=rowIn[i]; j<rowIn[i+1]; j++)
            {
                sir[j] = colIn[j];
                sr[j]  = elem_p[j];
            }
        }
        sjc[i] = nnz;

        fem1d_XFLT( N, Q, phi_i, q);
        fem1d_XCSRGMM(p, N, q, rowIn, colIn, elem_p);
        for(i=0; i<dim; i++)
        {
            for(j=rowIn[i]; j<rowIn[i+1]; j++)
            {
                si[j]  = elem_p[j];
            }
        }
    } 
    else
    {
        plhs[0] = mxCreateSparse( dim, dim, nnz, mxREAL); 
    
        mwSize* sir = mxGetIr(plhs[0]);
        mwSize* sjc = mxGetJc(plhs[0]);

        matlib_real* sr = mxGetPr(plhs[0]);

        for(i=0; i<dim; i++)
        {
            sjc[i] = rowIn[i];
            for(j=rowIn[i]; j<rowIn[i+1]; j++)
            {
                sir[j] = colIn[j];
                sr[j]  = elem_p[j];
            }
        }
        sjc[i] = nnz;
    }

    matlib_free(q.elem_p); 
    matlib_free(elem_p); 
    matlib_free(rowIn); 
    matlib_free(colIn); 
}
