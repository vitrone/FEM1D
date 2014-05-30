/*============================================================================+/
 | mxfem1d_quadM.c
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

    if(nrhs!=2) 
    {
        mexErrMsgTxt("Two inputs required.");
    }
    if(nlhs!=1) 
    {
        mexErrMsgTxt("One output required.");
    } 
    
    if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) 
    {
        mexErrMsgTxt("Input vector must be type double.");
    }

    matlib_xv quadW = { .len = mxGetNumberOfElements(prhs[0]), 
                        .elem_p = mxGetPr(prhs[0]), 
                        .type   = MATLIB_COL_VECT};

    if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) 
    {
        mexErrMsgTxt("Input matrix must be type double.");
    }
    matlib_xm IM = { .lenc = mxGetM(prhs[1]), 
                     .lenr = mxGetN(prhs[1]), 
                     .elem_p = mxGetPr(prhs[1]),
                     .op     = MATLIB_NO_TRANS,
                     .order  = MATLIB_COL_MAJOR };

    /* Minimum degree of polynomial allowed is p = 2
     * */ 
    if((IM.lenr<3) || (IM.lenc<3))
    {
        mexErrMsgTxt("Size of the inverse transform matrix is incorrect");
    } 
    if(quadW.len != IM.lenc) 
    {
        mexErrMsgTxt("Dimension mismatch.");
    } 
    matlib_xm Q;
    fem1d_quadM(quadW, IM, &Q);

    plhs[0] = mxCreateDoubleMatrix( Q.lenc, Q.lenr, mxREAL);
    double* ptr1 = mxGetPr(plhs[0]);
    double* ptr2 = Q.elem_p;

    for(matlib_index i=0; i<Q.lenc; i++)
    {
        for(matlib_index j=0; j<Q.lenr; j++)
        {
            ptr1[i+Q.lenc*j] = *ptr2;
            ptr2++;
        }
    }
    matlib_free(Q.elem_p);

}

                 
                 
                 
                 
