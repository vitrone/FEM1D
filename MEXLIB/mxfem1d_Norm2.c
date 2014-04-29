/*============================================================================+/
 | mxfem1d_Norm2.c
 | Handles complex as well as real inut data.
 | Complex mxArray is not compatible with C complex types. Hence, a complex 
 | number is treated as a 2-tuple.
 | 
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
        mexErrMsgIdAndTxt("FEM1D:Norm2:nrhs","Three inputs required.");
    }
    if(nlhs!=1) 
    {
        mexErrMsgIdAndTxt("FEM1D:Norm2:nlhs","One output required.");
    } 
    
    matlib_index p = (matlib_index)mxGetScalar(prhs[0]);
    matlib_index N = (matlib_index)mxGetScalar(prhs[1]);
    
    matlib_index dim = N*(p+1);
    matlib_dv ur = { .len = dim, 
                     .elem_p = mxGetPr(prhs[2]),
                     .type   = MATLIB_COL_VECT };

    double rnorm2, inorm2 = 0;
    rnorm2 = fem1d_DNorm2(p, N, ur);

    if(mxIsComplex(prhs[2]))
    {

        matlib_dv ui = { .len = dim, 
                         .elem_p = mxGetPi(prhs[2]),
                         .type   = MATLIB_COL_VECT };
        inorm2 = fem1d_DNorm2(p, N, ui);
    }
    plhs[0] = mxCreateDoubleScalar(sqrt(rnorm2*rnorm2 + inorm2*inorm2));
}
    
