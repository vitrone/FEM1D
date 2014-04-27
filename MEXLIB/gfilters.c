#include "mex.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
    double *mxGr, alpha;
    int i,n;
    mxArray *mxG[1];
    
    alpha = mxGetScalar(prhs[0]);
    n = ceil(mxGetScalar(prhs[1]));
    
    mxG[0] = mxCreateDoubleMatrix(1,n,mxREAL);
    mxGr = mxGetPr(mxG[0]);

    *(mxGr+0) = 1;
    *(mxGr+1) = -2*alpha;
    if(n>2){
        for(i=2;i<n;i++){
            *(mxGr+i)=(-2*alpha*(*(mxGr+i-1))+(i-2)*(*(mxGr+i-2)))/(i);
        }
    }
    *(plhs+0)    =   *(mxG+0);
    
    /* Clean up allocated memory. */

}


