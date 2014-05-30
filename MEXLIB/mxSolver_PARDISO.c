/*============================================================================+/
 | mxSolver_PARDISO.c
 | Symmetric sparse matrix
 | CSC3 format with Lower triangluar part
/+============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mex.h"
#include "matrix.h"
#include "matlib.h"


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
    matlib_zm_sparse M = { .lenc   = mxGetM(prhs[0]),
                           .lenr   = mxGetN(prhs[0]), 
                           .elem_p = mxGetPr(prhs[0]),
                           .op     = MATLIB_NO_TRANS,
                           .order  = MATLIB_COL_MAJOR };

    matlib_zv Pvb = { .len    = mxGetNumberOfElements(prhs[1]), 
                      .elem_p = mxGetPr()}
    pardiso_solver_t data = { .nsparse  = 1, 
                              .mnum     = 1, 
                              .sol_enum = PARDISO_LHS, 
                              .mtype    = PARDISO_COMPLEX_SYM,
                              .smat_p   = (void*)&M,
                              .rhs_p    = (void*)&Pvb,
                              .sol_p    = (void*)&V_vb};

    data.phase_enum = PARDISO_INIT;
    matlib_pardiso(&data);

    data.phase_enum = PARDISO_ANALYSIS_AND_FACTOR;
    matlib_pardiso(&data);

}
