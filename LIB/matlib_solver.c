/*============================================================================+/
 | Name: matlib_solver.c
/+============================================================================*/
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#define NDEBUG
#define MATLIB_NTRACE_DATA

#include "mkl.h"
#include "mkl_pardiso.h"
#include "matlib.h"

/*============================================================================*/

void matlib_pardiso(pardiso_solver_t* data)
/* 
 * Handles complex as well as real matrices.
 *
 * */ 
{
    debug_enter("%s", "");
    if (data->phase_enum == PARDISO_INIT)
    {
        matlib_index i;
        /* SETUP PARDISO CONTROL PARAMETERS */
        for (i = 0; i < PARDISO_NIPARAM; i++)
        {
            data->iparam[i] = 0;
            data->ptr[i]    = 0;
        }
        data->iparam[0]  = 1; /* Don't use default values */ 
        data->iparam[1]  = 2; /* Fill-in reducing odering for input matrix */ 
        data->iparam[3]  = 0; /* Preconditioning */ 
        data->iparam[4]  = 0; /*  */ 

        if(data->sol_enum==PARDISO_RHS)
        {
            data->iparam[5]  = 1; /* Write solution into b */ 
        }
        else if(data->sol_enum==PARDISO_LHS)
        {
            data->iparam[5]  = 0; /* Write solution into x */
        }
        else
        {
          term_exec( "Incorrect option for solution vector (sol_enum:%d)", 
                     data->sol_enum);
        }

        data->iparam[7]  = 2; /* Maximum number of iterative refinement steps, output 
                                 reported in iparam[6] */ 
        data->iparam[9]  = 13; /* Perturbing pivot elements */ 
        data->iparam[10] = 0;  /* Disable scaling  */ 
        data->iparam[12] = 0;  /*  */ 
        data->iparam[17] = -1; /* Disable reporting of nnz  */
        data->iparam[18] =  1; /*  */ 
        data->iparam[34] =  1; /* Zero-based indexing  */ 
        debug_body("%s", "Initialized PARDISO control parameters");
    }
    else if (data->phase_enum == PARDISO_ANALYSIS_AND_FACTOR)
    {
        debug_body("%s", "Start testing PARDISO");
        int nrhs  = 1; /* Number of right hand sides  */ 

        int msglvl = 0; /* Print statistical information in file */
        int error  = 0; /* Initialize error flag */
        debug_body("nr. sparse matrices: %d", data->nsparse);
        matlib_index maxfct = 1;
        matlib_index mnum = 1;

        if(data->nsparse>1)
        {
            if(data->mtype == PARDISO_COMPLEX_SYM)
            {
                debug_body("%s", "Factoring a complex symmetric matrix");
                matlib_znsparsem* smat_p = (matlib_znsparsem*) data->smat_p;
                PARDISO ( data->ptr, &maxfct, &mnum,
                          &data->mtype, &data->phase_enum,
                          &smat_p->lenc,
                          smat_p->elem_p[data->mnum-1],
                          smat_p->rowIn,
                          smat_p->colIn,
                          NULL, &nrhs,
                          data->iparam,
                          &msglvl, NULL, NULL,
                          &error);
            }
            else
            {
                debug_body("%s", "Factoring a real symmetric matrix");
                matlib_dnsparsem* smat_p = (matlib_dnsparsem*) data->smat_p;
                PARDISO ( data->ptr, &maxfct, &mnum,
                          &data->mtype, &data->phase_enum,
                          &smat_p->lenc,
                          smat_p->elem_p[data->mnum-1],
                          smat_p->rowIn,
                          smat_p->colIn,
                          NULL, &nrhs,
                          data->iparam,
                          &msglvl, NULL, NULL,
                          &error);
            }
        }
        else
        {
            if(data->mtype == PARDISO_COMPLEX_SYM)
            {
                debug_body("%s", "Factoring a complex symmetric matrix");
                matlib_zsparsem* smat_p = (matlib_zsparsem*) data->smat_p;
                PARDISO ( data->ptr, &maxfct, &mnum,
                          &data->mtype, &data->phase_enum,
                          &smat_p->lenc,
                          smat_p->elem_p,
                          smat_p->rowIn,
                          smat_p->colIn,
                          NULL, &nrhs,
                          data->iparam,
                          &msglvl, NULL, NULL,
                          &error);
            }
            else
            {
                debug_body("%s", "Factoring a real symmetric matrix");
                matlib_dsparsem* smat_p = (matlib_dsparsem*) data->smat_p;
                PARDISO ( data->ptr, &maxfct, &mnum,
                          &data->mtype, &data->phase_enum,
                          &smat_p->lenc,
                          smat_p->elem_p,
                          smat_p->rowIn,
                          smat_p->colIn,
                          NULL, &nrhs,
                          data->iparam,
                          &msglvl, NULL, NULL,
                          &error);
            }
        }
        if (error != 0)
        {
          term_exec("Analysis and factorization failed (error code: %d)", error);
        }
        debug_body("%s", "Analysis and factorization completed");
    }
    else if(data->phase_enum == PARDISO_SOLVE_AND_REFINE)
    {
        int nrhs  = 1; /* Number of right hand sides  */ 

        int msglvl = 0; /* Print statistical information in file */
        int error  = 0; /* Initialize error flag */
        matlib_index maxfct = 1;
        matlib_index mnum = 1;

        debug_body("nr. sparse matrices: %d", data->nsparse);
        if(data->nsparse>1)
        {
            if(data->mtype == PARDISO_COMPLEX_SYM)
            {

                matlib_znsparsem* smat_p = (matlib_znsparsem*) data->smat_p;
                matlib_zv* rhs_p  = (matlib_zv*) data->rhs_p;
                matlib_zv* sol_p  = (matlib_zv*) data->sol_p;

                BEGIN_DTRACE
                    debug_print("dimension of the sparse square matrix: %d", smat_p->lenc);
                    matlib_complex* ptr = smat_p->elem_p[data->mnum-1];
                    for(matlib_index i=0; i<smat_p->lenc; i++)
                    {
                        for(matlib_index j=smat_p->rowIn[i]; j<smat_p->rowIn[i+1]; j++)
                        {
                            debug_print( "M(%d,%d): % 0.16f %+0.16fi", 
                                         i, smat_p->colIn[j], ptr[j]);
                        }
                    }
                END_DTRACE

                PARDISO ( data->ptr, &maxfct, &(mnum),
                          &data->mtype, &data->phase_enum,
                          &smat_p->lenc,
                          smat_p->elem_p[data->mnum-1],
                          smat_p->rowIn,
                          smat_p->colIn,
                          NULL, &nrhs,
                          data->iparam,
                          &msglvl, 
                          rhs_p->elem_p, 
                          sol_p->elem_p,
                          &error);
            }
            else
            {
                matlib_dnsparsem* smat_p = (matlib_dnsparsem*) data->smat_p;
                matlib_dv* rhs_p  = (matlib_dv*) data->rhs_p;
                matlib_dv* sol_p  = (matlib_dv*) data->sol_p;

                PARDISO ( data->ptr, &maxfct, &mnum,
                          &data->mtype, &data->phase_enum,
                          &smat_p->lenc,
                          smat_p->elem_p[data->mnum-1],
                          smat_p->rowIn,
                          smat_p->colIn,
                          NULL, &nrhs,
                          data->iparam,
                          &msglvl, 
                          rhs_p->elem_p, 
                          sol_p->elem_p,
                          &error);
            }
        }
        else
        {
            if(data->mtype == PARDISO_COMPLEX_SYM)
            {
                matlib_zsparsem* smat_p = (matlib_zsparsem*) data->smat_p;
                matlib_zv* rhs_p  = (matlib_zv*) data->rhs_p;
                matlib_zv* sol_p  = (matlib_zv*) data->sol_p;

                debug_body("%s", "Solving a complex symmetric system.");
                PARDISO ( data->ptr, &maxfct, &mnum,
                          &data->mtype, &data->phase_enum,
                          &smat_p->lenc,
                          smat_p->elem_p,
                          smat_p->rowIn,
                          smat_p->colIn,
                          NULL, &nrhs,
                          data->iparam,
                          &msglvl, 
                          rhs_p->elem_p, 
                          sol_p->elem_p,
                          &error);
            }
            else
            {
                matlib_dsparsem* smat_p  = (matlib_dsparsem*) data->smat_p;
                matlib_dv* rhs_p  = (matlib_dv*) data->rhs_p;
                matlib_dv* sol_p  = (matlib_dv*) data->sol_p;

                PARDISO ( data->ptr, &maxfct, &mnum,
                          &data->mtype, &data->phase_enum,
                          &smat_p->lenc,
                          smat_p->elem_p,
                          smat_p->rowIn,
                          smat_p->colIn,
                          NULL, &nrhs,
                          data->iparam,
                          &msglvl, 
                          rhs_p->elem_p, 
                          sol_p->elem_p,
                          &error);
            }
        }

        if (error != 0)
        {
            term_exec("ERROR during solution and refinement: %d", error);
        }
    
        debug_body("%s", "Solution and refinement completed");
    }
    else if(data->phase_enum == PARDISO_FREE)
    {

        int msglvl = 0; /* Print statistical information in file */
        int error  = 0; /* Initialize error flag */
        matlib_index maxfct = 1;
        matlib_index mnum = 1;

        if(data->mtype == PARDISO_COMPLEX_SYM)
        {
            matlib_znsparsem* smat_p = (matlib_znsparsem*) data->smat_p;

            PARDISO ( data->ptr, &maxfct, &mnum,
                      &data->mtype, &data->phase_enum,
                      &smat_p->lenc,
                      NULL,
                      smat_p->rowIn,
                      smat_p->colIn,
                      NULL, NULL,
                      data->iparam,
                      &msglvl, 
                      NULL, 
                      NULL,
                      &error);
        }
        else
        {
            matlib_dnsparsem* smat_p = (matlib_dnsparsem*) data->smat_p;

            PARDISO ( data->ptr, &maxfct, &mnum,
                      &data->mtype, &data->phase_enum,
                      &smat_p->lenc,
                      NULL,
                      smat_p->rowIn,
                      smat_p->colIn,
                      NULL, NULL,
                      data->iparam,
                      &msglvl, 
                      NULL, NULL,
                      &error);
        }
    }
    
    debug_exit("%s", "");
}

