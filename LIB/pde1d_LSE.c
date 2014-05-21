/*============================================================================+/
 | File: pde1d_LSE.c 
 |
/+============================================================================*/
#include <pthread.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>

/* MKL */ 
#include "mkl.h"
#include "mkl_pardiso.h"

//#define NDEBUG
#define MATLIB_NTRACE_DATA

#include "pde1d_solver.h"
#include "assert.h"

/*============================================================================*/
/* Some default constants */
#define p_DEFAULT  4
#define N_DEFAULT  400
#define Nt_DEFAULT 200

#define nr_LGL_DEFAULT 9 /* (2*p_DEFAULT + 1) */

#define TOL_DEFAULT 1e-9
#define dt_DEFAULT  1.0e-3/2.0

#define nsparse_DEFAULT  20

#define x_l_DEFAULT -15.0
#define x_r_DEFAULT  15.0

/*============================================================================+/
 |SEQUENTIAL ROUTINES
/+============================================================================*/

void fem1d_linespace
(
    matlib_real  t0, 
    matlib_real  dt, 
    matlib_index Nt, 
    matlib_xv    *t
)
{
    matlib_create_xv( Nt+1, t, MATLIB_COL_VECT);

    t->elem_p[0] = t0;

    for(matlib_index i=1; i<Nt+1; i++)
    {
        t->elem_p[i] = t->elem_p[i-1] + dt;
    }
}

void pde1d_zm_sparse_GSM
/* Global Stiffness Matrix */ 
(
    matlib_index     N,
    matlib_complex   coeff,
    matlib_zm_sparse M
)
{
    debug_enter("%s", "");
    matlib_index i;

    M.elem_p[0] =  coeff/2.0 + M.elem_p[0];
    M.elem_p[1] = -coeff/2.0 + M.elem_p[1];
    for(i=1; i<N; i++)
    {
        M.elem_p[M.rowIn[i]]   =  coeff     + M.elem_p[M.rowIn[i]  ];
        M.elem_p[M.rowIn[i]+1] = -coeff/2.0 + M.elem_p[M.rowIn[i]+1];
    }
    M.elem_p[M.rowIn[i]] = coeff/2.0 + M.elem_p[M.rowIn[i]];
    
    for(i=N+1; i<M.lenc; i++)
    {
        M.elem_p[M.rowIn[i]] = coeff + M.elem_p[M.rowIn[i]];
    }
    debug_exit("%s", "");
}

void pde1d_zm_nsparse_GSM
/* Global Stiffness Matrix */ 
(
    matlib_index      N,
    matlib_complex    coeff,
    matlib_zm_nsparse M
)
{
    debug_enter("%s", "");
    matlib_complex* ptr;
    matlib_index i, k;
    for( k=0; k<M.nsparse; k++)
    {
        ptr = M.elem_p[k];
        ptr[0] =  coeff/2.0 + ptr[0];
        ptr[1] = -coeff/2.0 + ptr[1];
        for(i=1; i<N; i++)
        {
            ptr[M.rowIn[i]]   =  coeff     + ptr[M.rowIn[i]  ];
            ptr[M.rowIn[i]+1] = -coeff/2.0 + ptr[M.rowIn[i]+1];
        }
        ptr[M.rowIn[i]] = coeff/2.0 + ptr[M.rowIn[i]];
        
        for(i=N+1; i<M.lenc; i++)
        {
            ptr[M.rowIn[i]] = coeff + ptr[M.rowIn[i]];
        }
    }
    debug_exit("%s", "");
}

/*============================================================================*/

void pde1d_LSE_set_defaultsIVP(pde1d_LSE_data_t* input)
{
    debug_enter("%s", "");

    input->p      = p_DEFAULT;
    input->nr_LGL = nr_LGL_DEFAULT;
    input->N      = N_DEFAULT;

    input->domain[0] = x_l_DEFAULT;
    input->domain[1] = x_r_DEFAULT;
    
    input->dt = dt_DEFAULT;
    input->Nt = Nt_DEFAULT;
    input->nsparse = nsparse_DEFAULT;

    input->tol = TOL_DEFAULT;

    input->alpha = 1.0;

    input->sol_mode = PDE1D_LSE_EVOLVE_ONLY;

    input->u_analytic = NULL;
    input->params  = NULL;
    input->phix_p  = NULL;
    input->phixt_p = NULL;

    debug_exit("%s", "");
}
/*============================================================================*/

void pde1d_LSE_set_potential
(
    pde1d_LSE_data_t*   input,
    PDE1D_LSE_POTENTIAL phi_type,
    void* phi_p
)
{
    input->phi_type = phi_type;
    switch (input->phi_type)
    {
        case PDE1D_LSE_STATIC :
            if(phi_p != NULL)
            {
                input->phix_p  = phi_p;
            }
            else
            {
                term_exec("%s", "potential function is a NULL pointer ");
            }
            input->phixt_p = NULL;
            break;

        case PDE1D_LSE_DYNAMIC :
            if(phi_p != NULL)
            {
                input->phixt_p = phi_p;
            }
            else
            {
                term_exec("%s", "potential function is a NULL pointer ");
            }
            input->phix_p  = NULL;
            break;
    }
}
/*============================================================================*/

void pde1d_analytic_evol
(
    void** params,
    void (*u_analytic)(),
    matlib_xv x,
    matlib_xv t,
    matlib_zm u
)
{
    debug_enter("%s", "");
    matlib_zv u_tmp = {.len = u.lenc, .elem_p = u.elem_p};
    if(t.len != u.lenr)
    {
        term_exec("%s", "dimension mismatch!");
    
    }
    for(matlib_index i=0; i<t.len; i++)
    {
        (*u_analytic)(params, x, t.elem_p[i], u_tmp);
        (u_tmp.elem_p) += u.lenc;
    }

    debug_exit("%s", "");

}

/*============================================================================*/

void pde1d_LSE_init_solverIVP
(
    pde1d_LSE_data_t*   input, 
    pde1d_LSE_solver_t* data 
)
{
    debug_enter( "polynomial degree: %d, "
                 "nr. of LGL points: %d, "
                 "nr. of FEM-elements: %d, " 
                 "nr. of time-steps: %d", 
                 input->p, input->nr_LGL, input->N, input->Nt );
    

    debug_body( "domain: [%0.4f, %0.4f]", 
                input->domain[0], input->domain[1]);

    if(input->domain[0] >= input->domain[1])
    {
        term_exec("%s", "The computational domain is ill-defined.");
    }
    else
    {
        /* Jacobian of transformation 
         * */ 
        data->J = 0.5*(input->domain[1]-input->domain[0])/input->N;
    }

    if(input->dt > 0)
    {
        data->rho  = 2.0/input->dt;
        data->irho = input->dt/2.0;
    }
    else
    {
        term_exec("%s", "The time-step must be a positive non-zero quantity.");
    }

    /* Stiffness Matrix Coeff.
     * */
    data->s_coeff = I*(data->irho)*(input->alpha)/((data->J)*(data->J));     
    debug_body("s_coeff: %0.16f%+0.16fi", data->s_coeff);
    /* Mass Matrix Coeff.
     * */
    //data->m_coeff =  (matlib_complex*)calloc(2, sizeof(matlib_complex));
    data->m_coeff[0] =  1.0;
    data->m_coeff[1] = -I*data->irho;

    /* LGL-points on the refrence domain: [-1, 1]
     * vector: xi
     * Quadrature weights on LGL-points
     * vector: quadW 
     * */ 
    legendre_LGLdataLT1( (input->nr_LGL)-1, 
                          input->tol, 
                          &(data->xi), 
                          &(data->quadW));
    
    /* Forward Transform matrix : FM
     * Backward Transform matrix: IM
     * */ 
    matlib_create_xm( (input->p)+1, 
                      (data->xi).len, 
                      &(data->FM), 
                      MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);

    matlib_create_xm( (data->xi).len, 
                      (input->p)+1, 
                      &(data->IM), 
                      MATLIB_COL_MAJOR, MATLIB_NO_TRANS);    

    legendre_LGLdataFM( data->xi, data->FM);
    legendre_LGLdataIM( data->xi, data->IM);

    /* Quadrature matrix needed for assembling Global Mass Matrix 
     * */ 
    fem1d_quadM( data->quadW, data->IM, &(data->Q));

    /* generate the grid: x */ 
    fem1d_ref2mesh( data->xi, 
                    input->N, 
                    (input->domain)[0], 
                    (input->domain)[1], 
                    &input->x);

    /* generate the grids: x and t
     * */ 
    matlib_real t0 = 0;
    fem1d_linespace( t0, input->dt, input->Nt, &(input->t));
    debug_body("length of x: %d", input->x.len);
    debug_body("length of t: %d", input->t.len);
    
    /* Initial data vector
     * */ 
    matlib_create_zv( input->x.len, &(input->u_init), MATLIB_COL_VECT);

    /* Initialize all temporary variables using this array of pointers 
     * */ 
    matlib_index i;
    data->nr_vars = 5;
    matlib_zv** var_p_ = (matlib_zv**)calloc(data->nr_vars, sizeof(matlib_zv*));
    for(i=0; i<data->nr_vars; i++)
    {
        var_p_[i] = (matlib_zv*)calloc(1, sizeof(matlib_zv));
    }
    BEGIN_DEBUG
        for(i=0; i<data->nr_vars; i++)
        {
            debug_print( "initialize temporary var pointers[%i]: %p", 
                         i, var_p_[i]);
        }
    END_DEBUG


    data->var_p = (void**)var_p_;
    /* Vectors in FEM-basis
     * var_p[1] : V_vb 
     * Projection onto FEM-basis
     * var_p[0] : Pvb  
     * */ 
    matlib_index dim = (input->N)*(input->p)+1;
    debug_body("FEM-basis dimension: %d", dim);
    matlib_create_zv( dim, var_p_[0], MATLIB_COL_VECT);
    matlib_create_zv( dim, var_p_[1], MATLIB_COL_VECT);
    
    dim = dim + (input->N)-1;
    debug_body("Legendre-basis dimension: %d", dim);
    /* vector in Legendre basis
     * var_p[2] : V_tmp 
     * var_p[3] : U_tmp 
     * */ 
    matlib_create_zv( dim, var_p_[2], MATLIB_COL_VECT);
    matlib_create_zv( dim, var_p_[3], MATLIB_COL_VECT);
    
    /* Vector containing values of the time-independent potential
     * var_p[3] : phi
     * */ 
    matlib_create_zv( (input->x).len, 
                      var_p_[4], 
                      MATLIB_COL_VECT);

    if(input->sol_mode==PDE1D_LSE_ERROR_ONLY)
    {
        matlib_create_xv( (input->Nt)+1,
                          &(input->e_rel), MATLIB_COL_VECT);
        matlib_create_xv( (input->Nt)+1,
                          &(input->e_abs), MATLIB_COL_VECT);
    }
    else
    {
        /* Initialize the evolution matrix in column major format 
         * Stores the field in Legendre basis
         * matrix: U_evol
         * */ 
        matlib_create_zm( dim, 
                          (input->t).len, 
                          &(input->U_evol), 
                          MATLIB_COL_MAJOR, MATLIB_NO_TRANS);
    

    }
    debug_exit("%s", "");
}

/*============================================================================*/

void pde1d_LSE_solve_IVP
(
    pde1d_LSE_data_t*  input,
    pde1d_LSE_solver_t* data
)
{
    debug_enter("%s", "");
    
    void (*LSE_solver_IVP_p[2])();
 
    switch(input->sol_mode)
    {
        case PDE1D_LSE_EVOLVE_ONLY:
            LSE_solver_IVP_p[0] = pde1d_LSE_solve_IVP_evol;
            LSE_solver_IVP_p[1] = pde1d_LSE_solve_IVP2_evol;
            break;
        case PDE1D_LSE_ERROR_ONLY:
            LSE_solver_IVP_p[0] = pde1d_LSE_solve_IVP_error;
            LSE_solver_IVP_p[1] = pde1d_LSE_solve_IVP2_error;
            break;
    }


    switch(input->phi_type)
    {
        case PDE1D_LSE_STATIC:
            (*LSE_solver_IVP_p[0])(input, data);
            break;

        case PDE1D_LSE_DYNAMIC:
            (*LSE_solver_IVP_p[1])(input, data);
            break;
    }

}

void pde1d_LSE_destroy_solverIVP
(
    pde1d_LSE_data_t*   input,
    pde1d_LSE_solver_t* data
)
{
    debug_enter("%s", "Freeing all initialized memory.");

    matlib_free(input->x.elem_p);
    debug_body("Freed: %s", "x");

    matlib_free(input->t.elem_p);
    debug_body("Freed: %s", "t");
    
    matlib_free(input->u_init.elem_p);
    debug_body("Freed: %s", "u_init");
    
    matlib_free(data->xi.elem_p);
    debug_body("Freed: %s", "xi");

    matlib_free(data->quadW.elem_p);
    debug_body("Freed: %s", "quadW");

    matlib_free(data->FM.elem_p);
    debug_body("Freed: %s", "FM");

    matlib_free(data->IM.elem_p);
    debug_body("Freed: %s", "IM");

    matlib_free(data->Q.elem_p);
    debug_body("Freed: %s", "Q");

    switch(input->phi_type)
    {
        case PDE1D_LSE_STATIC:
            matlib_free(data->M.elem_p);
            matlib_free(data->M.colIn);
            matlib_free(data->M.rowIn);
            debug_body("Freed: %s", "M");
            break;
        case PDE1D_LSE_DYNAMIC:
            break;
    }
    if(input->sol_mode==PDE1D_LSE_ERROR_ONLY)
    {
        matlib_free(input->e_rel.elem_p);
        debug_body("Freed: %s", "e_rel");

        matlib_free(input->e_abs.elem_p);
        debug_body("Freed: %s", "e_abs");
    }
    
    matlib_index i;
    void* ptr = NULL;
    for(i=0; i<data->nr_vars; i++)
    {
        ptr = ((matlib_zv*)data->var_p[i])->elem_p;
        matlib_free(ptr);
    }
    matlib_free(data->var_p);
    debug_body("Freed: %s", "var_p");

    debug_exit("%s", "");
}

/*============================================================================*/
/* Static potential case
 * */ 
void pde1d_LSE_solve_IVP_evol
(
    pde1d_LSE_data_t*  input,
    pde1d_LSE_solver_t* data
)
{
    debug_enter("%s", "");
    debug_enter( "polynomial degree: %d, "
                 "nr. of LGL points: %d", 
                 input->p, input->nr_LGL );

    matlib_index i, j;

    /* Other temporary variables 
     * */ 
    matlib_zv Pvb   = *(matlib_zv*)(data->var_p[0]);
    matlib_zv V_vb  = *(matlib_zv*)(data->var_p[1]);
    matlib_zv V_tmp = *(matlib_zv*)(data->var_p[2]);
    matlib_zv U_tmp = *(matlib_zv*)(data->var_p[3]);
    matlib_zv phi   = *(matlib_zv*)(data->var_p[4]);

    debug_body("%s", "declared temporary vars");
    /* Filling the potential vector
     * */ 
    void (*phi_p)() = input->phix_p;
    (*phi_p)(input->params, data->m_coeff, input->x, phi);
    debug_body("%s", "potential computed");

    /* Initialize the sparse matrix in order to store the 
     * Global Stiffness-Mass Matrix: M
     * */ 
    //matlib_zm_sparse M;
    fem1d_zm_sparse_GMM(input->p, data->Q, phi, &(data->M));
    pde1d_zm_sparse_GSM(input->N, data->s_coeff, data->M);

    BEGIN_DTRACE
        debug_print( "dimension of the sparse square matrix: %d",
                     data->M.lenc);
        for(i=0; i<data->M.lenc; i++)
        {
            for(j=data->M.rowIn[i]; j<data->M.rowIn[i+1]; j++)
            {
                debug_print( "M(%d,%d): % 0.16f %+0.16fi", 
                             i, data->M.colIn[j], data->M.elem_p[j]);
            }
        }
    END_DTRACE

    /* Solve the equation while marching in time 
     * Linear system M * V_vb = Pvb 
     * */

    matlib_zv U_tmp1 = { .len    = input->U_evol.lenc, 
                         .elem_p = input->U_evol.elem_p, 
                         .type   = MATLIB_COL_VECT};

    fem1d_ZFLT(input->N, data->FM, input->u_init, U_tmp);
    matlib_zcopy(U_tmp, U_tmp1 );
    (U_tmp1.elem_p) += (U_tmp1.len); 

    pardiso_solver_t eq_data = { .nsparse  = 1, 
                                 .mnum     = 1, 
                                 .sol_enum = PARDISO_LHS, 
                                 .mtype    = PARDISO_COMPLEX_SYM,
                                 .smat_p   = (void*) &(data->M),
                                 .rhs_p    = (void*) &Pvb,
                                 .sol_p    = (void*) &V_vb};

    eq_data.phase_enum = PARDISO_INIT;
    matlib_pardiso(&eq_data);

    eq_data.phase_enum = PARDISO_ANALYSIS_AND_FACTOR;
    matlib_pardiso(&eq_data);

    for (i=0; i<input->Nt; i++)
    {
        debug_body("begin iteration: %d", i);
        fem1d_ZPrjL2F(input->p, U_tmp, Pvb);

        eq_data.phase_enum = PARDISO_SOLVE_AND_REFINE;
        matlib_pardiso(&eq_data);
        
        fem1d_ZF2L(input->p, V_vb, V_tmp);

        /* 2.0 * V_tmp -U_tmp --> U_tmp
         * */ 
        matlib_zaxpby(2.0, V_tmp, -1.0, U_tmp );
        matlib_zcopy(U_tmp, U_tmp1 );
        (U_tmp1.elem_p) += (U_tmp1.len); 
    }

    eq_data.phase_enum = PARDISO_FREE;
    matlib_pardiso(&eq_data);

    debug_exit("%s", "");
}
/*============================================================================*/

void pde1d_LSE_solve_IVP_error
(
    pde1d_LSE_data_t*  input,
    pde1d_LSE_solver_t* data
)
{
    debug_enter( "polynomial degree: %d, "
                 "nr. of LGL points: %d", 
                 input->p, input->nr_LGL );

    matlib_index i, j;

    /* Other temporary variables 
     * */ 
    matlib_zv Pvb   = *(matlib_zv*)(data->var_p[0]);
    matlib_zv V_vb  = *(matlib_zv*)(data->var_p[1]);
    matlib_zv V_tmp = *(matlib_zv*)(data->var_p[2]);
    matlib_zv U_tmp = *(matlib_zv*)(data->var_p[3]);
    matlib_zv phi   = *(matlib_zv*)(data->var_p[4]);

    debug_body("%s", "declared temporary vars");
    /* Filling the potential vector
     * */ 
    void (*phi_p)() = input->phix_p;
    (*phi_p)(input->params, data->m_coeff, input->x, phi);
    debug_body("%s", "potential computed");

    /* Analytic solution
     * */ 
    void (*u_analytic)() = input->u_analytic;
    /* Provide the initial condition 
     * */ 
    (*u_analytic)(input->params, input->x, (input->t).elem_p[0], input->u_init);
    debug_body("%s", "provided initial data");
    

    /* Initialize the sparse matrix in order to store the 
     * Global Stiffness-Mass Matrix: M
     * */ 
    fem1d_zm_sparse_GMM(input->p, data->Q, phi, &data->M);
    pde1d_zm_sparse_GSM(input->N, data->s_coeff, data->M);

    BEGIN_DTRACE
        debug_print("dimension of the sparse square matrix: %d", data->M.lenc);
        for(i=0; i<data->M.lenc; i++)
        {
            for(j=data->M.rowIn[i]; j<data->M.rowIn[i+1]; j++)
            {
                debug_print( "M(%d,%d): % 0.16f %+0.16fi", 
                             i, data->M.colIn[j], data->M.elem_p[j]);
            }
        }
    END_DTRACE

    /* Solve the equation while marching in time 
     * Linear system M * V_vb = Pvb 
     * */

    fem1d_ZFLT(input->N, data->FM, input->u_init, U_tmp);

    pardiso_solver_t eq_data = { .nsparse  = 1, 
                                 .mnum     = 1, 
                                 .sol_enum = PARDISO_LHS, 
                                 .mtype    = PARDISO_COMPLEX_SYM,
                                 .smat_p   = (void*) &(data->M),
                                 .rhs_p    = (void*) &Pvb,
                                 .sol_p    = (void*) &V_vb};

    eq_data.phase_enum = PARDISO_INIT;
    matlib_pardiso(&eq_data);

    eq_data.phase_enum = PARDISO_ANALYSIS_AND_FACTOR;
    matlib_pardiso(&eq_data);

    matlib_real norm_actual;

    (input->e_abs).elem_p[0] = 0;
    (input->e_rel).elem_p[0] = 0;

    for (i=0; i<input->Nt; i++)
    {
        debug_body("begin iteration: %d", i);
        fem1d_ZPrjL2F(input->p, U_tmp, Pvb);

        eq_data.phase_enum = PARDISO_SOLVE_AND_REFINE;
        matlib_pardiso(&eq_data);
        
        fem1d_ZF2L(input->p, V_vb, V_tmp);

        /* 2.0 * V_tmp -U_tmp --> U_tmp
         * */ 
        matlib_zaxpby(2.0, V_tmp, -1.0, U_tmp );

        /* Using phi to store the analytic solution
         * */ 
        (*u_analytic)(input->params, input->x, input->t.elem_p[i+1], phi);
        fem1d_ZFLT(input->N, data->FM, phi, V_tmp);

        /* Error analysis: 
         * */ 
        norm_actual = fem1d_ZNorm2(input->p, input->N, V_tmp);
        matlib_zaxpy(-1.0, U_tmp, V_tmp );

        (input->e_abs).elem_p[i+1] = fem1d_ZNorm2(input->p, input->N, V_tmp);
        (input->e_rel).elem_p[i+1] = (input->e_abs).elem_p[i+1]/fmax(norm_actual, input->tol);
        
        debug_body("Absolute Error: %0.16f", (input->e_abs).elem_p[i+1]);
        debug_body("Relative Error: %0.16f", (input->e_rel).elem_p[i+1]);
    }

    eq_data.phase_enum = PARDISO_FREE;
    matlib_pardiso(&eq_data);

    /* Free the sparse matrix 
     * */ 
    debug_exit("%s", "");
}

/*============================================================================*/

void pde1d_LSE_solve_IVP2_evol
(
    pde1d_LSE_data_t*   input,
    pde1d_LSE_solver_t* data
)
{

    debug_enter( "polynomial degree: %d, "
                 "nr. of LGL points: %d",
                 input->p, input->nr_LGL );

    matlib_index i, j;
    
    /* Other temporary variables 
     * */ 
    matlib_zv Pvb   = *(matlib_zv*)(data->var_p[0]);
    matlib_zv V_vb  = *(matlib_zv*)(data->var_p[1]);
    matlib_zv V_tmp = *(matlib_zv*)(data->var_p[2]);
    matlib_zv U_tmp = *(matlib_zv*)(data->var_p[3]);
    matlib_zv u_exact = *(matlib_zv*)(data->var_p[4]);
    debug_body("%s", "declared temporary vars");

    /* declare potential function
     * */ 
    void (*phi_p)() = input->phixt_p;
    debug_body("%s", "potential function declared");

    /* Assemble the global mass matrix 
     * */ 
    matlib_zm phi, q;
    matlib_index nsparse = input->nsparse;
    fem1d_zm_nsparse_GMM( input->p, input->N, nsparse, 
                          data->Q, &phi, &q, &data->nM, FEM1D_GMM_INIT);

    fem1d_zm_nsparse_GMM( input->p, input->N, nsparse, data->Q, 
                          NULL, NULL, &data->nM, FEM1D_GET_SPARSITY_ONLY);
    
    debug_body("phi: %d-by-%d", phi.lenc, phi.lenr);
    /* Setup the sparse linear system */
    /* Solve the equation while marching in time 
     *
     * Linear system M * V_vb = Pvb 
     * */
    matlib_zv U_tmp1 = { .len    = input->U_evol.lenc, 
                         .elem_p = input->U_evol.elem_p, 
                         .type   = MATLIB_COL_VECT};

    fem1d_ZFLT(input->N, data->FM, input->u_init, U_tmp);
    matlib_zcopy(U_tmp, U_tmp1 );
    (U_tmp1.elem_p) += (U_tmp1.len); 

    pardiso_solver_t eq_data = { .nsparse  = nsparse, 
                                 .mnum     = 1, 
                                 .sol_enum = PARDISO_LHS, 
                                 .mtype    = PARDISO_COMPLEX_SYM,
                                 .smat_p   = (void*)&(data->nM),
                                 .rhs_p    = (void*)&Pvb,
                                 .sol_p    = (void*)&V_vb};

    debug_body("%s", "Solver data initialized");

    BEGIN_DTRACE
        matlib_complex *ptr;
        for(matlib_index k = 0; k< nsparse; k++)
        {
            debug_print( "dimension of the sparse square matrix: %d", 
                         data->nM.lenc);
            ptr = data->nM.elem_p[k];
            for(i=0; i<data->nM.lenc; i++)
            {
                for(j=data->nM.rowIn[i]; j<data->nM.rowIn[i+1]; j++)
                {
                    debug_print( "M(%d,%d): % 0.16f %+0.16fi",
                                 i, data->nM.colIn[j], ptr[j]);
                }
            }
        }
    END_DTRACE

    eq_data.phase_enum = PARDISO_INIT;
    matlib_pardiso(&eq_data);
    debug_body("%s", "PARDISO initialized");

    matlib_index Nt_ = input->Nt/nsparse;
    matlib_xv t_tmp  = {.len = (nsparse + 1), .elem_p = input->t.elem_p}; 
    
    for (i=0; i<Nt_; i++)
    {
        (*phi_p)(input->params, data->m_coeff, input->x, t_tmp, phi);
        fem1d_zm_nsparse_GMM( input->p, input->N, nsparse, data->Q, 
                              &phi, &q, &data->nM, FEM1D_GET_NZE_ONLY);
        pde1d_zm_nsparse_GSM(input->N, data->s_coeff, data->nM);

        for(j=0; j<nsparse; j++)
        {
            debug_body("begin iteration: %d", i*nsparse+j);
            fem1d_ZPrjL2F(input->p, U_tmp, Pvb);

            eq_data.mnum = j+1;
            eq_data.phase_enum = PARDISO_ANALYSIS_AND_FACTOR;
            matlib_pardiso(&eq_data);

            eq_data.phase_enum = PARDISO_SOLVE_AND_REFINE;
            matlib_pardiso(&eq_data);
            
            fem1d_ZF2L(input->p, V_vb, V_tmp);

            /* 2.0 * V_tmp -U_tmp --> U_tmp*/ 
            matlib_zaxpby(2.0, V_tmp, -1.0, U_tmp );
            matlib_zcopy(U_tmp, U_tmp1 );
            (U_tmp1.elem_p) += (U_tmp1.len);
        }
        t_tmp.elem_p += nsparse;
    }

    eq_data.phase_enum = PARDISO_FREE;
    matlib_pardiso(&eq_data);
    fem1d_zm_nsparse_GMM( input->p, input->N, nsparse, 
                          data->Q, &phi, &q, &data->nM, FEM1D_GMM_FREE);

    debug_exit("%s", "");
}


void pde1d_LSE_solve_IVP2_error
(
    pde1d_LSE_data_t*   input,
    pde1d_LSE_solver_t* data
)
{

    debug_enter( "polynomial degree: %d, "
                 "nr. of LGL points: %d",
                 input->p, input->nr_LGL );

    matlib_index i, j;
    
    /* Other temporary variables 
     * */ 
    matlib_zv Pvb   = *(matlib_zv*)(data->var_p[0]);
    matlib_zv V_vb  = *(matlib_zv*)(data->var_p[1]);
    matlib_zv V_tmp = *(matlib_zv*)(data->var_p[2]);
    matlib_zv U_tmp = *(matlib_zv*)(data->var_p[3]);
    matlib_zv u_exact = *(matlib_zv*)(data->var_p[4]);
    debug_body("%s", "declared temporary vars");

    /* declare potential function
     * */ 
    void (*phi_p)() = input->phixt_p;
    debug_body("%s", "potential function declared");

    /* Analytic solution
     * */ 
    void (*u_analytic)() = input->u_analytic;
    /* Provide the initial condition 
     * */ 
    (*u_analytic)(input->params, input->x, (input->t).elem_p[0], input->u_init);
    debug_body("%s", "provided initial data");

    /* Assemble the global mass matrix 
     * */ 
    matlib_zm phi, q;
    matlib_index nsparse = input->nsparse;
    fem1d_zm_nsparse_GMM( input->p, input->N, nsparse, 
                          data->Q, &phi, &q, &data->nM, FEM1D_GMM_INIT);

    fem1d_zm_nsparse_GMM( input->p, input->N, nsparse, data->Q, 
                          NULL, NULL, &data->nM, FEM1D_GET_SPARSITY_ONLY);
    
    debug_body("phi: %d-by-%d", phi.lenc, phi.lenr);
    /* Setup the sparse linear system */
    /* Solve the equation while marching in time 
     *
     * Linear system M * V_vb = Pvb 
     * */

    fem1d_ZFLT(input->N, data->FM, input->u_init, U_tmp);


    pardiso_solver_t eq_data = { .nsparse  = nsparse, 
                                 .mnum     = 1, 
                                 .sol_enum = PARDISO_LHS, 
                                 .mtype    = PARDISO_COMPLEX_SYM,
                                 .smat_p   = (void*)&(data->nM),
                                 .rhs_p    = (void*)&Pvb,
                                 .sol_p    = (void*)&V_vb};

    debug_body("%s", "Solver data initialized");

    BEGIN_DTRACE
        matlib_complex *ptr;
        for(matlib_index k = 0; k< nsparse; k++)
        {
            debug_print( "dimension of the sparse square matrix: %d", 
                         data->nM.lenc);
            ptr = data->nM.elem_p[k];
            for(i=0; i<data->nM.lenc; i++)
            {
                for(j=data->nM.rowIn[i]; j<data->nM.rowIn[i+1]; j++)
                {
                    debug_print( "M(%d,%d): % 0.16f %+0.16fi",
                                 i, data->nM.colIn[j], ptr[j]);
                }
            }
        }
    END_DTRACE

    eq_data.phase_enum = PARDISO_INIT;
    matlib_pardiso(&eq_data);
    debug_body("%s", "PARDISO initialized");

    matlib_index Nt_ = input->Nt/nsparse;
    matlib_xv t_tmp  = {.len = (nsparse + 1), .elem_p = input->t.elem_p}; 

    matlib_real norm_actual;
    (input->e_abs).elem_p[0] = 0;
    (input->e_rel).elem_p[0] = 0;
    
    for (i=0; i<Nt_; i++)
    {
        (*phi_p)(input->params, data->m_coeff, input->x, t_tmp, phi);
        fem1d_zm_nsparse_GMM( input->p, input->N, nsparse, data->Q, 
                              &phi, &q, &data->nM, FEM1D_GET_NZE_ONLY);
        pde1d_zm_nsparse_GSM(input->N, data->s_coeff, data->nM);

        for(j=0; j<nsparse; j++)
        {
            debug_body("begin iteration: %d", i*nsparse+j);
            fem1d_ZPrjL2F(input->p, U_tmp, Pvb);

            eq_data.mnum = j+1;
            eq_data.phase_enum = PARDISO_ANALYSIS_AND_FACTOR;
            matlib_pardiso(&eq_data);

            eq_data.phase_enum = PARDISO_SOLVE_AND_REFINE;
            matlib_pardiso(&eq_data);
            
            fem1d_ZF2L(input->p, V_vb, V_tmp);

            /* 2.0 * V_tmp -U_tmp --> U_tmp*/ 
            matlib_zaxpby(2.0, V_tmp, -1.0, U_tmp );

            /* Using u_exact to store the analytic solution
             * */ 
            (*u_analytic)(input->params, input->x, t_tmp.elem_p[j+1], u_exact);
            fem1d_ZFLT(input->N, data->FM, u_exact, V_tmp);

            /* Error analysis: 
             * */ 
            norm_actual = fem1d_ZNorm2(input->p, input->N, V_tmp);
            matlib_zaxpy(-1.0, U_tmp, V_tmp );

            (input->e_abs).elem_p[j+nsparse*i+1] = fem1d_ZNorm2(input->p, input->N, V_tmp);
            (input->e_rel).elem_p[j+nsparse*i+1] = 
                   (input->e_abs).elem_p[j+nsparse*i+1]/fmax(norm_actual, input->tol);
            
            debug_body("Absolute Error: %0.16f", (input->e_abs).elem_p[j+nsparse*i+1]);
            debug_body("Relative Error: %0.16f", (input->e_rel).elem_p[j+nsparse*i+1]);
        }
        t_tmp.elem_p += nsparse;
    }

    eq_data.phase_enum = PARDISO_FREE;
    matlib_pardiso(&eq_data);
    fem1d_zm_nsparse_GMM( input->p, input->N, nsparse, 
                          data->Q, &phi, &q, &data->nM, FEM1D_GMM_FREE);

    debug_exit("%s", "");
}

/*============================================================================*/
void pde1d_LSE_Gaussian_WP_constant_potential
(
    void** params,
    matlib_xv   x,
    matlib_real t,
    matlib_zv   u
)
{
    debug_enter("%s", "");
    
    matlib_complex A_0 = *(matlib_complex*) params[0];
    matlib_complex a   = *(matlib_complex*) params[1];

    if(creal(a)<0)
    {
        term_exec( "%s", "real part of params[1] < 0");
    }

    matlib_real c = *(matlib_real*) params[2];

    matlib_complex phi_0 = *(matlib_complex*) params[3];

    debug_body( "Gaussian Wave Packet: "
                "apmlitude = %0.4f%+0.4fi, "
                "speed = %0.4f",
                creal(A_0), cimag(A_0), c);

    debug_body( "Gaussian Wave Packet: "
                "constant potential = %0.4f%+0.4fi, "
                "time = %0.4f",
                creal(phi_0), cimag(phi_0), t);

    matlib_complex C1 = 1.0/(4.0*a*t*I+1.0);
    matlib_complex C2 = A_0*csqrt(C1)*cexp(-0.25*I*c*c*t + I*phi_0*t);
    matlib_real C3 = c*t;
    matlib_real C4 = 0.5*c;


    matlib_real x_ = 0;
    matlib_real *xptr;
    for (xptr = x.elem_p; xptr<(x.elem_p+x.len); xptr++)
    {
        x_ = *xptr-C3;
        *(u.elem_p) = C2*cexp(-a*C1*x_*x_ + I*C4**xptr);
        u.elem_p++;
    }
    debug_exit("%s", "");
}

void pde1d_LSE_constant_potential
( 
    void** params,
    matlib_complex m_coeff[2], 
    matlib_xv      x, 
    matlib_zv      y
)
{
    matlib_index i;
    matlib_complex a = 0, b = 0;

    matlib_complex phi_0 = *(matlib_complex*) params[3];
    
    if(m_coeff != NULL)
    {
        a = m_coeff[0];
        b = m_coeff[1];
    }
    if(y.len >= x.len)
    {
        for (i=0; i<x.len; i++)
        {
            *(y.elem_p) = a + b*phi_0;
            y.elem_p++;
        }
    }
    else
    {
        term_exec("%s", "size mismatch for vectors");
    }
}

void pde1d_LSE_Gaussian_CW_constant_potential
(
    void** params,
    matlib_xv   x,
    matlib_real t,
    matlib_zv   u
)
{
    debug_enter("%s", "");
    
    matlib_complex A_0 = *(matlib_complex*) params[0];
    matlib_complex a   = *(matlib_complex*) params[1];

    if(creal(a)<0)
    {
        term_exec( "%s", "real part of params[1] < 0");
    }

    matlib_real c = *(matlib_real*) params[2];

    matlib_complex phi_0 = *(matlib_complex*) params[3];
    matlib_complex A_CW  = *(matlib_complex*) params[4];

    debug_body( "Gaussian Wave Packet: "
                "apmlitude = %0.4f%+0.4fi, "
                "speed = %0.4f",
                creal(A_0), cimag(A_0), c);

    debug_body( "Gaussian Wave Packet: "
                "constant potential = %0.4f%+0.4fi, "
                "time = %0.4f",
                creal(phi_0), cimag(phi_0), t);

    matlib_complex C1 = 1.0/(4.0*a*t*I+1.0);
    matlib_complex E  = cexp(-0.25*I*c*c*t + I*phi_0*t);
    matlib_complex C2 = A_0*csqrt(C1)*E;
    matlib_real C3 = c*t;
    matlib_real C4 = 0.5*c;
    matlib_complex C5 = A_CW*E;


    matlib_real x_ = 0;
    matlib_real *xptr;
    for (xptr = x.elem_p; xptr<(x.elem_p+x.len); xptr++)
    {
        x_ = *xptr-C3;
        *(u.elem_p) = (C5 + C2*cexp(-a*C1*x_*x_))*cexp(I*C4**xptr);
        u.elem_p++;
    }
    debug_exit("%s", "");
}
/*============================================================================*/

void pde1d_LSE_HermiteGaussian_WP_harmonic_potential
(
    void** params,
    matlib_xv   x,
    matlib_real t,
    matlib_zv   u
)
/* 
 * phi = - x^2
 *
 * u = psi_0(x)exp(-it)+psi_1(x)exp(-i3*t)
 * psi_0 = (1/pi)^1/4 * exp(-x^2/2)
 * psi_1 = (1/pi)^1/4 * sqrt(2) x exp(-x^2/2)
 *
 *
 * */
{
    debug_enter("%s", "");
    
    matlib_real coeff[2] = {1.0/pow(M_PI, 0.25), sqrt(2.0)/pow(M_PI, 0.25)};
    matlib_complex e[2];
    e[0] = cexp(-I*t); 
    e[1] = cexp(-I*3*t);

    matlib_real *xptr;

    for (xptr = x.elem_p; xptr<(x.elem_p+x.len); xptr++)
    {
        *(u.elem_p) = (e[0]*coeff[0]+e[1]*coeff[1]**xptr)*exp(-0.5**xptr**xptr);
        u.elem_p++;
    }
    debug_exit("%s", "");
}

void pde1d_LSE_harmonic_potential
( 
    void** params,
    matlib_complex m_coeff[2], 
    matlib_xv      x, 
    matlib_zv      y
)
{
    matlib_index i;
    if(y.len >= x.len)
    {
        for (i=0; i<x.len; i++)
        {
            *(y.elem_p) = m_coeff[0]
                          + m_coeff[1]*(-*(x.elem_p+i)**(x.elem_p+i));
            y.elem_p++;
        }
    }
}


/*============================================================================*/
void pde1d_LSE_Gaussian_WP_timedependent_linear_potential
(
    void** params,
    matlib_xv   x,
    matlib_real t,
    matlib_zv   u
)
/* 
 * 
 * phi    = g_0 cos(mu*t)
 * beta   = g_0 *sin(mu*t)/mu 
 * nu     = (2*g_0*cos(mu*t))/mu^2 - (2*g_0)/mu^2;
 * ibeta2 = (1/mu^3)(mu*(g_0^2*t)/2 - (g_0^2*sin(2*mu*t))/4;
 * Xi     = beta*x-ibeta2;
 *
 * */ 
{
    debug_enter("%s", "");
    
    matlib_complex A_0 = *(matlib_complex*) params[0];
    matlib_complex a   = *(matlib_complex*) params[1];

    if(creal(a)<0)
    {
        term_exec( "%s", "real part of params[1] < 0");
    }

    matlib_real c = *(matlib_real*) params[2];

    matlib_complex phi_0 = *(matlib_complex*) params[3];

    matlib_real g_0 = *(matlib_real*) params[4];
    matlib_real mu  = *(matlib_real*) params[5];

    debug_body( "Gaussian Wave Packet: "
                "apmlitude = %0.4f%+0.4fi, "
                "speed = %0.4f",
                creal(A_0), cimag(A_0), c);

    debug_body( "Gaussian Wave Packet: "
                "constant potential = %0.4f%+0.4fi, "
                "time = %0.4f",
                creal(phi_0), cimag(phi_0), t);

    matlib_complex C1 = 1.0/(4.0*a*t*I+1.0);
    matlib_complex C2 = A_0*csqrt(C1)*cexp(0.25*I*c*c*t + I*phi_0*t);
    matlib_real C3 = c*t;
    matlib_real C4 = 0.5*c;

    matlib_real *xptr;
    matlib_real Xi, x1;
    
    matlib_real theta  = mu*t;
    matlib_real beta   = g_0*sin(theta)/mu;
    matlib_real nu     = 2.0*g_0*(cos(theta)-1.0)/(mu*mu);
    matlib_real ibeta2 = g_0*g_0*(theta/2.0 - sin(2*theta)/4.0)/(mu*mu*mu);

    for (xptr = x.elem_p; xptr<(x.elem_p+x.len); xptr++)
    {
        Xi = beta**xptr - ibeta2;
        x1 = *xptr + nu - C3;
        *(u.elem_p) = C2*cexp(-a*C1*x1*x1 + I*C4*x1 + I*Xi);
        u.elem_p++;
    }
    debug_exit("%s", "");
}

void pde1d_LSE_timedependent_linear_potential
( 
    void**         params,
    matlib_complex m_coeff[2], 
    matlib_xv      x, 
    matlib_xv      t, 
    matlib_zm      y
)
{
    debug_enter("%s", "");
    matlib_index i, j;

    matlib_complex phi_0 = *(matlib_complex*) params[3];
    matlib_real g_0 = *(matlib_real*) params[4];
    matlib_real mu  = *(matlib_real*) params[5];

    debug_body("g_0: %0.4f, mu: %0.4f", g_0, mu);
    assert(t.len == (y.lenr+1));
    
    matlib_complex a = 0, b = 0;
    if(m_coeff != NULL)
    {
        a = m_coeff[0];
        b = m_coeff[1];
    }

    for(j=0; j<y.lenr; j++)
    {
        for (i=0; i<x.len; i++)
        {
            *(y.elem_p) =  a + b*phi_0 + 0.5*b**(x.elem_p+i)*g_0*
                           (cos(mu*(t.elem_p[j]))+cos(mu*(t.elem_p[j+1])));
            y.elem_p++;
        }
    }
    debug_exit("%s", "");
}
