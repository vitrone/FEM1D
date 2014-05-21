#ifndef PDE1D_SOLVER_H
#define  PDE1D_SOLVER_H

/*============================================================================+/
 | Include all the dependencies
/+============================================================================*/
#include "legendre.h"
#include "fem1d.h"

/*============================================================================+/
 | Linear Schroedinger Equation (LSE)
 | iu_t + alpha u_xx + phi(x,t)u = 0
 | Domain : [x_l, x_r]
 | 
 | Boundary Conditions
 |
 | DIRICHLET: u(x_l,t) = u_L, u(x_r, t) = u_R.
 |
 | NEUMANN  : u_x(x_l,t) = u1_L, u_x(x_r, t) = u1_R.
 |
 | PERIODIC : u(x_l,t) = u(x_r,t), u_x(x_l,t) = u_x(x_r,t).
/+============================================================================*/

typedef enum
{
    PDE1D_LSE_DIRICHLET,
    PDE1D_LSE_ROBIN,

} PDE1D_LSE_BC;


typedef enum
{
    PDE1D_LSE_STATIC,
    PDE1D_LSE_DYNAMIC

} PDE1D_LSE_POTENTIAL;

typedef enum
{
    PDE1D_LSE_EVOLVE_ONLY,
    PDE1D_LSE_EVOLVE_ERROR,
    PDE1D_LSE_ERROR_ONLY

} PDE1D_LSE_SOLVE;

typedef struct
{
    matlib_index p;         /* solution space : H^p */
    matlib_index nr_LGL;    /* LGL-quadrature on reference domain [-1, 1] */ 
    matlib_index N;         /* number of finite elements */ 
    matlib_real  domain[2]; /* computational domain */ 
    matlib_real  dt;        /* size of time-step */
    matlib_index Nt;        /* number of time-steps */
    matlib_index nsparse;   /* number of sparse matrices for dynamic problem */
    void**       params;    /* Could be anything! */ 
    matlib_complex alpha;   /* Coeff of u_xx */

    matlib_zv u_init;       /* Initial Condition */
    matlib_zm U_evol;       /* Computed soltuion in Legendre basis */
    void*     u_analytic;   /* analytic solution */
    matlib_xv x;
    matlib_xv t;

    PDE1D_LSE_POTENTIAL phi_type;
    void* phix_p;  /* potential function x-dependent part */ 
    void* phixt_p; /* potential x- and t-parts non-separable */

    matlib_real  tol;

    matlib_xv e_abs;
    matlib_xv e_rel;

    PDE1D_LSE_SOLVE sol_mode;

} pde1d_LSE_data_t;

typedef struct
{
    matlib_xv xi;
    matlib_xm FM;
    matlib_xm IM;
    matlib_xv quadW;
    matlib_xm Q;
    matlib_zm q;

    /* Sparse matrics 
     * */ 
    matlib_zm_sparse  M;  /* Stiffness + Mass excluding all time-dependent part */ 

    /* Sparse Matrics - Global Mass Matrices
     * */ 
    matlib_zm_nsparse nM; /* mixed potentials */ 

    matlib_real rho;
    matlib_real irho;
    matlib_real J;
    matlib_complex s_coeff;
    matlib_complex m_coeff[2];

    matlib_index nr_vars;
    void**    var_p; /* for all temporary variables needed for solution */ 

} pde1d_LSE_solver_t;

/*============================================================================*/

void pde1d_LSE_set_defaultsIVP(pde1d_LSE_data_t* input);

void pde1d_LSE_set_potential
(
    pde1d_LSE_data_t*  input,
    PDE1D_LSE_POTENTIAL phi_type,
    void* phi_p
);

void pde1d_analytic_evol
(
    void** params,
    void (*u_analytic)(),
    matlib_xv x,
    matlib_xv t,
    matlib_zm u
);

void pde1d_LSE_init_solverIVP
(
    pde1d_LSE_data_t* input, 
    pde1d_LSE_solver_t* data 
);

void fem1d_linespace
(
    matlib_real  t0, 
    matlib_real  dt, 
    matlib_index Nt, 
    matlib_xv    *t
);

void fem1d_zm_sparse_GSM
/* Global Stiffness Matrix */ 
(
    matlib_index     N,
    matlib_complex   coeff,
    matlib_zm_sparse M
);

void pde1d_zm_nsparse_GSM
/* Global Stiffness Matrix */ 
(
    matlib_index      N,
    matlib_complex    coeff,
    matlib_zm_nsparse M
);

void pde1d_LSE_solve_IVP
(
    pde1d_LSE_data_t*  input,
    pde1d_LSE_solver_t* data
);
void pde1d_LSE_solve_IVP_evol
(
    pde1d_LSE_data_t*  input,
    pde1d_LSE_solver_t* data
);
void pde1d_LSE_solve_IVP_error
(
    pde1d_LSE_data_t*  input,
    pde1d_LSE_solver_t* data
);

void pde1d_LSE_solve_IVP2_evol
(
    pde1d_LSE_data_t*   input,
    pde1d_LSE_solver_t* data
);
void pde1d_LSE_solve_IVP2_error
(
    pde1d_LSE_data_t*   input,
    pde1d_LSE_solver_t* data
);
void pde1d_LSE_destroy_solverIVP
(
    pde1d_LSE_data_t*   input,
    pde1d_LSE_solver_t* data
);

/*============================================================================*/
void pde1d_LSE_Gaussian_WP_constant_potential
(
    void** params,
    matlib_xv   x,
    matlib_real t,
    matlib_zv   u
);

void pde1d_LSE_Gaussian_CW_constant_potential
(
    void** params,
    matlib_xv   x,
    matlib_real t,
    matlib_zv   u
);

void pde1d_LSE_constant_potential
( 
    void** params,
    matlib_complex m_coeff[2], 
    matlib_xv      x, 
    matlib_zv      y
);
/*============================================================================*/

void pde1d_LSE_HermiteGaussian_WP_harmonic_potential
(
    void** params,
    matlib_xv   x,
    matlib_real t,
    matlib_zv   u
);

void pde1d_LSE_harmonic_potential
( 
    void** params,
    matlib_complex m_coeff[2], 
    matlib_xv      x, 
    matlib_zv      y
);
/*============================================================================*/

void pde1d_LSE_Gaussian_WP_timedependent_linear_potential
(
    void** params,
    matlib_xv   x,
    matlib_real t,
    matlib_zv   u
);

void pde1d_LSE_timedependent_linear_potential
( 
    void**         params,
    matlib_complex m_coeff[2], 
    matlib_xv      x, 
    matlib_xv      t, 
    matlib_zm      y
);

#endif /* PDE1D_SOLVER_H */
