#ifndef MATLIB_H
#define MATLIB_H

/*============================================================================+/
 | Include all the dependencies
/+============================================================================*/
#include "complex.h"
#include "basic.h"
#include "debug.h"
#include "ehandler.h"


#ifdef MATLIB_NTRACE_DATA
    #define _TRACE_DATA 0
#else
    #define _TRACE_DATA 1
#endif


/* Writing data tracing blocks */
#define BEGIN_DTRACE if(_TRACE_DATA){
#define END_DTRACE }



/*============================================================================*/

typedef unsigned int   matlib_index;
typedef double complex matlib_complex;

#define _CBLAS_ORDER     CBLAS_ORDER
#define _CBLAS_TRANSPOSE CBLAS_TRANSPOSE


typedef enum
{
    MATLIB_VECT_NONE,
    MATLIB_COL_VECT,
    MATLIB_ROW_VECT

} MATLIB_VECT_T;

typedef enum
{
    MATLIB_NO_TRANS,
    MATLIB_TRANS,
    MATLIB_CONJ_TRANS

} MATLIB_TRANSPOSE;


/* ORDER: 
 * COL_MAJOR: elements within a column are contiguous in memory
 * ROW_MAJOR: elements  with in a row are contiguous in memory
 * 
 * Matrix elements are stored in a array as 
 * A[i][j]--> COL_MAJOR: *(pA+i+stride*j), ROW_MAJOR: *(pA+i*stride+j)
 * "stride" here represents the constant spacing two row/column elements
 * For COL_MAJOR, stride>= length of columns (lenc)
 * For ROW_MAJOR, stride>= length of rows    (lenr)
 * 
 * stride corresponds to leading dimension argument of arrays in FOTRAN
 * All matrix data types declared in this file assume stride=lenc or lenr unless
 * otherwise stated.
 *
 * */ 
typedef enum
{
    MATLIB_ORDER_NONE,
    MATLIB_COL_MAJOR, 
    MATLIB_ROW_MAJOR 

} MATLIB_ORDER;


/* Data types for matrix library 
 *
 * DV : Double Vector
 * DM : Double Matrix
 *
 * */ 
typedef struct
{
    matlib_index  len;
    MATLIB_VECT_T type;
    double*       elem_p;

} matlib_dv;

typedef struct
{
    matlib_index    len;
    MATLIB_VECT_T   type;
    matlib_complex* elem_p;

} matlib_zv;

typedef struct
{
    matlib_index         len;
    MATLIB_VECT_T type;
    double*       elem_pr;
    double*       elem_pi;

} MATLIB_TV;

/* The operation field is to tell if the meaningful entries belong to the current
 * matrix or the transposed version of it. An actual taranspose operation is
 * carried out only if strictly needed.
 * */ 
typedef struct
{
    matlib_index        lenc; /* length of columns */ 
    matlib_index        lenr; /* length of rows    */ 
    MATLIB_ORDER order;
    MATLIB_TRANSPOSE op; /* operation */ 
    double*      elem_p;

} matlib_dm;

typedef struct
{
    matlib_index           lenc; /* length of columns */ 
    matlib_index           lenr; /* length of rows    */ 
    MATLIB_ORDER    order;
    MATLIB_TRANSPOSE op; /* operation */ 
    matlib_complex* elem_p;

} matlib_zm;

/* A Tuple of real and imaginary parts */ 
typedef struct
{
    matlib_index        lenc; /* length of columns */ 
    matlib_index        lenr; /* length of rows    */ 
    MATLIB_ORDER order;
    MATLIB_TRANSPOSE op; /* operation */ 
    double*      elem_pr;
    double*      elem_pi;

} MATLIB_TM;
/*============================================================================*/
/* SPARSE FORMATS */
typedef enum
{
    MATLIB_CSR3,
    MATLIB_CSC3, 
    MATLIB_COO,
    MATLIB_DIA

} MATLIB_SPARSE;

typedef enum
{
    MATLIB_UPPER,
    MATLIB_LOWER

} MATLIB_UPLO;

typedef struct
{
    matlib_index         lenc; /* length of columns */ 
    matlib_index         lenr; /* length of rows    */
    matlib_index*        rowIn;
    matlib_index*        colIn;
    MATLIB_SPARSE format;
    double*       elem_p;

} matlib_dsparsem;

typedef struct
{
    matlib_index           lenc; /* length of columns */ 
    matlib_index           lenr; /* length of rows    */
    matlib_index*          rowIn;
    matlib_index*          colIn;
    MATLIB_SPARSE   format;
    matlib_complex* elem_p;

} matlib_zsparsem;

/* This data structure is meant to store several sparse matrices which have the
 * same sparsity structure. This proves useful for solving initial value
 * problems where the potential is a time dependent function. The mass matrix
 * for several time steps can be computed and stored at once.*/ 
typedef struct
{
    matlib_index         lenc; /* length of columns */ 
    matlib_index         lenr; /* length of rows    */
    matlib_index*        rowIn;
    matlib_index*        colIn;
    MATLIB_SPARSE format;
    matlib_index         nsparse;
    double**      elem_p;

} matlib_dnsparsem;

typedef struct
{
    matlib_index           lenc; /* length of columns */ 
    matlib_index           lenr; /* length of rows    */
    matlib_index*          rowIn;
    matlib_index*          colIn;
    MATLIB_SPARSE   format;
    matlib_index           nsparse;
    matlib_complex** elem_p;

} matlib_znsparsem;


/*============================================================================*/
/* Define MACROS */ 

#define matlib_free(ptr)    \
    do{ if(ptr!=NULL)       \
        free((void*) ptr);  \
      } while (0)

#define DEBUG_PRINT_DV(NAME, fmt,...)                                    \
    do {                                                                 \
         if(_TRACE_DATA){                                          \
            matlib_index NAME ## _i;                                     \
            for((NAME ## _i)=0;                                          \
                (NAME ## _i)<(NAME.len);                                 \
                (NAME ##_i)++){                                          \
                fprintf( stderr,                                         \
                         "%s:%d:%s: -(" fmt #NAME "[%d] = % 0.16f)\n",   \
                         __FILE__,                                       \
                         __LINE__,                                       \
                         __func__,                                       \
                        __VA_ARGS__,                                     \
                        NAME ## _i,                                      \
                        *((NAME.elem_p)+(NAME ## _i)));                  \
            }                                                            \
         }                                                               \
    } while (0)           

/* For complex Vectors */ 
#define DEBUG_PRINT_ZV(NAME, fmt,...)                                    \
    do {                                                                 \
         if(_TRACE_DATA){                                          \
            matlib_index NAME ## _i;                                     \
            for((NAME ## _i)=0;                                          \
                (NAME ## _i)<(NAME.len);                                 \
                (NAME ##_i)++){                                          \
                fprintf( stderr,                                         \
                    "%s:%d:%s: -(" fmt #NAME "[%d] = % 0.16f %+0.16fi)\n",\
                         __FILE__,                                       \
                         __LINE__,                                       \
                         __func__,                                       \
                        __VA_ARGS__,                                     \
                        NAME ## _i,                                      \
                        creal(*((NAME.elem_p)+(NAME ## _i))),            \
                        cimag(*((NAME.elem_p)+(NAME ## _i))));           \
            }                                                            \
         }                                                               \
    } while (0)                                                          

/* Debuging a double matrix */ 
#define DEBUG_PRINT_DM(NAME, fmt,...)                                    \
    do {                                                                 \
         if(_TRACE_DATA){                                          \
            matlib_index NAME ## _i, NAME ## _j;                         \
            matlib_index NAME ## _col_st, NAME ## _row_st;               \
            if ((NAME.order) == MATLIB_COL_MAJOR){                       \
                (NAME ## _col_st) = 1;                                   \
                (NAME ## _row_st) = (NAME.lenc);                         \
            }                                                            \
            else{                                                        \
                (NAME ## _col_st) = (NAME.lenr);                         \
                (NAME ## _row_st) = 1;                                   \
            }                                                            \
            for((NAME ## _i)=0;                                          \
                (NAME ## _i)<(NAME.lenc);                                \
                (NAME ##_i)++){                                          \
                for((NAME ## _j)=0;                                      \
                    (NAME ## _j)<(NAME.lenr);                            \
                    (NAME ## _j)++){                                     \
                    fprintf( stderr,                                     \
                         "%s:%d:%s: -(" fmt #NAME "[%d][%d] = % 0.16f)\n",\
                             __FILE__,                                   \
                             __LINE__,                                   \
                             __func__,                                   \
                             __VA_ARGS__,                                \
                            NAME ## _i, NAME ## _j,                      \
                            *((NAME.elem_p)+                             \
                            (NAME ## _i)*(NAME ## _col_st)+              \
                            (NAME ## _j)*(NAME ## _row_st)));            \
                }                                                        \
            }                                                            \
         }                                                               \
    } while (0)                                                          


/* Debuging a complex matrix */ 
#define DEBUG_PRINT_ZM(NAME, fmt,...)                                    \
    do {                                                                 \
         if(_TRACE_DATA){                                          \
            matlib_index NAME ## _i, NAME ## _j;                         \
            matlib_index NAME ## _col_st, NAME ## _row_st;               \
            if ((NAME.order) == MATLIB_COL_MAJOR){                       \
                (NAME ## _col_st) = 1;                                   \
                (NAME ## _row_st) = (NAME.lenc);                         \
            }                                                            \
            else{                                                        \
                (NAME ## _col_st) = (NAME.lenr);                         \
                (NAME ## _row_st) = 1;                                   \
            }                                                            \
            for((NAME ## _i)=0;                                          \
                (NAME ## _i)<(NAME.lenc);                                \
                (NAME ##_i)++){                                          \
                for((NAME ## _j)=0;                                      \
                    (NAME ## _j)<(NAME.lenr);                            \
                    (NAME ## _j)++){                                     \
                    fprintf( stderr,                                     \
                             "%s:%d:%s: -(" fmt                          \
                             #NAME "[%d][%d] = % 0.16f %+0.16fi)\n",     \
                             __FILE__,                                   \
                             __LINE__,                                   \
                             __func__,                                   \
                             __VA_ARGS__,                                \
                            NAME ## _i, NAME ## _j,                      \
                            creal(*((NAME.elem_p)+                       \
                            (NAME ## _i)*(NAME ## _col_st)+              \
                            (NAME ## _j)*(NAME ## _row_st))),            \
                            cimag(*((NAME.elem_p)+                       \
                            (NAME ## _i)*(NAME ## _col_st)+              \
                            (NAME ## _j)*(NAME ## _row_st))));           \
                }                                                        \
            }                                                            \
         }                                                               \
    } while (0)                                                          

/* Convert a matrix into vector 
 *
 * V = MK_VM(B);
 *
 * */ 
#define MK_VM(B) { .len  = B.lenc * B.lenr, \
                   .type = MATLIB_COL_VECT, \
                   .elem_p = B.elem_p }


#define GET_CBLASORDER_AND_STRIDE(order_enum, stride, M)      \
    do {                                                      \
        if (M.order == MATLIB_COL_MAJOR)                      \
        {                                                     \
            stride = M.lenc;                                  \
            order_enum = CblasColMajor;                       \
        }                                                     \
        else if (M.order == MATLIB_ROW_MAJOR)                 \
        {                                                     \
            stride = M.lenr;                                  \
            order_enum = CblasRowMajor;                       \
        }                                                     \
        else                                                  \
        {                                                     \
            term_execb( "Storage order of the matrix " #M     \
                        " is unknown: (order enum: %d)",      \
                        M.order);                             \
        }                                                     \
                                                              \
    } while (0)

 
#define MATLIB_TRANSPOSE(M)                               \
do {                                                      \
        matlib_index M ##_dim = M.lenc;                   \
        M.lenc = M.lenr;                                  \
        M.lenr = M ## _dim;                               \
        if (M.order == MATLIB_COL_MAJOR)                  \
        {                                                 \
            M.order = MATLIB_ROW_MAJOR;                   \
        }                                                 \
        else if (M.order == MATLIB_ROW_MAJOR)             \
        {                                                 \
            M.order = MATLIB_COL_MAJOR;                   \
        }                                                 \
        else                                              \
        {                                                 \
            term_execb( "Storage order of the matrix " #M \
                        " is unknown: (order enum: %d)",  \
                        M.order);                         \
        }                                                 \
                                                          \
} while (0)                                                        
/*============================================================================*/
void matlib_create_zv
(
    matlib_index         length,
    matlib_zv*    v,
    MATLIB_VECT_T type_enum
);
void matlib_create_tv
(
    matlib_index         length,
    MATLIB_TV*    v,
    MATLIB_VECT_T type_enum
);
void matlib_create_dv
(
    matlib_index           length,
    matlib_dv*    v,
    MATLIB_VECT_T type_enum
);

void matlib_create_zm
( 
    matlib_index        lenc,
    matlib_index        lenr,
    matlib_zm*   M,
    MATLIB_ORDER order_enum,
    MATLIB_TRANSPOSE trans_enum
);
void matlib_create_tm
( 
    matlib_index        lenc,
    matlib_index        lenr,
    MATLIB_TM*   M,
    MATLIB_ORDER order_enum,
    MATLIB_TRANSPOSE trans_enum
);
void matlib_create_dm
( 
    matlib_index        lenc,
    matlib_index        lenr,
    matlib_dm*   M,
    MATLIB_ORDER order_enum,
    MATLIB_TRANSPOSE trans_enum
);
/*============================================================================*/

double matlib_dnrm2(matlib_dv x);
double matlib_znrm2(matlib_zv x);
double matlib_tnrm2(MATLIB_TV x);

void matlib_daxpy
(
    const double    alpha,
    const matlib_dv x,
          matlib_dv y
);
void matlib_taxpy
(
    const matlib_complex alpha,
    const MATLIB_TV      x,
          MATLIB_TV      y
);
void matlib_zaxpy
(
    const matlib_complex alpha,
    const matlib_zv      x,
          matlib_zv      y
);
void matlib_zaxpby
(
    const matlib_complex alpha,
    const matlib_zv      x,
    const matlib_complex beta,
          matlib_zv      y
);
double matlib_ddot
(
    const matlib_dv x,
    const matlib_dv y
);

void matlib_dgemv
(
    const double    alpha,
    const matlib_dm A, 
    const matlib_dv u,
    const double    beta,
          matlib_dv v
);
void matlib_zgemv
(
    const matlib_complex alpha,
    const matlib_zm      A, 
    const matlib_zv      u,
    const matlib_complex beta,
          matlib_zv      v
);

void matlib_dcsrsymv
/* Double CSR Symmetric Matrix-Vector */ 
(
    const MATLIB_UPLO     uplo_enum,
          matlib_dsparsem A, 
    const matlib_dv       u,
          matlib_dv       v
);

void matlib_zcsrsymv
/* Double CSR Symmetric Matrix-Vector */ 
(
    const MATLIB_UPLO     uplo_enum,
          matlib_zsparsem A, 
    const matlib_zv       u,
          matlib_zv       v
);
/*============================================================================*/

void matlib_dgemm
(
    const double    alpha,
    const matlib_dm A, 
    const matlib_dm B, 
    const double    beta,
          matlib_dm C
);
void matlib_zgemm
(
    const matlib_complex alpha,
    const matlib_zm      A, 
    const matlib_zm      B, 
    const matlib_complex beta,
          matlib_zm      C
);

/*============================================================================+/
 | IO functions for vectors and matrices
/+============================================================================*/

void matlib_dmwrite_csv(char* file_name, matlib_dm M);
void matlib_zmwrite_csv(char* file_name, matlib_zm M);

void matlib_dvwrite_csv
(
    char*        file_name, 
    matlib_index n,
    matlib_dv    v[n]
);

void matlib_zvwrite_csv
(
    char*        file_name, 
    matlib_index n,
    matlib_zv    v[n]
);

void matlib_dzvwrite_csv
(
    char*        file_name, 
    matlib_index m,
    matlib_dv    u[m],
    matlib_index n,
    matlib_zv    v[n]
);
/*============================================================================+/
 | Linear Solver Interface
/+============================================================================*/

/*=====================[Linear solver based on PARDISO]=======================*/
#define PARDISO_NIPARAM (64)
typedef enum
{
    PARDISO_RHS,
    PARDISO_LHS

} PARDISO_SOLVEC;

typedef enum
{
    PARDISO_INIT                = -2,
    PARDISO_FREE                = -1, /* release all memory */ 
    PARDISO_FREE_LU             = 0,
    PARDISO_ANALYSIS_AND_FACTOR = 12,
    PARDISO_SOLVE_AND_REFINE    = 33

} PARDISO_PHASE;

typedef enum
{
    PARDISO_REAL_SYM_INDEF = -2,
    PARDISO_REAL_SYM_PDEF  =  2,
    PARDISO_COMPLEX_SYM    =  6

} PARDISO_MTYPE;

typedef struct
{
    void*          ptr[PARDISO_NIPARAM];
    int            iparam[PARDISO_NIPARAM];
    PARDISO_PHASE  phase_enum;
    PARDISO_MTYPE  mtype;
    PARDISO_SOLVEC sol_enum;
    matlib_index   nsparse; /* number of sparse matrices with same 
                                sparsity structure */ 
    matlib_index   mnum;   /* matrix number to be used for solution */
    void*          smat_p; /* Sparse matrix struct */ 
    void*          rhs_p;  /* vector struct        */ 
    void*          sol_p;  /* vector struct        */ 

} pardiso_solver_t;

void matlib_pardiso(pardiso_solver_t* data);

/*============================================================================*/


#endif

