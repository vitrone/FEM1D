/*============================================================================+/
 | pfem1d.c
 | Library defines pthreaded fem1d functions. It uses pthpool.c library.
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
#include <unistd.h>
#include "mkl.h"

#define NDEBUG
#define MATLIB_NTRACE_DATA
#include "assert.h"
#include "pfem1d.h"
/*============================================================================*/

static void* pfem1d_thfunc_XFLT(void* mp);
static void* pfem1d_thfunc_ZFLT(void* mp);
static void* pfem1d_thfunc_XILT(void* mp);
static void* pfem1d_thfunc_ZILT(void* mp);
static void* pfem1d_thfunc_XFLT2(void* mp);
static void* pfem1d_thfunc_ZFLT2(void* mp);
static void* thfunc_dshapefunc2lp(void* mp);
static void* thfunc_zshapefunc2lp(void* mp);


/* Printed with precision 0.20f.*/
#ifndef FEM1D_H
    #define _A01 -0.40824829046386307274
    #define _A02 -0.31622776601683794118
    #define _A03 0.40824829046386307274
    #define _B03 -0.26726124191242439654
    #define _A04 0.31622776601683794118
    #define _B04 -0.23570226039551586683
    #define _A05 0.26726124191242439654
    #define _B05 -0.21320071635561041457
    #define _A06 0.23570226039551586683
    #define _B06 -0.19611613513818404453
    #define _A07 0.21320071635561041457
    #define _B07 -0.18257418583505535814
    #define _A08 0.19611613513818404453
    #define _B08 -0.17149858514250881925
    #define _A09 0.18257418583505535814
    #define _B09 -0.16222142113076254422
    #define _A10 0.17149858514250881925
    #define _A11 0.16222142113076254422

    #define _E00 0.16329931618554521799
    #define _F00 -0.81649658092772614548
    #define _E01 0.09035079029052513200
    #define _F01 -0.21081851067789195153
    #define _E02 0.05939138709164986513
    #define _F02 -0.10690449676496975584
    #define _E03 0.04285495643554834005
    #define _F03 -0.06734350297014739251
    #define _E04 0.03280011020855545106
    #define _F04 -0.04737793696791342546
    #define _E05 0.02614881801842454043
    #define _F05 -0.03565747911603345949
    #define _E06 0.02147931598059474659
    #define _F06 -0.02808833628231621055
    #define _E07 0.01805248264657987461
    #define _F07 -0.02286647801900117474
    #define _E08 0.01544965915531071841
    #define _F08 -0.01908487307420735773

    #define _N0 2.00000000000000000000
    #define _N1 0.66666666666666662966
    #define _N2 0.40000000000000002220
    #define _N3 0.28571428571428569843
    #define _N4 0.22222222222222220989
    #define _N5 0.18181818181818182323
    #define _N6 0.15384615384615385469
    #define _N7 0.13333333333333333148
    #define _N8 0.11764705882352941013
    #define _N9 0.10526315789473683626
    #define _N10 0.09523809523809523281

#endif

static void* thfunc_dshapefunc2lp_2 ( void* mp);
static void* thfunc_dshapefunc2lp_3 ( void* mp);
static void* thfunc_dshapefunc2lp_4 ( void* mp);
static void* thfunc_dshapefunc2lp_5 ( void* mp);
static void* thfunc_dshapefunc2lp_6 ( void* mp);
static void* thfunc_dshapefunc2lp_7 ( void* mp);
static void* thfunc_dshapefunc2lp_8 ( void* mp);
static void* thfunc_dshapefunc2lp_9 ( void* mp);
static void* thfunc_dshapefunc2lp_10( void* mp);

static void* thfunc_zshapefunc2lp_2 ( void* mp);
static void* thfunc_zshapefunc2lp_3 ( void* mp);
static void* thfunc_zshapefunc2lp_4 ( void* mp);
static void* thfunc_zshapefunc2lp_5 ( void* mp);
static void* thfunc_zshapefunc2lp_6 ( void* mp);
static void* thfunc_zshapefunc2lp_7 ( void* mp);
static void* thfunc_zshapefunc2lp_8 ( void* mp);
static void* thfunc_zshapefunc2lp_9 ( void* mp);
static void* thfunc_zshapefunc2lp_10( void* mp);

static void* thfunc_dprjLP2FEM_ShapeFunc(void* mp);
static void* thfunc_dprjLP2FEM_ShapeFunc(void* mp);

static void* thfunc_dprjLP2FEM_ShapeFunc_2 (void* mp);
static void* thfunc_dprjLP2FEM_ShapeFunc_3 (void* mp);
static void* thfunc_dprjLP2FEM_ShapeFunc_4 (void* mp);
static void* thfunc_dprjLP2FEM_ShapeFunc_5 (void* mp);
static void* thfunc_dprjLP2FEM_ShapeFunc_6 (void* mp);
static void* thfunc_dprjLP2FEM_ShapeFunc_7 (void* mp);
static void* thfunc_dprjLP2FEM_ShapeFunc_8 (void* mp);
static void* thfunc_dprjLP2FEM_ShapeFunc_9 (void* mp);
static void* thfunc_dprjLP2FEM_ShapeFunc_10(void* mp);

static void* thfunc_zprjLP2FEM_ShapeFunc_2 (void* mp);
static void* thfunc_zprjLP2FEM_ShapeFunc_3 (void* mp);
static void* thfunc_zprjLP2FEM_ShapeFunc_4 (void* mp);
static void* thfunc_zprjLP2FEM_ShapeFunc_5 (void* mp);
static void* thfunc_zprjLP2FEM_ShapeFunc_6 (void* mp);
static void* thfunc_zprjLP2FEM_ShapeFunc_7 (void* mp);
static void* thfunc_zprjLP2FEM_ShapeFunc_8 (void* mp);
static void* thfunc_zprjLP2FEM_ShapeFunc_9 (void* mp);
static void* thfunc_zprjLP2FEM_ShapeFunc_10(void* mp);

static void* thfunc_dlp_snorm2_d(void* mp);
static void* thfunc_zlp_snorm2_d(void* mp);

static void* thfunc_dlp_snorm2_d_2(void* mp);
static void* thfunc_dlp_snorm2_d_3(void* mp);
static void* thfunc_dlp_snorm2_d_4(void* mp);
static void* thfunc_dlp_snorm2_d_5(void* mp);
static void* thfunc_dlp_snorm2_d_6(void* mp);
static void* thfunc_dlp_snorm2_d_7(void* mp);
static void* thfunc_dlp_snorm2_d_8(void* mp);
static void* thfunc_dlp_snorm2_d_9(void* mp);
static void* thfunc_dlp_snorm2_d_10(void* mp);

static void* thfunc_zlp_snorm2_d_2(void* mp);
static void* thfunc_zlp_snorm2_d_3(void* mp);
static void* thfunc_zlp_snorm2_d_4(void* mp);
static void* thfunc_zlp_snorm2_d_5(void* mp);
static void* thfunc_zlp_snorm2_d_6(void* mp);
static void* thfunc_zlp_snorm2_d_7(void* mp);
static void* thfunc_zlp_snorm2_d_8(void* mp);
static void* thfunc_zlp_snorm2_d_9(void* mp);
static void* thfunc_zlp_snorm2_d_10(void* mp);

static void* pfem1d_thfunc_XCSRGMM2(void* mp);
static void* pfem1d_thfunc_ZCSRGMM2(void* mp);

/*============================================================================*/
static void* pfem1d_thfunc_XFLT(void* mp)
{

    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;
    matlib_index N = *((matlib_index*) (ptr->shared_data[0]));
    matlib_xm FM   =    *((matlib_xm*) (ptr->shared_data[1]));
    matlib_xv u    =    *((matlib_xv*) (ptr->shared_data[2]));
    matlib_xv U    =    *((matlib_xv*) (ptr->shared_data[3]));
    
    debug_enter( "thread index: %d, "
                 "nr. finite-elements: %d, "
                 "matrix FM: %d-by-%d, "
                 "vectors u: %d, U:%d", 
                 ptr->thread_index, N, FM.lenc, FM.lenr, u.len, U.len );

    /* highest degree of Legendre polynomials retained 
     * */
    matlib_index p = FM.lenc-1; 
    
    /* nr. LGL-points -1 used for sampling in FEM-element 
     * */ 
    matlib_index P = FM.lenr-1; 
    
    debug_body( "thread index: %d, "
                "degree of polynomial: %d,"
                "nr. sampling point: %d", 
                ptr->thread_index, 
                p, FM.lenr );

    assert(((FM.elem_p!=NULL) && (u.elem_p !=NULL)) && (U.elem_p!=NULL));

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);
    if(P>1)
    {
        debug_body("nr. computed finite-elements: %d", (u.len-1)/P);
        assert(N == (u.len-1)/P);
        if (U.len == N*(p+1))
        {
            matlib_index strideFM;
            _CBLAS_ORDER order_enum;
            GET_CBLASORDER_AND_STRIDE(order_enum, strideFM, FM);

            matlib_real one = 1.0, zero = 0;
            matlib_index incu = 1;
            matlib_index i;

            (u.elem_p) += (P*start_end_index[0]);
            (U.elem_p) += (FM.lenc*start_end_index[0]);
            for(i=start_end_index[0]; i<start_end_index[1]; i++)
            {
                cblas_dgemv( order_enum, 
                             CblasNoTrans, 
                             FM.lenc, 
                             FM.lenr, 
                             one, 
                             FM.elem_p, 
                             strideFM, 
                             u.elem_p,
                             incu,
                             zero, 
                             U.elem_p, 
                             incu);
                (u.elem_p) += P;
                (U.elem_p) += (FM.lenc);
            }
        }
        else
        {
            term_execb( "size of vectors/matrices incorrect: matrix "
                        "FM: %d-by-%d, vectors u: %d, U:%d", 
                        FM.lenc, FM.lenr, u.len, U.len );
        }
    }
    else
    {
        term_exec( "incorrect nr. of sampling points in ref. domain: %d",
                    FM.lenr );
    }
    debug_exit("Thread id: %d", ptr->thread_index);
}

void pfem1d_XFLT
(
    const matlib_index    N,
    const matlib_xm       FM,
          matlib_xv       u,
          matlib_xv       U,
          matlib_index    num_threads,
          pthpool_data_t* mp
)
{
    matlib_index i; 
    /* define the shared data */ 
    void* shared_data[4] = { (void*) &N,
                             (void*) &FM,
                             (void*) &u,
                             (void*) &U };

    matlib_index nsdata[num_threads][2];


    pthpool_arg_t   arg[num_threads];
    pthpool_task_t  task[num_threads];

    /* define the block of data per thread */ 
    matlib_index Np = N/(num_threads);

    for(i=0; i<num_threads-1; i++)
    {
        nsdata[i][0] = i*Np;
        nsdata[i][1] = (i+1)*Np;
        arg[i].shared_data    = shared_data; 
        arg[i].nonshared_data = (void**)&nsdata[i];
        arg[i].thread_index   = i;
        /* Define the task */ 
        task[i].function  = (void*)pfem1d_thfunc_XFLT;
        task[i].argument  = &arg[i];
        mp[i].task = &task[i];
    }
    nsdata[i][0] = i*Np;
    nsdata[i][1] = N;
    arg[i].shared_data    = shared_data; 
    arg[i].nonshared_data = (void**)&nsdata[i];
    arg[i].thread_index   = i;
    /* Define the task */ 
    task[i].function  = (void*)pfem1d_thfunc_XFLT;
    task[i].argument  = &arg[i];
    mp[i].task = &task[i];

    debug_body("%s", "created task");

    pthpool_exec_task(num_threads, mp, task);
    
    debug_exit("%s", "");

}
/*============================================================================*/
static void* pfem1d_thfunc_ZFLT(void* mp)
/* 
 * N : can be deduced from the length of the matrices provided but
 *     this leads to untraceble bugs. Hence this is made as an argument.
 * */ 
{
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;
    matlib_index N = *((matlib_index*) (ptr->shared_data[0]));
    matlib_xm FM   =    *((matlib_xm*) (ptr->shared_data[1]));
    matlib_zv u    =    *((matlib_zv*) (ptr->shared_data[2]));
    matlib_zv U    =    *((matlib_zv*) (ptr->shared_data[3]));

    debug_enter( "thread index: %d, "
                 "nr. finite-elements: %d, "
                 "matrix FM: %d-by-%d, "
                 "vectors u: %d, U:%d", 
                 ptr->thread_index, N, FM.lenc, FM.lenr, u.len, U.len );

    /* highest degree of Legendre polynomials retained 
     * */
    matlib_index p = FM.lenc-1; 
    
    /* nr. LGL-points -1 used for sampling in FEM-element 
     * */ 
    matlib_index P = FM.lenr-1; 
    
    debug_body( "degree of polynomial: %d,"
                "nr. sampling point: %d", 
                p, FM.lenr );

    assert(((FM.elem_p!=NULL) && (u.elem_p !=NULL)) && (U.elem_p!=NULL));

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);
    if(P>1)
    {
        debug_body("nr. computed finite-elements: %d", (u.len-1)/P);
        assert(N == (u.len-1)/P);
        if (U.len == N*(p+1))
        {
            matlib_index strideFM;
            _CBLAS_ORDER order_enum;
            GET_CBLASORDER_AND_STRIDE(order_enum, strideFM, FM);

            matlib_real one = 1.0, zero = 0;
            matlib_index incu = 2;
            matlib_index i;

            (u.elem_p) += (P*start_end_index[0]);
            (U.elem_p) += (FM.lenc*start_end_index[0]);
            for(i=start_end_index[0]; i<start_end_index[1]; i++)
            {
                cblas_dgemv( order_enum, 
                             CblasNoTrans, 
                             FM.lenc, 
                             FM.lenr, 
                             one, 
                             FM.elem_p, 
                             strideFM, 
                             (matlib_real*)u.elem_p,
                             incu,
                             zero, 
                             (matlib_real*)U.elem_p, 
                             incu);
                cblas_dgemv( order_enum, 
                             CblasNoTrans, 
                             FM.lenc, 
                             FM.lenr, 
                             one, 
                             FM.elem_p, 
                             strideFM, 
                             ((matlib_real*)u.elem_p+1),
                             incu,
                             zero, 
                             ((matlib_real*)U.elem_p+1), 
                             incu);
                (u.elem_p) += P;
                (U.elem_p) += (FM.lenc);
            }
        }
        else
        {
            term_execb( "size of vectors/matrices incorrect: matrix "
                        "FM: %d-by-%d, vectors u: %d, U:%d", 
                        FM.lenc, FM.lenr, u.len, U.len );
        }
    }
    else
    {
        term_exec( "incorrect nr. of sampling points in ref. domain: %d",
                    FM.lenr );
    }
    debug_exit("%s","");
}

void pfem1d_ZFLT
(
    const matlib_index    N,
    const matlib_xm       FM,
          matlib_zv       u,
          matlib_zv       U,
          matlib_index    num_threads,
          pthpool_data_t* mp
)
{
    matlib_index i; 
    /* define the shared data */ 
    void* shared_data[4] = { (void*) &N,
                             (void*) &FM,
                             (void*) &u,
                             (void*) &U };

    matlib_index nsdata[num_threads][2];


    pthpool_arg_t   arg[num_threads];
    pthpool_task_t  task[num_threads];

    /* define the block of data per thread */ 
    matlib_index Np = N/(num_threads);

    for(i=0; i<num_threads-1; i++)
    {
        nsdata[i][0] = i*Np;
        nsdata[i][1] = (i+1)*Np;
        arg[i].shared_data    = shared_data; 
        arg[i].nonshared_data = (void**)&nsdata[i];
        arg[i].thread_index   = i;
        /* Define the task */ 
        task[i].function  = (void*)pfem1d_thfunc_ZFLT;
        task[i].argument  = &arg[i];
        mp[i].task = &task[i];
    }
    nsdata[i][0] = i*Np;
    nsdata[i][1] = N;
    arg[i].shared_data    = shared_data; 
    arg[i].nonshared_data = (void**)&nsdata[i];
    arg[i].thread_index   = i;
    /* Define the task */ 
    task[i].function  = (void*)pfem1d_thfunc_ZFLT;
    task[i].argument  = &arg[i];
    mp[i].task = &task[i];

    debug_body("%s", "created task");

    pthpool_exec_task(num_threads, mp, task);
    
    debug_exit("%s", "");

}
/*============================================================================*/
static void* pfem1d_thfunc_XILT(void* mp)
/* 
 * N : can be deduced from the length of the matrices provided but
 *     this leads to untraceble bugs. Hence this is made as an argument.
 * */ 
{
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;
    matlib_index N = *((matlib_index*) (ptr->shared_data[0]));
    matlib_xm IM   = *((matlib_xm*) (ptr->shared_data[1]));
    matlib_xv U    = *((matlib_xv*) (ptr->shared_data[2]));
    matlib_xv u    = *((matlib_xv*) (ptr->shared_data[3]));
    
    pthread_mutex_t* lock_common = ((pthread_mutex_t*) (ptr->shared_data[4]));
    
    debug_enter( "thread index: %d, "
                 "nr. finite-elements: %d, "
                 "matrix IM: %d-by-%d, "
                 "vectors U: %d, u:%d", 
                 ptr->thread_index, N, IM.lenc, IM.lenr, U.len, u.len );

    /* nr. LGL-points-1 used for sampling in FEM-element 
     * */
    matlib_index P = IM.lenc-1; 
    
    /* Highest degree of Legendre polynomials
     * */
    matlib_index p = IM.lenr-1; 

    debug_body( "Thread id: %d, degree of polynomial: %d, "
                "nr. sampling point: %d", 
                ptr->thread_index, p, IM.lenc );

    assert(((IM.elem_p!=NULL) && (u.elem_p !=NULL)) && (U.elem_p!=NULL));

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data[0]);
    
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);
    bool not_first = !(start_end_index[0]==0);
    bool not_last  = !(start_end_index[1]==N);

    if(p>1)
    {
        debug_body( "Thread id: %d, nr. computed finite-elements: %d",
                    ptr->thread_index, U.len/(p+1));
        assert(N == U.len/(p+1));
        if (u.len == N*P+1)
        {
            matlib_index strideIM;
            _CBLAS_ORDER order_enum;
            
            GET_CBLASORDER_AND_STRIDE(order_enum, strideIM, IM);
            
            matlib_index incU = 1;
            matlib_index incu = 1;
            matlib_index i;

            (u.elem_p) += (P*start_end_index[0]);
            (U.elem_p) += (IM.lenr*start_end_index[0]);
            if(not_first)
            {
                pthread_mutex_lock(&lock_common[ptr->thread_index-1]);
            }
            cblas_dgemv( order_enum, 
                         CblasNoTrans, 
                         IM.lenc, 
                         IM.lenr,
                         1.0, 
                         IM.elem_p, 
                         strideIM, 
                         U.elem_p, 
                         incU, 
                         0,
                         u.elem_p,
                         incu);
            if(not_first)
            {
                matlib_real* common_u = (matlib_real*) (ptr->nonshared_data[1]);
                *common_u = *(u.elem_p);
                pthread_mutex_unlock(&lock_common[ptr->thread_index-1]);
            }
            (U.elem_p) += (IM.lenr);
            (u.elem_p) += P;
            for(i=start_end_index[0]+1; i<start_end_index[1]-1; i++)
            {
                cblas_dgemv( order_enum, 
                             CblasNoTrans, 
                             IM.lenc, 
                             IM.lenr,
                             1.0, 
                             IM.elem_p, 
                             strideIM, 
                             U.elem_p, 
                             incU, 
                             0,
                             u.elem_p,
                             incu);
                (U.elem_p) += (IM.lenr);
                (u.elem_p) += P;
            }
            if(not_last)
            {
                pthread_mutex_lock(&lock_common[ptr->thread_index]);
            }
            cblas_dgemv( order_enum, 
                         CblasNoTrans, 
                         IM.lenc, 
                         IM.lenr,
                         1.0, 
                         IM.elem_p, 
                         strideIM, 
                         U.elem_p, 
                         incU, 
                         0,
                         u.elem_p,
                         incu);
            if(not_last)
            {
                pthread_mutex_unlock(&lock_common[ptr->thread_index]);
            }
        }
        else
        {
            term_execb( "size of vectors/matrices incorrect:- IM: %d-by-%d, "
                        "U:%d, u:%d",
                        IM.lenc, IM.lenr, U.len, u.len );
        }
    }
    else
    {
        term_exec( "highest degree of polynomials incorrect: %d",
                    IM.lenr );
    }
    debug_exit("Thread id: %d", ptr->thread_index);
}
void pfem1d_XILT
(
    const matlib_index    N,
    const matlib_xm       IM,
          matlib_xv       U,
          matlib_xv       u,
          matlib_index    num_threads,
          pthpool_data_t* mp
)
{

    matlib_index i; 
    /* define mutexes for adjacent elements of u w.r.t. threads */
    pthread_mutex_t lock_common[num_threads-1];
    for(i=0; i<num_threads-1; i++)
    {
       pthread_mutex_init(&(lock_common[i]), NULL);
    }
    /* define the shared data */ 
    void* shared_data[5] = { (void*) &N,
                             (void*) &IM,
                             (void*) &U,
                             (void*) &u,
                             (void*) &(lock_common)};
    void* nonshared_data[num_threads][2];
    matlib_real common_u[num_threads-1];

    matlib_index nsdata[num_threads][2];


    pthpool_arg_t   arg[num_threads];
    pthpool_task_t  task[num_threads];

    /* define the block of data per thread */ 
    matlib_index Np = N/(num_threads);

    i = 0;
    nsdata[i][0] = i*Np;
    nsdata[i][1] = (i+1)*Np;
    arg[i].shared_data    = shared_data; 

    nonshared_data[i][0]  = (void*)nsdata[i];
    nonshared_data[i][1]  = NULL;
    arg[i].nonshared_data = (void**)&nonshared_data[i];
    arg[i].thread_index   = i;
    /* Define the task */ 
    task[i].function  = (void*)pfem1d_thfunc_XILT;
    task[i].argument  = &arg[i];
    mp[i].task = &task[i];
    for(i=1; i<num_threads-1; i++)
    {
        nsdata[i][0] = i*Np;
        nsdata[i][1] = (i+1)*Np;
        arg[i].shared_data    = shared_data; 
        nonshared_data[i][0]  = (void*)nsdata[i];
        nonshared_data[i][1]  = (void*)&common_u[i-1];
        arg[i].nonshared_data = (void**)&nonshared_data[i];
        arg[i].thread_index   = i;
        /* Define the task */ 
        task[i].function  = (void*)pfem1d_thfunc_XILT;
        task[i].argument  = &arg[i];
        mp[i].task = &task[i];
    }
    if(i<num_threads)
    {
        nsdata[i][0] = i*Np;
        nsdata[i][1] = N;
        arg[i].shared_data    = shared_data; 
        nonshared_data[i][0]  = (void*)nsdata[i];
        nonshared_data[i][1]  = (void*)&common_u[i-1];
        arg[i].nonshared_data = (void**)&nonshared_data[i];
        arg[i].thread_index   = i;
        /* Define the task */ 
        task[i].function  = (void*)pfem1d_thfunc_XILT;
        task[i].argument  = &arg[i];
        mp[i].task = &task[i];
    }
    debug_body("%s", "created task");

    pthpool_exec_task(num_threads, mp, task);
    matlib_index P = IM.lenc-1; 

    for(i=0; i<num_threads-1; i++)
    {
        u.elem_p[(i+1)*Np*P] = common_u[i];
    }

    debug_body("%s", "destroying mutexes");
    for(i=0; i<num_threads-1; i++)
    {
       pthread_mutex_destroy(&(lock_common[i]));
    }
    
    debug_exit("%s", "");
}
/*============================================================================*/
static void* pfem1d_thfunc_ZILT(void* mp)
/* 
 * N : can be deduced from the length of the matrices provided but
 *     this leads to untraceble bugs. Hence this is made as an argument.
 * */ 
{
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;
    matlib_index N = *((matlib_index*) (ptr->shared_data[0]));
    matlib_xm IM   =    *((matlib_xm*) (ptr->shared_data[1]));
    matlib_zv U    =    *((matlib_zv*) (ptr->shared_data[2]));
    matlib_zv u    =    *((matlib_zv*) (ptr->shared_data[3]));
    
    pthread_mutex_t* lock_common = ((pthread_mutex_t*) (ptr->shared_data[4]));

    debug_enter( "thread index: %d, "
                 "nr. finite-elements: %d, "
                 "matrix IM: %d-by-%d, "
                 "vectors U: %d, u:%d", 
                 ptr->thread_index, N, IM.lenc, IM.lenr, U.len, u.len );

    /* nr. LGL-points-1 used for sampling in FEM-element 
     * */
    matlib_index P = IM.lenc-1; 
    
    /* Highest degree of Legendre polynomials
     * */
    matlib_index p = IM.lenr-1; 

    debug_body( "degree of polynomial: %d, "
                "nr. sampling point: %d", 
                p, IM.lenc );

    assert(((IM.elem_p!=NULL) && (u.elem_p !=NULL)) && (U.elem_p!=NULL));

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data[0]);
    
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);
    bool not_first = !(start_end_index[0]==0);
    bool not_last  = !(start_end_index[1]==N);
    if(p>1)
    {
        debug_body("nr. computed finite-elements: %d", U.len/(p+1));
        assert(N == U.len/(p+1));
        if (u.len == N*P+1)
        {
            matlib_index strideIM;
            _CBLAS_ORDER order_enum;
            
            GET_CBLASORDER_AND_STRIDE(order_enum, strideIM, IM);
            
            matlib_index incu = 2;
            matlib_index i;
            (u.elem_p) += (P*start_end_index[0]);
            (U.elem_p) += (IM.lenr*start_end_index[0]);
            if(not_first)
            {
                pthread_mutex_lock(&lock_common[ptr->thread_index-1]);
            }
            cblas_dgemv( order_enum, 
                         CblasNoTrans, 
                         IM.lenc, 
                         IM.lenr,
                         1.0, 
                         IM.elem_p, 
                         strideIM, 
                         (matlib_real*)U.elem_p, 
                         incu, 
                         0,
                         (matlib_real*)u.elem_p,
                         incu);
            cblas_dgemv( order_enum, 
                         CblasNoTrans, 
                         IM.lenc, 
                         IM.lenr,
                         1.0, 
                         IM.elem_p, 
                         strideIM, 
                         ((matlib_real*)U.elem_p+1), 
                         incu, 
                         0,
                         ((matlib_real*)u.elem_p+1),
                         incu);
            if(not_first)
            {
                matlib_complex* common_u = (matlib_complex*) (ptr->nonshared_data[1]);
                *common_u = *(u.elem_p);
                pthread_mutex_unlock(&lock_common[ptr->thread_index-1]);
            }
            (U.elem_p) += (IM.lenr);
            (u.elem_p) += P;
            for(i=start_end_index[0]+1; i<start_end_index[1]-1; i++)
            {
                cblas_dgemv( order_enum, 
                             CblasNoTrans, 
                             IM.lenc, 
                             IM.lenr,
                             1.0, 
                             IM.elem_p, 
                             strideIM, 
                             (matlib_real*)U.elem_p, 
                             incu, 
                             0,
                             (matlib_real*)u.elem_p,
                             incu);
                cblas_dgemv( order_enum, 
                             CblasNoTrans, 
                             IM.lenc, 
                             IM.lenr,
                             1.0, 
                             IM.elem_p, 
                             strideIM, 
                             ((matlib_real*)U.elem_p+1), 
                             incu, 
                             0,
                             ((matlib_real*)u.elem_p+1),
                             incu);
                (U.elem_p) += (IM.lenr);
                (u.elem_p) += P;
            }
            if(not_last)
            {
                pthread_mutex_lock(&lock_common[ptr->thread_index]);
            }
            cblas_dgemv( order_enum, 
                         CblasNoTrans, 
                         IM.lenc, 
                         IM.lenr,
                         1.0, 
                         IM.elem_p, 
                         strideIM, 
                         (matlib_real*)U.elem_p, 
                         incu, 
                         0,
                         (matlib_real*)u.elem_p,
                         incu);
            cblas_dgemv( order_enum, 
                         CblasNoTrans, 
                         IM.lenc, 
                         IM.lenr,
                         1.0, 
                         IM.elem_p, 
                         strideIM, 
                         ((matlib_real*)U.elem_p+1), 
                         incu, 
                         0,
                         ((matlib_real*)u.elem_p+1),
                         incu);
            if(not_last)
            {
                pthread_mutex_unlock(&lock_common[ptr->thread_index]);
            }
        }
        else
        {
            term_execb( "size of vectors/matrices incorrect: ",
                        IM.lenc, IM.lenr, U.len, u.len );
        }
    }
    else
    {
        term_exec( "highest degree of polynomials incorrect: %d",
                    IM.lenr );
    }
    debug_exit("Thread id: %d", ptr->thread_index);
}

void pfem1d_ZILT
(
    const matlib_index    N,
    const matlib_xm       IM,
          matlib_zv       U,
          matlib_zv       u,
          matlib_index    num_threads,
          pthpool_data_t* mp
)
{

    matlib_index i; 
    /* define mutexes for adjacent elements of u w.r.t. threads */
    pthread_mutex_t lock_common[num_threads-1];
    for(i=0; i<num_threads-1; i++)
    {
       pthread_mutex_init(&(lock_common[i]), NULL);
    }
    /* define the shared data */ 
    void* shared_data[5] = { (void*) &N,
                             (void*) &IM,
                             (void*) &U,
                             (void*) &u,
                             (void*) &(lock_common)};
    void* nonshared_data[num_threads][2];
    matlib_complex common_u[num_threads-1];

    matlib_index nsdata[num_threads][2];


    pthpool_arg_t   arg[num_threads];
    pthpool_task_t  task[num_threads];

    /* define the block of data per thread */ 
    matlib_index Np = N/(num_threads);

    i = 0;
    nsdata[i][0] = i*Np;
    nsdata[i][1] = (i+1)*Np;
    arg[i].shared_data    = shared_data; 

    nonshared_data[i][0]  = (void*)nsdata[i];
    nonshared_data[i][1]  = NULL;
    arg[i].nonshared_data = (void**)&nonshared_data[i];
    arg[i].thread_index   = i;
    /* Define the task */ 
    task[i].function  = (void*)pfem1d_thfunc_ZILT;
    task[i].argument  = &arg[i];
    mp[i].task = &task[i];
    for(i=1; i<num_threads-1; i++)
    {
        nsdata[i][0] = i*Np;
        nsdata[i][1] = (i+1)*Np;
        arg[i].shared_data    = shared_data; 
        nonshared_data[i][0]  = (void*)nsdata[i];
        nonshared_data[i][1]  = (void*)&common_u[i-1];
        arg[i].nonshared_data = (void**)&nonshared_data[i];
        arg[i].thread_index   = i;
        /* Define the task */ 
        task[i].function  = (void*)pfem1d_thfunc_ZILT;
        task[i].argument  = &arg[i];
        mp[i].task = &task[i];
    }
    if(i<num_threads)
    {
        nsdata[i][0] = i*Np;
        nsdata[i][1] = N;
        arg[i].shared_data    = shared_data; 
        nonshared_data[i][0]  = (void*)nsdata[i];
        nonshared_data[i][1]  = (void*)&common_u[i-1];
        arg[i].nonshared_data = (void**)&nonshared_data[i];
        arg[i].thread_index   = i;
        /* Define the task */ 
        task[i].function  = (void*)pfem1d_thfunc_ZILT;
        task[i].argument  = &arg[i];
        mp[i].task = &task[i];
    }
    debug_body("%s", "created task");

    pthpool_exec_task(num_threads, mp, task);
    matlib_index P = IM.lenc-1; 

    for(i=0; i<num_threads-1; i++)
    {
        u.elem_p[(i+1)*Np*P] = common_u[i];
    }

    debug_body("%s", "destroying mutexes");
    for(i=0; i<num_threads-1; i++)
    {
       pthread_mutex_destroy(&(lock_common[i]));
    }
    
    debug_exit("%s", "");
}
/*============================================================================*/
static void* pfem1d_thfunc_XFLT2(void* mp)
/* 
 * N : can be deduced from the length of the matrices provided but
 *     this leads to untraceble bugs. Hence this is made as an argument.
 * */ 
{
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;
    matlib_index N = *((matlib_index*) (ptr->shared_data[0]));
    matlib_xm FM   =    *((matlib_xm*) (ptr->shared_data[1]));
    matlib_xm u    =    *((matlib_xm*) (ptr->shared_data[2]));
    matlib_xm U    =    *((matlib_xm*) (ptr->shared_data[3]));

    debug_enter( "Thread index: %d, "
                 "nr. finite-elements: %d, matrices FM: %d-by-%d, "
                 "u: %d-by-%d, U: %d-by-%d", N, ptr->thread_index,
                 FM.lenc, FM.lenr, 
                 u.lenc , u.lenr,
                 U.lenc , U.lenr );

   /* highest degree of Legendre polynomials retained 
    * */ 
    matlib_index p = FM.lenc-1; 
    
    /* nr. LGL-points -1 used for sampling in FEM-element 
     * */
    matlib_index P = FM.lenr-1;  

    bool order_OK = (u.order == MATLIB_COL_MAJOR) && 
                    (U.order == MATLIB_COL_MAJOR);

    if(!order_OK)
    {
        term_exec( "Storage order of matrices incorrect: %s", "u and U");
    }

    assert(((FM.elem_p!=NULL) && (u.elem_p !=NULL)) && (U.elem_p!=NULL));
    assert(P>1);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    if(P>1)
    {
        debug_body("nr. computed finite-elements: %d", (u.lenc-1)/P);
        assert(N == (u.lenc-1)/P);

        if ((U.lenc == N*(p+1)) && (u.lenr == U.lenr))
        {
            debug_body( "highest degree: %d, "
                         "nr. samples on ref. domain: %d, "
                         "nr. of transform vectors: %d",
                         p, FM.lenr, u.lenr);

            matlib_index strideFM;
            _CBLAS_ORDER order_enum;

            GET_CBLASORDER_AND_STRIDE(order_enum, strideFM, FM);

            matlib_index incu = 1;
            matlib_index incU = 1;
            matlib_index i, j;
            
            (u.elem_p) += (u.lenc*start_end_index[0]);
            (U.elem_p) += (U.lenc*start_end_index[0]);

            for (j=start_end_index[0]; j<start_end_index[1]; j++)
            {
                for(i=0; i<N; i++)
                {
                    cblas_dgemv( order_enum, 
                                 CblasNoTrans, 
                                 FM.lenc, 
                                 FM.lenr, 
                                 1.0, 
                                 FM.elem_p, 
                                 strideFM,
                                 u.elem_p, 
                                 incu, 
                                 0, 
                                 U.elem_p, 
                                 incU);

                    (u.elem_p) += P;
                    (U.elem_p) += FM.lenc;
                }
                (u.elem_p)++;
            }
        }
        else
        {
            term_exec( "size of matrices incorrect: "
                       "FM: %d-by-%d, u: %d-by-%d, U: %d-by-%d", N, 
                       FM.lenc, FM.lenr, 
                       u.lenc , u.lenr,
                       U.lenc , U.lenr );
        }
    }
    else
    {
        term_exec( "incorrect nr. of sampling points in ref. domain: %d",
                    FM.lenr );
    }

    debug_exit("%s","");
}

void pfem1d_XFLT2
(
    const matlib_index    N,
    const matlib_xm       FM,
          matlib_xm       u,
          matlib_xm       U,
          matlib_index    num_threads,
          pthpool_data_t* mp
)
{
    matlib_index i; 
    /* define the shared data */ 
    void* shared_data[4] = { (void*) &N,
                             (void*) &FM,
                             (void*) &u,
                             (void*) &U };

    matlib_index nsdata[num_threads][2];


    pthpool_arg_t   arg[num_threads];
    pthpool_task_t  task[num_threads];

    /* define the block of data per thread */ 
    matlib_index Np = u.lenr/(num_threads);

    for(i=0; i<num_threads-1; i++)
    {
        nsdata[i][0] = i*Np;
        nsdata[i][1] = (i+1)*Np;
        arg[i].shared_data    = shared_data; 
        arg[i].nonshared_data = (void**)&nsdata[i];
        arg[i].thread_index   = i;
        /* Define the task */ 
        task[i].function  = (void*)pfem1d_thfunc_XFLT2;
        task[i].argument  = &arg[i];
        mp[i].task = &task[i];
    }
    nsdata[i][0] = i*Np;
    nsdata[i][1] = u.lenr;
    arg[i].shared_data    = shared_data; 
    arg[i].nonshared_data = (void**)&nsdata[i];
    arg[i].thread_index   = i;
    /* Define the task */ 
    task[i].function  = (void*)pfem1d_thfunc_XFLT2;
    task[i].argument  = &arg[i];
    mp[i].task = &task[i];

    debug_body("%s", "created task");

    pthpool_exec_task(num_threads, mp, task);
    
    debug_exit("%s", "");

}
/*============================================================================*/
static void* pfem1d_thfunc_ZFLT2(void* mp)
/* 
 * N : can be deduced from the length of the matrices provided but
 *     this leads to untraceble bugs. Hence this is made as an argument.
 * */ 
{
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;
    matlib_index N = *((matlib_index*) (ptr->shared_data[0]));
    matlib_xm FM   =    *((matlib_xm*) (ptr->shared_data[1]));
    matlib_zm u    =    *((matlib_zm*) (ptr->shared_data[2]));
    matlib_zm U    =    *((matlib_zm*) (ptr->shared_data[3]));
    debug_enter( "nr. finite-elements: %d, matrices FM: %d-by-%d, "
                 "u: %d-by-%d, U: %d-by-%d", N, 
                 FM.lenc, FM.lenr, 
                 u.lenc , u.lenr,
                 U.lenc , U.lenr );

   /* highest degree of Legendre polynomials retained 
    * */ 
    matlib_index p = FM.lenc-1; 
    
    /* nr. LGL-points -1 used for sampling in FEM-element 
     * */
    matlib_index P = FM.lenr-1;  


    bool order_OK = (u.order == MATLIB_COL_MAJOR) && 
                    (U.order == MATLIB_COL_MAJOR);

    if(!order_OK)
    {
        term_exec( "Storage order of matrices incorrect: %s", "u and U");
    }

    assert(((FM.elem_p!=NULL) && (u.elem_p !=NULL)) && (U.elem_p!=NULL));
    assert(P>1);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    if(P>1)
    {
        debug_body("nr. computed finite-elements: %d", (u.lenc-1)/P);
        assert(N == (u.lenc-1)/P);

        if ((U.lenc == N*(p+1)) && (u.lenr == U.lenr))
        {
            debug_body( "highest degree: %d, "
                         "nr. samples on ref. domain: %d, "
                         "nr. of transform vectors: %d",
                         p, FM.lenr, u.lenr);

            matlib_index strideFM;
            _CBLAS_ORDER order_enum;

            GET_CBLASORDER_AND_STRIDE(order_enum, strideFM, FM);

            matlib_index incu = 2;
            matlib_index incU = 2;
            matlib_index i, j;
            
            (u.elem_p) += (u.lenc*start_end_index[0]);
            (U.elem_p) += (U.lenc*start_end_index[0]);

            for (j=start_end_index[0]; j<start_end_index[1]; j++)
            {
                for(i=0; i<N; i++)
                {
                    cblas_dgemv( order_enum, 
                                 CblasNoTrans, 
                                 FM.lenc, 
                                 FM.lenr, 
                                 1.0, 
                                 FM.elem_p, 
                                 strideFM,
                                 (matlib_real*)u.elem_p, 
                                 incu, 
                                 0, 
                                 (matlib_real*)U.elem_p, 
                                 incU);
                    cblas_dgemv( order_enum, 
                                 CblasNoTrans, 
                                 FM.lenc, 
                                 FM.lenr, 
                                 1.0, 
                                 FM.elem_p, 
                                 strideFM,
                                 (matlib_real*)u.elem_p+1, 
                                 incu, 
                                 0, 
                                 (matlib_real*)U.elem_p+1, 
                                 incU);

                    (u.elem_p) += P;
                    (U.elem_p) += FM.lenc;
                }
                (u.elem_p)++;
            }
        }
        else
        {
            term_exec( "size of matrices incorrect: "
                       "FM: %d-by-%d, u: %d-by-%d, U: %d-by-%d", N, 
                       FM.lenc, FM.lenr, 
                       u.lenc , u.lenr,
                       U.lenc , U.lenr );
        }
    }
    else
    {
        term_exec( "incorrect nr. of sampling points in ref. domain: %d",
                    FM.lenr );
    }
    debug_exit("%s","");
}

void pfem1d_ZFLT2
(
    const matlib_index    N,
    const matlib_xm       FM,
          matlib_zm       u,
          matlib_zm       U,
          matlib_index    num_threads,
          pthpool_data_t* mp
)
{
    matlib_index i; 
    /* define the shared data */ 
    void* shared_data[4] = { (void*) &N,
                             (void*) &FM,
                             (void*) &u,
                             (void*) &U };

    matlib_index nsdata[num_threads][2];


    pthpool_arg_t   arg[num_threads];
    pthpool_task_t  task[num_threads];

    /* define the block of data per thread */ 
    matlib_index Np = u.lenr/(num_threads);

    for(i=0; i<num_threads-1; i++)
    {
        nsdata[i][0] = i*Np;
        nsdata[i][1] = (i+1)*Np;
        arg[i].shared_data    = shared_data; 
        arg[i].nonshared_data = (void**)&nsdata[i];
        arg[i].thread_index   = i;
        /* Define the task */ 
        task[i].function  = (void*)pfem1d_thfunc_ZFLT2;
        task[i].argument  = &arg[i];
        mp[i].task = &task[i];
    }
    nsdata[i][0] = i*Np;
    nsdata[i][1] = u.lenr;
    arg[i].shared_data    = shared_data; 
    arg[i].nonshared_data = (void**)&nsdata[i];
    arg[i].thread_index   = i;
    /* Define the task */ 
    task[i].function  = (void*)pfem1d_thfunc_ZFLT2;
    task[i].argument  = &arg[i];
    mp[i].task = &task[i];

    debug_body("%s", "created task");

    pthpool_exec_task(num_threads, mp, task);
    
    debug_exit("%s", "");

}

/*============================================================================+/
 | Transformation from Legendre basis to FEM-basis and vice-versa
/+============================================================================*/

static void* thfunc_dshapefunc2lp(void* mp)
/* (elem_n+1)-by-1 vector, vertex function basis */
/* (p-1)-by-1 vector, bubble function basis      */
/* (p+1)*elem_n-by-1 vector in FEM-basis         */

{
    /* p must be greater than 3 */ 
    matlib_index i, j;

    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_index p = *(matlib_index*) (ptr->shared_data[0]);

    matlib_real *v = (matlib_real*) (ptr->shared_data[1]);
    matlib_real *b = (matlib_real*) (ptr->shared_data[2]);
    matlib_real *u = (matlib_real*) (ptr->shared_data[3]);
    matlib_real *B = (matlib_real*) (ptr->shared_data[5]);
    matlib_real *C = (matlib_real*) (ptr->shared_data[4]);
    matlib_real *A = (matlib_real*) (ptr->shared_data[6]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    v += start_end_index[0];

    b += ((p-1) * start_end_index[0]);
    u += ((p+1) * start_end_index[0]);

    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v)+A[0]**b    ; u++;
        *u = 0.5*(-*(v-1) + *v)+A[1]**(b+1); u++;

        //pB = (matlib_real*) (ptr->shared_data[5]);
        //pC = (matlib_real*) (ptr->shared_data[6]);
        *u = B[0]**b + C[0]**(b+2); b++, u++;
        *u = B[1]**b + C[1]**(b+2); b++, u++;
        *u = B[2]**b + C[2]**(b+2); b++, u++;
        *u = B[3]**b + C[3]**(b+2); b++, u++;
        *u = B[4]**b + C[4]**(b+2); b++, u++;
        *u = B[5]**b + C[5]**(b+2); b++, u++;
        *u = B[6]**b + C[6]**(b+2); b++, u++;

        for(j=7; j<(p-3); j++, b++, u++)
        {
            *u = B[j]**b + C[j]**(b+2);
        }
        *u = A[2]**b, u++, b++; 
        *u = A[3]**b, u++, b++; 
    }
    debug_exit("%s", "");
}
 
void pfem1d_XF2L
(
    const matlib_index    p, 
    const matlib_xv       vb,
          matlib_xv       u,
          matlib_index    num_threads,
          pthpool_data_t* mp
)
{
    debug_enter( "highest polynomial degree: %d "
                 "length of vectors vb: %d, u: %d", 
                 p, vb.len, u.len );

    /* vb.len = N*p+1, u.len = N*(p+1)
     * vb.len-1+N = u.len
     * */ 
    matlib_index N = (vb.len-1)/p;
    debug_body( "nr finite elements: %d ", N);

    assert((vb.elem_p != NULL) && (u.elem_p != NULL));


    void* (*fp[9])(void*) = { thfunc_dshapefunc2lp_2,
                              thfunc_dshapefunc2lp_3, 
                              thfunc_dshapefunc2lp_4, 
                              thfunc_dshapefunc2lp_5, 
                              thfunc_dshapefunc2lp_6, 
                              thfunc_dshapefunc2lp_7, 
                              thfunc_dshapefunc2lp_8, 
                              thfunc_dshapefunc2lp_9, 
                              thfunc_dshapefunc2lp_10};

    if(u.len == vb.len+(N-1))
    {
        matlib_index i,j;
        matlib_index nsdata[num_threads][2];

        pthpool_arg_t   arg[num_threads];
        pthpool_task_t  task[num_threads];

        /* define the block of data per thread */ 
        matlib_index Np = N/(num_threads);

        if(p>10)
        {

            /* define the shared data */ 
            matlib_real A[4];

            matlib_xv B, C;
            matlib_create_xv(p-3, &B, MATLIB_COL_VECT);
            matlib_create_xv(p-3, &C, MATLIB_COL_VECT);
            
            A[0] = -1.0/sqrt(6);
            A[1] = -1.0/sqrt(10);
            A[2] =  1.0/sqrt(2*(2*p-3));
            A[3] =  1.0/sqrt(2*(2*p-1));

            for( j=0; j<(p-3); j++)
            {
                B.elem_p[j] =  1.0/sqrt(2*(2*j+3)); 
                C.elem_p[j] = -1.0/sqrt(2*(2*j+7)); 
            }

            void* shared_data[7] = { (void*) &p,
                                     (void*) (vb.elem_p),
                                     (void*) (vb.elem_p+N+1),
                                     (void*) (u.elem_p),
                                     (void*) C.elem_p,
                                     (void*) B.elem_p,
                                     (void*) A
                                    };
            

            for(i=0; i<num_threads-1; i++)
            {
                nsdata[i][0] = i*Np;
                nsdata[i][1] = (i+1)*Np;
                arg[i].shared_data    = shared_data; 
                arg[i].nonshared_data = (void**)&nsdata[i];
                arg[i].thread_index   = i;
                /* Define the task */ 
                task[i].function  = (void*)thfunc_dshapefunc2lp;
                task[i].argument  = &arg[i];
                mp[i].task = &task[i];
            }
            nsdata[i][0] = i*Np;
            nsdata[i][1] = N;
            arg[i].shared_data    = shared_data; 
            arg[i].nonshared_data = (void**)&nsdata[i];
            arg[i].thread_index   = i;
            /* Define the task */ 
            task[i].function  = (void*)thfunc_dshapefunc2lp;
            task[i].argument  = &arg[i];
            mp[i].task = &task[i];

            debug_body("%s", "created task");

            pthpool_exec_task(num_threads, mp, task);

            matlib_free(B.elem_p);
            matlib_free(C.elem_p);
        }
        else
        {
            matlib_index func_index = p-2;
            void* shared_data[3] = { (void*) (vb.elem_p),
                                     (void*) (vb.elem_p+N+1),
                                     (void*) (u.elem_p)
                                    };
            
            for(i=0; i<num_threads-1; i++)
            {
                nsdata[i][0] = i*Np;
                nsdata[i][1] = (i+1)*Np;
                arg[i].shared_data    = shared_data; 
                arg[i].nonshared_data = (void*)&nsdata[i];
                arg[i].thread_index   = i;
                /* Define the task */ 
                task[i].function  = (void*)fp[func_index];
                task[i].argument  = &arg[i];
                mp[i].task = &task[i];
            }
            nsdata[i][0] = i*Np;
            nsdata[i][1] = N;
            arg[i].shared_data    = shared_data; 
            arg[i].nonshared_data = (void*)&nsdata[i];
            arg[i].thread_index   = i;
            /* Define the task */ 
            task[i].function  = (void*)fp[func_index];
            task[i].argument  = &arg[i];
            mp[i].task = &task[i];

            debug_body("%s", "created task");
            pthpool_exec_task(num_threads, mp, task);

        }        
    }
    else
    {
        term_execb( "length of vectors incorrect: "
                    "vb: %d, u: %d",
                    vb.len, u.len);
    }

    debug_exit("%s", "");
}
/*============================================================================*/

static void* thfunc_zshapefunc2lp(void* mp)
/* (elem_n+1)-by-1 vector, vertex function basis */
/* (p-1)-by-1 vector, bubble function basis      */
/* (p+1)*elem_n-by-1 vector in FEM-basis         */

{
    /* p must be greater than 3 */ 
    matlib_index i;

    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_index p = *(matlib_index*) (ptr->shared_data[0]);

    matlib_complex *v = (matlib_complex*) (ptr->shared_data[1]);
    matlib_complex *b = (matlib_complex*) (ptr->shared_data[2]);
    matlib_complex *u = (matlib_complex*) (ptr->shared_data[3]);
    matlib_real *A = (matlib_real*) (ptr->shared_data[4]);
    matlib_real *pB, *pC;

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    v += start_end_index[0];

    b += ((p-1) * start_end_index[0]);
    u += ((p+1) * start_end_index[0]);

    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v)+A[0]**b    ; u++;
        *u = 0.5*(-*(v-1) + *v)+A[1]**(b+1); u++;

        pB = (matlib_real*) (ptr->shared_data[5]);
        pC = (matlib_real*) (ptr->shared_data[6]);
        *u = *pB**b + *pC**(b+2); pB++, pC++, b++, u++;
        *u = *pB**b + *pC**(b+2); pB++, pC++, b++, u++;
        *u = *pB**b + *pC**(b+2); pB++, pC++, b++, u++;
        *u = *pB**b + *pC**(b+2); pB++, pC++, b++, u++;
        *u = *pB**b + *pC**(b+2); pB++, pC++, b++, u++;
        *u = *pB**b + *pC**(b+2); pB++, pC++, b++, u++;
        *u = *pB**b + *pC**(b+2); pB++, pC++, b++, u++;

        for 
        (
            ; (pB < (matlib_real*) (ptr->shared_data[5]) + (p-3));
            pB++, pC++, b++, u++
        )
        {
            *u = *pB**b + *pC**(b+2);
        }
        *u = A[2]**b, u++, b++; 
        *u = A[3]**b, u++, b++; 
    }
    debug_exit("%s", "");
}
 
void pfem1d_ZF2L
(
    const matlib_index    p, 
    const matlib_zv       vb,
          matlib_zv       u,
          matlib_index    num_threads,
          pthpool_data_t* mp
)
{
    debug_enter( "highest polynomial degree: %d "
                 "length of vectors vb: %d, u: %d", 
                 p, vb.len, u.len );

    /* vb.len = N*p+1, u.len = N*(p+1)
     * vb.len-1+N = u.len
     * */ 
    matlib_index N = (vb.len-1)/p;
    debug_body( "nr finite elements: %d ", N);
    void* (*fp[9])(void*) = { thfunc_zshapefunc2lp_2,
                              thfunc_zshapefunc2lp_3, 
                              thfunc_zshapefunc2lp_4, 
                              thfunc_zshapefunc2lp_5, 
                              thfunc_zshapefunc2lp_6, 
                              thfunc_zshapefunc2lp_7, 
                              thfunc_zshapefunc2lp_8, 
                              thfunc_zshapefunc2lp_9, 
                              thfunc_zshapefunc2lp_10};


    assert((vb.elem_p != NULL) && (u.elem_p != NULL));

    if(u.len == vb.len+(N-1))
    {
        matlib_index i,j;
        matlib_index nsdata[num_threads][2];


        pthpool_arg_t   arg[num_threads];
        pthpool_task_t  task[num_threads];

        /* define the block of data per thread */ 
        matlib_index Np = N/(num_threads);
        if(p>10)
        {

            /* define the shared data */ 
            matlib_real A[4];

            matlib_xv B, C;
            matlib_create_xv(p-3, &B, MATLIB_COL_VECT);
            matlib_create_xv(p-3, &C, MATLIB_COL_VECT);
            
            A[0] = -1.0/sqrt(6);
            A[1] = -1.0/sqrt(10);
            A[2] =  1.0/sqrt(2*(2*p-3));
            A[3] =  1.0/sqrt(2*(2*p-1));

            for( j=0; j<(p-3); j++)
            {
                B.elem_p[j] =  1.0/sqrt(2*(2*j+3.0)); 
                C.elem_p[j] = -1.0/sqrt(2*(2*j+7.0)); 
            }
            void* shared_data[7] = { (void*) &p,
                                     (void*) (vb.elem_p),
                                     (void*) (vb.elem_p+N+1),
                                     (void*) (u.elem_p),
                                     (void*) A,
                                     (void*) B.elem_p,
                                     (void*) C.elem_p
                                    };
            

            for(i=0; i<num_threads-1; i++)
            {
                nsdata[i][0] = i*Np;
                nsdata[i][1] = (i+1)*Np;
                arg[i].shared_data    = shared_data; 
                arg[i].nonshared_data = (void*)&nsdata[i];
                arg[i].thread_index   = i;
                /* Define the task */ 
                task[i].function  = (void*)thfunc_zshapefunc2lp;
                task[i].argument  = &arg[i];
                mp[i].task = &task[i];
            }
            nsdata[i][0] = i*Np;
            nsdata[i][1] = N;
            arg[i].shared_data    = shared_data; 
            arg[i].nonshared_data = (void*)&nsdata[i];
            arg[i].thread_index   = i;
            /* Define the task */ 
            task[i].function  = (void*)thfunc_zshapefunc2lp;
            task[i].argument  = &arg[i];
            mp[i].task = &task[i];

            debug_body("%s", "created task");
            pthpool_exec_task(num_threads, mp, task);

            matlib_free(B.elem_p);
            matlib_free(C.elem_p);

        }
        else
        {
            matlib_index func_index = p-2;

            void* shared_data[3] = { (void*) (vb.elem_p),
                                     (void*) (vb.elem_p+N+1),
                                     (void*) (u.elem_p)
                                    };
            

            for(i=0; i<num_threads-1; i++)
            {
                nsdata[i][0] = i*Np;
                nsdata[i][1] = (i+1)*Np;
                arg[i].shared_data    = shared_data; 
                arg[i].nonshared_data = (void*)&nsdata[i];
                arg[i].thread_index   = i;
                /* Define the task */ 
                task[i].function  = (void*)fp[func_index];
                task[i].argument  = &arg[i];
                mp[i].task = &task[i];
            }
            nsdata[i][0] = i*Np;
            nsdata[i][1] = N;
            arg[i].shared_data    = shared_data; 
            arg[i].nonshared_data = (void*)&nsdata[i];
            arg[i].thread_index   = i;
            /* Define the task */ 
            task[i].function  = (void*)fp[func_index];
            task[i].argument  = &arg[i];
            mp[i].task = &task[i];

            debug_body("%s", "created task");
            pthpool_exec_task(num_threads, mp, task);
        
        }
        

    }
    else
    {
        term_execb( "length of vectors incorrect: "
                    "vb: %d, u: %d",
                    vb.len, u.len);
    }

    debug_exit("%s", "");
}
/*============================================================================+/
 | Unrolled version 
/+============================================================================*/

 

static void* thfunc_dshapefunc2lp_2(void* mp)
{

    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *v = (matlib_real*) (ptr->shared_data[0]);
    matlib_real *b = (matlib_real*) (ptr->shared_data[1]);
    matlib_real *u = (matlib_real*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    v += start_end_index[0];

    b += (start_end_index[0]);
    u += (3 * start_end_index[0]);

    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b; u++;
        *u = 0.5*(-*(v-1) + *v)          ; u++;

        *u = _A03**b, u++, b++;
    }
}
 

 

static void* thfunc_dshapefunc2lp_3(void* mp)
{

    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *v = (matlib_real*) (ptr->shared_data[0]);
    matlib_real *b = (matlib_real*) (ptr->shared_data[1]);
    matlib_real *u = (matlib_real*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    v += start_end_index[0];

    b += (2 * start_end_index[0]);
    u += (4 * start_end_index[0]);

    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;

        *u = _A03**b, u++, b++;
        *u = _A04**b, u++, b++;
    }
}
 

 

static void* thfunc_dshapefunc2lp_4(void* mp)
{

    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *v = (matlib_real*) (ptr->shared_data[0]);
    matlib_real *b = (matlib_real*) (ptr->shared_data[1]);
    matlib_real *u = (matlib_real*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    v += start_end_index[0];

    b += (3 * start_end_index[0]);
    u += (5 * start_end_index[0]);

    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;

        *u = _A04**b, u++, b++;
        *u = _A05**b, u++, b++;
    }
}
 

 

static void* thfunc_dshapefunc2lp_5(void* mp)
{

    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *v = (matlib_real*) (ptr->shared_data[0]);
    matlib_real *b = (matlib_real*) (ptr->shared_data[1]);
    matlib_real *u = (matlib_real*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    v += start_end_index[0];

    b += (4 * start_end_index[0]);
    u += (6 * start_end_index[0]);

    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;

        *u = _A05**b, u++, b++;
        *u = _A06**b, u++, b++;
    }
}
 

 

static void* thfunc_dshapefunc2lp_6(void* mp)
{

    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *v = (matlib_real*) (ptr->shared_data[0]);
    matlib_real *b = (matlib_real*) (ptr->shared_data[1]);
    matlib_real *u = (matlib_real*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    v += start_end_index[0];

    b += (5 * start_end_index[0]);
    u += (7 * start_end_index[0]);

    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;
        *u = _A05**b + _B05**(b+2), u++, b++;

        *u = _A06**b, u++, b++;
        *u = _A07**b, u++, b++;
    }
}
 

 

static void* thfunc_dshapefunc2lp_7(void* mp)
{

    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *v = (matlib_real*) (ptr->shared_data[0]);
    matlib_real *b = (matlib_real*) (ptr->shared_data[1]);
    matlib_real *u = (matlib_real*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    v += start_end_index[0];

    b += ( 6 * start_end_index[0]);
    u += ( 8 * start_end_index[0]);

    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;
        *u = _A05**b + _B05**(b+2), u++, b++;
        *u = _A06**b + _B06**(b+2), u++, b++;

        *u = _A07**b, u++, b++;
        *u = _A08**b, u++, b++;
    }
}
 

 

static void* thfunc_dshapefunc2lp_8(void* mp)
{

    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *v = (matlib_real*) (ptr->shared_data[0]);
    matlib_real *b = (matlib_real*) (ptr->shared_data[1]);
    matlib_real *u = (matlib_real*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    v += start_end_index[0];

    b += ( 7 * start_end_index[0]);
    u += ( 9 * start_end_index[0]);

    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;
        *u = _A05**b + _B05**(b+2), u++, b++;
        *u = _A06**b + _B06**(b+2), u++, b++;
        *u = _A07**b + _B07**(b+2), u++, b++;

        *u = _A08**b, u++, b++;
        *u = _A09**b, u++, b++;
    }
}
 

 

static void* thfunc_dshapefunc2lp_9(void* mp)
{

    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *v = (matlib_real*) (ptr->shared_data[0]);
    matlib_real *b = (matlib_real*) (ptr->shared_data[1]);
    matlib_real *u = (matlib_real*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    v += start_end_index[0];

    b += ( 8 * start_end_index[0]);
    u += (10 * start_end_index[0]);

    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;
        *u = _A05**b + _B05**(b+2), u++, b++;
        *u = _A06**b + _B06**(b+2), u++, b++;
        *u = _A07**b + _B07**(b+2), u++, b++;
        *u = _A08**b + _B08**(b+2), u++, b++;

        *u = _A09**b, u++, b++;
        *u = _A10**b, u++, b++;
    }
}
 

 

static void* thfunc_dshapefunc2lp_10(void* mp)
{

    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *v = (matlib_real*) (ptr->shared_data[0]);
    matlib_real *b = (matlib_real*) (ptr->shared_data[1]);
    matlib_real *u = (matlib_real*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    v += start_end_index[0];

    b += ( 9 * start_end_index[0]);
    u += (11 * start_end_index[0]);

    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;
        *u = _A05**b + _B05**(b+2), u++, b++;
        *u = _A06**b + _B06**(b+2), u++, b++;
        *u = _A07**b + _B07**(b+2), u++, b++;
        *u = _A08**b + _B08**(b+2), u++, b++;
        *u = _A09**b + _B09**(b+2), u++, b++;

        *u = _A10**b, u++, b++;
        *u = _A11**b, u++, b++;
    }
}
 

/* COMPLEX VERSION */ 

 

static void* thfunc_zshapefunc2lp_2(void* mp)
{

    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *v = (matlib_complex*) (ptr->shared_data[0]);
    matlib_complex *b = (matlib_complex*) (ptr->shared_data[1]);
    matlib_complex *u = (matlib_complex*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    v += start_end_index[0];

    b += (start_end_index[0]);
    u += (3 * start_end_index[0]);

    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b; u++;
        *u = 0.5*(-*(v-1) + *v)          ; u++;

        *u = _A03**b, u++, b++;
    }
}
 

 

static void* thfunc_zshapefunc2lp_3(void* mp)
{

    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *v = (matlib_complex*) (ptr->shared_data[0]);
    matlib_complex *b = (matlib_complex*) (ptr->shared_data[1]);
    matlib_complex *u = (matlib_complex*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    v += start_end_index[0];

    b += (2 * start_end_index[0]);
    u += (4 * start_end_index[0]);

    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;

        *u = _A03**b, u++, b++;
        *u = _A04**b, u++, b++;
    }
}
 

 

static void* thfunc_zshapefunc2lp_4(void* mp)
{

    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *v = (matlib_complex*) (ptr->shared_data[0]);
    matlib_complex *b = (matlib_complex*) (ptr->shared_data[1]);
    matlib_complex *u = (matlib_complex*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    v += start_end_index[0];

    b += (3 * start_end_index[0]);
    u += (5 * start_end_index[0]);

    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;

        *u = _A04**b, u++, b++;
        *u = _A05**b, u++, b++;
    }
}
 

 

static void* thfunc_zshapefunc2lp_5(void* mp)
{

    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *v = (matlib_complex*) (ptr->shared_data[0]);
    matlib_complex *b = (matlib_complex*) (ptr->shared_data[1]);
    matlib_complex *u = (matlib_complex*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    v += start_end_index[0];

    b += (4 * start_end_index[0]);
    u += (6 * start_end_index[0]);

    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;

        *u = _A05**b, u++, b++;
        *u = _A06**b, u++, b++;
    }
}
 

 

static void* thfunc_zshapefunc2lp_6(void* mp)
{

    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *v = (matlib_complex*) (ptr->shared_data[0]);
    matlib_complex *b = (matlib_complex*) (ptr->shared_data[1]);
    matlib_complex *u = (matlib_complex*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    v += start_end_index[0];

    b += (5 * start_end_index[0]);
    u += (7 * start_end_index[0]);

    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;
        *u = _A05**b + _B05**(b+2), u++, b++;

        *u = _A06**b, u++, b++;
        *u = _A07**b, u++, b++;
    }
}
 

 

static void* thfunc_zshapefunc2lp_7(void* mp)
{

    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *v = (matlib_complex*) (ptr->shared_data[0]);
    matlib_complex *b = (matlib_complex*) (ptr->shared_data[1]);
    matlib_complex *u = (matlib_complex*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    v += start_end_index[0];

    b += ( 6 * start_end_index[0]);
    u += ( 8 * start_end_index[0]);

    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;
        *u = _A05**b + _B05**(b+2), u++, b++;
        *u = _A06**b + _B06**(b+2), u++, b++;

        *u = _A07**b, u++, b++;
        *u = _A08**b, u++, b++;
    }
}
 

 

static void* thfunc_zshapefunc2lp_8(void* mp)
{

    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *v = (matlib_complex*) (ptr->shared_data[0]);
    matlib_complex *b = (matlib_complex*) (ptr->shared_data[1]);
    matlib_complex *u = (matlib_complex*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    v += start_end_index[0];

    b += ( 7 * start_end_index[0]);
    u += ( 9 * start_end_index[0]);

    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;
        *u = _A05**b + _B05**(b+2), u++, b++;
        *u = _A06**b + _B06**(b+2), u++, b++;
        *u = _A07**b + _B07**(b+2), u++, b++;

        *u = _A08**b, u++, b++;
        *u = _A09**b, u++, b++;
    }
}
 

 

static void* thfunc_zshapefunc2lp_9(void* mp)
{

    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *v = (matlib_complex*) (ptr->shared_data[0]);
    matlib_complex *b = (matlib_complex*) (ptr->shared_data[1]);
    matlib_complex *u = (matlib_complex*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    v += start_end_index[0];

    b += ( 8 * start_end_index[0]);
    u += (10 * start_end_index[0]);

    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;
        *u = _A05**b + _B05**(b+2), u++, b++;
        *u = _A06**b + _B06**(b+2), u++, b++;
        *u = _A07**b + _B07**(b+2), u++, b++;
        *u = _A08**b + _B08**(b+2), u++, b++;

        *u = _A09**b, u++, b++;
        *u = _A10**b, u++, b++;
    }
}
 

 

static void* thfunc_zshapefunc2lp_10(void* mp)
{

    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *v = (matlib_complex*) (ptr->shared_data[0]);
    matlib_complex *b = (matlib_complex*) (ptr->shared_data[1]);
    matlib_complex *u = (matlib_complex*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    v += start_end_index[0];

    b += ( 9 * start_end_index[0]);
    u += (11 * start_end_index[0]);

    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        v++;
        *u = 0.5*( *(v-1) + *v) + _A01**b    ; u++;
        *u = 0.5*(-*(v-1) + *v) + _A02**(b+1); u++;
        *u = _A03**b + _B03**(b+2), u++, b++;
        *u = _A04**b + _B04**(b+2), u++, b++;
        *u = _A05**b + _B05**(b+2), u++, b++;
        *u = _A06**b + _B06**(b+2), u++, b++;
        *u = _A07**b + _B07**(b+2), u++, b++;
        *u = _A08**b + _B08**(b+2), u++, b++;
        *u = _A09**b + _B09**(b+2), u++, b++;

        *u = _A10**b, u++, b++;
        *u = _A11**b, u++, b++;
    }
}
 

/*============================================================================+/
 | Projection from LP basis representation to FEM-basis
/+============================================================================*/

static void* thfunc_dprjLP2FEM_ShapeFunc(void* mp)
/* 
 * Pv : vector of size (elem_n+1)     
 * Pb : vector of size elem_n*(p-1)   
 *
 * */ 
{
    matlib_index i, j;

    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_index p = *(matlib_index*) (ptr->shared_data[0]);

    matlib_real *u  = (matlib_real*) (ptr->shared_data[1]);
    matlib_real *Pv = (matlib_real*) (ptr->shared_data[2]);
    matlib_real *Pb = (matlib_real*) (ptr->shared_data[3]);
    matlib_real *B  = (matlib_real*) (ptr->shared_data[4]);
    matlib_real *C  = (matlib_real*) (ptr->shared_data[5]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real tmp;
    u += ((p+1) * start_end_index[0]);
    Pv += start_end_index[0];
    Pb += ((p-1) * start_end_index[0]);
    if(start_end_index[0]==0)
    {
        tmp = 0;
    }
    else
    {
        tmp = (*(u-p-1) + *(u-p)/3);
    }

    for (i=start_end_index[0]; i<start_end_index[1]; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        
        *Pb = B[0]**(u+2) + C[0]**u; Pb++, u++;
        *Pb = B[1]**(u+2) + C[1]**u; Pb++, u++;
        *Pb = B[2]**(u+2) + C[2]**u; Pb++, u++;
        *Pb = B[3]**(u+2) + C[3]**u; Pb++, u++;
        *Pb = B[4]**(u+2) + C[4]**u; Pb++, u++;
        *Pb = B[5]**(u+2) + C[5]**u; Pb++, u++;
        *Pb = B[6]**(u+2) + C[6]**u; Pb++, u++;
        *Pb = B[7]**(u+2) + C[7]**u; Pb++, u++;
        *Pb = B[8]**(u+2) + C[8]**u; Pb++, u++;

        for( j=9; j < (p-1); Pb++, u++, j++)
        {
            *Pb = B[j]**(u+2) + C[j]**u;
        }
        u += 2;
    }
    debug_exit("%s", "");

}

void pfem1d_XPrjL2F
(
    matlib_index    p,
    matlib_xv       u,
    matlib_xv       Pvb,
    matlib_index    num_threads,
    pthpool_data_t* mp
)
{
    debug_enter( "highest polynomial degree: %d "
                 "length of vectors u: %d, Pbv: %d", 
                 p, u.len, Pvb.len );

    matlib_index N = u.len/(p+1);
    debug_body( "nr finite elements: %d ", N);

    assert((u.elem_p != NULL) && (Pvb.elem_p != NULL));

    if(u.len == Pvb.len+N-1)
    {
        void* (*fp[9])(void*) = { thfunc_dprjLP2FEM_ShapeFunc_2,
                                  thfunc_dprjLP2FEM_ShapeFunc_3, 
                                  thfunc_dprjLP2FEM_ShapeFunc_4, 
                                  thfunc_dprjLP2FEM_ShapeFunc_5, 
                                  thfunc_dprjLP2FEM_ShapeFunc_6, 
                                  thfunc_dprjLP2FEM_ShapeFunc_7, 
                                  thfunc_dprjLP2FEM_ShapeFunc_8, 
                                  thfunc_dprjLP2FEM_ShapeFunc_9, 
                                  thfunc_dprjLP2FEM_ShapeFunc_10};



        matlib_index i;
        matlib_index nsdata[num_threads][2];

        pthpool_arg_t   arg[num_threads];
        pthpool_task_t  task[num_threads];

        /* define the block of data per thread */ 
        matlib_index Np = N/(num_threads);
        matlib_real tmp = 0;

        if(p>10)
        {
            matlib_xv B, C;
            matlib_create_xv(p-1, &B, MATLIB_COL_VECT);
            matlib_create_xv(p-1, &C, MATLIB_COL_VECT);
            
            for(i=0; i<p-1; i++)
            {
                tmp =  1.0/sqrt(4*i+6);
                B.elem_p[i] =  tmp/(i+2.5);
                C.elem_p[i] = -tmp/(i+0.5);
            }

            void* shared_data[6] = { (void*) &p,
                                     (void*) (u.elem_p),
                                     (void*) (Pvb.elem_p),
                                     (void*) (Pvb.elem_p+N+1),
                                     (void*) B.elem_p,
                                     (void*) C.elem_p
                                    };
            

            for(i=0; i<num_threads-1; i++)
            {
                nsdata[i][0] = i*Np;
                nsdata[i][1] = (i+1)*Np;
                arg[i].shared_data    = shared_data; 
                arg[i].nonshared_data = (void**)&nsdata[i];
                arg[i].thread_index   = i;
                /* Define the task */ 
                task[i].function  = (void*)thfunc_dprjLP2FEM_ShapeFunc;
                task[i].argument  = &arg[i];
                mp[i].task = &task[i];
            }
            nsdata[i][0] = i*Np;
            nsdata[i][1] = N;
            arg[i].shared_data    = shared_data; 
            arg[i].nonshared_data = (void**)&nsdata[i];
            arg[i].thread_index   = i;
            /* Define the task */ 
            task[i].function  = (void*)thfunc_dprjLP2FEM_ShapeFunc;
            task[i].argument  = &arg[i];
            mp[i].task = &task[i];

            debug_body("%s", "created task");
            pthpool_exec_task(num_threads, mp, task);
            
            i = (p+1)*(N-1);
            *(Pvb.elem_p+N) = *(u.elem_p+i) + *(u.elem_p+i+1)/3;

            matlib_free(B.elem_p);
            matlib_free(C.elem_p);
        }
        else
        {
            matlib_index func_index = p-2;
            void* shared_data[3] = { (void*) (u.elem_p),
                                     (void*) (Pvb.elem_p),
                                     (void*) (Pvb.elem_p+N+1)
                                    };

            for(i=0; i<num_threads-1; i++)
            {
                nsdata[i][0] = i*Np;
                nsdata[i][1] = (i+1)*Np;
                arg[i].shared_data    = shared_data; 
                arg[i].nonshared_data = (void**)&nsdata[i];
                arg[i].thread_index   = i;
                /* Define the task */ 
                task[i].function  = (void*)fp[func_index];
                task[i].argument  = &arg[i];
                mp[i].task = &task[i];
            }
            nsdata[i][0] = i*Np;
            nsdata[i][1] = N;
            arg[i].shared_data    = shared_data; 
            arg[i].nonshared_data = (void**)&nsdata[i];
            arg[i].thread_index   = i;
            /* Define the task */ 
            task[i].function  = (void*)fp[func_index];
            task[i].argument  = &arg[i];
            mp[i].task = &task[i];

            debug_body("%s", "created task");
            pthpool_exec_task(num_threads, mp, task);
            
            i = (p+1)*(N-1);
            *(Pvb.elem_p+N) = *(u.elem_p+i) + *(u.elem_p+i+1)/3;
        }
    }
    else
    {
        term_execb( "length of vectors incorrect: "
                    "u: %d, Pvb: %d",
                    u.len, Pvb.len);
    }

    debug_exit("%s", "");

}

static void* thfunc_zprjLP2FEM_ShapeFunc(void* mp)
/* 
 * Pv : vector of size (elem_n+1)     
 * Pb : vector of size elem_n*(p-1)   
 *
 * */ 
{
    matlib_index i, j;

    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_index p = *(matlib_index*) (ptr->shared_data[0]);

    matlib_complex *u  = (matlib_complex*) (ptr->shared_data[1]);
    matlib_complex *Pv = (matlib_complex*) (ptr->shared_data[2]);
    matlib_complex *Pb = (matlib_complex*) (ptr->shared_data[3]);
    matlib_real *B  = (matlib_real*) (ptr->shared_data[4]);
    matlib_real *C  = (matlib_real*) (ptr->shared_data[5]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_complex tmp;
    u += ((p+1) * start_end_index[0]);
    Pv += start_end_index[0];
    Pb += ((p-1) * start_end_index[0]);
    if(start_end_index[0]==0)
    {
        tmp = 0;
    }
    else
    {
        tmp = (*(u-p-1) + *(u-p)/3);
    }

    for (i=start_end_index[0]; i<start_end_index[1]; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        
        *Pb = B[0]**(u+2) + C[0]**u; Pb++, u++;
        *Pb = B[1]**(u+2) + C[1]**u; Pb++, u++;
        *Pb = B[2]**(u+2) + C[2]**u; Pb++, u++;
        *Pb = B[3]**(u+2) + C[3]**u; Pb++, u++;
        *Pb = B[4]**(u+2) + C[4]**u; Pb++, u++;
        *Pb = B[5]**(u+2) + C[5]**u; Pb++, u++;
        *Pb = B[6]**(u+2) + C[6]**u; Pb++, u++;
        *Pb = B[7]**(u+2) + C[7]**u; Pb++, u++;
        *Pb = B[8]**(u+2) + C[8]**u; Pb++, u++;

        for( j=9; j < (p-1); Pb++, u++, j++)
        {
            *Pb = B[j]**(u+2) + C[j]**u;
        }
        u += 2;
    }
    debug_exit("%s", "");

}

void pfem1d_ZPrjL2F
(
    matlib_index    p,
    matlib_zv       u,
    matlib_zv       Pvb,
    matlib_index    num_threads,
    pthpool_data_t* mp
)
{
    debug_enter( "highest polynomial degree: %d "
                 "length of vectors u: %d, Pbv: %d", 
                 p, u.len, Pvb.len );

    matlib_index N = u.len/(p+1);
    debug_body( "nr finite elements: %d ", N);

    assert((u.elem_p != NULL) && (Pvb.elem_p != NULL));

    if(u.len == Pvb.len+N-1)
    {
        void* (*fp[9])(void*) = { thfunc_zprjLP2FEM_ShapeFunc_2,
                                  thfunc_zprjLP2FEM_ShapeFunc_3, 
                                  thfunc_zprjLP2FEM_ShapeFunc_4, 
                                  thfunc_zprjLP2FEM_ShapeFunc_5, 
                                  thfunc_zprjLP2FEM_ShapeFunc_6, 
                                  thfunc_zprjLP2FEM_ShapeFunc_7, 
                                  thfunc_zprjLP2FEM_ShapeFunc_8, 
                                  thfunc_zprjLP2FEM_ShapeFunc_9, 
                                  thfunc_zprjLP2FEM_ShapeFunc_10};



        matlib_index i;
        matlib_index nsdata[num_threads][2];

        pthpool_arg_t   arg[num_threads];
        pthpool_task_t  task[num_threads];

        /* define the block of data per thread */ 
        matlib_index Np = N/(num_threads);
        matlib_real tmp = 0;

        if(p>10)
        {
            matlib_xv B, C;
            matlib_create_xv(p-1, &B, MATLIB_COL_VECT);
            matlib_create_xv(p-1, &C, MATLIB_COL_VECT);
            
            for(i=0; i<p-1; i++)
            {
                tmp =  1.0/sqrt(4*i+6);
                B.elem_p[i] =  tmp/(i+2.5);
                C.elem_p[i] = -tmp/(i+0.5);
            }

            void* shared_data[6] = { (void*) &p,
                                     (void*) (u.elem_p),
                                     (void*) (Pvb.elem_p),
                                     (void*) (Pvb.elem_p+N+1),
                                     (void*) B.elem_p,
                                     (void*) C.elem_p
                                    };
            

            for(i=0; i<num_threads-1; i++)
            {
                nsdata[i][0] = i*Np;
                nsdata[i][1] = (i+1)*Np;
                arg[i].shared_data    = shared_data; 
                arg[i].nonshared_data = (void**)&nsdata[i];
                arg[i].thread_index   = i;
                /* Define the task */ 
                task[i].function  = (void*)thfunc_zprjLP2FEM_ShapeFunc;
                task[i].argument  = &arg[i];
                mp[i].task = &task[i];
            }
            nsdata[i][0] = i*Np;
            nsdata[i][1] = N;
            arg[i].shared_data    = shared_data; 
            arg[i].nonshared_data = (void**)&nsdata[i];
            arg[i].thread_index   = i;
            /* Define the task */ 
            task[i].function  = (void*)thfunc_zprjLP2FEM_ShapeFunc;
            task[i].argument  = &arg[i];
            mp[i].task = &task[i];

            debug_body("%s", "created task");
            pthpool_exec_task(num_threads, mp, task);
            
            i = (p+1)*(N-1);
            *(Pvb.elem_p+N) = *(u.elem_p+i) + *(u.elem_p+i+1)/3;

            matlib_free(B.elem_p);
            matlib_free(C.elem_p);
        }
        else
        {
            matlib_index func_index = p-2;
            void* shared_data[3] = { (void*) (u.elem_p),
                                     (void*) (Pvb.elem_p),
                                     (void*) (Pvb.elem_p+N+1)
                                    };

            for(i=0; i<num_threads-1; i++)
            {
                nsdata[i][0] = i*Np;
                nsdata[i][1] = (i+1)*Np;
                arg[i].shared_data    = shared_data; 
                arg[i].nonshared_data = (void**)&nsdata[i];
                arg[i].thread_index   = i;
                /* Define the task */ 
                task[i].function  = (void*)fp[func_index];
                task[i].argument  = &arg[i];
                mp[i].task = &task[i];
            }
            nsdata[i][0] = i*Np;
            nsdata[i][1] = N;
            arg[i].shared_data    = shared_data; 
            arg[i].nonshared_data = (void**)&nsdata[i];
            arg[i].thread_index   = i;
            /* Define the task */ 
            task[i].function  = (void*)fp[func_index];
            task[i].argument  = &arg[i];
            mp[i].task = &task[i];

            debug_body("%s", "created task");
            pthpool_exec_task(num_threads, mp, task);
            
            i = (p+1)*(N-1);
            *(Pvb.elem_p+N) = *(u.elem_p+i) + *(u.elem_p+i+1)/3;
        }
    }
    else
    {
        term_execb( "length of vectors incorrect: "
                    "u: %d, Pvb: %d",
                    u.len, Pvb.len);
    }

    debug_exit("%s", "");

}
/*============================================================================*/

 

static void* thfunc_dprjLP2FEM_ShapeFunc_2(void* mp)
{
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *u  = (matlib_real*) (ptr->shared_data[0]);
    matlib_real *Pv = (matlib_real*) (ptr->shared_data[1]);
    matlib_real *Pb = (matlib_real*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real tmp;
    u += (3 * start_end_index[0]);
    Pv += start_end_index[0];
    Pb += (start_end_index[0]);
    if(start_end_index[0]==0)
    {
        tmp = 0;
    }
    else
    {
        tmp = (*(u-3) + *(u-2)/3);
    }
    for (i=start_end_index[0]; i<start_end_index[1]; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        u += 2;

    }
}
 

 

static void* thfunc_dprjLP2FEM_ShapeFunc_3(void* mp)
{
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *u  = (matlib_real*) (ptr->shared_data[0]);
    matlib_real *Pv = (matlib_real*) (ptr->shared_data[1]);
    matlib_real *Pb = (matlib_real*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real tmp;
    u += (4 * start_end_index[0]);
    Pv += start_end_index[0];
    Pb += (2 * start_end_index[0]);
    if(start_end_index[0]==0)
    {
        tmp = 0;
    }
    else
    {
        tmp = (*(u-4) + *(u-3)/3);
    }
    for (i=start_end_index[0]; i<start_end_index[1]; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        u += 2;

    }
}
 

 

static void* thfunc_dprjLP2FEM_ShapeFunc_4(void* mp)
{
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *u  = (matlib_real*) (ptr->shared_data[0]);
    matlib_real *Pv = (matlib_real*) (ptr->shared_data[1]);
    matlib_real *Pb = (matlib_real*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real tmp;
    u += (5 * start_end_index[0]);
    Pv += start_end_index[0];
    Pb += (3 * start_end_index[0]);
    if(start_end_index[0]==0)
    {
        tmp = 0;
    }
    else
    {
        tmp = (*(u-5) + *(u-4)/3);
    }
    for (i=start_end_index[0]; i<start_end_index[1]; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        u += 2;

    }
}
 

 

static void* thfunc_dprjLP2FEM_ShapeFunc_5(void* mp)
{
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *u  = (matlib_real*) (ptr->shared_data[0]);
    matlib_real *Pv = (matlib_real*) (ptr->shared_data[1]);
    matlib_real *Pb = (matlib_real*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real tmp;
    u  += (6 * start_end_index[0]);
    Pv += start_end_index[0];
    Pb += (4 * start_end_index[0]);
    if(start_end_index[0]==0)
    {
        tmp = 0;
    }
    else
    {
        tmp = (*(u-6) + *(u-5)/3);
    }

    for (i=start_end_index[0]; i<start_end_index[1]; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        u += 2;
    }
}
 

 

static void* thfunc_dprjLP2FEM_ShapeFunc_6(void* mp)
{
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *u  = (matlib_real*) (ptr->shared_data[0]);
    matlib_real *Pv = (matlib_real*) (ptr->shared_data[1]);
    matlib_real *Pb = (matlib_real*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real tmp;
    u += (7 * start_end_index[0]);
    Pv += start_end_index[0];
    Pb += (5 * start_end_index[0]);
    if(start_end_index[0]==0)
    {
        tmp = 0;
    }
    else
    {
        tmp = (*(u-7) + *(u-6)/3);
    }

    for (i=start_end_index[0]; i<start_end_index[1]; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        *Pb = _E04**(u+2) + _F04**u; Pb++; u++;
        u += 2;
    }
}
 

 

static void* thfunc_dprjLP2FEM_ShapeFunc_7(void* mp)
{
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *u  = (matlib_real*) (ptr->shared_data[0]);
    matlib_real *Pv = (matlib_real*) (ptr->shared_data[1]);
    matlib_real *Pb = (matlib_real*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real tmp;
    u += (8 * start_end_index[0]);
    Pv += start_end_index[0];
    Pb += (6 * start_end_index[0]);
    if(start_end_index[0]==0)
    {
        tmp = 0;
    }
    else
    {
        tmp = (*(u-8) + *(u-7)/3);
    }
    for (i=start_end_index[0]; i<start_end_index[1]; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        *Pb = _E04**(u+2) + _F04**u; Pb++; u++;
        *Pb = _E05**(u+2) + _F05**u; Pb++; u++;
        u += 2;
    }
}
 

 

static void* thfunc_dprjLP2FEM_ShapeFunc_8(void* mp)
{
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *u  = (matlib_real*) (ptr->shared_data[0]);
    matlib_real *Pv = (matlib_real*) (ptr->shared_data[1]);
    matlib_real *Pb = (matlib_real*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real tmp;
    u += (9 * start_end_index[0]);
    Pv += start_end_index[0];
    Pb += (7 * start_end_index[0]);
    if(start_end_index[0]==0)
    {
        tmp = 0;
    }
    else
    {
        tmp = (*(u-9) + *(u-8)/3);
    }
    for (i=start_end_index[0]; i<start_end_index[1]; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        *Pb = _E04**(u+2) + _F04**u; Pb++; u++;
        *Pb = _E05**(u+2) + _F05**u; Pb++; u++;
        *Pb = _E06**(u+2) + _F06**u; Pb++; u++;
        u += 2;
    }
}
 

 

static void* thfunc_dprjLP2FEM_ShapeFunc_9(void* mp)
{
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *u  = (matlib_real*) (ptr->shared_data[0]);
    matlib_real *Pv = (matlib_real*) (ptr->shared_data[1]);
    matlib_real *Pb = (matlib_real*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real tmp;
    u += (10 * start_end_index[0]);
    Pv += start_end_index[0];
    Pb += (8 * start_end_index[0]);
    if(start_end_index[0]==0)
    {
        tmp = 0;
    }
    else
    {
        tmp = (*(u-10) + *(u-9)/3);
    }
    for (i=start_end_index[0]; i<start_end_index[1]; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        *Pb = _E04**(u+2) + _F04**u; Pb++; u++;
        *Pb = _E05**(u+2) + _F05**u; Pb++; u++;
        *Pb = _E06**(u+2) + _F06**u; Pb++; u++;
        *Pb = _E07**(u+2) + _F07**u; Pb++; u++;
        u += 2;
    }
}
 

 

static void* thfunc_dprjLP2FEM_ShapeFunc_10(void* mp)
{
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *u  = (matlib_real*) (ptr->shared_data[0]);
    matlib_real *Pv = (matlib_real*) (ptr->shared_data[1]);
    matlib_real *Pb = (matlib_real*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real tmp;
    
    u += (11 * start_end_index[0]);
    Pv += start_end_index[0];
    Pb += (9 * start_end_index[0]);

    if(start_end_index[0]==0)
    {
        tmp = 0;
    }
    else
    {
        tmp = (*(u-11) + *(u-10)/3);
    }

    for (i=start_end_index[0]; i<start_end_index[1]; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        *Pb = _E04**(u+2) + _F04**u; Pb++; u++;
        *Pb = _E05**(u+2) + _F05**u; Pb++; u++;
        *Pb = _E06**(u+2) + _F06**u; Pb++; u++;
        *Pb = _E07**(u+2) + _F07**u; Pb++; u++;
        *Pb = _E08**(u+2) + _F08**u; Pb++; u++;
        u += 2;
    }
}
 

/* COMPLEX VERSION */ 

 

static void* thfunc_zprjLP2FEM_ShapeFunc_2(void* mp)
{
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *u  = (matlib_complex*) (ptr->shared_data[0]);
    matlib_complex *Pv = (matlib_complex*) (ptr->shared_data[1]);
    matlib_complex *Pb = (matlib_complex*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_complex tmp;

    u += (3 * start_end_index[0]);
    Pv += start_end_index[0];
    Pb += (start_end_index[0]);
    if(start_end_index[0]==0)
    {
        tmp = 0;
    }
    else
    {
        tmp = (*(u-3) + *(u-2)/3);
    }
    for (i=start_end_index[0]; i<start_end_index[1]; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        u += 2;

    }
}
 

 

static void* thfunc_zprjLP2FEM_ShapeFunc_3(void* mp)
{
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *u  = (matlib_complex*) (ptr->shared_data[0]);
    matlib_complex *Pv = (matlib_complex*) (ptr->shared_data[1]);
    matlib_complex *Pb = (matlib_complex*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_complex tmp;

    u += (4 * start_end_index[0]);
    Pv += start_end_index[0];
    Pb += (2 * start_end_index[0]);
    if(start_end_index[0]==0)
    {
        tmp = 0;
    }
    else
    {
        tmp = (*(u-4) + *(u-3)/3);
    }
    for (i=start_end_index[0]; i<start_end_index[1]; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        u += 2;

    }
}
 

 

static void* thfunc_zprjLP2FEM_ShapeFunc_4(void* mp)
{
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *u  = (matlib_complex*) (ptr->shared_data[0]);
    matlib_complex *Pv = (matlib_complex*) (ptr->shared_data[1]);
    matlib_complex *Pb = (matlib_complex*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_complex tmp;

    u += (5 * start_end_index[0]);
    Pv += start_end_index[0];
    Pb += (3 * start_end_index[0]);
    if(start_end_index[0]==0)
    {
        tmp = 0;
    }
    else
    {
        tmp = (*(u-5) + *(u-4)/3);
    }
    for (i=start_end_index[0]; i<start_end_index[1]; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        u += 2;

    }
}
 

 

static void* thfunc_zprjLP2FEM_ShapeFunc_5(void* mp)
{
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *u  = (matlib_complex*) (ptr->shared_data[0]);
    matlib_complex *Pv = (matlib_complex*) (ptr->shared_data[1]);
    matlib_complex *Pb = (matlib_complex*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_complex tmp;

    u  += (6 * start_end_index[0]);
    Pv += start_end_index[0];
    Pb += (4 * start_end_index[0]);
    if(start_end_index[0]==0)
    {
        tmp = 0;
    }
    else
    {
        tmp = (*(u-6) + *(u-5)/3);
    }

    for (i=start_end_index[0]; i<start_end_index[1]; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        u += 2;
    }
}
 

 

static void* thfunc_zprjLP2FEM_ShapeFunc_6(void* mp)
{
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *u  = (matlib_complex*) (ptr->shared_data[0]);
    matlib_complex *Pv = (matlib_complex*) (ptr->shared_data[1]);
    matlib_complex *Pb = (matlib_complex*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_complex tmp;

    u += (7 * start_end_index[0]);
    Pv += start_end_index[0];
    Pb += (5 * start_end_index[0]);
    if(start_end_index[0]==0)
    {
        tmp = 0;
    }
    else
    {
        tmp = (*(u-7) + *(u-6)/3);
    }

    for (i=start_end_index[0]; i<start_end_index[1]; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        *Pb = _E04**(u+2) + _F04**u; Pb++; u++;
        u += 2;
    }
}
 

 

static void* thfunc_zprjLP2FEM_ShapeFunc_7(void* mp)
{
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *u  = (matlib_complex*) (ptr->shared_data[0]);
    matlib_complex *Pv = (matlib_complex*) (ptr->shared_data[1]);
    matlib_complex *Pb = (matlib_complex*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_complex tmp;

    u += (8 * start_end_index[0]);
    Pv += start_end_index[0];
    Pb += (6 * start_end_index[0]);
    if(start_end_index[0]==0)
    {
        tmp = 0;
    }
    else
    {
        tmp = (*(u-8) + *(u-7)/3);
    }
    for (i=start_end_index[0]; i<start_end_index[1]; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        *Pb = _E04**(u+2) + _F04**u; Pb++; u++;
        *Pb = _E05**(u+2) + _F05**u; Pb++; u++;
        u += 2;
    }
}
 

 

static void* thfunc_zprjLP2FEM_ShapeFunc_8(void* mp)
{
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *u  = (matlib_complex*) (ptr->shared_data[0]);
    matlib_complex *Pv = (matlib_complex*) (ptr->shared_data[1]);
    matlib_complex *Pb = (matlib_complex*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_complex tmp;

    u += (9 * start_end_index[0]);
    Pv += start_end_index[0];
    Pb += (7 * start_end_index[0]);
    if(start_end_index[0]==0)
    {
        tmp = 0;
    }
    else
    {
        tmp = (*(u-9) + *(u-8)/3);
    }
    for (i=start_end_index[0]; i<start_end_index[1]; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        *Pb = _E04**(u+2) + _F04**u; Pb++; u++;
        *Pb = _E05**(u+2) + _F05**u; Pb++; u++;
        *Pb = _E06**(u+2) + _F06**u; Pb++; u++;
        u += 2;
    }
}
 

 

static void* thfunc_zprjLP2FEM_ShapeFunc_9(void* mp)
{
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *u  = (matlib_complex*) (ptr->shared_data[0]);
    matlib_complex *Pv = (matlib_complex*) (ptr->shared_data[1]);
    matlib_complex *Pb = (matlib_complex*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_complex tmp;

    u += (10 * start_end_index[0]);
    Pv += start_end_index[0];
    Pb += (8 * start_end_index[0]);
    if(start_end_index[0]==0)
    {
        tmp = 0;
    }
    else
    {
        tmp = (*(u-10) + *(u-9)/3);
    }
    for (i=start_end_index[0]; i<start_end_index[1]; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        *Pb = _E04**(u+2) + _F04**u; Pb++; u++;
        *Pb = _E05**(u+2) + _F05**u; Pb++; u++;
        *Pb = _E06**(u+2) + _F06**u; Pb++; u++;
        *Pb = _E07**(u+2) + _F07**u; Pb++; u++;
        u += 2;
    }
}
 

 

static void* thfunc_zprjLP2FEM_ShapeFunc_10(void* mp)
{
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *u  = (matlib_complex*) (ptr->shared_data[0]);
    matlib_complex *Pv = (matlib_complex*) (ptr->shared_data[1]);
    matlib_complex *Pb = (matlib_complex*) (ptr->shared_data[2]);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_complex tmp;
    
    u += (11 * start_end_index[0]);
    Pv += start_end_index[0];
    Pb += (9 * start_end_index[0]);

    if(start_end_index[0]==0)
    {
        tmp = 0;
    }
    else
    {
        tmp = (*(u-11) + *(u-10)/3);
    }

    for (i=start_end_index[0]; i<start_end_index[1]; i++, Pv++)
    {
        *Pv = (*u - *(u+1)/3) + tmp;
        tmp = (*u + *(u+1)/3);
        *Pb = _E00**(u+2) + _F00**u; Pb++; u++;
        *Pb = _E01**(u+2) + _F01**u; Pb++; u++;
        *Pb = _E02**(u+2) + _F02**u; Pb++; u++;
        *Pb = _E03**(u+2) + _F03**u; Pb++; u++;
        *Pb = _E04**(u+2) + _F04**u; Pb++; u++;
        *Pb = _E05**(u+2) + _F05**u; Pb++; u++;
        *Pb = _E06**(u+2) + _F06**u; Pb++; u++;
        *Pb = _E07**(u+2) + _F07**u; Pb++; u++;
        *Pb = _E08**(u+2) + _F08**u; Pb++; u++;
        u += 2;
    }
}
 
/*============================================================================+/
 | Norm using Parseval's theorem
/+============================================================================*/
static void* thfunc_dlp_snorm2_d(void* mp)
{
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *u   = (matlib_real*) (ptr->shared_data[0]);
    matlib_xv A = *(matlib_xv*) (ptr->shared_data[1]);

    matlib_index i, j;

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data[0]);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real* snorm   = (matlib_real*) (ptr->nonshared_data[1]);
    u += (A.len * start_end_index[0]);

    *snorm = 0;
    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    { 
        (*snorm) += (*u**u*A.elem_p[0]); u++;
        (*snorm) += (*u**u*A.elem_p[1]); u++;
        (*snorm) += (*u**u*A.elem_p[2]); u++;
        (*snorm) += (*u**u*A.elem_p[3]); u++;
        (*snorm) += (*u**u*A.elem_p[4]); u++;
        (*snorm) += (*u**u*A.elem_p[5]); u++;
        (*snorm) += (*u**u*A.elem_p[6]); u++;
        (*snorm) += (*u**u*A.elem_p[7]); u++;
        (*snorm) += (*u**u*A.elem_p[8]); u++;
        (*snorm) += (*u**u*A.elem_p[9]); u++;
        (*snorm) += (*u**u*A.elem_p[10]); u++;
        for( j=11; j < A.len; u++, j++)
        {
            (*snorm) += (*u**u*A.elem_p[j]);
        }
    }
    debug_exit("Thread id: %d, snorm2: %0.16f", ptr->thread_index, *snorm);
}

matlib_real pfem1d_XNorm2
(
    matlib_index    p,
    matlib_index    N,
    matlib_xv       u,
    matlib_index    num_threads,
    pthpool_data_t* mp
)
{
    debug_enter( "Highest degree of polynomial: %d "
                 "nr. finite-elements : %d "
                 "length of u: %d", p, N, u.len);

    void* (*fp[9])(void*) = { thfunc_dlp_snorm2_d_2, 
                              thfunc_dlp_snorm2_d_3, 
                              thfunc_dlp_snorm2_d_4, 
                              thfunc_dlp_snorm2_d_5, 
                              thfunc_dlp_snorm2_d_6, 
                              thfunc_dlp_snorm2_d_7, 
                              thfunc_dlp_snorm2_d_8, 
                              thfunc_dlp_snorm2_d_9, 
                              thfunc_dlp_snorm2_d_10};
 
    matlib_real norm2 = 0;
    assert(u.elem_p!=NULL);
    if(u.len == (p+1)*N)
    {
        matlib_index i,j;
        matlib_index nsdata[num_threads][2];
        matlib_real snorm2[num_threads];
        void* nonshared_data[num_threads][2];

        pthpool_arg_t   arg[num_threads];
        pthpool_task_t  task[num_threads];

        /* define the block of data per thread */ 
        matlib_index Np = N/(num_threads);
        if(p>10)
        {
            matlib_xv A;
            matlib_create_xv(p+1, &A, MATLIB_COL_VECT);

            for( j=0; j<A.len; j++)
            {
                A.elem_p[j] =  1.0/(j+0.5); 
            }

            void* shared_data[2] = { (void*) u.elem_p,
                                     (void*) &A };
            

            for(i=0; i<num_threads-1; i++)
            {
                nsdata[i][0] = i*Np;
                nsdata[i][1] = (i+1)*Np;
                arg[i].shared_data = shared_data; 

                nonshared_data[i][0]  = (void*)nsdata[i];
                nonshared_data[i][1]  = (void*)&snorm2[i];
                arg[i].nonshared_data = (void**)&nonshared_data[i];
                arg[i].thread_index   = i;
                /* Define the task */ 
                task[i].function  = (void*)thfunc_dlp_snorm2_d;
                task[i].argument  = &arg[i];
                mp[i].task = &task[i];
            }
            nsdata[i][0] = i*Np;
            nsdata[i][1] = N;
            arg[i].shared_data    = shared_data; 
            nonshared_data[i][0]  = (void*)nsdata[i];
            nonshared_data[i][1]  = (void*)&snorm2[i];
            arg[i].nonshared_data = (void**)&nonshared_data[i];
            arg[i].thread_index   = i;
            /* Define the task */ 
            task[i].function = (void*)thfunc_dlp_snorm2_d;
            task[i].argument = &arg[i];
            mp[i].task = &task[i];

            debug_body("%s", "created task");

            pthpool_exec_task(num_threads, mp, task);

            matlib_free((void*)A.elem_p);
        }
        else
        {
            matlib_index func_index = p-2;
            void* shared_data = (void*) u.elem_p;
            

            for(i=0; i<num_threads-1; i++)
            {
                nsdata[i][0] = i*Np;
                nsdata[i][1] = (i+1)*Np;
                arg[i].shared_data = shared_data; 

                nonshared_data[i][0]  = (void*)nsdata[i];
                nonshared_data[i][1]  = (void*)&snorm2[i];
                arg[i].nonshared_data = (void**)&nonshared_data[i];
                arg[i].thread_index   = i;
                /* Define the task */ 
                task[i].function  = (void*)fp[func_index];
                task[i].argument  = &arg[i];
                mp[i].task = &task[i];
            }
            nsdata[i][0] = i*Np;
            nsdata[i][1] = N;
            arg[i].shared_data    = shared_data; 
            nonshared_data[i][0]  = (void*)nsdata[i];
            nonshared_data[i][1]  = (void*)&snorm2[i];
            arg[i].nonshared_data = (void**)&nonshared_data[i];
            arg[i].thread_index   = i;
            /* Define the task */ 
            task[i].function = (void*)fp[func_index];
            task[i].argument = &arg[i];
            mp[i].task = &task[i];

            debug_body("%s", "created task");

            pthpool_exec_task(num_threads, mp, task);

        }
        for(i=0; i<num_threads; i++)
        {
            norm2 += (snorm2[i]);
        }
        norm2 = sqrt(norm2);
        
    }
    else
    {
        term_execb( "length of vector incorrect: "
                    "u: %d (degree : %d, nr. fem-elements: %d)",
                    u.len, p, N);
    }

    debug_exit("L2 norm: %0.16f", (norm2));
    return(norm2);

}
/*============================================================================*/

static void* thfunc_zlp_snorm2_d(void* mp)
{
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *u   = (matlib_complex*) (ptr->shared_data[0]);
    matlib_xv A = *(matlib_xv*) (ptr->shared_data[1]);

    matlib_index i, j;

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data[0]);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real* snorm   = (matlib_real*) (ptr->nonshared_data[1]);
    u += (A.len * start_end_index[0]);

    *snorm = 0;
    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    { 
        (*snorm) += (*u*conj(*u)*A.elem_p[0]); u++;
        (*snorm) += (*u*conj(*u)*A.elem_p[1]); u++;
        (*snorm) += (*u*conj(*u)*A.elem_p[2]); u++;
        (*snorm) += (*u*conj(*u)*A.elem_p[3]); u++;
        (*snorm) += (*u*conj(*u)*A.elem_p[4]); u++;
        (*snorm) += (*u*conj(*u)*A.elem_p[5]); u++;
        (*snorm) += (*u*conj(*u)*A.elem_p[6]); u++;
        (*snorm) += (*u*conj(*u)*A.elem_p[7]); u++;
        (*snorm) += (*u*conj(*u)*A.elem_p[8]); u++;
        (*snorm) += (*u*conj(*u)*A.elem_p[9]); u++;
        (*snorm) += (*u*conj(*u)*A.elem_p[10]); u++;
        for( j=11; j < A.len; u++, j++)
        {
            (*snorm) += (*u*conj(*u)*A.elem_p[j]);
        }
    }
    debug_exit("Thread id: %d, snorm2: %0.16f", ptr->thread_index, *snorm);
}

matlib_real pfem1d_ZNorm2
(
    matlib_index    p,
    matlib_index    N,
    matlib_zv       u,
    matlib_index    num_threads,
    pthpool_data_t* mp
)
{
    debug_enter( "Highest degree of polynomial: %d "
                 "nr. finite-elements : %d "
                 "length of u: %d", p, N, u.len);

    void* (*fp[9])(void*) = { thfunc_zlp_snorm2_d_2, 
                              thfunc_zlp_snorm2_d_3, 
                              thfunc_zlp_snorm2_d_4, 
                              thfunc_zlp_snorm2_d_5, 
                              thfunc_zlp_snorm2_d_6, 
                              thfunc_zlp_snorm2_d_7, 
                              thfunc_zlp_snorm2_d_8, 
                              thfunc_zlp_snorm2_d_9, 
                              thfunc_zlp_snorm2_d_10};
 
    matlib_real norm2 = 0;
    assert(u.elem_p!=NULL);
    if(u.len == (p+1)*N)
    {
        matlib_index i,j;
        matlib_index nsdata[num_threads][2];
        matlib_real snorm2[num_threads];
        void* nonshared_data[num_threads][2];

        pthpool_arg_t   arg[num_threads];
        pthpool_task_t  task[num_threads];

        /* define the block of data per thread */ 
        matlib_index Np = N/(num_threads);
        if(p>10)
        {
            matlib_xv A;
            matlib_create_xv(p+1, &A, MATLIB_COL_VECT);

            for( j=0; j<A.len; j++)
            {
                A.elem_p[j] =  1.0/(j+0.5); 
            }

            void* shared_data[2] = { (void*) u.elem_p,
                                     (void*) &A };
            

            for(i=0; i<num_threads-1; i++)
            {
                nsdata[i][0] = i*Np;
                nsdata[i][1] = (i+1)*Np;
                arg[i].shared_data = shared_data; 

                nonshared_data[i][0]  = (void*)nsdata[i];
                nonshared_data[i][1]  = (void*)&snorm2[i];
                arg[i].nonshared_data = (void**)&nonshared_data[i];
                arg[i].thread_index   = i;
                /* Define the task */ 
                task[i].function  = (void*)thfunc_zlp_snorm2_d;
                task[i].argument  = &arg[i];
                mp[i].task = &task[i];
            }
            nsdata[i][0] = i*Np;
            nsdata[i][1] = N;
            arg[i].shared_data    = shared_data; 
            nonshared_data[i][0]  = (void*)nsdata[i];
            nonshared_data[i][1]  = (void*)&snorm2[i];
            arg[i].nonshared_data = (void**)&nonshared_data[i];
            arg[i].thread_index   = i;
            /* Define the task */ 
            task[i].function = (void*)thfunc_zlp_snorm2_d;
            task[i].argument = &arg[i];
            mp[i].task = &task[i];

            debug_body("%s", "created task");

            pthpool_exec_task(num_threads, mp, task);

            matlib_free((void*)A.elem_p);

        }
        else
        {
            matlib_index func_index = p-2;
            void* shared_data = (void*) u.elem_p;
            

            for(i=0; i<num_threads-1; i++)
            {
                nsdata[i][0] = i*Np;
                nsdata[i][1] = (i+1)*Np;
                arg[i].shared_data = shared_data; 

                nonshared_data[i][0]  = (void*)nsdata[i];
                nonshared_data[i][1]  = (void*)&snorm2[i];
                arg[i].nonshared_data = (void**)&nonshared_data[i];
                arg[i].thread_index   = i;
                /* Define the task */ 
                task[i].function  = (void*)fp[func_index];
                task[i].argument  = &arg[i];
                mp[i].task = &task[i];
            }
            nsdata[i][0] = i*Np;
            nsdata[i][1] = N;
            arg[i].shared_data    = shared_data; 
            nonshared_data[i][0]  = (void*)nsdata[i];
            nonshared_data[i][1]  = (void*)&snorm2[i];
            arg[i].nonshared_data = (void**)&nonshared_data[i];
            arg[i].thread_index   = i;
            /* Define the task */ 
            task[i].function = (void*)fp[func_index];
            task[i].argument = &arg[i];
            mp[i].task = &task[i];

            debug_body("%s", "created task");

            pthpool_exec_task(num_threads, mp, task);
        
        }
        for(i=0; i<num_threads; i++)
        {
            norm2 += (snorm2[i]);
        }
        norm2 = sqrt(norm2);
    }
    else
    {
        term_execb( "length of vector incorrect: "
                    "u: %d (degree : %d, nr. fem-elements: %d)",
                    u.len, p, N);
    }

    debug_exit("L2 norm: %0.16f", (norm2));
    return(norm2);

}
/*============================================================================*/

 
static void* thfunc_dlp_snorm2_d_2(void* mp)
{
    matlib_real tmp;
    matlib_index i;

    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *u   = (matlib_real*) (ptr->shared_data);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data[0]);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real* snorm   = (matlib_real*) (ptr->nonshared_data[1]);
    u += (3 * start_end_index[0]);

    *snorm = 0;
    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        tmp = 	  _N0 * *(u+0) * *(u+0)
        	+ _N1 * *(u+1) * *(u+1)
        	+ _N2 * *(u+2) * *(u+2);
        u += 3;
        *snorm += tmp;
    }
}
 
 
static void* thfunc_dlp_snorm2_d_3(void* mp)
{
    matlib_real tmp;
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *u   = (matlib_real*) (ptr->shared_data);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data[0]);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real* snorm   = (matlib_real*) (ptr->nonshared_data[1]);
    u += (4 * start_end_index[0]);

    *snorm = 0;
    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        tmp = 	  _N0 * *(u+0) * *(u+0)
        	+ _N1 * *(u+1) * *(u+1)
        	+ _N2 * *(u+2) * *(u+2)
        	+ _N3 * *(u+3) * *(u+3);
        u += 4;
        *snorm += tmp;
    }
}
 
 
static void* thfunc_dlp_snorm2_d_4(void* mp)
{
    matlib_real tmp;
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *u   = (matlib_real*) (ptr->shared_data);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data[0]);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real* snorm   = (matlib_real*) (ptr->nonshared_data[1]);
    u += (5 * start_end_index[0]);

    *snorm = 0;
    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        tmp = 	  _N0 * *(u+0) * *(u+0)
        	+ _N1 * *(u+1) * *(u+1)
        	+ _N2 * *(u+2) * *(u+2)
        	+ _N3 * *(u+3) * *(u+3)
        	+ _N4 * *(u+4) * *(u+4);
        u += 5;
        *snorm += tmp;
    }
}
 
 
static void* thfunc_dlp_snorm2_d_5(void* mp)
{
    matlib_real tmp;
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *u   = (matlib_real*) (ptr->shared_data);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data[0]);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real* snorm   = (matlib_real*) (ptr->nonshared_data[1]);
    u += (6 * start_end_index[0]);

    *snorm = 0;
    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        tmp = 	  _N0 * *(u+0) * *(u+0)
        	+ _N1 * *(u+1) * *(u+1)
        	+ _N2 * *(u+2) * *(u+2)
        	+ _N3 * *(u+3) * *(u+3)
        	+ _N4 * *(u+4) * *(u+4)
        	+ _N5 * *(u+5) * *(u+5);
        u += 6;
        *snorm += tmp;
    }
}
 
 
static void* thfunc_dlp_snorm2_d_6(void* mp)
{
    matlib_real tmp;
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *u   = (matlib_real*) (ptr->shared_data);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data[0]);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real* snorm   = (matlib_real*) (ptr->nonshared_data[1]);
    u += (7 * start_end_index[0]);

    *snorm = 0;
    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        tmp = 	  _N0 * *(u+0) * *(u+0)
        	+ _N1 * *(u+1) * *(u+1)
        	+ _N2 * *(u+2) * *(u+2)
        	+ _N3 * *(u+3) * *(u+3)
        	+ _N4 * *(u+4) * *(u+4)
        	+ _N5 * *(u+5) * *(u+5)
        	+ _N6 * *(u+6) * *(u+6);
        u += 7;
        *snorm += tmp;
    }
}
 
 
static void* thfunc_dlp_snorm2_d_7(void* mp)
{
    matlib_real tmp;
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *u   = (matlib_real*) (ptr->shared_data);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data[0]);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real* snorm   = (matlib_real*) (ptr->nonshared_data[1]);
    u += (8 * start_end_index[0]);

    *snorm = 0;
    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        tmp = 	  _N0 * *(u+0) * *(u+0)
        	+ _N1 * *(u+1) * *(u+1)
        	+ _N2 * *(u+2) * *(u+2)
        	+ _N3 * *(u+3) * *(u+3)
        	+ _N4 * *(u+4) * *(u+4)
        	+ _N5 * *(u+5) * *(u+5)
        	+ _N6 * *(u+6) * *(u+6)
        	+ _N7 * *(u+7) * *(u+7);
        u += 8;
        *snorm += tmp;
    }
}
 
 
static void* thfunc_dlp_snorm2_d_8(void* mp)
{
    matlib_real tmp;
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *u   = (matlib_real*) (ptr->shared_data);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data[0]);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real* snorm   = (matlib_real*) (ptr->nonshared_data[1]);
    u += (9 * start_end_index[0]);

    *snorm = 0;
    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        tmp = 	  _N0 * *(u+0) * *(u+0)
        	+ _N1 * *(u+1) * *(u+1)
        	+ _N2 * *(u+2) * *(u+2)
        	+ _N3 * *(u+3) * *(u+3)
        	+ _N4 * *(u+4) * *(u+4)
        	+ _N5 * *(u+5) * *(u+5)
        	+ _N6 * *(u+6) * *(u+6)
        	+ _N7 * *(u+7) * *(u+7)
        	+ _N8 * *(u+8) * *(u+8);
        u += 9;
        *snorm += tmp;
    }
}
 
 
static void* thfunc_dlp_snorm2_d_9(void* mp)
{
    matlib_real tmp;
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *u   = (matlib_real*) (ptr->shared_data);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data[0]);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real* snorm   = (matlib_real*) (ptr->nonshared_data[1]);
    u += (10 * start_end_index[0]);

    *snorm = 0;
    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        tmp = 	  _N0 * *(u+0) * *(u+0)
        	+ _N1 * *(u+1) * *(u+1)
        	+ _N2 * *(u+2) * *(u+2)
        	+ _N3 * *(u+3) * *(u+3)
        	+ _N4 * *(u+4) * *(u+4)
        	+ _N5 * *(u+5) * *(u+5)
        	+ _N6 * *(u+6) * *(u+6)
        	+ _N7 * *(u+7) * *(u+7)
        	+ _N8 * *(u+8) * *(u+8)
        	+ _N9 * *(u+9) * *(u+9);
        u += 10;
        *snorm += tmp;
    }
}
 
 
static void* thfunc_dlp_snorm2_d_10(void* mp)
{
    matlib_real tmp;
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_real *u   = (matlib_real*) (ptr->shared_data);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data[0]);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real* snorm   = (matlib_real*) (ptr->nonshared_data[1]);
    u += (11 * start_end_index[0]);

    *snorm = 0;
    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        tmp = 	  _N0 * *(u+0) * *(u+0)
        	+ _N1 * *(u+1) * *(u+1)
        	+ _N2 * *(u+2) * *(u+2)
        	+ _N3 * *(u+3) * *(u+3)
        	+ _N4 * *(u+4) * *(u+4)
        	+ _N5 * *(u+5) * *(u+5)
        	+ _N6 * *(u+6) * *(u+6)
        	+ _N7 * *(u+7) * *(u+7)
        	+ _N8 * *(u+8) * *(u+8)
        	+ _N9 * *(u+9) * *(u+9)
        	+ _N10 * *(u+10) * *(u+10);
        u += 11;
        *snorm += tmp;
    }
}
 

/* COMPLEX VERSION */ 

 
static void* thfunc_zlp_snorm2_d_2(void* mp)
{
    matlib_real tmp;
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *u = (matlib_complex*) (ptr->shared_data);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data[0]);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real* snorm   = (matlib_real*) (ptr->nonshared_data[1]);
    u += (3 * start_end_index[0]);

    *snorm = 0;
    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        tmp = 	  _N0 * *(u+0) * conj(*(u+0))
        	+ _N1 * *(u+1) * conj(*(u+1))
        	+ _N2 * *(u+2) * conj(*(u+2));
        u += 3;
        *snorm += tmp;
    }
}
 
 
static void* thfunc_zlp_snorm2_d_3(void* mp)
{
    matlib_real tmp;
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *u = (matlib_complex*) (ptr->shared_data);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data[0]);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real* snorm   = (matlib_real*) (ptr->nonshared_data[1]);
    u += (4 * start_end_index[0]);

    *snorm = 0;
    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        tmp = 	  _N0 * *(u+0) * conj(*(u+0))
        	+ _N1 * *(u+1) * conj(*(u+1))
        	+ _N2 * *(u+2) * conj(*(u+2))
        	+ _N3 * *(u+3) * conj(*(u+3));
        u += 4;
        *snorm += tmp;
    }
}
 
 
static void* thfunc_zlp_snorm2_d_4(void* mp)
{
    matlib_real tmp;
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *u = (matlib_complex*) (ptr->shared_data);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data[0]);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real* snorm   = (matlib_real*) (ptr->nonshared_data[1]);
    u += (5 * start_end_index[0]);

    *snorm = 0;
    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        tmp = 	  _N0 * *(u+0) * conj(*(u+0))
        	+ _N1 * *(u+1) * conj(*(u+1))
        	+ _N2 * *(u+2) * conj(*(u+2))
        	+ _N3 * *(u+3) * conj(*(u+3))
        	+ _N4 * *(u+4) * conj(*(u+4));
        u += 5;
        *snorm += tmp;
    }
}
 
 
static void* thfunc_zlp_snorm2_d_5(void* mp)
{
    matlib_real tmp;
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *u = (matlib_complex*) (ptr->shared_data);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data[0]);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real* snorm   = (matlib_real*) (ptr->nonshared_data[1]);
    u += (6 * start_end_index[0]);

    *snorm = 0;
    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        tmp = 	  _N0 * *(u+0) * conj(*(u+0))
        	+ _N1 * *(u+1) * conj(*(u+1))
        	+ _N2 * *(u+2) * conj(*(u+2))
        	+ _N3 * *(u+3) * conj(*(u+3))
        	+ _N4 * *(u+4) * conj(*(u+4))
        	+ _N5 * *(u+5) * conj(*(u+5));
        u += 6;
        *snorm += tmp;
    }
}
 
 
static void* thfunc_zlp_snorm2_d_6(void* mp)
{
    matlib_real tmp;
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *u = (matlib_complex*) (ptr->shared_data);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data[0]);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real* snorm   = (matlib_real*) (ptr->nonshared_data[1]);
    u += (7 * start_end_index[0]);

    *snorm = 0;
    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        tmp = 	  _N0 * *(u+0) * conj(*(u+0))
        	+ _N1 * *(u+1) * conj(*(u+1))
        	+ _N2 * *(u+2) * conj(*(u+2))
        	+ _N3 * *(u+3) * conj(*(u+3))
        	+ _N4 * *(u+4) * conj(*(u+4))
        	+ _N5 * *(u+5) * conj(*(u+5))
        	+ _N6 * *(u+6) * conj(*(u+6));
        u += 7;
        *snorm += tmp;
    }
}
 
 
static void* thfunc_zlp_snorm2_d_7(void* mp)
{
    matlib_real tmp;
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *u = (matlib_complex*) (ptr->shared_data);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data[0]);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real* snorm   = (matlib_real*) (ptr->nonshared_data[1]);
    u += (8 * start_end_index[0]);

    *snorm = 0;
    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        tmp = 	  _N0 * *(u+0) * conj(*(u+0))
        	+ _N1 * *(u+1) * conj(*(u+1))
        	+ _N2 * *(u+2) * conj(*(u+2))
        	+ _N3 * *(u+3) * conj(*(u+3))
        	+ _N4 * *(u+4) * conj(*(u+4))
        	+ _N5 * *(u+5) * conj(*(u+5))
        	+ _N6 * *(u+6) * conj(*(u+6))
        	+ _N7 * *(u+7) * conj(*(u+7));
        u += 8;
        *snorm += tmp;
    }
}
 
 
static void* thfunc_zlp_snorm2_d_8(void* mp)
{
    matlib_real tmp;
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *u = (matlib_complex*) (ptr->shared_data);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data[0]);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real* snorm   = (matlib_real*) (ptr->nonshared_data[1]);
    u += (9 * start_end_index[0]);

    *snorm = 0;
    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        tmp = 	  _N0 * *(u+0) * conj(*(u+0))
        	+ _N1 * *(u+1) * conj(*(u+1))
        	+ _N2 * *(u+2) * conj(*(u+2))
        	+ _N3 * *(u+3) * conj(*(u+3))
        	+ _N4 * *(u+4) * conj(*(u+4))
        	+ _N5 * *(u+5) * conj(*(u+5))
        	+ _N6 * *(u+6) * conj(*(u+6))
        	+ _N7 * *(u+7) * conj(*(u+7))
        	+ _N8 * *(u+8) * conj(*(u+8));
        u += 9;
        *snorm += tmp;
    }
}
 
 
static void* thfunc_zlp_snorm2_d_9(void* mp)
{
    matlib_real tmp;
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *u = (matlib_complex*) (ptr->shared_data);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data[0]);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real* snorm   = (matlib_real*) (ptr->nonshared_data[1]);
    u += (10 * start_end_index[0]);

    *snorm = 0;
    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        tmp = 	  _N0 * *(u+0) * conj(*(u+0))
        	+ _N1 * *(u+1) * conj(*(u+1))
        	+ _N2 * *(u+2) * conj(*(u+2))
        	+ _N3 * *(u+3) * conj(*(u+3))
        	+ _N4 * *(u+4) * conj(*(u+4))
        	+ _N5 * *(u+5) * conj(*(u+5))
        	+ _N6 * *(u+6) * conj(*(u+6))
        	+ _N7 * *(u+7) * conj(*(u+7))
        	+ _N8 * *(u+8) * conj(*(u+8))
        	+ _N9 * *(u+9) * conj(*(u+9));
        u += 10;
        *snorm += tmp;
    }
}
 
 
static void* thfunc_zlp_snorm2_d_10(void* mp)
{
    matlib_real tmp;
    matlib_index i;
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;

    matlib_complex *u = (matlib_complex*) (ptr->shared_data);

    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data[0]);
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    matlib_real* snorm   = (matlib_real*) (ptr->nonshared_data[1]);
    u += (11 * start_end_index[0]);

    *snorm = 0;
    for (i=start_end_index[0]; i<start_end_index[1]; i++)
    {
        tmp = 	  _N0 * *(u+0) * conj(*(u+0))
        	+ _N1 * *(u+1) * conj(*(u+1))
        	+ _N2 * *(u+2) * conj(*(u+2))
        	+ _N3 * *(u+3) * conj(*(u+3))
        	+ _N4 * *(u+4) * conj(*(u+4))
        	+ _N5 * *(u+5) * conj(*(u+5))
        	+ _N6 * *(u+6) * conj(*(u+6))
        	+ _N7 * *(u+7) * conj(*(u+7))
        	+ _N8 * *(u+8) * conj(*(u+8))
        	+ _N9 * *(u+9) * conj(*(u+9))
        	+ _N10 * *(u+10) * conj(*(u+10));
        u += 11;
        *snorm += tmp;
    }
}
 
/*============================================================================+/
 | Building Global Mass Matrix 
/+============================================================================*/

static void* pfem1d_thfunc_XCSRGMM2(void* mp)
/* Real CSR - Assemble Global Mass Matrix*/ 
{

    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;
    matlib_index p = *((matlib_index*) (ptr->shared_data[0]));
    matlib_index N = *((matlib_index*) (ptr->shared_data[1]));
    matlib_xm q    =    *((matlib_xm*) (ptr->shared_data[2]));

    matlib_xm_nsparse M = *((matlib_xm_nsparse*) (ptr->shared_data[3]));

    debug_enter( "Thread id: %d, degree of polynomial: %d, "
                 "Number of finite-elements: %d, "
                 "length of q: %d-by-%d", 
                 ptr->thread_index,
                 p, N, q.lenc, q.lenr);

   matlib_index nr_combi = q.lenc/N; 
   debug_body("nr. of combinations of basis functions: %d", nr_combi);
   

    /* number of nonzero elements nnz = N*p*(p+3)/2+1;
     * */
    matlib_index s, i, l, l0, m, st;

    i = 0;
    s = 0;
    
    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    q.elem_p += (q.lenc*start_end_index[0]);

    matlib_real **ugpmm;

    for( ugpmm = M.elem_p+start_end_index[0]; 
         ugpmm < M.elem_p+start_end_index[1]; ugpmm++)
    {
        (*ugpmm)[s] = *(q.elem_p);
        s++;
    
        (*ugpmm)[s] = *(q.elem_p+1);
        s++;

        matlib_index shiftq = 3;
        
        for (l=N+1; l<(N+p); l++)
        {
            (*ugpmm)[s] = (*(q.elem_p + shiftq));
            shiftq ++;
            s++;
        }

        l0 = N+1;
        st = 0;
        for( i=1; i<N; i++)
        {
            /* i-th row: start from the diagonal 
             * (v_i,v_i)_{K_i} + (v_i,v_i)_{K_{i+1}}
             * (v_i,v_{i+1})_{K_{i+1}}
             * */ 
            (*ugpmm)[s] = *(q.elem_p+st+2) + *(q.elem_p+st+nr_combi);
            s++;
            (*ugpmm)[s] = *(q.elem_p+st+nr_combi+1); 
            s++;
            /* v_i interacts with [b_i?], [b_{i+1,?}]
             * (p-1) bubble functions in each element
             */
            shiftq = st+3+p-1;
            for (l=l0; l<(l0+p-1); l++)
            {
                (*ugpmm)[s] = *(q.elem_p + shiftq);
                shiftq ++;
                s++;
            }
            l0 = l;

            shiftq = st + nr_combi+3;
            for (l=l0; l<(l0+p-1); l++)
            {
                (*ugpmm)[s] = *(q.elem_p+shiftq);
                shiftq++;
                s++;
            
            }
            st = st+nr_combi;
        }
        debug_body("s=%d", s);
        /* last row: only the diagonal 
         * (v_N,v_N)_{K_N}
         * */ 
        (*ugpmm)[s] = *(q.elem_p+st+2);
        s++;
        shiftq = st+3+(p-1);
        for (l=l0; l<(l0+p-1); l++)
        {
            (*ugpmm)[s] = *(q.elem_p+shiftq);
            shiftq++;
            s++;
        }       
        i++;

        /* vertex-vertex and vertex-bubble combinations are done
         * do the bubble-bubble combination
         * */ 
        matlib_index shiftq0 = 3+2*(p-1);

        for( st=0; st<N*nr_combi; st +=nr_combi)
        {
            shiftq = shiftq0;
            for( l=0; l<p-1; l++)
            {
                for (m=i; m<(i+p-1-l); m++)
                {
                    (*ugpmm)[s] = *(q.elem_p+st+shiftq);
                    shiftq++;
                    s++;
                }
                i++;
            }
        }
        q.elem_p += q.lenc;
    }

    debug_exit("%s", "");
}

void pfem1d_xm_nsparse_GMM
/* Double - Assemble Global Mass Matrix*/ 
(
    matlib_index       p,
    matlib_index       N,
    matlib_xm          Q,
    matlib_xm*         phi,
    matlib_xm*         q,
    matlib_xm_nsparse* M,
    matlib_index       num_threads,
    pthpool_data_t*    mp
    
)
{
    debug_enter( "poynomial degree: %d, "
                 "nr. of fem-elements: %d, "
                 "nr. of sparse matrices: %d, "
                 "size of Q: %d-by-%d), ",
                 p, N, M->nsparse, Q.lenc, Q.lenr);

    assert(Q.lenc==(p-1)*(p+4)/2+3);
    assert(q.lenr==M->nsparse);

    matlib_index i; 
    pfem1d_XFLT2( N, Q, *phi, *q, num_threads, mp);
    /* define the shared data */ 
    void* shared_data[4] = { (void*) &p,
                             (void*) &N,
                             (void*) q,
                             (void*) M };

    matlib_index nsdata[num_threads][2];

    pthpool_arg_t   arg[num_threads];
    pthpool_task_t  task[num_threads];

    /* define the block of data per thread */ 
    matlib_index Np = M->nsparse/(num_threads);

    for(i=0; i<num_threads-1; i++)
    {
        nsdata[i][0] = i*Np;
        nsdata[i][1] = (i+1)*Np;
        arg[i].shared_data    = shared_data; 
        arg[i].nonshared_data = (void**)&nsdata[i];
        arg[i].thread_index   = i;
        /* Define the task */ 
        task[i].function  = (void*)pfem1d_thfunc_XCSRGMM2;
        task[i].argument  = &arg[i];
        mp[i].task = &task[i];
    }
    nsdata[i][0] = i*Np;
    nsdata[i][1] = M->nsparse;
    arg[i].shared_data    = shared_data; 
    arg[i].nonshared_data = (void**)&nsdata[i];
    arg[i].thread_index   = i;
    /* Define the task */ 
    task[i].function  = (void*)pfem1d_thfunc_XCSRGMM2;
    task[i].argument  = &arg[i];
    mp[i].task = &task[i];

    debug_body("%s", "created task");

    pthpool_exec_task(num_threads, mp, task);
    
    debug_exit("%s", "");
}

/*============================================================================*/
static void* pfem1d_thfunc_ZCSRGMM2(void* mp)
{
    pthpool_arg_t *ptr = (pthpool_arg_t*) mp;
    matlib_index p = *((matlib_index*) (ptr->shared_data[0]));
    matlib_index N = *((matlib_index*) (ptr->shared_data[1]));
    matlib_zm q    =    *((matlib_zm*) (ptr->shared_data[2]));

    matlib_zm_nsparse M = *((matlib_zm_nsparse*) (ptr->shared_data[3]));

    debug_enter( "Thread id: %d, degree of polynomial: %d, "
                 "Number of finite-elements: %d, "
                 "length of q: %d-by-%d", 
                 ptr->thread_index,
                 p, N, q.lenc, q.lenr);


   matlib_index nr_combi = q.lenc/N; 
   debug_body("nr. of combinations of basis functions: %d", nr_combi);
   

    /* number of nonzero elements nnz = N*p*(p+3)/2+1;
     * */
    matlib_index s, i, l, l0, m, st;

    /* zero-th row 
     *
     * (v_0,v_0)_{K_1}
     * (v_0,v_1)_{K_1}
     * */ 
    matlib_index* start_end_index = (matlib_index*) (ptr->nonshared_data);
    
    debug_body( "Thread id: %d, start_index: %d, end_index: %d",
                ptr->thread_index, 
                start_end_index[0], start_end_index[1]);

    q.elem_p += (q.lenc*start_end_index[0]);

    matlib_complex **ugpmm;

    for( ugpmm = M.elem_p+start_end_index[0]; 
         ugpmm < M.elem_p+start_end_index[1]; ugpmm++)
    {
        i = 0;
        s = 0;
        (*ugpmm)[s] = *(q.elem_p);
        s++;
    
        (*ugpmm)[s] = *(q.elem_p+1);
        s++;

        matlib_index shiftq = 3;
        
        for (l=N+1; l<(N+p); l++)
        {
            (*ugpmm)[s] = (*(q.elem_p + shiftq));
            shiftq ++;
            s++;
        }

        l0 = N+1;
        st = 0;
        for( i=1; i<N; i++)
        {
            /* i-th row: start from the diagonal 
             * (v_i,v_i)_{K_i} + (v_i,v_i)_{K_{i+1}}
             * (v_i,v_{i+1})_{K_{i+1}}
             * */ 
            (*ugpmm)[s] = *(q.elem_p+st+2) + *(q.elem_p+st+nr_combi);
            s++;
            (*ugpmm)[s] = *(q.elem_p+st+nr_combi+1); 
            s++;
            /* v_i interacts with [b_i?], [b_{i+1,?}]
             * (p-1) bubble functions in each element
             */
            shiftq = st+3+p-1;
            for (l=l0; l<(l0+p-1); l++)
            {
                (*ugpmm)[s] = *(q.elem_p + shiftq);
                shiftq ++;
                s++;
            }
            l0 = l;

            shiftq = st + nr_combi+3;
            for (l=l0; l<(l0+p-1); l++)
            {
                (*ugpmm)[s] = *(q.elem_p+shiftq);
                shiftq++;
                s++;
            
            }
            st = st+nr_combi;
        }
        debug_body("s=%d", s);
        /* last row: only the diagonal 
         * (v_N,v_N)_{K_N}
         * */ 
        (*ugpmm)[s] = *(q.elem_p+st+2);
        s++;
        shiftq = st+3+(p-1);
        for (l=l0; l<(l0+p-1); l++)
        {
            (*ugpmm)[s] = *(q.elem_p+shiftq);
            shiftq++;
            s++;
        }
        i++;

        /* vertex-vertex and vertex-bubble combinations are done
         * do the bubble-bubble combination
         * */ 
        matlib_index shiftq0 = 3+2*(p-1);

        for( st=0; st<N*nr_combi; st +=nr_combi)
        {
            shiftq = shiftq0;
            for( l=0; l<p-1; l++)
            {
                for (m=i; m<(i+p-1-l); m++)
                {
                    (*ugpmm)[s] = *(q.elem_p+st+shiftq);
                    shiftq++;
                    s++;
                }
                i++;
            }
        }
        q.elem_p += q.lenc;
    }

    debug_exit("%s", "");
}

void pfem1d_zm_nsparse_GMM
/* Double - Assemble Global Mass Matrix*/ 
(
    matlib_index       p,
    matlib_index       N,
    matlib_xm          Q,
    matlib_zm*         phi,
    matlib_zm*         q,
    matlib_zm_nsparse* M,
    matlib_index       num_threads,
    pthpool_data_t*    mp
    
)
{
    debug_enter( "poynomial degree: %d, "
                 "nr. of fem-elements: %d, "
                 "nr. of sparse matrices: %d, "
                 "size of Q: %d-by-%d), ",
                 p, N, M->nsparse, Q.lenc, Q.lenr);

    assert(Q.lenc==(p-1)*(p+4)/2+3);
    assert(q.lenr==M->nsparse);

    matlib_index i; 
    pfem1d_ZFLT2( N, Q, *phi, *q, num_threads, mp);
    /* define the shared data */ 
    void* shared_data[4] = { (void*) &p,
                             (void*) &N,
                             (void*) q,
                             (void*) M };

    matlib_index nsdata[num_threads][2];

    pthpool_arg_t   arg[num_threads];
    pthpool_task_t  task[num_threads];

    /* define the block of data per thread */ 
    matlib_index Np = M->nsparse/(num_threads);

    for(i=0; i<num_threads-1; i++)
    {
        nsdata[i][0] = i*Np;
        nsdata[i][1] = (i+1)*Np;
        arg[i].shared_data    = shared_data; 
        arg[i].nonshared_data = (void**)&nsdata[i];
        arg[i].thread_index   = i;
        /* Define the task */ 
        task[i].function  = (void*)pfem1d_thfunc_ZCSRGMM2;
        task[i].argument  = &arg[i];
        mp[i].task = &task[i];
    }
    nsdata[i][0] = i*Np;
    nsdata[i][1] = M->nsparse;
    arg[i].shared_data    = shared_data; 
    arg[i].nonshared_data = (void**)&nsdata[i];
    arg[i].thread_index   = i;
    /* Define the task */ 
    task[i].function  = (void*)pfem1d_thfunc_ZCSRGMM2;
    task[i].argument  = &arg[i];
    mp[i].task = &task[i];

    debug_body("%s", "created task");

    pthpool_exec_task(num_threads, mp, task);

    debug_exit("%s", "");
}

