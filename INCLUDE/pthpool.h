#ifndef PTHPOOL_H
#define PTHPOOL_H

/*============================================================================+/
 | Include all the dependencies
/+============================================================================*/
#include "basic.h"
#include "matlib.h"
#include "debug.h"
#include "ehandler.h"
/*============================================================================+/
 |DATA STRUCTURES AND ENUMS
/+============================================================================*/
#define MAX_NUM_CPU 4

typedef enum
{
    PTHPOOL_WAIT,
    PTHPOOL_PERFORM,
    PTHPOOL_NOTIFY,
    PTHPOOL_EXIT

} PTHPOOL_ACTION;

/* Argument of the function to be executed in any thread */ 
typedef struct
{
    void**       shared_data; /* array of pointers */ 
    void**       nonshared_data;
    matlib_index thread_index;

} pthpool_arg_t;
/* A task consist of executing a function.
 * This can be a generic function.
 *
 * */ 
typedef struct
{
    void          (*function) (void *);
    pthpool_arg_t *argument;

} pthpool_task_t;

typedef struct
{
    pthread_t       thread;
    cpu_set_t       cpu;
    pthread_mutex_t lock;
    pthread_cond_t  notify;
    pthpool_task_t* task;
    matlib_index    thread_index;
    PTHPOOL_ACTION  action;

} pthpool_data_t;
/*============================================================================*/

void pthpool_create_threads
( 
    matlib_index    num_threads, 
    pthpool_data_t* mp
);

void pthpool_exec_task
( 
    matlib_index    num_threads, 
    pthpool_data_t* mp, 
    pthpool_task_t* task
);

void pthpool_destroy_threads
( 
    matlib_index    num_threads, 
    pthpool_data_t* mp 
);
/* Parallelize evaluation of functions defined for vectors */ 
void pthpool_func
(
    matlib_index*   Np,
    void**          shared_data,
    void*           thfunc,
    matlib_index    num_threads,
    pthpool_data_t* mp

);
#endif

