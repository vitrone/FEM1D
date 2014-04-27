#ifndef PFEM1D_H
#define PFEM1D_H

/*============================================================================+/
 | Include all the dependencies
/+============================================================================*/
#include "basic.h"
#include "matlib.h"
#include "pthpool.h"
#include "debug.h"
#include "ehandler.h"
/*============================================================================+/
 |DATA STRUCTURES AND ENUMS
/+============================================================================*/
void pfem1d_DFLT
(
    const matlib_index    N,
    const matlib_dm       FM,
          matlib_dv       u,
          matlib_dv       U,
          matlib_index    num_threads,
          pthpool_data_t* mp
);

void pfem1d_ZFLT
(
    const matlib_index    N,
    const matlib_dm       FM,
          matlib_zv       u,
          matlib_zv       U,
          matlib_index    num_threads,
          pthpool_data_t* mp
);

void pfem1d_DILT
(
    const matlib_index    N,
    const matlib_dm       IM,
          matlib_dv       U,
          matlib_dv       u,
          matlib_index    num_threads,
          pthpool_data_t* mp
);

void pfem1d_ZILT
(
    const matlib_index    N,
    const matlib_dm       IM,
          matlib_zv       U,
          matlib_zv       u,
          matlib_index    num_threads,
          pthpool_data_t* mp
);

void pfem1d_DFLT2
(
    const matlib_index    N,
    const matlib_dm       FM,
          matlib_dm       u,
          matlib_dm       U,
          matlib_index    num_threads,
          pthpool_data_t* mp
);
void pfem1d_ZFLT2
(
    const matlib_index    N,
    const matlib_dm       FM,
          matlib_zm       u,
          matlib_zm       U,
          matlib_index    num_threads,
          pthpool_data_t* mp
);
void pfem1d_DF2L
(
    const matlib_index p, 
    const matlib_dv    vb,
          matlib_dv    u,
          matlib_index    num_threads,
          pthpool_data_t* mp
);

void pfem1d_ZF2L
(
    const matlib_index    p, 
    const matlib_zv       vb,
          matlib_zv       u,
          matlib_index    num_threads,
          pthpool_data_t* mp
);

void pfem1d_DPrjL2F
(
    matlib_index    p,
    matlib_dv       u,
    matlib_dv       Pvb,
    matlib_index    num_threads,
    pthpool_data_t* mp
);

void pfem1d_ZPrjL2F
(
    matlib_index    p,
    matlib_zv       u,
    matlib_zv       Pvb,
    matlib_index    num_threads,
    pthpool_data_t* mp
);

double pfem1d_DNorm2
(
    matlib_index    p,
    matlib_index    N,
    matlib_dv       u,
    matlib_index    num_threads,
    pthpool_data_t* mp
);

double pfem1d_ZNorm2
(
    matlib_index    p,
    matlib_index    N,
    matlib_zv       u,
    matlib_index    num_threads,
    pthpool_data_t* mp
);









#endif
