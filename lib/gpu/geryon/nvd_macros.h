#ifndef NVD_MACROS_H
#define NVD_MACROS_H

#include <stdio.h>
#include <cassert>
#include <cuda.h>

#if CUDA_VERSION >= 3020
#define CUDA_INT_TYPE size_t
#else
#define CUDA_INT_TYPE unsigned
#endif

#ifdef MPI_GERYON
#include "mpi.h"
#define NVD_GERYON_EXIT do {                                               \
  int is_final;                                                            \
  MPI_Finalized(&is_final);                                                \
  if (!is_final)                                                           \
    MPI_Abort(MPI_COMM_WORLD,-1);                                          \
  } while(0)
#else
#define NVD_GERYON_EXIT assert(0==1)
#endif

#ifndef UCL_GERYON_EXIT
#define UCL_GERYON_EXIT NVD_GERYON_EXIT
#endif

#ifdef UCL_DEBUG
#define UCL_SYNC_DEBUG
#define UCL_DESTRUCT_CHECK
#endif

#ifndef UCL_NO_API_CHECK

#define CU_SAFE_CALL_NS( call ) do {                                         \
    CUresult err = call;                                                     \
    if( CUDA_SUCCESS != err) {                                               \
        fprintf(stderr, "Cuda driver error %d in call at file '%s' in line %i.\n",   \
                err, __FILE__, __LINE__ );                                   \
        NVD_GERYON_EXIT;                                                     \
    } } while (0)

#ifdef UCL_SYNC_DEBUG

#define CU_SAFE_CALL( call ) do {                                            \
    CU_SAFE_CALL_NS( call );                                                 \
    CUresult err=cuCtxSynchronize();                                                  \
    if( CUDA_SUCCESS != err) {                                               \
        fprintf(stderr, "Cuda driver error %d in file '%s' in line %i.\n",   \
                err, __FILE__, __LINE__ );                                   \
        NVD_GERYON_EXIT;                                                     \
    } } while (0)

#else

#define CU_SAFE_CALL( call ) CU_SAFE_CALL_NS( call )

#endif

#else  // not DEBUG

// void macros for performance reasons
#define CU_SAFE_CALL_NS( call ) call
#define CU_SAFE_CALL( call) call

#endif

#ifdef UCL_DESTRUCT_CHECK

#define CU_DESTRUCT_CALL( call) CU_SAFE_CALL( call)
#define CU_DESTRUCT_CALL_NS( call) CU_SAFE_CALL_NS( call)

#else

#define CU_DESTRUCT_CALL( call) call
#define CU_DESTRUCT_CALL_NS( call) call

#endif

#endif

