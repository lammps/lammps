#ifndef HIP_MACROS_H
#define HIP_MACROS_H

#include <cstdio>
#include <cassert>
#include <hip/hip_runtime.h>

#define CUDA_INT_TYPE size_t

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
    hipError_t err = call;                                                     \
    if( hipSuccess != err) {                                               \
        fprintf(stderr, "HIP runtime error %d in call at file '%s' in line %i.\n",   \
                err, __FILE__, __LINE__ );                                   \
        NVD_GERYON_EXIT;                                                     \
    } } while (0)

#ifdef UCL_SYNC_DEBUG

#define CU_SAFE_CALL( call ) do {                                            \
    CU_SAFE_CALL_NS( call );                                                 \
    hipError_t err=hipCtxSynchronize();                                                  \
    if( hipSuccess != err) {                                               \
        fprintf(stderr, "HIP runtime error %d in file '%s' in line %i.\n",   \
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

