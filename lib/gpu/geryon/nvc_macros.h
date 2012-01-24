#ifndef NVC_MACROS_H
#define NVC_MACROS_H

#if defined(__APPLE__)
#if _GLIBCXX_ATOMIC_BUILTINS == 1
#undef _GLIBCXX_ATOMIC_BUILTINS
#endif // _GLIBCXX_ATOMIC_BUILTINS
#endif // __APPLE__

#include <stdio.h>
#include <cassert>
#include <cuda_runtime.h>

#ifdef MPI_GERYON
#include "mpi.h"
#define NVC_GERYON_EXIT MPI_Abort(MPI_COMM_WORLD,-1)
#else
#define NVC_GERYON_EXIT assert(0==1)
#endif

#ifndef UCL_GERYON_EXIT
#define UCL_GERYON_EXIT NVC_GERYON_EXIT
#endif

#ifdef UCL_DEBUG
#define UCL_SYNC_DEBUG
#define UCL_DESTRUCT_CHECK
#endif

#ifndef UCL_NO_API_CHECK

#define CUDA_SAFE_CALL_NS( call) do {                                        \
    cudaError err = call;                                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in call at file '%s' in line %i : %s.\n", \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        NVC_GERYON_EXIT;                                                     \
    } } while (0)

#ifdef UCL_SYNC_DEBUG

#define CUDA_SAFE_CALL( call) do {                                           \
    CUDA_SAFE_CALL_NS( call);                                                \
    cudaError err=cudaThreadSynchronize();                                   \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        NVC_GERYON_EXIT;                                                     \
    } } while (0)

#else

#define CUDA_SAFE_CALL( call) CUDA_SAFE_CALL_NS( call)

#endif

#else  // not DEBUG

// void macros for performance reasons
#define CUDA_SAFE_CALL( call) call
#define CUDA_SAFE_CALL_NS( call) call

#endif

#ifdef UCL_DESTRUCT_CHECK

#define CUDA_DESTRUCT_CALL( call) CUDA_SAFE_CALL( call)
#define CUDA_DESTRUCT_CALL_NS( call) CUDA_SAFE_CALL_NS( call)

#else

#define CUDA_DESTRUCT_CALL( call) call
#define CUDA_DESTRUCT_CALL_NS( call) call

#endif

#endif

