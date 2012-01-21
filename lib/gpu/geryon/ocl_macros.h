#ifndef OCL_MACROS_H
#define OCL_MACROS_H

#include <stdio.h>
#include <cassert>

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#ifdef MPI_GERYON
#include "mpi.h"
#define OCL_GERYON_EXIT MPI_Abort(MPI_COMM_WORLD,-1)
#else
#define OCL_GERYON_EXIT assert(0==1)
#endif

#ifndef UCL_GERYON_EXIT
#define UCL_GERYON_EXIT OCL_GERYON_EXIT
#endif

#ifdef UCL_DEBUG
#define UCL_SYNC_DEBUG
#define UCL_DESTRUCT_CHECK
#endif

#ifndef UCL_NO_API_CHECK

#  define CL_SAFE_CALL( call) do {                                         \
    cl_int err = call;                                                     \
    if( err != CL_SUCCESS) {                                               \
        fprintf(stderr, "OpenCL error in file '%s' in line %i : %d.\n",    \
                __FILE__, __LINE__, err );                                 \
        OCL_GERYON_EXIT;                                                           \
    } } while (0)

#  define CL_CHECK_ERR( val) do {                                        \
    if( val != CL_SUCCESS) {                                               \
        fprintf(stderr, "OpenCL error in file '%s' in line %i : %d.\n",    \
                __FILE__, __LINE__, val );                                 \
        OCL_GERYON_EXIT;                                                           \
    } } while (0)

#else  // not DEBUG

// void macros for performance reasons
#  define CL_SAFE_CALL( call) call
#  define CL_CHECK_ERR( val)

#endif

#ifdef UCL_DESTRUCT_CHECK

#define CL_DESTRUCT_CALL( call) CL_SAFE_CALL( call)

#else

#define CL_DESTRUCT_CALL( call) call

#endif

#endif
