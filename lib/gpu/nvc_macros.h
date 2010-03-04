/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Mike Brown (SNL), wmbrown@sandia.gov
                         Peng Wang (Nvidia), penwang@nvidia.com
                         Paul Crozier (SNL), pscrozi@sandia.gov
------------------------------------------------------------------------- */

#ifndef NVC_MACROS_H
#define NVC_MACROS_H

#if defined(__APPLE__)
#if _GLIBCXX_ATOMIC_BUILTINS == 1
#undef _GLIBCXX_ATOMIC_BUILTINS
#endif // _GLIBCXX_ATOMIC_BUILTINS
#endif // __APPLE__

#include <stdio.h>
#include "math_constants.h"
#define INT_MUL(x,y)	(__mul24(x,y))
//#define INT_MUL(x,y) ((x)*(y))

template <class numbr>
static __inline__ __device__ numbr cuda_inf() { return CUDART_INF_F; }

#ifdef CUDA_DOUBLE
template <>
static __inline__ __device__ double cuda_inf<double>() { return CUDART_INF; }
#endif

template <class numbr>
static __inline__ __device__ numbr cuda_zero() { return 0.0; }

template <>
static __inline__ __device__ float cuda_zero<float>() { return 0.0f; }

#ifndef NO_DEBUG

#  define CU_SAFE_CALL_NO_SYNC( call ) do {                                  \
    CUresult err = call;                                                     \
    if( CUDA_SUCCESS != err) {                                               \
        fprintf(stderr, "Cuda driver error %x in file '%s' in line %i.\n",   \
                err, __FILE__, __LINE__ );                                   \
        exit(EXIT_FAILURE);                                                  \
    } } while (0)

#  define CU_SAFE_CALL( call ) do {                                          \
    CU_SAFE_CALL_NO_SYNC(call);                                              \
    CUresult err = cuCtxSynchronize();                                       \
    if( CUDA_SUCCESS != err) {                                               \
        fprintf(stderr, "Cuda driver error %x in file '%s' in line %i.\n",   \
                err, __FILE__, __LINE__ );                                   \
        exit(EXIT_FAILURE);                                                  \
    } } while (0)

#  define CUDA_SAFE_CALL_NO_SYNC( call) do {                                 \
    cudaError err = call;                                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        exit(EXIT_FAILURE);                                                  \
    } } while (0)

#  define CUDA_SAFE_CALL( call) do {                                         \
    CUDA_SAFE_CALL_NO_SYNC(call);                                            \
    cudaError err = cudaThreadSynchronize();                                 \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        exit(EXIT_FAILURE);                                                  \
    } } while (0)

#  define CUFFT_SAFE_CALL( call) do {                                        \
    cufftResult err = call;                                                  \
    if( CUFFT_SUCCESS != err) {                                              \
        fprintf(stderr, "CUFFT error in file '%s' in line %i.\n",            \
                __FILE__, __LINE__);                                         \
        exit(EXIT_FAILURE);                                                  \
    } } while (0)

#  define CUT_SAFE_CALL( call)                                               \
    if( CUTTrue != call) {                                                   \
        fprintf(stderr, "Cut error in file '%s' in line %i.\n",              \
                __FILE__, __LINE__);                                         \
        exit(EXIT_FAILURE);                                                  \
    } 

    //! Check for CUDA error
#  define CUT_CHECK_ERROR(errorMessage) do {                                 \
    cudaError_t err = cudaGetLastError();                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
                errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    err = cudaThreadSynchronize();                                           \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
                errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
        exit(EXIT_FAILURE);                                                  \
    } } while (0)

    //! Check for malloc error
#  define CUT_SAFE_MALLOC( mallocCall ) do{                                  \
    if( !(mallocCall)) {                                                     \
        fprintf(stderr, "Host malloc failure in file '%s' in line %i\n",     \
                __FILE__, __LINE__);                                         \
        exit(EXIT_FAILURE);                                                  \
    } } while(0);

    //! Check if conditon is true (flexible assert)
#  define CUT_CONDITION( val)                                                \
    if( CUTFalse == cutCheckCondition( val, __FILE__, __LINE__)) {           \
        exit(EXIT_FAILURE);                                                  \
    }

#else  // not DEBUG

#define CUT_BANK_CHECKER( array, index)  array[index]

    // void macros for performance reasons
#  define CUT_CHECK_ERROR(errorMessage)
#  define CUT_CHECK_ERROR_GL()
#  define CUT_CONDITION( val) 
#  define CU_SAFE_CALL_NO_SYNC( call) call
#  define CU_SAFE_CALL( call) call
#  define CUDA_SAFE_CALL_NO_SYNC( call) call
#  define CUDA_SAFE_CALL( call) call
#  define CUT_SAFE_CALL( call) call
#  define CUFFT_SAFE_CALL( call) call
#  define CUT_SAFE_MALLOC( mallocCall ) mallocCall

#endif

#endif
