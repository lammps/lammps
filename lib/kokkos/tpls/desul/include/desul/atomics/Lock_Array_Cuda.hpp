/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_LOCK_ARRAY_CUDA_HPP_
#define DESUL_ATOMICS_LOCK_ARRAY_CUDA_HPP_

#include "desul/atomics/Common.hpp"
#include "desul/atomics/Macros.hpp"

#ifdef DESUL_HAVE_CUDA_ATOMICS

#include <cstdint>

namespace desul {
namespace Impl {

#ifdef __CUDA_ARCH__
#define DESUL_IMPL_BALLOT_MASK(m, x) __ballot_sync(m, x)
#define DESUL_IMPL_ACTIVEMASK __activemask()
#else
#define DESUL_IMPL_BALLOT_MASK(m, x) m == 0 ? 0 : 1
#define DESUL_IMPL_ACTIVEMASK 0
#endif

/// \brief This global variable in Host space is the central definition
///        of these arrays.
extern int32_t* CUDA_SPACE_ATOMIC_LOCKS_DEVICE_h;
extern int32_t* CUDA_SPACE_ATOMIC_LOCKS_NODE_h;

/// \brief After this call, the g_host_cuda_lock_arrays variable has
///        valid, initialized arrays.
///
/// This call is idempotent.
/// The function is templated to make it a weak symbol to deal with Kokkos/RAJA
///   snapshotted version while also linking against pure Desul
template <typename /*AlwaysInt*/ = int>
void init_lock_arrays_cuda();

/// \brief After this call, the g_host_cuda_lock_arrays variable has
///        all null pointers, and all array memory has been freed.
///
/// This call is idempotent.
/// The function is templated to make it a weak symbol to deal with Kokkos/RAJA
///   snapshotted version while also linking against pure Desul
template <typename /*AlwaysInt*/ = int>
void finalize_lock_arrays_cuda();

}  // namespace Impl
}  // namespace desul

#if defined(__CUDACC__)

namespace desul {
namespace Impl {

/// \brief This global variable in CUDA space is what kernels use
///        to get access to the lock arrays.
///
/// When relocatable device code is enabled, there can be one single
/// instance of this global variable for the entire executable,
/// whose definition will be in Kokkos_Cuda_Locks.cpp (and whose declaration
/// here must then be extern.
/// This one instance will be initialized by initialize_host_cuda_lock_arrays
/// and need not be modified afterwards.
///
/// When relocatable device code is disabled, an instance of this variable
/// will be created in every translation unit that sees this header file
/// (we make this clear by marking it static, meaning no other translation
///  unit can link to it).
/// Since the Kokkos_Cuda_Locks.cpp translation unit cannot initialize the
/// instances in other translation units, we must update this CUDA global
/// variable based on the Host global variable prior to running any kernels
/// that will use it.
/// That is the purpose of the KOKKOS_ENSURE_CUDA_LOCK_ARRAYS_ON_DEVICE macro.
__device__
#ifdef __CUDACC_RDC__
    __constant__ extern
#endif
    int32_t* CUDA_SPACE_ATOMIC_LOCKS_DEVICE;

__device__
#ifdef __CUDACC_RDC__
    __constant__ extern
#endif
    int32_t* CUDA_SPACE_ATOMIC_LOCKS_NODE;

#define CUDA_SPACE_ATOMIC_MASK 0x1FFFF

/// \brief Acquire a lock for the address
///
/// This function tries to acquire the lock for the hash value derived
/// from the provided ptr. If the lock is successfully acquired the
/// function returns true. Otherwise it returns false.
__device__ inline bool lock_address_cuda(void* ptr, desul::MemoryScopeDevice) {
  size_t offset = size_t(ptr);
  offset = offset >> 2;
  offset = offset & CUDA_SPACE_ATOMIC_MASK;
  return (0 == atomicExch(&desul::Impl::CUDA_SPACE_ATOMIC_LOCKS_DEVICE[offset], 1));
}
__device__ inline bool lock_address_cuda(void* ptr, desul::MemoryScopeNode) {
  size_t offset = size_t(ptr);
  offset = offset >> 2;
  offset = offset & CUDA_SPACE_ATOMIC_MASK;
  return (0 == atomicExch(&desul::Impl::CUDA_SPACE_ATOMIC_LOCKS_NODE[offset], 1));
}

/// \brief Release lock for the address
///
/// This function releases the lock for the hash value derived
/// from the provided ptr. This function should only be called
/// after previously successfully acquiring a lock with
/// lock_address.
__device__ inline void unlock_address_cuda(void* ptr, desul::MemoryScopeDevice) {
  size_t offset = size_t(ptr);
  offset = offset >> 2;
  offset = offset & CUDA_SPACE_ATOMIC_MASK;
  atomicExch(&desul::Impl::CUDA_SPACE_ATOMIC_LOCKS_DEVICE[offset], 0);
}
__device__ inline void unlock_address_cuda(void* ptr, desul::MemoryScopeNode) {
  size_t offset = size_t(ptr);
  offset = offset >> 2;
  offset = offset & CUDA_SPACE_ATOMIC_MASK;
  atomicExch(&desul::Impl::CUDA_SPACE_ATOMIC_LOCKS_NODE[offset], 0);
}

}  // namespace Impl
}  // namespace desul

// Make lock_array_copied an explicit translation unit scope thingy
namespace desul {
namespace Impl {
namespace {
static int lock_array_copied = 0;
inline int eliminate_warning_for_lock_array() { return lock_array_copied; }
}  // namespace
}  // namespace Impl
}  // namespace desul
/* It is critical that this code be a macro, so that it will
   capture the right address for desul::Impl::CUDA_SPACE_ATOMIC_LOCKS_DEVICE
   putting this in an inline function will NOT do the right thing! */
#define DESUL_IMPL_COPY_CUDA_LOCK_ARRAYS_TO_DEVICE()                       \
  {                                                                        \
    if (::desul::Impl::lock_array_copied == 0) {                           \
      cudaMemcpyToSymbol(::desul::Impl::CUDA_SPACE_ATOMIC_LOCKS_DEVICE,    \
                         &::desul::Impl::CUDA_SPACE_ATOMIC_LOCKS_DEVICE_h, \
                         sizeof(int32_t*));                                \
      cudaMemcpyToSymbol(::desul::Impl::CUDA_SPACE_ATOMIC_LOCKS_NODE,      \
                         &::desul::Impl::CUDA_SPACE_ATOMIC_LOCKS_NODE_h,   \
                         sizeof(int32_t*));                                \
    }                                                                      \
    ::desul::Impl::lock_array_copied = 1;                                  \
  }

#endif /* defined( __CUDACC__ ) */

#endif /* defined( DESUL_HAVE_CUDA_ATOMICS ) */

#if defined(__CUDACC_RDC__) || (!defined(__CUDACC__))
#define DESUL_ENSURE_CUDA_LOCK_ARRAYS_ON_DEVICE()
#else
#define DESUL_ENSURE_CUDA_LOCK_ARRAYS_ON_DEVICE() \
  DESUL_IMPL_COPY_CUDA_LOCK_ARRAYS_TO_DEVICE()
#endif

#endif /* #ifndef KOKKOS_CUDA_LOCKS_HPP_ */
