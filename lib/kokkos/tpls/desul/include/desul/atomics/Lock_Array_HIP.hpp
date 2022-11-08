/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_LOCK_ARRAY_HIP_HPP_
#define DESUL_ATOMICS_LOCK_ARRAY_HIP_HPP_

#include "desul/atomics/Common.hpp"
#include "desul/atomics/Macros.hpp"

#ifdef DESUL_HAVE_HIP_ATOMICS

#include <hip/hip_runtime.h>

#include <cstdint>

namespace desul {
namespace Impl {

#ifdef __HIP_DEVICE_COMPILE__
#define DESUL_IMPL_BALLOT_MASK(x) __ballot(x)
#else
#define DESUL_IMPL_BALLOT_MASK(x) 0
#endif

/**
 * \brief This global variable in Host space is the central definition of these
 * arrays.
 */
extern int32_t* HIP_SPACE_ATOMIC_LOCKS_DEVICE_h;
extern int32_t* HIP_SPACE_ATOMIC_LOCKS_NODE_h;

/// \brief After this call, the g_host_cuda_lock_arrays variable has
///        valid, initialized arrays.
///
/// This call is idempotent.
/// The function is templated to make it a weak symbol to deal with Kokkos/RAJA
///   snapshotted version while also linking against pure Desul
template <typename /*AlwaysInt*/ = int>
void init_lock_arrays_hip();

/// \brief After this call, the g_host_cuda_lock_arrays variable has
///        all null pointers, and all array memory has been freed.
///
/// This call is idempotent.
/// The function is templated to make it a weak symbol to deal with Kokkos/RAJA
///   snapshotted version while also linking against pure Desul
template <typename /*AlwaysInt*/ = int>
void finalize_lock_arrays_hip();
}  // namespace Impl
}  // namespace desul

#ifdef __HIPCC__
namespace desul {
namespace Impl {

/**
 * \brief This global variable in HIP space is what kernels use to get access
 * to the lock arrays.
 *
 * When relocatable device code is enabled, there can be one single instance of
 * this global variable for the entire executable, whose definition will be in
 * Kokkos_HIP_Locks.cpp (and whose declaration here must then be extern.  This
 * one instance will be initialized by initialize_host_hip_lock_arrays and need
 * not be modified afterwards.
 *
 * When relocatable device code is disabled, an instance of this variable will
 * be created in every translation unit that sees this header file (we make this
 * clear by marking it static, meaning no other translation unit can link to
 * it). Since the Kokkos_HIP_Locks.cpp translation unit cannot initialize the
 * instances in other translation units, we must update this CUDA global
 * variable based on the Host global variable prior to running any kernels that
 * will use it.  That is the purpose of the
 * KOKKOS_ENSURE_HIP_LOCK_ARRAYS_ON_DEVICE macro.
 */
__device__
#ifdef DESUL_HIP_RDC
    __constant__ extern
#endif
    int32_t* HIP_SPACE_ATOMIC_LOCKS_DEVICE;

__device__
#ifdef DESUL_HIP_RDC
    __constant__ extern
#endif
    int32_t* HIP_SPACE_ATOMIC_LOCKS_NODE;

#define HIP_SPACE_ATOMIC_MASK 0x1FFFF

/// \brief Acquire a lock for the address
///
/// This function tries to acquire the lock for the hash value derived
/// from the provided ptr. If the lock is successfully acquired the
/// function returns true. Otherwise it returns false.
__device__ inline bool lock_address_hip(void* ptr, desul::MemoryScopeDevice) {
  size_t offset = size_t(ptr);
  offset = offset >> 2;
  offset = offset & HIP_SPACE_ATOMIC_MASK;
  return (0 == atomicExch(&desul::Impl::HIP_SPACE_ATOMIC_LOCKS_DEVICE[offset], 1));
}

__device__ inline bool lock_address_hip(void* ptr, desul::MemoryScopeNode) {
  size_t offset = size_t(ptr);
  offset = offset >> 2;
  offset = offset & HIP_SPACE_ATOMIC_MASK;
  return (0 == atomicExch(&desul::Impl::HIP_SPACE_ATOMIC_LOCKS_NODE[offset], 1));
}

/**
 * \brief Release lock for the address
 *
 * This function releases the lock for the hash value derived from the provided
 * ptr. This function should only be called after previously successfully
 * acquiring a lock with lock_address.
 */
__device__ inline void unlock_address_hip(void* ptr, desul::MemoryScopeDevice) {
  size_t offset = size_t(ptr);
  offset = offset >> 2;
  offset = offset & HIP_SPACE_ATOMIC_MASK;
  atomicExch(&desul::Impl::HIP_SPACE_ATOMIC_LOCKS_DEVICE[offset], 0);
}

__device__ inline void unlock_address_hip(void* ptr, desul::MemoryScopeNode) {
  size_t offset = size_t(ptr);
  offset = offset >> 2;
  offset = offset & HIP_SPACE_ATOMIC_MASK;
  atomicExch(&desul::Impl::HIP_SPACE_ATOMIC_LOCKS_NODE[offset], 0);
}
#endif
}  // namespace Impl
}  // namespace desul

// Make lock_array_copied an explicit translation unit scope thing
namespace desul {
namespace Impl {
namespace {
static int lock_array_copied = 0;
inline int eliminate_warning_for_lock_array() { return lock_array_copied; }
}  // namespace
}  // namespace Impl
}  // namespace desul

/* It is critical that this code be a macro, so that it will
   capture the right address for g_device_hip_lock_arrays!
   putting this in an inline function will NOT do the right thing! */
#define DESUL_IMPL_COPY_HIP_LOCK_ARRAYS_TO_DEVICE()                                   \
  {                                                                                   \
    if (::desul::Impl::lock_array_copied == 0) {                                      \
      (void)hipMemcpyToSymbol(                                                        \
          HIP_SYMBOL(::desul::Impl::HIP_SPACE_ATOMIC_LOCKS_DEVICE),                   \
          &::desul::Impl::HIP_SPACE_ATOMIC_LOCKS_DEVICE_h,                            \
          sizeof(int32_t*));                                                          \
      (void)hipMemcpyToSymbol(HIP_SYMBOL(::desul::Impl::HIP_SPACE_ATOMIC_LOCKS_NODE), \
                              &::desul::Impl::HIP_SPACE_ATOMIC_LOCKS_NODE_h,          \
                              sizeof(int32_t*));                                      \
    }                                                                                 \
    ::desul::Impl::lock_array_copied = 1;                                             \
  }

#endif

#if defined(DESUL_HIP_RDC) || (!defined(__HIPCC__))
#define DESUL_ENSURE_HIP_LOCK_ARRAYS_ON_DEVICE()
#else
#define DESUL_ENSURE_HIP_LOCK_ARRAYS_ON_DEVICE() \
  DESUL_IMPL_COPY_HIP_LOCK_ARRAYS_TO_DEVICE()
#endif

#endif
