/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_LOCK_ARRAY_CUDA_HPP_
#define DESUL_ATOMICS_LOCK_ARRAY_CUDA_HPP_

#include <cstdint>

#include "desul/atomics/Common.hpp"
#include "desul/atomics/Macros.hpp"

namespace desul {
namespace Impl {

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
/// That is the purpose of the ensure_cuda_lock_arrays_on_device function.
#ifdef DESUL_ATOMICS_ENABLE_CUDA_SEPARABLE_COMPILATION
extern
#endif
    __device__ __constant__ int32_t* CUDA_SPACE_ATOMIC_LOCKS_DEVICE;

#ifdef DESUL_ATOMICS_ENABLE_CUDA_SEPARABLE_COMPILATION
extern
#endif
    __device__ __constant__ int32_t* CUDA_SPACE_ATOMIC_LOCKS_NODE;

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

#ifdef DESUL_ATOMICS_ENABLE_CUDA_SEPARABLE_COMPILATION
inline
#else
inline static
#endif
    void
    copy_cuda_lock_arrays_to_device() {
  static bool once = []() {
    cudaMemcpyToSymbol(CUDA_SPACE_ATOMIC_LOCKS_DEVICE,
                       &CUDA_SPACE_ATOMIC_LOCKS_DEVICE_h,
                       sizeof(int32_t*));
    cudaMemcpyToSymbol(CUDA_SPACE_ATOMIC_LOCKS_NODE,
                       &CUDA_SPACE_ATOMIC_LOCKS_NODE_h,
                       sizeof(int32_t*));
    return true;
  }();
  (void)once;
}

}  // namespace Impl
}  // namespace desul

namespace desul {

#ifdef DESUL_ATOMICS_ENABLE_CUDA_SEPARABLE_COMPILATION
inline void ensure_cuda_lock_arrays_on_device() {}
#else
static inline void ensure_cuda_lock_arrays_on_device() {
  Impl::copy_cuda_lock_arrays_to_device();
}
#endif

}  // namespace desul

#endif /* #ifndef DESUL_ATOMICS_LOCK_ARRAY_CUDA_HPP_ */
