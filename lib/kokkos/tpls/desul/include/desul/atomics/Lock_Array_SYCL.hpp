/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_LOCK_ARRAY_SYCL_HPP_
#define DESUL_ATOMICS_LOCK_ARRAY_SYCL_HPP_

#include <cstdint>

#include "desul/atomics/Adapt_SYCL.hpp"
#include "desul/atomics/Common.hpp"
#include "desul/atomics/Macros.hpp"

// FIXME_SYCL
#if __has_include(<sycl/sycl.hpp>)
#include <sycl/sycl.hpp>
#else
#include <CL/sycl.hpp>
#endif

namespace desul {
namespace Impl {

// FIXME_SYCL Use SYCL_EXT_ONEAPI_DEVICE_GLOBAL when available instead
#ifdef DESUL_SYCL_DEVICE_GLOBAL_SUPPORTED

/**
 * \brief This global variable in Host space is the central definition of these
 * arrays.
 */
extern int32_t* SYCL_SPACE_ATOMIC_LOCKS_DEVICE_h;
extern int32_t* SYCL_SPACE_ATOMIC_LOCKS_NODE_h;

/// \brief After this call, the lock arrays used in [un]lock_address_sycl
///        are initialized and ready to be used.
///
/// This call is idempotent.
/// The function is templated to make it a weak symbol to deal with Kokkos/RAJA
///   snapshotted version while also linking against pure Desul
template <typename /*AlwaysInt*/ = int>
void init_lock_arrays_sycl(sycl::queue q);

/// \brief After this call, the lock arrays used in [un]lock_address_sycl
///        are freed and can't be used anymore.
///
/// This call is idempotent.
/// The function is templated to make it a weak symbol to deal with Kokkos/RAJA
///   snapshotted version while also linking against pure Desul
template <typename /*AlwaysInt*/ = int>
void finalize_lock_arrays_sycl(sycl::queue q);

/**
 * \brief This global variable in SYCL space is what kernels use to get access
 * to the lock arrays.
 *
 * There is only one single instance of this global variable for the entire
 * executable, whose definition will be in Kokkos_SYCL_Locks.cpp (and whose
 * declaration here must be extern). This one instance will be initialized
 * by initialize_host_sycl_lock_arrays and need not be modified afterwards.
 */
SYCL_EXTERNAL extern sycl_device_global<int32_t*> SYCL_SPACE_ATOMIC_LOCKS_DEVICE;

SYCL_EXTERNAL extern sycl_device_global<int32_t*> SYCL_SPACE_ATOMIC_LOCKS_NODE;

#define SYCL_SPACE_ATOMIC_MASK 0x1FFFF

/// \brief Acquire a lock for the address
///
/// This function tries to acquire the lock for the hash value derived
/// from the provided ptr. If the lock is successfully acquired the
/// function returns true. Otherwise it returns false.
inline bool lock_address_sycl(void* ptr, MemoryScopeDevice) {
  size_t offset = size_t(ptr);
  offset = offset >> 2;
  offset = offset & SYCL_SPACE_ATOMIC_MASK;
  sycl::atomic_ref<int32_t,
                   sycl::memory_order::relaxed,
                   sycl::memory_scope::device,
                   sycl::access::address_space::global_space>
      lock_device_ref(SYCL_SPACE_ATOMIC_LOCKS_DEVICE[offset]);
  return (0 == lock_device_ref.exchange(1));
}

inline bool lock_address_sycl(void* ptr, MemoryScopeNode) {
  size_t offset = size_t(ptr);
  offset = offset >> 2;
  offset = offset & SYCL_SPACE_ATOMIC_MASK;
  sycl::atomic_ref<int32_t,
                   sycl::memory_order::relaxed,
                   sycl::memory_scope::system,
                   sycl::access::address_space::global_space>
      lock_node_ref(SYCL_SPACE_ATOMIC_LOCKS_NODE[offset]);
  return (0 == lock_node_ref.exchange(1));
}

/**
 * \brief Release lock for the address
 *
 * This function releases the lock for the hash value derived from the provided
 * ptr. This function should only be called after previously successfully
 * acquiring a lock with lock_address.
 */
inline void unlock_address_sycl(void* ptr, MemoryScopeDevice) {
  size_t offset = size_t(ptr);
  offset = offset >> 2;
  offset = offset & SYCL_SPACE_ATOMIC_MASK;
  sycl::atomic_ref<int32_t,
                   sycl::memory_order::relaxed,
                   sycl::memory_scope::device,
                   sycl::access::address_space::global_space>
      lock_device_ref(SYCL_SPACE_ATOMIC_LOCKS_DEVICE[offset]);
  lock_device_ref.exchange(0);
}

inline void unlock_address_sycl(void* ptr, MemoryScopeNode) {
  size_t offset = size_t(ptr);
  offset = offset >> 2;
  offset = offset & SYCL_SPACE_ATOMIC_MASK;
  sycl::atomic_ref<int32_t,
                   sycl::memory_order::relaxed,
                   sycl::memory_scope::system,
                   sycl::access::address_space::global_space>
      lock_node_ref(SYCL_SPACE_ATOMIC_LOCKS_NODE[offset]);
  lock_node_ref.exchange(0);
}

#else  // not supported

template <typename /*AlwaysInt*/ = int>
void init_lock_arrays_sycl(sycl::queue) {
  assert(false);
}

template <typename /*AlwaysInt*/ = int>
void finalize_lock_arrays_sycl(sycl::queue) {
  assert(false);
}

inline bool lock_address_sycl(void*, MemoryScopeDevice) {
  assert(false);
  // return true so that the CAS loops don't hang.
  return true;
}

inline bool lock_address_sycl(void*, MemoryScopeNode) {
  assert(false);
  // return true so that the CAS loops don't hang.
  return true;
}

inline void unlock_address_sycl(void*, MemoryScopeDevice) { assert(false); }

inline void unlock_address_sycl(void*, MemoryScopeNode) { assert(false); }
#endif
}  // namespace Impl
}  // namespace desul
#endif
