/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_LOAD_AND_STORE_HIP_HPP_
#define DESUL_ATOMICS_LOAD_AND_STORE_HIP_HPP_

#include <desul/atomics/Adapt_HIP.hpp>
#include <desul/atomics/Lock_Free_Types_HIP.hpp>

namespace desul {
namespace Impl {

template <class T,
          class MemoryOrder,
          class MemoryScope,
          std::enable_if_t<device_atomic_always_lock_free<T>, int> = 0>
__device__ T device_atomic_load(T const* ptr, MemoryOrder, MemoryScope) {
  return __hip_atomic_load(
      ptr, HIPMemoryOrder<MemoryOrder>::value, HIPMemoryScope<MemoryScope>::value);
}

template <class T,
          class MemoryOrder,
          class MemoryScope,
          std::enable_if_t<device_atomic_always_lock_free<T>, int> = 0>
__device__ void device_atomic_store(T* ptr, T val, MemoryOrder, MemoryScope) {
  __hip_atomic_store(
      ptr, val, HIPMemoryOrder<MemoryOrder>::value, HIPMemoryScope<MemoryScope>::value);
}

template <class T,
          class MemoryOrder,
          class MemoryScope,
          std::enable_if_t<!device_atomic_always_lock_free<T>, int> = 0>
__device__ T device_atomic_load(T const* ptr,
                                MemoryOrder /*order*/,
                                MemoryScope scope) {
  // This is a way to avoid deadlock in a warp or wave front
  T ret;
  int done = 0;
  unsigned long long int active = __ballot(1);
  unsigned long long int done_active = 0;
  while (active != done_active) {
    if (!done) {
      if (lock_address_hip((void*)ptr, scope)) {
        device_atomic_thread_fence(MemoryOrderAcquire(), scope);
        ret = *ptr;
        device_atomic_thread_fence(MemoryOrderRelease(), scope);
        unlock_address_hip((void*)ptr, scope);
        done = 1;
      }
    }
    done_active = __ballot(done);
  }
  return ret;
}

template <class T,
          class MemoryOrder,
          class MemoryScope,
          std::enable_if_t<!device_atomic_always_lock_free<T>, int> = 0>
__device__ void device_atomic_store(T* ptr,
                                    T val,
                                    MemoryOrder /*order*/,
                                    MemoryScope scope) {
  // This is a way to avoid deadlock in a warp or wave front
  int done = 0;
  unsigned long long int active = __ballot(1);
  unsigned long long int done_active = 0;
  while (active != done_active) {
    if (!done) {
      if (lock_address_hip((void*)ptr, scope)) {
        device_atomic_thread_fence(MemoryOrderAcquire(), scope);
        *ptr = val;
        device_atomic_thread_fence(MemoryOrderRelease(), scope);
        unlock_address_hip((void*)ptr, scope);
        done = 1;
      }
    }
    done_active = __ballot(done);
  }
}

}  // namespace Impl
}  // namespace desul

#endif
