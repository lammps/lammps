/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_COMPARE_EXCHANGE_HIP_HPP_
#define DESUL_ATOMICS_COMPARE_EXCHANGE_HIP_HPP_

#include <desul/atomics/Adapt_HIP.hpp>
#include <desul/atomics/Common.hpp>
#include <desul/atomics/Lock_Array_HIP.hpp>
#include <desul/atomics/Thread_Fence_HIP.hpp>
#include <type_traits>

namespace desul {
namespace Impl {

template <class T>
struct atomic_exchange_available_hip {
  constexpr static bool value =
      ((sizeof(T) == 1 && alignof(T) == 1) || (sizeof(T) == 4 && alignof(T) == 4) ||
       (sizeof(T) == 8 && alignof(T) == 8)) &&
      std::is_trivially_copyable<T>::value;
};

template <class T, class MemoryOrder, class MemoryScope>
__device__ std::enable_if_t<atomic_exchange_available_hip<T>::value, T>
device_atomic_compare_exchange(
    T* const dest, T compare, T value, MemoryOrder, MemoryScope) {
  (void)__hip_atomic_compare_exchange_strong(
      dest,
      &compare,
      value,
      HIPMemoryOrder<MemoryOrder>::value,
      HIPMemoryOrder<cmpexch_failure_memory_order<MemoryOrder>>::value,
      HIPMemoryScope<MemoryScope>::value);
  return compare;
}

template <class T, class MemoryOrder, class MemoryScope>
__device__ std::enable_if_t<atomic_exchange_available_hip<T>::value, T>
device_atomic_exchange(T* const dest, T value, MemoryOrder, MemoryScope) {
  T return_val = __hip_atomic_exchange(dest,
                                       value,
                                       HIPMemoryOrder<MemoryOrder>::value,
                                       HIPMemoryScope<MemoryScope>::value);
  return return_val;
}

template <class T, class MemoryOrder, class MemoryScope>
__device__ std::enable_if_t<!atomic_exchange_available_hip<T>::value, T>
device_atomic_compare_exchange(
    T* const dest, T compare, T value, MemoryOrder, MemoryScope scope) {
  // This is a way to avoid deadlock in a warp or wave front
  T return_val;
  int done = 0;
  unsigned long long int active = __ballot(1);
  unsigned long long int done_active = 0;
  while (active != done_active) {
    if (!done) {
      if (lock_address_hip((void*)dest, scope)) {
        if (std::is_same<MemoryOrder, MemoryOrderSeqCst>::value)
          atomic_thread_fence(MemoryOrderRelease(), scope);
        atomic_thread_fence(MemoryOrderAcquire(), scope);
        return_val = *dest;
        if (return_val == compare) {
          *dest = value;
          device_atomic_thread_fence(MemoryOrderRelease(), scope);
        }
        unlock_address_hip((void*)dest, scope);
        done = 1;
      }
    }
    done_active = __ballot(done);
  }
  return return_val;
}

template <class T, class MemoryOrder, class MemoryScope>
__device__ std::enable_if_t<!atomic_exchange_available_hip<T>::value, T>
device_atomic_exchange(T* const dest, T value, MemoryOrder, MemoryScope scope) {
  // This is a way to avoid deadlock in a warp or wave front
  T return_val;
  int done = 0;
  unsigned long long int active = __ballot(1);
  unsigned long long int done_active = 0;
  while (active != done_active) {
    if (!done) {
      if (lock_address_hip((void*)dest, scope)) {
        if (std::is_same<MemoryOrder, MemoryOrderSeqCst>::value)
          atomic_thread_fence(MemoryOrderRelease(), scope);
        device_atomic_thread_fence(MemoryOrderAcquire(), scope);
        return_val = *dest;
        *dest = value;
        device_atomic_thread_fence(MemoryOrderRelease(), scope);
        unlock_address_hip((void*)dest, scope);
        done = 1;
      }
    }
    done_active = __ballot(done);
  }
  return return_val;
}

}  // namespace Impl
}  // namespace desul

#endif
