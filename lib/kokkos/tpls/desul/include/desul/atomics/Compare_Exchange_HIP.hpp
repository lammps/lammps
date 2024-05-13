/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_COMPARE_EXCHANGE_HIP_HPP_
#define DESUL_ATOMICS_COMPARE_EXCHANGE_HIP_HPP_

#include <desul/atomics/Common.hpp>
#include <desul/atomics/Lock_Array_HIP.hpp>
#include <desul/atomics/Thread_Fence_HIP.hpp>
#include <type_traits>

namespace desul {
namespace Impl {

template <class T, class MemoryScope>
__device__ std::enable_if_t<sizeof(T) == 4, T> device_atomic_compare_exchange(
    T* const dest, T compare, T value, MemoryOrderRelaxed, MemoryScope) {
  static_assert(sizeof(unsigned int) == 4,
                "this function assumes an unsigned int is 32-bit");
  unsigned int return_val = atomicCAS(reinterpret_cast<unsigned int*>(dest),
                                      reinterpret_cast<unsigned int&>(compare),
                                      reinterpret_cast<unsigned int&>(value));
  return reinterpret_cast<T&>(return_val);
}
template <class T, class MemoryScope>
__device__ std::enable_if_t<sizeof(T) == 8, T> device_atomic_compare_exchange(
    T* const dest, T compare, T value, MemoryOrderRelaxed, MemoryScope) {
  static_assert(sizeof(unsigned long long int) == 8,
                "this function assumes an unsigned long long is 64-bit");
  unsigned long long int return_val =
      atomicCAS(reinterpret_cast<unsigned long long int*>(dest),
                reinterpret_cast<unsigned long long int&>(compare),
                reinterpret_cast<unsigned long long int&>(value));
  return reinterpret_cast<T&>(return_val);
}

template <class T, class MemoryScope>
__device__ std::enable_if_t<sizeof(T) == 4 || sizeof(T) == 8, T>
device_atomic_compare_exchange(
    T* const dest, T compare, T value, MemoryOrderRelease, MemoryScope) {
  T return_val = atomic_compare_exchange(
      dest, compare, value, MemoryOrderRelaxed(), MemoryScope());
  atomic_thread_fence(MemoryOrderRelease(), MemoryScope());
  return return_val;
}

template <class T, class MemoryScope>
__device__ std::enable_if_t<sizeof(T) == 4 || sizeof(T) == 8, T>
device_atomic_compare_exchange(
    T* const dest, T compare, T value, MemoryOrderAcquire, MemoryScope) {
  atomic_thread_fence(MemoryOrderAcquire(), MemoryScope());
  T return_val = atomic_compare_exchange(
      dest, compare, value, MemoryOrderRelaxed(), MemoryScope());
  return return_val;
}

template <class T, class MemoryScope>
__device__ std::enable_if_t<sizeof(T) == 4 || sizeof(T) == 8, T>
device_atomic_compare_exchange(
    T* const dest, T compare, T value, MemoryOrderAcqRel, MemoryScope) {
  atomic_thread_fence(MemoryOrderAcquire(), MemoryScope());
  T return_val = atomic_compare_exchange(
      dest, compare, value, MemoryOrderRelaxed(), MemoryScope());
  atomic_thread_fence(MemoryOrderRelease(), MemoryScope());
  return return_val;
}

template <class T, class MemoryScope>
__device__ std::enable_if_t<sizeof(T) == 4, T> device_atomic_exchange(
    T* const dest, T value, MemoryOrderRelaxed, MemoryScope) {
  static_assert(sizeof(unsigned int) == 4,
                "this function assumes an unsigned int is 32-bit");
  unsigned int return_val = atomicExch(reinterpret_cast<unsigned int*>(dest),
                                       reinterpret_cast<unsigned int&>(value));
  return reinterpret_cast<T&>(return_val);
}
template <class T, class MemoryScope>
__device__ std::enable_if_t<sizeof(T) == 8, T> device_atomic_exchange(
    T* const dest, T value, MemoryOrderRelaxed, MemoryScope) {
  static_assert(sizeof(unsigned long long int) == 8,
                "this function assumes an unsigned long long is 64-bit");
  unsigned long long int return_val =
      atomicExch(reinterpret_cast<unsigned long long int*>(dest),
                 reinterpret_cast<unsigned long long int&>(value));
  return reinterpret_cast<T&>(return_val);
}

template <class T, class MemoryScope>
__device__ std::enable_if_t<sizeof(T) == 4 || sizeof(T) == 8, T> device_atomic_exchange(
    T* const dest, T compare, T value, MemoryOrderRelease, MemoryScope) {
  T return_val = device_atomic_compare_exchange(
      dest, compare, value, MemoryOrderRelaxed(), MemoryScope());
  device_atomic_thread_fence(MemoryOrderRelease(), MemoryScope());
  return reinterpret_cast<T&>(return_val);
}

template <class T, class MemoryScope>
__device__ std::enable_if_t<sizeof(T) == 4 || sizeof(T) == 8, T> device_atomic_exchange(
    T* const dest, T /*compare*/, T value, MemoryOrderAcquire, MemoryScope) {
  device_atomic_thread_fence(MemoryOrderAcquire(), MemoryScope());
  T return_val =
      device_atomic_exchange(dest, value, MemoryOrderRelaxed(), MemoryScope());
  return reinterpret_cast<T&>(return_val);
}

template <class T, class MemoryScope>
__device__ std::enable_if_t<sizeof(T) == 4 || sizeof(T) == 8, T> device_atomic_exchange(
    T* const dest, T value, MemoryOrderAcqRel, MemoryScope) {
  device_atomic_thread_fence(MemoryOrderAcquire(), MemoryScope());
  T return_val =
      device_atomic_exchange(dest, value, MemoryOrderRelaxed(), MemoryScope());
  device_atomic_thread_fence(MemoryOrderRelease(), MemoryScope());
  return reinterpret_cast<T&>(return_val);
}

template <class T, class MemoryScope>
__device__ std::enable_if_t<sizeof(T) == 4 || sizeof(T) == 8, T> device_atomic_exchange(
    T* const dest, T value, MemoryOrderSeqCst, MemoryScope) {
  device_atomic_thread_fence(MemoryOrderAcquire(), MemoryScope());
  T return_val =
      device_atomic_exchange(dest, value, MemoryOrderRelaxed(), MemoryScope());
  device_atomic_thread_fence(MemoryOrderRelease(), MemoryScope());
  return reinterpret_cast<T&>(return_val);
}

template <class T, class MemoryScope>
__device__ std::enable_if_t<sizeof(T) == 4 || sizeof(T) == 8, T>
device_atomic_compare_exchange(
    T* const dest, T compare, T value, MemoryOrderSeqCst, MemoryScope) {
  device_atomic_thread_fence(MemoryOrderAcquire(), MemoryScope());
  T return_val = device_atomic_compare_exchange(
      dest, compare, value, MemoryOrderRelaxed(), MemoryScope());
  device_atomic_thread_fence(MemoryOrderRelease(), MemoryScope());
  return return_val;
}

template <class T, class MemoryOrder, class MemoryScope>
__device__ std::enable_if_t<(sizeof(T) != 8) && (sizeof(T) != 4), T>
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
__device__ std::enable_if_t<(sizeof(T) != 8) && (sizeof(T) != 4), T>
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
