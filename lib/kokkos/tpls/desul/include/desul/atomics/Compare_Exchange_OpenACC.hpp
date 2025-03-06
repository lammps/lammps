/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_COMPARE_EXCHANGE_OPENACC_HPP_
#define DESUL_ATOMICS_COMPARE_EXCHANGE_OPENACC_HPP_

#include <openacc.h>

#include <desul/atomics/Common.hpp>
#include <desul/atomics/Thread_Fence_OpenACC.hpp>
#include <type_traits>

namespace desul {
namespace Impl {

#ifdef __NVCOMPILER

#pragma acc routine seq
template <class T, class MemoryOrder, class MemoryScope>
T device_atomic_exchange(T* dest, T value, MemoryOrder, MemoryScope /*scope*/) {
  if constexpr (std::is_arithmetic_v<T> && ((sizeof(T) == 4) || (sizeof(T) == 8))) {
    T return_val;
#pragma acc atomic capture
    {
      return_val = *dest;
      *dest = value;
    }
    return return_val;
  } else {
    // FIXME_OPENACC
    if (acc_on_device(acc_device_not_host)) {
      printf(
          "DESUL error in device_atomic_exchange(): Not supported atomic operation in "
          "the OpenACC backend\n");
    }
    //  Acquire a lock for the address
    // while (!lock_address_openacc((void*)dest, scope)) {
    // }
    // device_atomic_thread_fence(MemoryOrderAcquire(), scope);
    T return_val = *dest;
    *dest = value;
    // device_atomic_thread_fence(MemoryOrderRelease(), scope);
    // unlock_address_openacc((void*)dest, scope);
    return return_val;
  }
}

#pragma acc routine seq
template <class T, class MemoryOrder, class MemoryScope>
T device_atomic_compare_exchange(
    T* dest, T compare, T value, MemoryOrder, MemoryScope scope) {
  // Floating point types treated separetely to work around compiler errors
  // "parse invalid cast opcode for cast from 'i32' to 'float'".
  // Also not just "forwarding" arguments to atomicCAS because it does not have an
  // overload that takes int64_t
  if constexpr (std::is_integral_v<T> && ((sizeof(T) == 4) || (sizeof(T) == 8))) {
    static_assert(sizeof(unsigned int) == 4);
    static_assert(sizeof(unsigned long long int) == 8);
    using cas_t =
        std::conditional_t<(sizeof(T) == 4), unsigned int, unsigned long long int>;
    cas_t return_val = atomicCAS(reinterpret_cast<cas_t*>(dest),
                                 reinterpret_cast<cas_t&>(compare),
                                 reinterpret_cast<cas_t&>(value));
    return reinterpret_cast<T&>(return_val);
#ifdef DESUL_CUDA_ARCH_IS_PRE_PASCAL
  } else if constexpr (std::is_same_v<T, float>) {
#else
  } else if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>) {
#endif
    return atomicCAS(dest, compare, value);
  } else {
    // FIXME_OPENACC
    if (acc_on_device(acc_device_not_host)) {
      printf(
          "DESUL error in device_atomic_compare_exchange(): Not supported atomic "
          "operation in the OpenACC backend\n");
    }
    T current_val = *dest;
    // Acquire a lock for the address
    // while (!lock_address_openacc((void*)dest, scope)) {
    //}
    // device_atomic_thread_fence(MemoryOrderAcquire(), scope);
    if (current_val == compare) {
      *dest = value;
      // device_atomic_thread_fence(MemoryOrderRelease(), scope);
    }
    // unlock_address_openacc((void*)dest, scope);
    return current_val;
  }
}

#else  // not NVHPC

#pragma acc routine seq
template <class T, class MemoryOrder, class MemoryScope>
T device_atomic_exchange(T* dest, T value, MemoryOrder, MemoryScope) {
  if constexpr (std::is_arithmetic_v<T>) {
    T return_val;
#pragma acc atomic capture
    {
      return_val = *dest;
      *dest = value;
    }
    return return_val;
  } else {
    // FIXME_OPENACC
    printf(
        "DESUL error in device_atomic_exchange(): Not supported atomic operation in "
        "the OpenACC backend\n");
    //  Acquire a lock for the address
    // while (!lock_address_openacc((void*)dest, scope)) {
    // }
    // device_atomic_thread_fence(MemoryOrderAcquire(), scope);
    T return_val = *dest;
    *dest = value;
    // device_atomic_thread_fence(MemoryOrderRelease(), scope);
    // unlock_address_openacc((void*)dest, scope);
    return return_val;
  }
}

#pragma acc routine seq
template <class T, class MemoryOrder, class MemoryScope>
T device_atomic_compare_exchange(
    T* dest, T compare, T value, MemoryOrder, MemoryScope scope) {
  // FIXME_OPENACC
  printf(
      "DESUL error in device_atomic_compare_exchange(): Not supported atomic operation "
      "in the OpenACC backend\n");
  T current_val = *dest;
  // Acquire a lock for the address
  // while (!lock_address_openacc((void*)dest, scope)) {
  //}
  // device_atomic_thread_fence(MemoryOrderAcquire(), scope);
  if (current_val == compare) {
    *dest = value;
    // device_atomic_thread_fence(MemoryOrderRelease(), scope);
  }
  // unlock_address_openacc((void*)dest, scope);
  return current_val;
}

#endif

}  // namespace Impl
}  // namespace desul

#endif
