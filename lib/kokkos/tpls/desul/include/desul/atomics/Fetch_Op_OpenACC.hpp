/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/
#ifndef DESUL_ATOMICS_FETCH_OP_OPENACC_HPP_
#define DESUL_ATOMICS_FETCH_OP_OPENACC_HPP_

#include <algorithm>  // min, max
#include <desul/atomics/Common.hpp>
#include <type_traits>

namespace desul {
namespace Impl {

#ifdef __NVCOMPILER

template <class T>
inline constexpr bool is_openacc_integral_type_v =
    std::is_same_v<T, int> || std::is_same_v<T, unsigned int> ||
    std::is_same_v<T, unsigned long long>;

template <class T>
inline constexpr bool is_openacc_arithmetic_type_v = std::is_same_v<T, float> ||
#ifndef DESUL_CUDA_ARCH_IS_PRE_PASCAL
                                                     std::is_same_v<T, double> ||
#endif
                                                     is_openacc_integral_type_v<T>;

#else

template <class T>
inline constexpr bool is_openacc_integral_type_v = std::is_integral_v<T>;

template <class T>
inline constexpr bool is_openacc_arithmetic_type_v = std::is_arithmetic_v<T>;

#endif

//<editor-fold
// desc="device_atomic_fetch_{add,sub,mul,div,lshift,rshift,mod,max,min,and,or,xor}">
#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_arithmetic_type_v<T>, T> device_atomic_fetch_add(
    T* ptr, const T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T old;
#pragma acc atomic capture
  {
    old = *ptr;
    *ptr += val;
  }
  return old;
}

#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_arithmetic_type_v<T>, T> device_atomic_fetch_inc(
    T* ptr, MemoryOrderRelaxed, MemoryScopeDevice) {
  T old;
#pragma acc atomic capture
  {
    old = *ptr;
    *ptr += T(1);
  }
  return old;
}

#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_arithmetic_type_v<T>, T> device_atomic_fetch_sub(
    T* ptr, const T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T old;
#pragma acc atomic capture
  {
    old = *ptr;
    *ptr -= val;
  }
  return old;
}

#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_arithmetic_type_v<T>, T> device_atomic_fetch_dec(
    T* ptr, MemoryOrderRelaxed, MemoryScopeDevice) {
  T old;
#pragma acc atomic capture
  {
    old = *ptr;
    *ptr -= T(1);
  }
  return old;
}

#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_arithmetic_type_v<T>, T> device_atomic_fetch_mul(
    T* ptr, const T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T old;
#pragma acc atomic capture
  {
    old = *ptr;
    *ptr *= val;
  }
  return old;
}

#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_arithmetic_type_v<T>, T> device_atomic_fetch_div(
    T* ptr, const T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T old;
#pragma acc atomic capture
  {
    old = *ptr;
    *ptr /= val;
  }
  return old;
}

#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_integral_type_v<T>, T> device_atomic_fetch_lshift(
    T* ptr, const unsigned int val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T old;
#pragma acc atomic capture
  {
    old = *ptr;
    *ptr = *ptr << val;
  }
  return old;
}

#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_integral_type_v<T>, T> device_atomic_fetch_rshift(
    T* ptr, const unsigned int val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T old;
#pragma acc atomic capture
  {
    old = *ptr;
    *ptr = *ptr >> val;
  }
  return old;
}

#ifdef __NVCOMPILER
#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_arithmetic_type_v<T>, T> device_atomic_fetch_max(
    T* ptr, const T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T old;
  old = atomicMax(ptr, val);
  return old;
}
#endif

#ifdef __NVCOMPILER
#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_arithmetic_type_v<T>, T> device_atomic_fetch_min(
    T* ptr, const T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  int old;
  old = atomicMin(ptr, val);
  return old;
}
#endif

#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_integral_type_v<T>, T> device_atomic_fetch_and(
    T* ptr, const T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T old;
#pragma acc atomic capture
  {
    old = *ptr;
    *ptr &= val;
  }
  return old;
}

#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_integral_type_v<T>, T> device_atomic_fetch_or(
    T* ptr, const T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T old;
#pragma acc atomic capture
  {
    old = *ptr;
    *ptr |= val;
  }
  return old;
}

#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_integral_type_v<T>, T> device_atomic_fetch_xor(
    T* ptr, const T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T old;
#pragma acc atomic capture
  {
    old = *ptr;
    *ptr ^= val;
  }
  return old;
}
//</editor-fold>

//<editor-fold
// desc="device_atomic_{add,sub,mul,div,lshift,rshift,mod,max,min,and,or,xor}_fetch">
#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_arithmetic_type_v<T>, T> device_atomic_add_fetch(
    T* ptr, const T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T tmp;
#pragma acc atomic capture
  {
    *ptr += val;
    tmp = *ptr;
  }
  return tmp;
}

#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_arithmetic_type_v<T>, T> device_atomic_inc_fetch(
    T* ptr, MemoryOrderRelaxed, MemoryScopeDevice) {
  T tmp;
#pragma acc atomic capture
  {
    *ptr += T(1);
    tmp = *ptr;
  }
  return tmp;
}

#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_arithmetic_type_v<T>, T> device_atomic_sub_fetch(
    T* ptr, const T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T tmp;
#pragma acc atomic capture
  {
    *ptr -= val;
    tmp = *ptr;
  }
  return tmp;
}

#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_arithmetic_type_v<T>, T> device_atomic_dec_fetch(
    T* ptr, MemoryOrderRelaxed, MemoryScopeDevice) {
  T tmp;
#pragma acc atomic capture
  {
    *ptr -= T(1);
    tmp = *ptr;
  }
  return tmp;
}

#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_arithmetic_type_v<T>, T> device_atomic_mul_fetch(
    T* ptr, const T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T tmp;
#pragma acc atomic capture
  {
    *ptr *= val;
    tmp = *ptr;
  }
  return tmp;
}

#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_arithmetic_type_v<T>, T> device_atomic_div_fetch(
    T* ptr, const T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T tmp;
#pragma acc atomic capture
  {
    *ptr /= val;
    tmp = *ptr;
  }
  return tmp;
}

#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_integral_type_v<T>, T> device_atomic_lshift_fetch(
    T* ptr, const unsigned int val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T tmp;
#pragma acc atomic capture
  {
    *ptr = *ptr << val;
    tmp = *ptr;
  }
  return tmp;
}

#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_integral_type_v<T>, T> device_atomic_rshift_fetch(
    T* ptr, const unsigned int val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T tmp;
#pragma acc atomic capture
  {
    *ptr = *ptr >> val;
    tmp = *ptr;
  }
  return tmp;
}

#ifdef __NVCOMPILER
#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_arithmetic_type_v<T>, T> device_atomic_max_fetch(
    T* ptr, const T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T tmp;
  tmp = atomicMax(ptr, val);
  tmp = std::max(tmp, val);
  return tmp;
}
#endif

#ifdef __NVCOMPILER
#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_arithmetic_type_v<T>, T> device_atomic_min_fetch(
    T* ptr, const T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T tmp;
  tmp = atomicMin(ptr, val);
  tmp = std::min(tmp, val);
  return tmp;
}
#endif

#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_integral_type_v<T>, T> device_atomic_and_fetch(
    T* ptr, const T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T tmp;
#pragma acc atomic capture
  {
    *ptr &= val;
    tmp = *ptr;
  }
  return tmp;
}

#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_integral_type_v<T>, T> device_atomic_or_fetch(
    T* ptr, const T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T tmp;
#pragma acc atomic capture
  {
    *ptr |= val;
    tmp = *ptr;
  }
  return tmp;
}

#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_integral_type_v<T>, T> device_atomic_xor_fetch(
    T* ptr, const T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T tmp;
#pragma acc atomic capture
  {
    *ptr ^= val;
    tmp = *ptr;
  }
  return tmp;
}
//</editor-fold>

//<editor-fold desc="device_atomic_{store,load}">
#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_arithmetic_type_v<T>, void> device_atomic_store(
    T* const ptr, const T val, MemoryOrderRelaxed, MemoryScopeDevice) {
#pragma acc atomic write
  *ptr = val;
}

#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_arithmetic_type_v<T>, void> device_atomic_store(
    T* const ptr, const T val, MemoryOrderRelease, MemoryScopeDevice) {
  if (acc_on_device(acc_device_not_host)) {
    printf(
        "DESUL error in device_atomic_store(MemoryOrderRelease): Not supported atomic "
        "operation in the OpenACC backend\n");
  }
#pragma acc atomic write
  *ptr = val;
}

#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_arithmetic_type_v<T>, T> device_atomic_load(
    const T* const ptr, MemoryOrderRelaxed, MemoryScopeDevice) {
  T retval;
#pragma acc atomic read
  retval = *ptr;
  return retval;
}

#pragma acc routine seq
template <class T>
std::enable_if_t<is_openacc_arithmetic_type_v<T>, T> device_atomic_load(
    const T* const ptr, MemoryOrderAcquire, MemoryScopeDevice) {
  if (acc_on_device(acc_device_not_host)) {
    printf(
        "DESUL error in device_atomic_load(MemoryOrderAcquire): Not supported atomic "
        "operation in the OpenACC backend\n");
  }
  T retval;
#pragma acc atomic read
  retval = *ptr;
  return retval;
}
//</editor-fold>

}  // namespace Impl
}  // namespace desul

#endif
