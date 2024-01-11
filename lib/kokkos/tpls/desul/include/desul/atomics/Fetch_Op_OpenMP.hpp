/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/
#ifndef DESUL_ATOMICS_FETCH_OP_OPENMP_HPP_
#define DESUL_ATOMICS_FETCH_OP_OPENMP_HPP_

#include <desul/atomics/Common.hpp>
#include <desul/atomics/openmp/OpenMP_40.hpp>

#if 0  // FIXME_OPENMP
namespace desul {
namespace Impl {

// clang-format off
//<editor-fold desc="atomic_fetch_{add,sub,and,or,xor}">
template <class T>
std::enable_if_t<std::is_arithmetic<T>::value, T> host_atomic_fetch_add(
T* ptr, T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T tmp;
#pragma omp atomic capture
  { tmp = *ptr; *ptr += val; }
  return tmp;
}

template <class T>
std::enable_if_t<std::is_arithmetic<T>::value, T> host_atomic_fetch_sub(
T* ptr, T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T tmp;
#pragma omp atomic capture
  { tmp = *ptr; *ptr -= val; }
  return tmp;
}

template <class T>
std::enable_if_t<std::is_arithmetic<T>::value, T> host_atomic_fetch_and(
T* ptr, T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T tmp;
#pragma omp atomic capture
  { tmp = *ptr; *ptr &= val; }
  return tmp;
}

template <class T>
std::enable_if_t<std::is_arithmetic<T>::value, T> host_atomic_fetch_or(
T* ptr, T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T tmp;
#pragma omp atomic capture
  { tmp = *ptr; *ptr |= val; }
  return tmp;
}

template <class T>
std::enable_if_t<std::is_arithmetic<T>::value, T> host_atomic_fetch_xor(
T* ptr, T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T tmp;
#pragma omp atomic capture
  { tmp = *ptr; *ptr ^= val; }
  return tmp;
}
//</editor-fold>

//<editor-fold desc="atomic_{add,sub,and,or,xor}_fetch">
template <class T>
std::enable_if_t<std::is_arithmetic<T>::value, T> host_atomic_add_fetch(
    T* ptr, T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T tmp;
#pragma omp atomic capture
  { *ptr += val; tmp = *ptr; }
  return tmp;
}

template <class T>
std::enable_if_t<std::is_arithmetic<T>::value, T> host_atomic_sub_fetch(
    T* ptr, T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T tmp;
#pragma omp atomic capture
  { *ptr -= val; tmp = *ptr; }
  return tmp;
}

template <class T>
std::enable_if_t<std::is_arithmetic<T>::value, T> host_atomic_and_fetch(
    T* ptr, T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T tmp;
#pragma omp atomic capture
  { *ptr &= val; tmp = *ptr; }
  return tmp;
}

template <class T>
std::enable_if_t<std::is_arithmetic<T>::value, T> host_atomic_or_fetch(
    T* ptr, T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T tmp;
#pragma omp atomic capture
  { *ptr |= val; tmp = *ptr; }
  return tmp;
}

template <class T>
std::enable_if_t<std::is_arithmetic<T>::value, T> host_atomic_xor_fetch(
    T* ptr, T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  T tmp;
#pragma omp atomic capture
  { *ptr ^= val; tmp = *ptr; }
  return tmp;
}
//</editor-fold>
// clang-format on

#define DESUL_IMPL_OPENMP_HOST_ATOMIC_FETCH_OP_ARITHMETIC(OP, MEMORY_ORDER)   \
  template <class T>                                                          \
  std::enable_if_t<std::is_arithmetic<T>::value, T> host_atomic_fetch_##OP(   \
      T* ptr, T val, MemoryOrderRelaxed, MemoryScopeCore) {                   \
    return host_atomic_fetch_##OP(                                            \
        ptr, val, MemoryOrderRelaxed(), MemoryScopeDevice());                 \
  }                                                                           \
  template <class T>                                                          \
  std::enable_if_t<std::is_arithmetic<T>::value, T> host_atomic_##OP##_fetch( \
      T* ptr, T val, MemoryOrderRelaxed, MemoryScopeCore) {                   \
    return host_atomic_##OP##_fetch(                                          \
        ptr, val, MemoryOrderRelaxed(), MemoryScopeDevice());                 \
  }

}  // namespace Impl
}  // namespace desul
#endif

#endif
