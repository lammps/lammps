//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_SPINWAIT_HPP
#define KOKKOS_SPINWAIT_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Atomic.hpp>

#include <cstdint>

#include <type_traits>

namespace Kokkos {
namespace Impl {

enum class WaitMode : int {
  ACTIVE  // Used for tight loops to keep threads active longest
  ,
  PASSIVE  // Used to quickly yield the thread to quite down the system
  ,
  ROOT  // Never sleep or yield the root thread
};

void host_thread_yield(const uint32_t i, const WaitMode mode);

template <typename T>
std::enable_if_t<std::is_integral<T>::value, void> root_spinwait_while_equal(
    T const volatile& flag, const T value) {
  Kokkos::store_fence();
  uint32_t i = 0;
  while (value == flag) {
    host_thread_yield(++i, WaitMode::ROOT);
  }
  Kokkos::load_fence();
}

template <typename T>
std::enable_if_t<std::is_integral<T>::value, void> root_spinwait_until_equal(
    T const volatile& flag, const T value) {
  Kokkos::store_fence();
  uint32_t i = 0;
  while (value != flag) {
    host_thread_yield(++i, WaitMode::ROOT);
  }
  Kokkos::load_fence();
}

template <typename T>
std::enable_if_t<std::is_integral<T>::value, void> spinwait_while_equal(
    T const volatile& flag, const T value) {
  Kokkos::store_fence();
  uint32_t i = 0;
  while (value == flag) {
    host_thread_yield(++i, WaitMode::ACTIVE);
  }
  Kokkos::load_fence();
}

template <typename T>
std::enable_if_t<std::is_integral<T>::value, void> yield_while_equal(
    T const volatile& flag, const T value) {
  Kokkos::store_fence();
  uint32_t i = 0;
  while (value == flag) {
    host_thread_yield(++i, WaitMode::PASSIVE);
  }
  Kokkos::load_fence();
}

template <typename T>
std::enable_if_t<std::is_integral<T>::value, void> spinwait_until_equal(
    T const volatile& flag, const T value) {
  Kokkos::store_fence();
  uint32_t i = 0;
  while (value != flag) {
    host_thread_yield(++i, WaitMode::ACTIVE);
  }
  Kokkos::load_fence();
}

template <typename T>
std::enable_if_t<std::is_integral<T>::value, void> yield_until_equal(
    T const volatile& flag, const T value) {
  Kokkos::store_fence();
  uint32_t i = 0;
  while (value != flag) {
    host_thread_yield(++i, WaitMode::PASSIVE);
  }
  Kokkos::load_fence();
}

} /* namespace Impl */
} /* namespace Kokkos */

#endif /* #ifndef KOKKOS_SPINWAIT_HPP */
