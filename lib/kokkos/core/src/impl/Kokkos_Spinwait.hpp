/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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
