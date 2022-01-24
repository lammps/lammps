/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//              Copyright (2019) Sandia Corporation
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

#ifndef KOKKOS_KOKKOS_ATOMIC_MEMORY_ORDER_HPP
#define KOKKOS_KOKKOS_ATOMIC_MEMORY_ORDER_HPP

#include <Kokkos_Macros.hpp>

#include <atomic>

namespace Kokkos {
namespace Impl {

/** @file
 * Provides strongly-typed analogs of the standard memory order enumerators.
 * In addition to (very slightly) reducing the constant propagation burden on
 * the compiler, this allows us to give compile-time errors for things that
 * don't make sense, like atomic_load with memory order release.
 */

struct memory_order_seq_cst_t {
  using memory_order = memory_order_seq_cst_t;
#if defined(KOKKOS_ENABLE_GNU_ATOMICS) ||   \
    defined(KOKKOS_ENABLE_INTEL_ATOMICS) || \
    defined(KOKKOS_ENABLE_CUDA_ASM_ATOMICS)
  static constexpr auto gnu_constant = __ATOMIC_SEQ_CST;
#endif
  static constexpr auto std_constant = std::memory_order_seq_cst;
};
constexpr memory_order_seq_cst_t memory_order_seq_cst = {};

struct memory_order_relaxed_t {
  using memory_order = memory_order_relaxed_t;
#if defined(KOKKOS_ENABLE_GNU_ATOMICS) ||   \
    defined(KOKKOS_ENABLE_INTEL_ATOMICS) || \
    defined(KOKKOS_ENABLE_CUDA_ASM_ATOMICS)
  static constexpr auto gnu_constant = __ATOMIC_RELAXED;
#endif
  static constexpr auto std_constant = std::memory_order_relaxed;
};
constexpr memory_order_relaxed_t memory_order_relaxed = {};

struct memory_order_acquire_t {
  using memory_order = memory_order_acquire_t;
#if defined(KOKKOS_ENABLE_GNU_ATOMICS) ||   \
    defined(KOKKOS_ENABLE_INTEL_ATOMICS) || \
    defined(KOKKOS_ENABLE_CUDA_ASM_ATOMICS)
  static constexpr auto gnu_constant = __ATOMIC_ACQUIRE;
#endif
  static constexpr auto std_constant = std::memory_order_acquire;
};
constexpr memory_order_acquire_t memory_order_acquire = {};

struct memory_order_release_t {
  using memory_order = memory_order_release_t;
#if defined(KOKKOS_ENABLE_GNU_ATOMICS) ||   \
    defined(KOKKOS_ENABLE_INTEL_ATOMICS) || \
    defined(KOKKOS_ENABLE_CUDA_ASM_ATOMICS)
  static constexpr auto gnu_constant = __ATOMIC_RELEASE;
#endif
  static constexpr auto std_constant = std::memory_order_release;
};
constexpr memory_order_release_t memory_order_release = {};

struct memory_order_acq_rel_t {
  using memory_order = memory_order_acq_rel_t;
#if defined(KOKKOS_ENABLE_GNU_ATOMICS) ||   \
    defined(KOKKOS_ENABLE_INTEL_ATOMICS) || \
    defined(KOKKOS_ENABLE_CUDA_ASM_ATOMICS)
  static constexpr auto gnu_constant = __ATOMIC_ACQ_REL;
#endif
  static constexpr auto std_constant = std::memory_order_acq_rel;
};
constexpr memory_order_acq_rel_t memory_order_acq_rel = {};

// Intentionally omit consume (for now)

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_ATOMIC_MEMORY_ORDER_HPP
