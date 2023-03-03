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
