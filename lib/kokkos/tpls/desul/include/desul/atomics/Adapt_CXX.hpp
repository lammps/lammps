/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_ADAPT_CXX_HPP_
#define DESUL_ATOMICS_ADAPT_CXX_HPP_

#include <atomic>
#include <desul/atomics/Common.hpp>

namespace desul {
namespace Impl {

template <class MemoryOrderDesul>
struct CXXMemoryOrder;

template <>
struct CXXMemoryOrder<MemoryOrderRelaxed> {
  static constexpr std::memory_order value = std::memory_order_relaxed;
};

template <>
struct CXXMemoryOrder<MemoryOrderAcquire> {
  static constexpr std::memory_order value = std::memory_order_acquire;
};

template <>
struct CXXMemoryOrder<MemoryOrderRelease> {
  static constexpr std::memory_order value = std::memory_order_release;
};

template <>
struct CXXMemoryOrder<MemoryOrderAcqRel> {
  static constexpr std::memory_order value = std::memory_order_acq_rel;
};

template <>
struct CXXMemoryOrder<MemoryOrderSeqCst> {
  static constexpr std::memory_order value = std::memory_order_seq_cst;
};

}  // namespace Impl
}  // namespace desul

#endif
