/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_ADAPT_HIP_HPP_
#define DESUL_ATOMICS_ADAPT_HIP_HPP_

#include <desul/atomics/Common.hpp>

namespace desul {
namespace Impl {

// FIXME same code as GCCMemoryOrder
template <class MemoryOrder>
struct HIPMemoryOrder;

template <>
struct HIPMemoryOrder<MemoryOrderRelaxed> {
  static constexpr int value = __ATOMIC_RELAXED;
};

template <>
struct HIPMemoryOrder<MemoryOrderAcquire> {
  static constexpr int value = __ATOMIC_ACQUIRE;
};

template <>
struct HIPMemoryOrder<MemoryOrderRelease> {
  static constexpr int value = __ATOMIC_RELEASE;
};

template <>
struct HIPMemoryOrder<MemoryOrderAcqRel> {
  static constexpr int value = __ATOMIC_ACQ_REL;
};

template <>
struct HIPMemoryOrder<MemoryOrderSeqCst> {
  static constexpr int value = __ATOMIC_SEQ_CST;
};

// __HIP_MEMORY_SCOPE_SYSTEM
// __HIP_MEMORY_SCOPE_AGENT
// __HIP_MEMORY_SCOPE_WORKGROUP
// __HIP_MEMORY_SCOPE_WAVEFRONT
// __HIP_MEMORY_SCOPE_SINGLETHREAD
template <class MemoryScope>
struct HIPMemoryScope;

template <>
struct HIPMemoryScope<MemoryScopeCore> {
  static constexpr int value = __HIP_MEMORY_SCOPE_WORKGROUP;
};

template <>
struct HIPMemoryScope<MemoryScopeDevice> {
  static constexpr int value = __HIP_MEMORY_SCOPE_AGENT;
};

template <>
struct HIPMemoryScope<MemoryScopeNode> {
  static constexpr int value = __HIP_MEMORY_SCOPE_SYSTEM;
};

template <>
struct HIPMemoryScope<MemoryScopeSystem> {
  static constexpr int value = __HIP_MEMORY_SCOPE_SYSTEM;
};

}  // namespace Impl
}  // namespace desul

#endif
