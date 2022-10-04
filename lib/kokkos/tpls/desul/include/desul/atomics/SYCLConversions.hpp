/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_SYCL_CONVERSIONS_HPP_
#define DESUL_ATOMICS_SYCL_CONVERSIONS_HPP_
#ifdef DESUL_HAVE_SYCL_ATOMICS

// clang-format off
#include "desul/atomics/Common.hpp"

#include <CL/sycl.hpp>
// clang-format on

namespace desul {
namespace Impl {

#ifdef __clang__
namespace sycl_sync_and_atomics = ::sycl::ext::oneapi;
#else
namespace sycl_sync_and_atomics = ::sycl;
#endif

template <bool extended_namespace>
using sycl_memory_order = std::conditional_t<extended_namespace,
                                             sycl_sync_and_atomics::memory_order,
                                             sycl::memory_order>;
template <bool extended_namespace>
using sycl_memory_scope = std::conditional_t<extended_namespace,
                                             sycl_sync_and_atomics::memory_scope,
                                             sycl::memory_scope>;

template <class MemoryOrder, bool extended_namespace = true>
struct DesulToSYCLMemoryOrder;
template <bool extended_namespace>
struct DesulToSYCLMemoryOrder<MemoryOrderSeqCst, extended_namespace> {
  static constexpr sycl_memory_order<extended_namespace> value =
      sycl_memory_order<extended_namespace>::seq_cst;
};
template <bool extended_namespace>
struct DesulToSYCLMemoryOrder<MemoryOrderAcquire, extended_namespace> {
  static constexpr sycl_memory_order<extended_namespace> value =
      sycl_memory_order<extended_namespace>::acquire;
};
template <bool extended_namespace>
struct DesulToSYCLMemoryOrder<MemoryOrderRelease, extended_namespace> {
  static constexpr sycl_memory_order<extended_namespace> value =
      sycl_memory_order<extended_namespace>::release;
};
template <bool extended_namespace>
struct DesulToSYCLMemoryOrder<MemoryOrderAcqRel, extended_namespace> {
  static constexpr sycl_memory_order<extended_namespace> value =
      sycl_memory_order<extended_namespace>::acq_rel;
};
template <bool extended_namespace>
struct DesulToSYCLMemoryOrder<MemoryOrderRelaxed, extended_namespace> {
  static constexpr sycl_memory_order<extended_namespace> value =
      sycl_memory_order<extended_namespace>::relaxed;
};

template <class MemoryScope, bool extended_namespace = true>
struct DesulToSYCLMemoryScope;
template <bool extended_namespace>
struct DesulToSYCLMemoryScope<MemoryScopeCore, extended_namespace> {
  static constexpr sycl_memory_scope<extended_namespace> value =
      sycl_memory_scope<extended_namespace>::work_group;
};
template <bool extended_namespace>
struct DesulToSYCLMemoryScope<MemoryScopeDevice, extended_namespace> {
  static constexpr sycl_memory_scope<extended_namespace> value =
      sycl_memory_scope<extended_namespace>::device;
};
template <bool extended_namespace>
struct DesulToSYCLMemoryScope<MemoryScopeSystem, extended_namespace> {
  static constexpr sycl_memory_scope<extended_namespace> value =
      sycl_memory_scope<extended_namespace>::system;
};

// FIXME_SYCL generic_space isn't available yet for CUDA.
#ifdef __NVPTX__
template <class T, class MemoryOrder, class MemoryScope>
using sycl_atomic_ref = sycl::atomic_ref<T,
                                         DesulToSYCLMemoryOrder<MemoryOrder>::value,
                                         DesulToSYCLMemoryScope<MemoryScope>::value,
                                         sycl::access::address_space::global_space>;
#else
template <class T, class MemoryOrder, class MemoryScope>
using sycl_atomic_ref = sycl::atomic_ref<T,
                                         DesulToSYCLMemoryOrder<MemoryOrder>::value,
                                         DesulToSYCLMemoryScope<MemoryScope>::value,
                                         sycl::access::address_space::generic_space>;
#endif
}  // namespace Impl
}  // namespace desul

#endif
#endif
