/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_ADAPT_SYCL_HPP_
#define DESUL_ATOMICS_ADAPT_SYCL_HPP_

#include <desul/atomics/Common.hpp>

// FIXME_SYCL SYCL2020 dictates that <sycl/sycl.hpp> is the header to include
// but icpx 2022.1.0 and earlier versions only provide <CL/sycl.hpp>
#if __has_include(<sycl/sycl.hpp>)
#include <sycl/sycl.hpp>
#else
#include <CL/sycl.hpp>
#endif

namespace desul {
namespace Impl {

//<editor-fold desc="SYCL memory order">
// The default memory order for sycl::atomic_ref
// can be seq_cst, acq_rel, or relaxed according to the
// "SYCL 2020 Specification (revision 6)", see
// https://registry.khronos.org/SYCL/specs/sycl-2020/html/sycl-2020.html#sec:atomic-references.
// Thus, we map MemoryOrderAcquire and MemoryOrderRelease to acq_rel.
template <class MemoryOrder>
struct SYCLMemoryOrder;

template <>
struct SYCLMemoryOrder<MemoryOrderSeqCst> {
  static constexpr sycl::memory_order value = sycl::memory_order::seq_cst;
};
template <>
struct SYCLMemoryOrder<MemoryOrderAcquire> {
  static constexpr sycl::memory_order value = sycl::memory_order::acq_rel;
};
template <>
struct SYCLMemoryOrder<MemoryOrderRelease> {
  static constexpr sycl::memory_order value = sycl::memory_order::acq_rel;
};
template <>
struct SYCLMemoryOrder<MemoryOrderAcqRel> {
  static constexpr sycl::memory_order value = sycl::memory_order::acq_rel;
};
template <>
struct SYCLMemoryOrder<MemoryOrderRelaxed> {
  static constexpr sycl::memory_order value = sycl::memory_order::relaxed;
};
//</editor-fold>

//<editor-fold desc="SYCL memory scope">
template <class MemoryScope>
struct SYCLMemoryScope;

template <>
struct SYCLMemoryScope<MemoryScopeCore> {
  static constexpr sycl::memory_scope value = sycl::memory_scope::work_group;
};

template <>
struct SYCLMemoryScope<MemoryScopeDevice> {
  static constexpr sycl::memory_scope value = sycl::memory_scope::device;
};

template <>
struct SYCLMemoryScope<MemoryScopeSystem> {
  static constexpr sycl::memory_scope value = sycl::memory_scope::system;
};
//</editor-fold>

// FIXME_SYCL generic_space isn't available yet for CUDA.
#ifdef __NVPTX__
template <class T, class MemoryOrder, class MemoryScope>
using sycl_atomic_ref = sycl::atomic_ref<T,
                                         SYCLMemoryOrder<MemoryOrder>::value,
                                         SYCLMemoryScope<MemoryScope>::value,
                                         sycl::access::address_space::global_space>;
#else
template <class T, class MemoryOrder, class MemoryScope>
using sycl_atomic_ref = sycl::atomic_ref<T,
                                         SYCLMemoryOrder<MemoryOrder>::value,
                                         SYCLMemoryScope<MemoryScope>::value,
                                         sycl::access::address_space::generic_space>;
#endif

#ifdef DESUL_SYCL_DEVICE_GLOBAL_SUPPORTED
#ifdef SYCL_EXT_ONEAPI_DEVICE_GLOBAL
template <class T>
using sycl_device_global = sycl::ext::oneapi::experimental::device_global<T>;
#else
template <class T>
using sycl_device_global = sycl::ext::oneapi::experimental::device_global<
    T,
    decltype(sycl::ext::oneapi::experimental::properties(
        sycl::ext::oneapi::experimental::device_image_scope))>;
#endif
#endif

}  // namespace Impl
}  // namespace desul

#endif
