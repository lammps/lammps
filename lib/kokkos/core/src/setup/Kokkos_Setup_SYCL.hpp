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

#ifndef KOKKOS_SETUP_SYCL_HPP_
#define KOKKOS_SETUP_SYCL_HPP_

// FIXME_SYCL Using in-order queues currently gives better performance on Intel
// GPUs and we run into correctness issues with out-of-order queues on NVIDIA
// GPUs.
#define KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES

// FIXME_SYCL the fallback assert is temporarily disabled by default in the
// compiler so we need to force it
#ifndef SYCL_ENABLE_FALLBACK_ASSERT
#define SYCL_ENABLE_FALLBACK_ASSERT
#endif
#ifndef SYCL_FALLBACK_ASSERT
#define SYCL_FALLBACK_ASSERT 1
#endif

// FIXME_SYCL
#if __has_include(<sycl/sycl.hpp>)
#include <sycl/sycl.hpp>
#else
#include <CL/sycl.hpp>
#endif

#if defined(__INTEL_LLVM_COMPILER) && __INTEL_LLVM_COMPILER >= 20230200
#define KOKKOS_IMPL_SYCL_GET_MULTI_PTR(accessor) \
  accessor.get_multi_ptr<sycl::access::decorated::yes>()
#else
#define KOKKOS_IMPL_SYCL_GET_MULTI_PTR(accessor) accessor.get_pointer()
#endif

// FIXME_SYCL Use type directly once it has stabilized in SYCL.
namespace Kokkos::Impl {
#ifndef SYCL_EXT_INTEL_USM_ADDRESS_SPACES
#error SYCL_EXT_INTEL_USM_ADDRESS_SPACES undefined!
#elif SYCL_EXT_INTEL_USM_ADDRESS_SPACES >= 2
template <typename T>
using sycl_device_ptr = sycl::ext::intel::device_ptr<T>;
template <typename T>
using sycl_host_ptr = sycl::ext::intel::host_ptr<T>;
#else
template <typename T>
using sycl_device_ptr = sycl::device_ptr<T>;
template <typename T>
using sycl_host_ptr = sycl::host_ptr<T>;
#endif
}  // namespace Kokkos::Impl

#endif
