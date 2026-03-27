// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_SETUP_SYCL_HPP_
#define KOKKOS_SETUP_SYCL_HPP_

// FIXME_SYCL
#if __has_include(<sycl/sycl.hpp>)
#include <sycl/sycl.hpp>
#else
#include <CL/sycl.hpp>
#endif

#ifndef KOKKOS_ENABLE_IMPL_SYCL_OUT_OF_ORDER_QUEUES
#define KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
#endif

#if defined(SYCL_EXT_ONEAPI_GRAPH) && \
    defined(KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES)
#define KOKKOS_IMPL_SYCL_GRAPH_SUPPORT
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

// clang-format off
#ifdef KOKKOS_ENABLE_SYCL_RELOCATABLE_DEVICE_CODE
#define KOKKOS_IMPL_RELOCATABLE_FUNCTION SYCL_EXTERNAL
#else
#define KOKKOS_IMPL_RELOCATABLE_FUNCTION @"KOKKOS_RELOCATABLE_FUNCTION requires Kokkos_ENABLE_SYCL_RELOCATABLE_DEVICE_CODE=ON"
#endif
// clang-format on

#endif
