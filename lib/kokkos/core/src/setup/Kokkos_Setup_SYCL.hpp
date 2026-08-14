// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_SETUP_SYCL_HPP_
#define KOKKOS_SETUP_SYCL_HPP_

#include <sycl/sycl.hpp>

#ifndef KOKKOS_ENABLE_IMPL_SYCL_OUT_OF_ORDER_QUEUES
#define KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
#endif

#if defined(SYCL_EXT_ONEAPI_GRAPH) && \
    defined(KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES)
#define KOKKOS_IMPL_SYCL_GRAPH_SUPPORT
#endif

// clang-format off
#ifdef KOKKOS_ENABLE_SYCL_RELOCATABLE_DEVICE_CODE
#define KOKKOS_IMPL_RELOCATABLE_FUNCTION SYCL_EXTERNAL
#else
#define KOKKOS_IMPL_RELOCATABLE_FUNCTION @"KOKKOS_RELOCATABLE_FUNCTION requires Kokkos_ENABLE_SYCL_RELOCATABLE_DEVICE_CODE=ON"
#endif
// clang-format on

#endif
