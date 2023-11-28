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

#ifdef __SYCL_DEVICE_ONLY__
#define KOKKOS_IMPL_DO_NOT_USE_PRINTF(format, ...)                \
  do {                                                            \
    const __attribute__((opencl_constant)) char fmt[] = (format); \
    sycl::ext::oneapi::experimental::printf(fmt, ##__VA_ARGS__);  \
  } while (0)
#endif

#endif
