/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_THREAD_FENCE_SYCL_HPP_
#define DESUL_ATOMICS_THREAD_FENCE_SYCL_HPP_

#include <desul/atomics/Adapt_SYCL.hpp>
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

template <class MemoryOrder, class MemoryScope>
void device_atomic_thread_fence(MemoryOrder, MemoryScope) {
  sycl::atomic_fence(SYCLMemoryOrder<MemoryOrder>::value,
                     SYCLMemoryScope<MemoryScope>::value);
}

}  // namespace Impl
}  // namespace desul

#endif
