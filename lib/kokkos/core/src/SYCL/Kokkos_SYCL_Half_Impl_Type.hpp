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

#ifndef KOKKOS_SYCL_HALF_IMPL_TYPE_HPP_
#define KOKKOS_SYCL_HALF_IMPL_TYPE_HPP_

#include <Kokkos_Macros.hpp>

// FIXME_SYCL
#if __has_include(<sycl/sycl.hpp>)
#include <sycl/sycl.hpp>
#else
#include <CL/sycl.hpp>
#endif

// Make sure no one else tries to define half_t
#ifndef KOKKOS_IMPL_HALF_TYPE_DEFINED
#define KOKKOS_IMPL_HALF_TYPE_DEFINED
#define KOKKOS_IMPL_SYCL_HALF_TYPE_DEFINED

namespace Kokkos::Impl {
struct half_impl_t {
  using type = sycl::half;
};
}  // namespace Kokkos::Impl
#endif  // KOKKOS_IMPL_HALF_TYPE_DEFINED

// Make sure no one else tries to define bhalf_t
#ifndef KOKKOS_IMPL_BHALF_TYPE_DEFINED
// FIXME_SYCL Evaluate when to drop the check
#if __has_include(<sycl/ext/oneapi/bfloat16.hpp>)
#define KOKKOS_IMPL_BHALF_TYPE_DEFINED
#define KOKKOS_IMPL_SYCL_BHALF_TYPE_DEFINED
namespace Kokkos::Impl {
struct bhalf_impl_t {
  using type = sycl::ext::oneapi::bfloat16;
};
}  // namespace Kokkos::Impl
#elif defined(SYCL_EXT_ONEAPI_BFLOAT16) && defined(KOKKOS_ARCH_INTEL_GPU)
// FIXME_SYCL bfloat16 is only supported for compute capability 8.0 or higher
// on Nvidia GPUs but SYCL_EXT_ONEAPI_BFLOAT16 is defined even for lower compute
// capability.
#define KOKKOS_IMPL_BHALF_TYPE_DEFINED
#define KOKKOS_IMPL_SYCL_BHALF_TYPE_DEFINED
namespace Kokkos::Impl {
struct bhalf_impl_t {
  using type = sycl::ext::oneapi::experimental::bfloat16;
};
}  // namespace Kokkos::Impl
#endif  // test for bfloat16 support
#endif  // KOKKOS_IMPL_BHALF_TYPE_DEFINED
#endif  // KOKKOS_SYCL_HALF_IMPL_TYPE_HPP_
