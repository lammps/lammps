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
#ifdef KOKKOS_ENABLE_SYCL

// FIXME_SYCL
#if __has_include(<sycl/sycl.hpp>)
#include <sycl/sycl.hpp>
#else
#include <CL/sycl.hpp>
#endif

#ifndef KOKKOS_IMPL_HALF_TYPE_DEFINED
// Make sure no one else tries to define half_t
#define KOKKOS_IMPL_HALF_TYPE_DEFINED
#define KOKKOS_IMPL_SYCL_HALF_TYPE_DEFINED

namespace Kokkos {
namespace Impl {
struct half_impl_t {
  using type = sycl::half;
};
}  // namespace Impl
}  // namespace Kokkos
#endif  // KOKKOS_IMPL_HALF_TYPE_DEFINED
#endif  // KOKKOS_ENABLE_SYCL
#endif
