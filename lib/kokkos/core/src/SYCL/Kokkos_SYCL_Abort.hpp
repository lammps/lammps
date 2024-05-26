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

#ifndef KOKKOS_SYCL_ABORT_HPP
#define KOKKOS_SYCL_ABORT_HPP

#include <Kokkos_Printf.hpp>
#if defined(KOKKOS_ENABLE_SYCL)
// FIXME_SYCL
#if __has_include(<sycl/sycl.hpp>)
#include <sycl/sycl.hpp>
#else
#include <CL/sycl.hpp>
#endif

namespace Kokkos {
namespace Impl {

inline void sycl_abort(char const* msg) {
#ifdef NDEBUG
  Kokkos::printf("Aborting with message %s.\n", msg);
#else
  // Choosing "" here causes problems but a single whitespace character works.
  const char* empty = " ";
  __assert_fail(msg, empty, 0, empty);
#endif
}

}  // namespace Impl
}  // namespace Kokkos

#endif
#endif
