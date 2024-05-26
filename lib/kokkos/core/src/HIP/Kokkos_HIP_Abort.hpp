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

#ifndef KOKKOS_HIP_ABORT_HPP
#define KOKKOS_HIP_ABORT_HPP

#include <Kokkos_Macros.hpp>

#include <hip/hip_runtime.h>

namespace Kokkos {
namespace Impl {

// The two keywords below are not contradictory. `noinline` is a
// directive to the optimizer.
[[noreturn]] __device__ __attribute__((noinline)) inline void hip_abort(
    char const *msg) {
  const char empty[] = "";
  __assert_fail(msg, empty, 0, empty);
  // This loop is never executed. It's intended to suppress warnings that the
  // function returns, even though it does not. This is necessary because
  // abort() is not marked as [[noreturn]], even though it does not return.
  while (true)
    ;
}

}  // namespace Impl
}  // namespace Kokkos

#endif
