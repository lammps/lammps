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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif

#ifndef KOKKOS_EXPERIMENTAL_MDSPAN_HPP
#define KOKKOS_EXPERIMENTAL_MDSPAN_HPP

// Look for the right mdspan
#if __cplusplus >= 202002L
#include <version>
#endif

// Only use standard library mdspan if we are not running Cuda or HIP.
// Likely these implementations won't be supported on device, so we should use
// our own device-compatible version for now.
#if (__cpp_lib_mdspan >= 202207L) && !defined(KOKKOS_ENABLE_CUDA) && \
    !defined(KOKKOS_ENABLE_HIP)
#include <mdspan>
namespace Kokkos {
using std::default_accessor;
using std::dextents;
using std::dynamic_extent;
using std::extents;
using std::layout_left;
using std::layout_right;
using std::layout_stride;
using std::mdspan;
}  // namespace Kokkos
#else
#include <mdspan/mdspan.hpp>
#endif

#endif  // KOKKOS_EXPERIMENTAL_MDSPAN_HPP
