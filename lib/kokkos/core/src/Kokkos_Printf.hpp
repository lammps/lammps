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

#ifndef KOKKOS_PRINTF_HPP
#define KOKKOS_PRINTF_HPP

#include <Kokkos_Macros.hpp>

#ifdef KOKKOS_ENABLE_SYCL
#include <sycl/sycl.hpp>
#else
#include <cstdio>
#endif

namespace Kokkos {

// In contrast to std::printf, return void to get a consistent behavior across
// backends. The GPU backends always return 1 and NVHPC only compiles if we
// don't ask for the return value.
template <typename... Args>
KOKKOS_FORCEINLINE_FUNCTION void printf(const char* format, Args... args) {
#ifdef KOKKOS_ENABLE_SYCL
  // Some compilers warn if "args" is empty and format is not a string literal
  if constexpr (sizeof...(Args) == 0)
    sycl::ext::oneapi::experimental::printf("%s", format);
  else
    sycl::ext::oneapi::experimental::printf(format, args...);
#else
  if constexpr (sizeof...(Args) == 0) ::printf("%s", format);
    // FIXME_OPENMPTARGET non-string-literal argument used in printf is not
    // supported for spir64
#if !(defined(KOKKOS_ENABLE_OPENMPTARGET) && defined(KOKKOS_ARCH_INTEL_GPU))
  else
    ::printf(format, args...);
#endif
#endif
}

}  // namespace Kokkos

#endif /* #ifndef KOKKOS_PRINTF_HPP */
