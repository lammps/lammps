// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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
  if constexpr (sizeof...(Args) == 0)
    ::printf("%s", format);
  else
    ::printf(format, args...);
#endif
}

}  // namespace Kokkos

#endif /* #ifndef KOKKOS_PRINTF_HPP */
