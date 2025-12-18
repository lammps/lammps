// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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
