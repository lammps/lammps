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
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <cstdlib>
#include <iostream>
#include <Kokkos_Abort.hpp>
#include <impl/Kokkos_Stacktrace.hpp>

namespace Kokkos {
namespace Impl {

void host_abort(const char *const message) {
  std::cerr << message;

#ifdef KOKKOS_IMPL_ENABLE_STACKTRACE
  std::cerr << "\nBacktrace:\n";
  save_stacktrace();
  print_demangled_saved_stacktrace(std::cerr);
#else
  std::cerr << "\nTraceback functionality not available\n";
#endif

  ::abort();
}

}  // namespace Impl
}  // namespace Kokkos
