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

#include <iostream>
#include "Kokkos_Core.hpp"

#include <impl/Kokkos_Stacktrace.hpp>

namespace Test {

void stacktrace_test_f0(std::ostream& out);

int stacktrace_test_f1(std::ostream& out) {
  out << "Top of f1" << std::endl;
  stacktrace_test_f0(out);
  Kokkos::Impl::save_stacktrace();
  stacktrace_test_f0(out);

  return 42;
}

}  // namespace Test
