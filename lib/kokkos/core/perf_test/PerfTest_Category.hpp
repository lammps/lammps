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

#ifndef KOKKOS_TEST_PERFTEST_CAT_HPP
#define KOKKOS_TEST_PERFTEST_CAT_HPP

namespace Test {

inline int command_line_num_args(int n = 0) {
  static int n_args = 0;
  if (n > 0) n_args = n;
  return n_args;
}

inline const char* command_line_arg(int k, char** input_args = nullptr) {
  static char** args;
  if (input_args != nullptr) args = input_args;
  if (command_line_num_args() > k)
    return args[k];
  else
    return nullptr;
}

}  // namespace Test

#define TEST_CATEGORY default_exec
#define TEST_EXECSPACE Kokkos::DefaultExecutionSpace

#endif
