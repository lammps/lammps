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

#include <cstdlib>
#include <iostream>
#include <exception>
#include <Kokkos_Core.hpp>

// If any of the finalize hooks given to Kokkos::push_finalize_hook
// throws but does not catch an exception, make sure that
// Kokkos::finalize calls std::terminate.

namespace {  // (anonymous)

// If you change this, change CMakeLists.txt in this directory too!
// I verified that changing this string makes the test fail.
const char my_terminate_str[] =
    "PASSED: I am the custom std::terminate handler.";

// Tell compilers not to complain that this function doesn't return.
[[noreturn]] void my_terminate_handler() {
  std::cerr << my_terminate_str << std::endl;
  std::abort();  // terminate handlers normally would end by calling this
}

}  // namespace

int main(int argc, char *argv[]) {
  // If std::terminate is called, it will call my_terminate_handler.
  std::set_terminate(my_terminate_handler);

  Kokkos::initialize(argc, argv);
  Kokkos::push_finalize_hook(
      [] { throw std::runtime_error("I am an uncaught exception!"); });

  // This should call std::terminate, which in turn will call
  // my_terminate_handler above.  That will print the message that
  // makes this test count as passed.
  Kokkos::finalize();

  // The test actually failed if we got to this point.
  std::cerr << "FAILED to call std::terminate!" << std::endl;
  return EXIT_FAILURE;
}
