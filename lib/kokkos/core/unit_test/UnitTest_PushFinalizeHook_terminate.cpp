/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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
