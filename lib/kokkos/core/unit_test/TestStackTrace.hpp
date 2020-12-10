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
#include <iostream>
#include <gtest/gtest.h>
#include "Kokkos_Core.hpp"

#include <impl/Kokkos_Stacktrace.hpp>

namespace Test {

void stacktrace_test_f0(std::ostream& out);

int stacktrace_test_f1(std::ostream& out);

void stacktrace_test_f2(std::ostream& out);

int stacktrace_test_f3(std::ostream& out, const int level);

void stacktrace_test_f4();

void my_fancy_handler();

void test_stacktrace(bool bTerminate, bool bCustom = true) {
  stacktrace_test_f1(std::cout);
  bool bDynamic = false;
  {
    std::stringstream sstream;
    Kokkos::Impl::print_saved_stacktrace(sstream);
    std::string foutput = sstream.str();

    bDynamic = std::string::npos != foutput.find("stacktrace");

    if (bDynamic) {
      printf("test_f1: %s \n", foutput.c_str());
      ASSERT_TRUE(std::string::npos != foutput.find("stacktrace_test_f1"));
      for (auto x : {"stacktrace_test_f0", "stacktrace_test_f2",
                     "stacktrace_test_f3", "stacktrace_test_f4"}) {
        ASSERT_TRUE(std::string::npos == foutput.find(x));
      }
    }
  }

  {
    std::stringstream sstream;
    Kokkos::Impl::print_demangled_saved_stacktrace(sstream);

    if (bDynamic) {
      std::string foutput = sstream.str();
      printf("demangled test_f1: %s \n", foutput.c_str());
      ASSERT_TRUE(std::string::npos !=
                  foutput.find("Test::stacktrace_test_f1"));
      for (auto x : {"stacktrace_test_f0", "stacktrace_test_f2",
                     "stacktrace_test_f3", "stacktrace_test_f4"}) {
        ASSERT_TRUE(std::string::npos == foutput.find(x));
      }
    }
  }

  int val = stacktrace_test_f3(std::cout, 4);

  // Don't remove this, otherwise the compiler will optimize away call sequences
  // via
  printf("StackTrace f3(std::cout, 4) returned: %i\n", val);

  // TODO test by making sure that f3 and f1, but no other functions,
  // appear in the stack trace, and that f3 appears 5 times.
  // Fix that f3 doesn't show up when compiling with -O3
  {
    std::stringstream sstream;
    Kokkos::Impl::print_saved_stacktrace(sstream);

    if (bDynamic) {
      std::string foutput = sstream.str();
      printf("test_f3: %s \n", foutput.c_str());
      for (auto x : {"stacktrace_test_f1", "stacktrace_test_f3"}) {
        ASSERT_TRUE(std::string::npos != foutput.find(x));
      }
    }
    // TODO make sure stacktrace_test_f2/4 don't show up
    // TODO make sure stacktrace_test_f3 shows up 5 times
  }

  {
    std::stringstream sstream;
    Kokkos::Impl::print_demangled_saved_stacktrace(sstream);

    if (bDynamic) {
      std::string foutput = sstream.str();
      printf("demangled test_f3: %s \n", foutput.c_str());
      for (auto x : {"stacktrace_test_f1", "stacktrace_test_f3"}) {
        ASSERT_TRUE(std::string::npos != foutput.find(x));
      }
    }

    // TODO make sure stacktrace_test_f2/4 don't show up
    // TODO make sure stacktrace_test_f3 shows up 5 times
  }
  std::cout << "Test setting std::terminate handler that prints "
               "the last saved stack trace"
            << std::endl;

  stacktrace_test_f4();

  if (bCustom) {
    Kokkos::Impl::set_kokkos_terminate_handler(my_fancy_handler);
  } else {
    Kokkos::Impl::set_kokkos_terminate_handler();
  }

  // TODO test that this prints "Oh noes!" and the correct stacktrace.
  if (bTerminate) {
    std::terminate();
  }
}

TEST(defaultdevicetype, stacktrace_normal) { test_stacktrace(false); }

TEST(defaultdevicetype_DeathTest, stacktrace_terminate) {
  ASSERT_DEATH({ test_stacktrace(true); },
               "I am the custom std::terminate handler.");
}

TEST(defaultdevicetype_DeathTest, stacktrace_generic_term) {
  ASSERT_DEATH({ test_stacktrace(true, false); },
               "Kokkos observes that std::terminate has been called");
}

}  // namespace Test
