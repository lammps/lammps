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

size_t find_first_non_whitespace(const std::string& s, const size_t start_pos) {
  constexpr size_t num_ws_chars = 3;
  const char ws_chars[]         = "\n\t ";
  return s.find_first_not_of(ws_chars, start_pos, num_ws_chars);
}

void test_stacktrace(bool bTerminate, bool bCustom = true) {
  stacktrace_test_f1(std::cout);
  bool bDynamic = false;
  {
    std::stringstream sstream;
    Kokkos::Impl::print_saved_stacktrace(sstream);
    std::string foutput = sstream.str();

    bDynamic = std::string::npos != foutput.find("stacktrace");

    if (bDynamic) {
      printf("test_f1:\n%s \n", foutput.c_str());
      ASSERT_NE(std::string::npos, foutput.find("stacktrace_test_f1"));
      for (auto x : {"stacktrace_test_f0", "stacktrace_test_f2",
                     "stacktrace_test_f3", "stacktrace_test_f4"}) {
        ASSERT_EQ(std::string::npos, foutput.find(x));
      }
    }
  }

  {
    std::stringstream sstream;
    Kokkos::Impl::print_demangled_saved_stacktrace(sstream);

    if (bDynamic) {
      std::string foutput = sstream.str();
      printf("demangled test_f1:\n%s \n", foutput.c_str());
      ASSERT_NE(std::string::npos, foutput.find("Test::stacktrace_test_f1"));
      for (auto x : {"stacktrace_test_f0", "stacktrace_test_f2",
                     "stacktrace_test_f3", "stacktrace_test_f4"}) {
        ASSERT_EQ(std::string::npos, foutput.find(x));
      }
      EXPECT_EQ(0u, find_first_non_whitespace(foutput, 0));
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
      printf("test_f3:\n%s \n", foutput.c_str());
      for (auto x : {"stacktrace_test_f1", "stacktrace_test_f3"}) {
        ASSERT_NE(std::string::npos, foutput.find(x));
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
      printf("demangled test_f3:\n%s \n", foutput.c_str());
      for (auto x : {"stacktrace_test_f1", "stacktrace_test_f3"}) {
        ASSERT_NE(std::string::npos, foutput.find(x));
      }
      EXPECT_EQ(0u, find_first_non_whitespace(foutput, 0));
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
