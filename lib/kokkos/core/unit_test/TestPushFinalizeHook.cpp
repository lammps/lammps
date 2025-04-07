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

#include <Kokkos_Core.hpp>

#include <gtest/gtest.h>

#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>

#include "KokkosExecutionEnvironmentNeverInitializedFixture.hpp"

namespace {

using PushFinalizeHook_DeathTest = KokkosExecutionEnvironmentNeverInitialized;

// Output for the finalize hooks.  Use this to make sure that all the hooks
// ran, and that they ran in the correct order.
std::ostringstream hookOutput;

const char hook1str[] = "Behold, I am Hook 1; first pushed, last to be called.";
const char hook2str[] = "Yea verily, I am Hook 2.";
const char hook3str[] = "Indeed, I am Hook 3.";
const char hook4str[] = "Last but not least, I am Hook 4.";

// Don't just have all the hooks print the same thing except for a number.
// Have them print different things, so we can detect interleaving.  The hooks
// need to run sequentially, in LIFO order.  Also, make sure that the function
// accepts at least the following kinds of hooks:
//
// 1. A plain old function that takes no arguments and returns nothing.
// 2. Lambda, that can be assigned to std::function<void()>
// 3. An actual std::function<void()>
// 4. A named object with operator().  This is what C++ programmers
//    unfortunately like to call "functor," even though this word means
//    something different in other languages.

void hook1() { hookOutput << hook1str << '\n'; }

struct Hook4 {
  void operator()() const { hookOutput << hook4str << '\n'; }
};

TEST_F(PushFinalizeHook_DeathTest, called_in_reverse_order) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  std::string const expectedOutput([] {
    std::ostringstream os;
    os << hook4str << '\n'
       << hook3str << '\n'
       << hook2str << '\n'
       << hook1str << '\n';
    return os.str();
  }());

  EXPECT_EXIT(
      {
        Kokkos::push_finalize_hook(hook1);  // plain old function
        Kokkos::push_finalize_hook(
            [] { hookOutput << hook2str << '\n'; });  // lambda
        Kokkos::initialize();
        std::function<void()> hook3 = [] { hookOutput << hook3str << '\n'; };
        Kokkos::push_finalize_hook(hook3);  // actual std::function
        Hook4 hook4;
        Kokkos::push_finalize_hook(hook4);  // function object instance

        Kokkos::finalize();
        std::exit(hookOutput.str() == expectedOutput ? EXIT_SUCCESS
                                                     : EXIT_FAILURE);
      },
      ::testing::ExitedWithCode(EXIT_SUCCESS), "");
}

char const my_terminate_handler_msg[] = "my terminate handler was called\n";
[[noreturn]] void my_terminate_handler() {
  std::cerr << my_terminate_handler_msg;
  std::abort();
}

TEST_F(PushFinalizeHook_DeathTest, terminate_on_throw) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  auto terminate_handler = std::get_terminate();

  std::set_terminate(my_terminate_handler);

  EXPECT_DEATH(
      {
        Kokkos::push_finalize_hook(
            [] { throw std::runtime_error("uncaught exception"); });
        Kokkos::initialize();
        Kokkos::finalize();
      },
      my_terminate_handler_msg);

  std::set_terminate(terminate_handler);
}

}  // namespace
