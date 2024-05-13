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
#include <exception>
#include <iostream>
#include <sstream>
#include <Kokkos_Core.hpp>

namespace {  // (anonymous)

// Output for the finalize hooks.  Use this to make sure that all the
// hooks ran, and that they ran in the correct order.
std::ostringstream hookOutput;

const char hook1str[] = "Behold, I am Hook 1; first pushed, last to be called.";
const char hook2str[] = "Yea verily, I am Hook 2.";
const char hook3str[] = "Indeed, I am Hook 3.";
const char hook4str[] = "Last but not least, I am Hook 4.";

}  // namespace

// Don't just have all the hooks print the same thing except for a
// number.  Have them print different things, so we can detect
// interleaving.  The hooks need to run sequentially, in LIFO order.
// Also, make sure that the function accepts at least the following
// kinds of hooks:
//
// 1. A plain old function that takes no arguments and returns nothing.
// 2. Lambda, that can be assigned to std::function<void()>
// 3. An actual std::function<void()>
// 4. A named object with operator().  This is what C++ programmers
//    unfortunately like to call "functor," even though this word
//    means something different in other languages.

void hook1() { hookOutput << hook1str << std::endl; }

struct Hook4 {
  void operator()() const { hookOutput << hook4str << std::endl; }
};

int main(int argc, char* argv[]) {
  using std::cout;
  using std::endl;

  const std::string expectedOutput([] {
    std::ostringstream os;
    os << hook4str << endl
       << hook3str << endl
       << hook2str << endl
       << hook1str << endl;
    return os.str();
  }());

  Kokkos::initialize(argc, argv);

  Kokkos::push_finalize_hook(hook1);  // plain old function
  Kokkos::push_finalize_hook([] { hookOutput << hook2str << endl; });  // lambda
  std::function<void()> hook3 = [] { hookOutput << hook3str << endl; };
  Kokkos::push_finalize_hook(hook3);  // actual std::function
  Hook4 hook4;
  Kokkos::push_finalize_hook(hook4);  // function object instance

  // This should invoke the finalize hooks in reverse order.
  // Furthermore, it should not throw an exception.
  try {
    Kokkos::finalize();
  } catch (std::exception& e) {
    cout << "FAILED: Kokkos::finalize threw an exception: " << e.what() << endl;
    return EXIT_FAILURE;
  } catch (...) {
    cout << "FAILED: Kokkos::finalize threw an exception whose base class "
            "is not std::exception."
         << endl;
    return EXIT_FAILURE;
  }

  const bool success = (hookOutput.str() == expectedOutput);
  if (success) {
    cout << "SUCCESS" << endl;
  } else {
    cout << "FAILED:" << endl
         << "  Expected output:" << endl
         << expectedOutput << endl
         << "  Actual output:" << endl
         << hookOutput.str() << endl;
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
