/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
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
#include <exception>
#include <iostream>
#include <sstream>
#include <Kokkos_Core.hpp>

namespace { // (anonymous)

// Output for the finalize hooks.  Use this to make sure that all the
// hooks ran, and that they ran in the correct order.
std::ostringstream hookOutput;

const char hook1str[] = "Behold, I am Hook 1; first pushed, last to be called.";
const char hook2str[] = "Yea verily, I am Hook 2.";
const char hook3str[] = "Indeed, I am Hook 3.";
const char hook4str[] = "Last but not least, I am Hook 4.";

} // namespace (anonymous)

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

void hook1 () {
  hookOutput << hook1str << std::endl;
}

struct Hook4 {
  void operator () () const {
    hookOutput << hook4str << std::endl;
  }
};

int main( int argc, char *argv[] ) {
  using std::cout;
  using std::endl;

  const std::string expectedOutput ([] {
      std::ostringstream os;
      os << hook4str << endl
         << hook3str << endl
         << hook2str << endl
         << hook1str << endl;
      return os.str();
    }());

  Kokkos::initialize(argc, argv);

  Kokkos::push_finalize_hook(hook1); // plain old function
  Kokkos::push_finalize_hook ([] {
      hookOutput << hook2str << endl;
    }); // lambda
  std::function<void()> hook3 = [] {
    hookOutput << hook3str << endl;
  };
  Kokkos::push_finalize_hook(hook3); // actual std::function
  Hook4 hook4;
  Kokkos::push_finalize_hook(hook4); // function object instance

  // This should invoke the finalize hooks in reverse order.
  // Furthermore, it should not throw an exception.
  try {
    Kokkos::finalize();
  }
  catch (std::exception& e) {
    cout << "FAILED: Kokkos::finalize threw an exception: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (...) {
    cout << "FAILED: Kokkos::finalize threw an exception whose base class "
      "is not std::exception." << endl;
    return EXIT_FAILURE;
  }

  const bool success = (hookOutput.str() == expectedOutput);
  if (success) {
    cout << "SUCCESS" << endl;
  }
  else {
    cout << "FAILED:" << endl
         << "  Expected output:" << endl
         << expectedOutput << endl
         << "  Actual output:" << endl
         << hookOutput.str() << endl;
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
