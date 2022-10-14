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

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

namespace {

template <class ExecutionSpace>
struct CheckClassWithExecutionSpaceAsDataMemberIsCopyable {
  Kokkos::DefaultExecutionSpace device;
  Kokkos::DefaultHostExecutionSpace host;

  KOKKOS_FUNCTION void operator()(int, int& e) const {
    auto copy = *this;
    // not actually doing anything useful, mostly checking that
    // ExecutionSpace::in_parallel() is callable
    if (static_cast<int>(copy.device.in_parallel()) < 0) {
      ++e;
    }
  }

  CheckClassWithExecutionSpaceAsDataMemberIsCopyable() {
    int errors;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecutionSpace>(0, 1), *this,
                            errors);
    EXPECT_EQ(errors, 0);
  }
};

// FIXME_OPENMPTARGET nvlink error: Undefined reference to
// '_ZSt25__throw_bad_function_callv' in
// '/tmp/TestOpenMPTarget_ExecutionSpace-434d81.cubin'
#ifndef KOKKOS_ENABLE_OPENMPTARGET
TEST(TEST_CATEGORY, execution_space_as_class_data_member) {
  CheckClassWithExecutionSpaceAsDataMemberIsCopyable<TEST_EXECSPACE>();
}
#endif

}  // namespace
