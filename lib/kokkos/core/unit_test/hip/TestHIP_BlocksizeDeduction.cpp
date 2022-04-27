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

#include <Kokkos_Core.hpp>
#include <TestHIP_Category.hpp>

namespace Test {

struct TestNone {
  Kokkos::View<size_t*, TEST_EXECSPACE> view;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const { view(i) = i; }

  TestNone() { view = Kokkos::View<size_t*, TEST_EXECSPACE>("dummy", 1); }
};

struct TestSpiller {
  Kokkos::View<size_t*, TEST_EXECSPACE> view;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    size_t array[1000] = {0};
    // and update flag
    size_t value = 0;
    for (int ii = i; ii < 1000; ++ii) {
      array[ii] = value;
      value += ii;
    }
    for (int ii = i; ii < 1000; ++ii) {
      value *= array[ii];
    }
    Kokkos::atomic_add(&view[0], value);
  }

  TestSpiller() { view = Kokkos::View<size_t*, TEST_EXECSPACE>("dummy", 1); }
};

TEST(hip, preferred_blocksize_deduction) {
  using execution_space =
      typename Kokkos::Impl::FunctorPolicyExecutionSpace<TestSpiller,
                                                         void>::execution_space;
  using policy = Kokkos::RangePolicy<execution_space>;

  {
    using DriverType = Kokkos::Impl::ParallelFor<TestNone, policy>;
    ASSERT_TRUE(Kokkos::Experimental::Impl::HIPParallelLaunch<
                    DriverType>::get_scratch_size() == 0);
  }

  {
    using DriverType = Kokkos::Impl::ParallelFor<TestSpiller, policy>;
    ASSERT_TRUE(Kokkos::Experimental::Impl::HIPParallelLaunch<
                    DriverType>::get_scratch_size() > 0);
  }
}

}  // namespace Test
