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

#include <cstdio>
#include <sstream>
#include <iostream>

#include <Kokkos_Core.hpp>

namespace Test {
namespace {
struct SumFunctor {
  KOKKOS_INLINE_FUNCTION
  void operator()(int i, int& lsum) const { lsum += i; }
};

template <class ExecSpace>
void check_distinctive(ExecSpace, ExecSpace) {}

#ifdef KOKKOS_ENABLE_CUDA
void check_distinctive(Kokkos::Cuda exec1, Kokkos::Cuda exec2) {
  ASSERT_NE(exec1.cuda_stream(), exec2.cuda_stream());
}
#endif
#ifdef KOKKOS_ENABLE_HIP
void check_distinctive(Kokkos::Experimental::HIP exec1,
                       Kokkos::Experimental::HIP exec2) {
  ASSERT_NE(exec1.hip_stream(), exec2.hip_stream());
}
#endif
#ifdef KOKKOS_ENABLE_SYCL
void check_distinctive(Kokkos::Experimental::SYCL exec1,
                       Kokkos::Experimental::SYCL exec2) {
  ASSERT_NE(*exec1.impl_internal_space_instance()->m_queue,
            *exec2.impl_internal_space_instance()->m_queue);
}
#endif
}  // namespace

void test_partitioning(std::vector<TEST_EXECSPACE>& instances) {
  check_distinctive(instances[0], instances[1]);
  int sum1, sum2;
  int N = 3910;
  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<TEST_EXECSPACE>(instances[0], 0, N), SumFunctor(),
      sum1);
  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<TEST_EXECSPACE>(instances[1], 0, N), SumFunctor(),
      sum2);
  ASSERT_EQ(sum1, sum2);
  ASSERT_EQ(sum1, N * (N - 1) / 2);

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || \
    defined(KOKKOS_ENABLE_SYCL)
  // Eliminate unused function warning
  // (i.e. when compiling for Serial and CUDA, during Serial compilation the
  // Cuda overload is unused ...)
  if (sum1 != sum2) {
#ifdef KOKKOS_ENABLE_CUDA
    check_distinctive(Kokkos::Cuda(), Kokkos::Cuda());
#endif
#ifdef KOKKOS_ENABLE_HIP
    check_distinctive(Kokkos::Experimental::HIP(), Kokkos::Experimental::HIP());
#endif
#ifdef KOKKOS_ENABLE_SYCL
    check_distinctive(Kokkos::Experimental::SYCL(),
                      Kokkos::Experimental::SYCL());
#endif
  }
#endif
}

TEST(TEST_CATEGORY, partitioning_by_args) {
  auto instances =
      Kokkos::Experimental::partition_space(TEST_EXECSPACE(), 1, 1.);
  ASSERT_EQ(int(instances.size()), 2);
  test_partitioning(instances);
}

TEST(TEST_CATEGORY, partitioning_by_vector) {
  std::vector<int> weights{1, 1};
  auto instances =
      Kokkos::Experimental::partition_space(TEST_EXECSPACE(), weights);
  ASSERT_EQ(int(instances.size()), 2);
  test_partitioning(instances);
}
}  // namespace Test
