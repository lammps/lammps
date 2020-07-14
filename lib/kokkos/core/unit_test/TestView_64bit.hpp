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

namespace Test {

template <class Device>
void test_64bit() {
#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
  int64_t N   = 5000000000;
  int64_t sum = 0;
  {
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<typename Device::execution_space,
                            Kokkos::IndexType<int64_t>>(0, N),
        KOKKOS_LAMBDA(const int64_t& /*i*/, int64_t& lsum) { lsum += 1; }, sum);
    ASSERT_EQ(N, sum);
  }
  {
    Kokkos::View<char*, Device> a("A", N);
    Kokkos::deep_copy(a, char(1));
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<typename Device::execution_space,
                            Kokkos::IndexType<int64_t>>(0, N),
        KOKKOS_LAMBDA(const int64_t& i, int64_t& lsum) {
          lsum += int64_t(a(i));
        },
        sum);
    ASSERT_EQ(N, sum);
    Kokkos::parallel_for(
        Kokkos::RangePolicy<typename Device::execution_space,
                            Kokkos::IndexType<int64_t>>(0, N),
        KOKKOS_LAMBDA(const int64_t& i) { a(i) = 3; });
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<typename Device::execution_space,
                            Kokkos::IndexType<int64_t>>(0, N),
        KOKKOS_LAMBDA(const int64_t& i, int64_t& lsum) {
          lsum += int64_t(a(i));
        },
        sum);
    ASSERT_EQ(N * 3, sum);
  }
  {
    int64_t N0 = 56925;
    int64_t N1 = 56927;

    Kokkos::View<char**, Device> m("Matrix", N0, N1);
    Kokkos::deep_copy(m, char(1));
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<typename Device::execution_space,
                            Kokkos::IndexType<int64_t>>(0, N0 * N1),
        KOKKOS_LAMBDA(const int64_t& i, int64_t& lsum) {
          lsum += int64_t(m(i % N0, i / N0));
        },
        sum);
    ASSERT_EQ(N0 * N1, sum);
    Kokkos::parallel_reduce(
        Kokkos::MDRangePolicy<typename Device::execution_space, Kokkos::Rank<2>,
                              Kokkos::IndexType<int64_t>>({0, 0}, {N0, N1}),
        KOKKOS_LAMBDA(const int64_t& i0, const int64_t& i1, int64_t& lsum) {
          lsum += int64_t(m(i0, i1));
        },
        sum);
    ASSERT_EQ(N0 * N1, sum);
  }
  {
    int N0    = 1024 * 1024 * 1500;
    int64_t P = 1713091;
    Kokkos::View<int*, Device> a("A", N0);
    Kokkos::parallel_for(
        "FillA",
        Kokkos::RangePolicy<typename Device::execution_space,
                            Kokkos::IndexType<int>>(0, N0),
        KOKKOS_LAMBDA(const int& i) { a(i) = i % P; });
    int64_t sum0 = 0;
    Kokkos::parallel_reduce(
        "FillA",
        Kokkos::RangePolicy<typename Device::execution_space,
                            Kokkos::IndexType<int>>(0, N0),
        KOKKOS_LAMBDA(const int& i, int64_t& lsum) { lsum += a(i); }, sum0);
    int64_t expected =
        (P * (P - 1) / 2) * int64_t(N0 / P) + (N0 % P) * (N0 % P - 1) / 2;
    ASSERT_EQ(expected, sum0);
  }
#endif
}

#ifdef KOKKOS_ENABLE_LARGE_MEM_TESTS
TEST(TEST_CATEGORY, view_64bit) { test_64bit<TEST_EXECSPACE>(); }
#endif

}  // namespace Test
