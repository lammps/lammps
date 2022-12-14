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
#include <cstdio>

namespace Test {

template <class Device>
struct TestScan {
  using execution_space = Device;
  using value_type      = int64_t;

  Kokkos::View<int, Device, Kokkos::MemoryTraits<Kokkos::Atomic> > errors;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int iwork, value_type& update,
                  const bool final_pass) const {
    const value_type n         = iwork + 1;
    const value_type imbalance = ((1000 <= n) && (0 == n % 1000)) ? 1000 : 0;

    // Insert an artificial load imbalance

    for (value_type i = 0; i < imbalance; ++i) {
      ++update;
    }

    update += n - imbalance;

    if (final_pass) {
      const value_type answer =
          n & 1 ? (n * ((n + 1) / 2)) : ((n / 2) * (n + 1));

      if (answer != update) {
        int fail = errors()++;

        if (fail < 20) {
          KOKKOS_IMPL_DO_NOT_USE_PRINTF("TestScan(%d,%ld) != %ld\n", iwork,
                                        static_cast<long>(update),
                                        static_cast<long>(answer));
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& update) const { update = 0; }

  KOKKOS_INLINE_FUNCTION
  void join(value_type& update, const value_type& input) const {
    update += input;
  }

  TestScan(const size_t N) {
    Kokkos::View<int, Device> errors_a("Errors");
    Kokkos::deep_copy(errors_a, 0);
    errors = errors_a;

    Kokkos::parallel_scan(N, *this);

    value_type total = 0;
    Kokkos::parallel_scan(N, *this, total);

    // We can't return a value in a constructor so use a lambda as wrapper to
    // ignore it.
    [&] { ASSERT_EQ(size_t((N + 1) * N / 2), size_t(total)); }();
    check_error();
  }

  TestScan(const size_t Start, const size_t N) {
    using exec_policy = Kokkos::RangePolicy<execution_space>;

    Kokkos::View<int, Device> errors_a("Errors");
    Kokkos::deep_copy(errors_a, 0);
    errors = errors_a;

    Kokkos::parallel_scan(exec_policy(Start, N), *this);
    Kokkos::fence();

    check_error();
  }

  void check_error() {
    int total_errors;
    Kokkos::deep_copy(total_errors, errors);
    ASSERT_EQ(total_errors, 0);
  }

  static void test_range(const size_t begin, const size_t end) {
    for (auto i = begin; i < end; ++i) {
      (void)TestScan(i);
    }
  }
};

TEST(TEST_CATEGORY, scan) {
  TestScan<TEST_EXECSPACE>::test_range(1, 1000);
  TestScan<TEST_EXECSPACE>(0);
  TestScan<TEST_EXECSPACE>(100000);
  TestScan<TEST_EXECSPACE>(10000000);
  TEST_EXECSPACE().fence();
}
}  // namespace Test
