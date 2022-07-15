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

#ifndef KOKKOS_TEST_TEAM_BASIC_HPP
#define KOKKOS_TEST_TEAM_BASIC_HPP
#include <TestTeam.hpp>

namespace Test {

TEST(TEST_CATEGORY, team_for) {
  TestTeamPolicy<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static> >::test_for(
      0);
  TestTeamPolicy<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic> >::test_for(
      0);

  TestTeamPolicy<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static> >::test_for(
      2);
  TestTeamPolicy<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic> >::test_for(
      2);

  TestTeamPolicy<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static> >::test_for(
      1000);
  TestTeamPolicy<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic> >::test_for(
      1000);
}

// FIXME_OPENMPTARGET wrong results
#ifndef KOKKOS_ENABLE_OPENMPTARGET
TEST(TEST_CATEGORY, team_reduce) {
  TestTeamPolicy<TEST_EXECSPACE,
                 Kokkos::Schedule<Kokkos::Static> >::test_reduce(0);
  TestTeamPolicy<TEST_EXECSPACE,
                 Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce(0);
  TestTeamPolicy<TEST_EXECSPACE,
                 Kokkos::Schedule<Kokkos::Static> >::test_reduce(2);
  TestTeamPolicy<TEST_EXECSPACE,
                 Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce(2);
  TestTeamPolicy<TEST_EXECSPACE,
                 Kokkos::Schedule<Kokkos::Static> >::test_reduce(1000);
  TestTeamPolicy<TEST_EXECSPACE,
                 Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce(1000);
}
#endif

template <typename ExecutionSpace>
struct TestTeamReduceLarge {
  using team_policy_t = Kokkos::TeamPolicy<ExecutionSpace>;
  using member_t      = typename team_policy_t::member_type;

  int m_range;

  TestTeamReduceLarge(const int range) : m_range(range) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_t& t, int& update) const {
    Kokkos::single(Kokkos::PerTeam(t), [&]() { update++; });
  }

  void run() {
    int result = 0;
    Kokkos::parallel_reduce(team_policy_t(m_range, Kokkos::AUTO), *this,
                            result);
    EXPECT_EQ(m_range, result);
  }
};

TEST(TEST_CATEGORY, team_reduce_large) {
  std::vector<int> ranges{(2LU << 23) - 1, 2LU << 23, (2LU << 24),
                          (2LU << 24) + 1, 1LU << 29};
  for (const auto range : ranges) {
    TestTeamReduceLarge<TEST_EXECSPACE> test(range);
    test.run();
  }
}

TEST(TEST_CATEGORY, team_broadcast_long) {
  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                    long>::test_teambroadcast(0, 1);
  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                    long>::test_teambroadcast(0, 1);

  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                    long>::test_teambroadcast(2, 1);
  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                    long>::test_teambroadcast(2, 1);

  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                    long>::test_teambroadcast(16, 1);
  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                    long>::test_teambroadcast(16, 1);

  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                    long>::test_teambroadcast(1000, 1);
  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                    long>::test_teambroadcast(1000, 1);
}

// FIXME_OPENMPTARGET CI fails with
// Libomptarget error: Copying data from device failed.
// Possibly, because long_wrapper is not trivially-copyable.
#ifndef KOKKOS_ENABLE_OPENMPTARGET
struct long_wrapper {
  long value;

  KOKKOS_FUNCTION
  long_wrapper() : value(0) {}

  KOKKOS_FUNCTION
  long_wrapper(long val) : value(val) {}

  KOKKOS_FUNCTION
  friend void operator+=(long_wrapper& lhs, const long_wrapper& rhs) {
    lhs.value += rhs.value;
  }

  KOKKOS_FUNCTION
  friend void operator+=(volatile long_wrapper& lhs,
                         const volatile long_wrapper& rhs) {
    lhs.value += rhs.value;
  }

  KOKKOS_FUNCTION
  void operator=(const long_wrapper& other) { value = other.value; }

  KOKKOS_FUNCTION
  void operator=(const volatile long_wrapper& other) volatile {
    value = other.value;
  }
  KOKKOS_FUNCTION
  operator long() const { return value; }
};
}  // namespace Test

namespace Kokkos {
template <>
struct reduction_identity<Test::long_wrapper>
    : public reduction_identity<long> {};
}  // namespace Kokkos

namespace Test {

// Test for non-arithmetic type
TEST(TEST_CATEGORY, team_broadcast_long_wrapper) {
  static_assert(!std::is_arithmetic<long_wrapper>::value, "");

  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                    long_wrapper>::test_teambroadcast(0, 1);
  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                    long_wrapper>::test_teambroadcast(0, 1);

  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                    long_wrapper>::test_teambroadcast(2, 1);
  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                    long_wrapper>::test_teambroadcast(2, 1);
  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                    long_wrapper>::test_teambroadcast(16, 1);
  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                    long_wrapper>::test_teambroadcast(16, 1);

  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                    long_wrapper>::test_teambroadcast(1000, 1);
  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                    long_wrapper>::test_teambroadcast(1000, 1);
}
#endif

TEST(TEST_CATEGORY, team_broadcast_char) {
  {
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                      unsigned char>::test_teambroadcast(0, 1);
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                      unsigned char>::test_teambroadcast(0, 1);

    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                      unsigned char>::test_teambroadcast(2, 1);
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                      unsigned char>::test_teambroadcast(2, 1);

    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                      unsigned char>::test_teambroadcast(16, 1);
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                      unsigned char>::test_teambroadcast(16, 1);

    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                      long>::test_teambroadcast(1000, 1);
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                      long>::test_teambroadcast(1000, 1);
  }
}

TEST(TEST_CATEGORY, team_broadcast_float) {
  {
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                      float>::test_teambroadcast(0, 1.3);
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                      float>::test_teambroadcast(0, 1.3);

    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                      float>::test_teambroadcast(2, 1.3);
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                      float>::test_teambroadcast(2, 1.3);

    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                      float>::test_teambroadcast(16, 1.3);
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                      float>::test_teambroadcast(16, 1.3);

    // FIXME_CUDA
#ifdef KOKKOS_ENABLE_CUDA
    if (!std::is_same<TEST_EXECSPACE, Kokkos::Cuda>::value)
#endif
    // FIXME_HIP
#ifdef KOKKOS_ENABLE_HIP
      if (!std::is_same<TEST_EXECSPACE, Kokkos::Experimental::HIP>::value)
#endif
      {
        TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                          float>::test_teambroadcast(1000, 1.3);
        TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                          float>::test_teambroadcast(1000, 1.3);
      }
  }
}

TEST(TEST_CATEGORY, team_broadcast_double) {
  {
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                      double>::test_teambroadcast(0, 1.3);
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                      double>::test_teambroadcast(0, 1.3);

    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                      double>::test_teambroadcast(2, 1.3);
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                      double>::test_teambroadcast(2, 1.3);

    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                      double>::test_teambroadcast(16, 1.3);
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                      double>::test_teambroadcast(16, 1.3);

    // FIXME_CUDA
#ifdef KOKKOS_ENABLE_CUDA
    if (!std::is_same<TEST_EXECSPACE, Kokkos::Cuda>::value)
#endif
    // FIXME_HIP
#ifdef KOKKOS_ENABLE_HIP
      if (!std::is_same<TEST_EXECSPACE, Kokkos::Experimental::HIP>::value)
#endif
      {
        TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                          double>::test_teambroadcast(1000, 1.3);
        TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,

                          double>::test_teambroadcast(1000, 1.3);
      }
  }
}

TEST(TEST_CATEGORY, team_handle_by_value) {
  { TestTeamPolicyHandleByValue<TEST_EXECSPACE>(); }
}

}  // namespace Test

#ifndef KOKKOS_ENABLE_OPENMPTARGET
#include <TestTeamVector.hpp>
#endif
#endif
