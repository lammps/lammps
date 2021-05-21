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

#include <iostream>

#include <Kokkos_Core.hpp>

namespace Test {

template <class Space, Kokkos::Experimental::UniqueTokenScope Scope>
class TestUniqueToken {
 public:
  using execution_space = typename Space::execution_space;
  using view_type       = Kokkos::View<int*, execution_space>;

  Kokkos::Experimental::UniqueToken<execution_space, Scope> tokens;

  view_type verify;
  view_type counts;
  view_type errors;

  struct count_test_start_tag {};
  struct count_test_check_tag {};

  KOKKOS_INLINE_FUNCTION
  void operator()(long) const {
    Kokkos::Experimental::AcquireUniqueToken<execution_space, Scope> token_val(
        tokens);
    const int32_t t = token_val.value();

    bool ok = true;

    ok = ok && 0 <= t;
    ok = ok && t < tokens.size();
    ok = ok && 0 == Kokkos::atomic_fetch_add(&verify(t), 1);

    Kokkos::atomic_fetch_add(&counts(t), 1);

    ok = ok && 1 == Kokkos::atomic_fetch_add(&verify(t), -1);

    if (!ok) {
      Kokkos::atomic_fetch_add(&errors(0), 1);
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(count_test_start_tag, long) const {
    constexpr int R = 10;
    int id          = tokens.acquire();
    for (int j = 0; j < R; j++) counts(id)++;
    tokens.release(id);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(count_test_check_tag, long i, int64_t& lsum) const {
    lsum += counts(i);
  }

  TestUniqueToken()
      : tokens(execution_space()),
        verify("TestUniqueTokenVerify", tokens.size()),
        counts("TestUniqueTokenCounts", tokens.size()),
        errors("TestUniqueTokenErrors", 1) {}

  static void run() {
    using policy = Kokkos::RangePolicy<execution_space>;

    TestUniqueToken self;

    {
      const int duplicate = 100;
      const long n        = duplicate * self.tokens.size();

      Kokkos::parallel_for(policy(0, n), self);
      Kokkos::parallel_for(policy(0, n), self);
      Kokkos::parallel_for(policy(0, n), self);
      Kokkos::fence();
    }

    typename view_type::HostMirror host_counts =
        Kokkos::create_mirror_view(self.counts);

    Kokkos::deep_copy(host_counts, self.counts);

    int32_t max = 0;

    {
      const long n = host_counts.extent(0);
      for (long i = 0; i < n; ++i) {
        if (max < host_counts[i]) max = host_counts[i];
      }
    }

    // Count test for pull request #3260
    {
      constexpr int N = 1000000;
      constexpr int R = 10;
      int num         = self.tokens.size();
      Kokkos::resize(self.counts, num);
      Kokkos::deep_copy(self.counts, 0);
      Kokkos::parallel_for(
          "Start", Kokkos::RangePolicy<Space, count_test_start_tag>(0, N),
          self);
      int64_t sum = 0;
      Kokkos::parallel_reduce(
          "Check", Kokkos::RangePolicy<Space, count_test_check_tag>(0, num),
          self, sum);
      ASSERT_EQ(sum, int64_t(N) * R);
    }

    std::cout << "TestUniqueToken max reuse = " << max << std::endl;

    typename view_type::HostMirror host_errors =
        Kokkos::create_mirror_view(self.errors);

    Kokkos::deep_copy(host_errors, self.errors);

    ASSERT_EQ(host_errors(0), 0);
  }
};

TEST(TEST_CATEGORY, unique_token_global) {
  TestUniqueToken<TEST_EXECSPACE,
                  Kokkos::Experimental::UniqueTokenScope::Global>::run();
}

TEST(TEST_CATEGORY, unique_token_instance) {
  TestUniqueToken<TEST_EXECSPACE,
                  Kokkos::Experimental::UniqueTokenScope::Instance>::run();
}

template <class Space>
class TestAcquireTeamUniqueToken {
 public:
  using execution_space = typename Space::execution_space;
  using view_type       = Kokkos::View<int*, execution_space>;
  using scratch_view =
      Kokkos::View<int, typename execution_space::scratch_memory_space,
                   Kokkos::MemoryUnmanaged>;
  using team_policy_type = Kokkos::TeamPolicy<execution_space>;
  using team_member_type = typename team_policy_type::member_type;
  using tokens_type      = Kokkos::Experimental::UniqueToken<execution_space>;

  tokens_type tokens;

  view_type verify;
  view_type counts;
  view_type errors;

  KOKKOS_INLINE_FUNCTION
  void operator()(team_member_type team) const {
    Kokkos::Experimental::AcquireTeamUniqueToken<team_policy_type> token_val(
        tokens, team);
    scratch_view team_rank_0_token_val(team.team_scratch(0));
    const int32_t t = token_val.value();

    bool ok = true;

    ok = ok && 0 <= t;
    ok = ok && t < tokens.size();

    Kokkos::single(Kokkos::PerTeam(team), [&]() {
      ok = ok && 0 == Kokkos::atomic_fetch_add(&verify(t), 1);

      Kokkos::atomic_fetch_add(&counts(t), 1);

      ok = ok && 1 == Kokkos::atomic_fetch_add(&verify(t), -1);
    });

    if (team.team_rank() == 0) {
      team_rank_0_token_val() = t;
    }
    team.team_barrier();
    ok = ok && team_rank_0_token_val() == t;

    if (!ok) {
      Kokkos::atomic_fetch_add(&errors(0), 1);
    }
  }

  TestAcquireTeamUniqueToken(int team_size)
      : tokens(execution_space::concurrency() / team_size, execution_space()),
        verify("TestAcquireTeamUniqueTokenVerify", tokens.size()),
        counts("TestAcquireTeamUniqueTokenCounts", tokens.size()),
        errors("TestAcquireTeamUniqueTokenErrors", 1) {}

  static void run() {
    const int max_team_size = team_policy_type(1, 1).team_size_max(
        TestAcquireTeamUniqueToken(1), Kokkos::ParallelForTag());
    const int team_size = std::min(2, max_team_size);
    TestAcquireTeamUniqueToken self(team_size);

    {
      const int duplicate = 100;
      const long n        = duplicate * self.tokens.size();

      team_policy_type team_policy(n, team_size);
      team_policy.set_scratch_size(
          0, Kokkos::PerTeam(Kokkos::Experimental::AcquireTeamUniqueToken<
                                 team_policy_type>::shmem_size() +
                             scratch_view::shmem_size()));

      Kokkos::parallel_for(team_policy, self);
      Kokkos::fence();
    }

    typename view_type::HostMirror host_counts =
        Kokkos::create_mirror_view(self.counts);

    Kokkos::deep_copy(host_counts, self.counts);

    int32_t max = 0;

    {
      const long n = host_counts.extent(0);
      for (long i = 0; i < n; ++i) {
        if (max < host_counts[i]) max = host_counts[i];
      }
    }

    std::cout << "TestAcquireTeamUniqueToken max reuse = " << max << std::endl;

    typename view_type::HostMirror host_errors =
        Kokkos::create_mirror_view(self.errors);

    Kokkos::deep_copy(host_errors, self.errors);

    ASSERT_EQ(host_errors(0), 0);
  }
};

TEST(TEST_CATEGORY, acquire_team_unique_token) {
  TestAcquireTeamUniqueToken<TEST_EXECSPACE>::run();
}

}  // namespace Test
