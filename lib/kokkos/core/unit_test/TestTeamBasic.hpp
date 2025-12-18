// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_TEAM_BASIC_HPP
#define KOKKOS_TEST_TEAM_BASIC_HPP
#include <TestTeam.hpp>

namespace Test {

TEST(TEST_CATEGORY, team_for) {
  TestTeamPolicy<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>>::test_for(0);
  TestTeamPolicy<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>>::test_for(
      0);

  TestTeamPolicy<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>>::test_for(2);
  TestTeamPolicy<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>>::test_for(
      2);

  TestTeamPolicy<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>>::test_for(
      1000);
  TestTeamPolicy<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>>::test_for(
      1000);
}

// FIXME_OPENMPTARGET wrong results
#ifndef KOKKOS_ENABLE_OPENMPTARGET
TEST(TEST_CATEGORY, team_reduce) {
  TestTeamPolicy<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>>::test_reduce(
      0);
  TestTeamPolicy<TEST_EXECSPACE,
                 Kokkos::Schedule<Kokkos::Dynamic>>::test_reduce(0);
  TestTeamPolicy<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>>::test_reduce(
      2);
  TestTeamPolicy<TEST_EXECSPACE,
                 Kokkos::Schedule<Kokkos::Dynamic>>::test_reduce(2);
  TestTeamPolicy<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>>::test_reduce(
      1000);
  TestTeamPolicy<TEST_EXECSPACE,
                 Kokkos::Schedule<Kokkos::Dynamic>>::test_reduce(1000);
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

/*! \brief Test passing an aggregate to Kokkos::single in a parallel_for with
           team policy
*/
template <typename ExecutionSpace>
struct TestTeamForAggregate {
  using range_policy_t = Kokkos::RangePolicy<ExecutionSpace>;
  using team_policy_t  = Kokkos::TeamPolicy<ExecutionSpace>;
  using member_t       = typename team_policy_t::member_type;
  using memory_space   = typename ExecutionSpace::memory_space;
  using results_type   = Kokkos::View<double*, memory_space>;

  static constexpr double INIT_VALUE   = -1.0;
  static constexpr double EXPECT_VALUE = 1.0;

  struct Agg {
    double d;
  };
  results_type results_;

  TestTeamForAggregate(const size_t size) : results_("results", size) {}
  TestTeamForAggregate() : TestTeamForAggregate(0) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_t& t) const {
    Agg lagg;
    lagg.d = INIT_VALUE;
    Kokkos::single(
        Kokkos::PerTeam(t), [&](Agg& myAgg) { myAgg.d = EXPECT_VALUE; }, lagg);
    size_t i = t.league_rank() * t.team_size() + t.team_rank();
    if (i < results_.size()) {
      results_(i) = lagg.d;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, int& lNumErrs) const {
    if (EXPECT_VALUE != results_(i)) {
      ++lNumErrs;
    }
  }

  static void run() {
    int minTeamSize = 1;
    /* OpenMPTarget hard-codes 32 as the minimum size
       FIXME OPENMPTARGET
    */
#ifdef KOKKOS_ENABLE_OPENMPTARGET
    if constexpr (std::is_same<ExecutionSpace,
                               Kokkos::Experimental::OpenMPTarget>::value) {
      minTeamSize = 32;
    }
#endif

    int maxTeamSize;
    {
      TestTeamForAggregate test;
      maxTeamSize = team_policy_t(1, minTeamSize)
                        .team_size_max(test, Kokkos::ParallelForTag());
    }

    for (int teamSize = minTeamSize; teamSize <= maxTeamSize; teamSize *= 2) {
      for (int problemSize : {1, 100, 10'000, 1'000'000}) {
        const int leagueSize = (problemSize + teamSize - 1) / teamSize;
        TestTeamForAggregate test(problemSize);
        Kokkos::parallel_for(team_policy_t(leagueSize, teamSize), test);
        int numErrs = 0;
        Kokkos::parallel_reduce(range_policy_t(0, problemSize), test, numErrs);
        EXPECT_EQ(numErrs, 0)
            << " teamSize=" << teamSize << " problemSize=" << problemSize;
      }
    }
  }
};

TEST(TEST_CATEGORY, team_parallel_single) {
  TestTeamForAggregate<TEST_EXECSPACE>::run();
}

template <typename ExecutionSpace>
struct LargeTeamScratchFunctor {
  using team_member = typename Kokkos::TeamPolicy<ExecutionSpace>::member_type;
  const size_t m_per_team_bytes;

  KOKKOS_FUNCTION void operator()(const team_member& member) const {
    double* team_shared = static_cast<double*>(
        member.team_scratch(/*level*/ 1).get_shmem(m_per_team_bytes));
    if (team_shared == nullptr)
      Kokkos::abort("Couldn't allocate required size!\n");
    double* team_shared_1 = static_cast<double*>(
        member.team_scratch(/*level*/ 1).get_shmem(sizeof(double)));
    if (team_shared_1 != nullptr)
      Kokkos::abort("Allocated more memory than requested!\n");
  }
};

TEST(TEST_CATEGORY, large_team_scratch_size) {
#ifdef KOKKOS_IMPL_32BIT
  GTEST_SKIP() << "Fails on 32-bit";  // FIXME_32BIT
#endif
  const int level   = 1;
  const int n_teams = 1;

#if defined(KOKKOS_ENABLE_LARGE_MEM_TESTS) || defined(KOKKOS_ENABLE_CUDA) || \
    defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_SYCL)
  //  The GPU backends use unsigned as size_type and we need to check that we
  //  didn't screw the scratch space calculation up, CPU backends anyway use
  //  size_t so less likely to screw up
  const size_t per_team_extent = 502795560;
#else
  const size_t per_team_extent = 268435460;
#endif

  const size_t per_team_bytes = per_team_extent * sizeof(double);

#ifdef KOKKOS_ENABLE_OPENMPTARGET
  Kokkos::TeamPolicy<TEST_EXECSPACE> policy(
      n_teams,
      std::is_same<TEST_EXECSPACE, Kokkos::Experimental::OpenMPTarget>::value
          ? 32
          : 1);
#else
  Kokkos::TeamPolicy<TEST_EXECSPACE> policy(n_teams, 1);
#endif
  policy.set_scratch_size(level, Kokkos::PerTeam(per_team_bytes));

  Kokkos::parallel_for(policy,
                       LargeTeamScratchFunctor<TEST_EXECSPACE>{per_team_bytes});
  Kokkos::fence();
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
  long_wrapper(const long_wrapper& val) : value(val.value) {}

  KOKKOS_FUNCTION
  friend void operator+=(long_wrapper& lhs, const long_wrapper& rhs) {
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
  static_assert(!std::is_arithmetic_v<long_wrapper>);

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
    if (!std::is_same_v<TEST_EXECSPACE, Kokkos::Cuda>)
#endif
    // FIXME_HIP
#ifdef KOKKOS_ENABLE_HIP
      if (!std::is_same_v<TEST_EXECSPACE, Kokkos::HIP>)
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
    if (!std::is_same_v<TEST_EXECSPACE, Kokkos::Cuda>)
#endif
    // FIXME_HIP
#ifdef KOKKOS_ENABLE_HIP
      if (!std::is_same_v<TEST_EXECSPACE, Kokkos::HIP>)
#endif
      {
        TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                          double>::test_teambroadcast(1000, 1.3);
        TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,

                          double>::test_teambroadcast(1000, 1.3);
      }
  }
}

struct TeamBroadcastIntPtrFunctor {
  using TeamMember = typename Kokkos::TeamPolicy<TEST_EXECSPACE>::member_type;

  TeamBroadcastIntPtrFunctor() : view("view") {}

  void run() {
    auto team_size = Kokkos::TeamPolicy<TEST_EXECSPACE>(1, Kokkos::AUTO)
                         .team_size_max(*this, Kokkos::ParallelReduceTag());
    int errors;
    Kokkos::parallel_reduce(Kokkos::TeamPolicy<TEST_EXECSPACE>(1, team_size),
                            *this, errors);
    ASSERT_EQ(errors, 0);
  }

  KOKKOS_FUNCTION void operator()(const TeamMember& team_member,
                                  int& error_count) const {
    int* ptr = team_member.team_rank() == 0 ? view.data() : nullptr;
    team_member.team_broadcast(ptr, 0);
    if (ptr != view.data()) error_count++;
  }

  Kokkos::View<int, TEST_EXECSPACE> view;
};

TEST(TEST_CATEGORY, team_broadcast_int_ptr) {
  TeamBroadcastIntPtrFunctor team_broadcast_functor;
  team_broadcast_functor.run();
}

struct TeamSingleThreadIntPtrFunctor {
  using TeamMember = typename Kokkos::TeamPolicy<TEST_EXECSPACE>::member_type;

  TeamSingleThreadIntPtrFunctor() : view("view") {}

  void run() {
    auto vector_length =
        Kokkos::TeamPolicy<TEST_EXECSPACE>::vector_length_max();
    int errors;
    Kokkos::parallel_reduce(
        Kokkos::TeamPolicy<TEST_EXECSPACE>(1, 1, vector_length), *this, errors);
    ASSERT_EQ(errors, 0);
  }

  KOKKOS_FUNCTION void operator()(const TeamMember& team_member,
                                  int& error_count) const {
    int* ptr = nullptr;
    Kokkos::single(
        Kokkos::PerThread(team_member),
        [&](int*& my_ptr) { my_ptr = view.data(); }, ptr);
    if (ptr != view.data()) error_count++;
  }

  Kokkos::View<int, TEST_EXECSPACE> view;
};

TEST(TEST_CATEGORY, team_single_thread_int_ptr) {
  TeamSingleThreadIntPtrFunctor team_single_thread_functor;
  team_single_thread_functor.run();
}

struct TeamSingleTeamIntPtrFunctor {
  using TeamMember = typename Kokkos::TeamPolicy<TEST_EXECSPACE>::member_type;

  TeamSingleTeamIntPtrFunctor() : view("view") {}

  void run() {
    auto team_size = Kokkos::TeamPolicy<TEST_EXECSPACE>(1, Kokkos::AUTO)
                         .team_size_max(*this, Kokkos::ParallelReduceTag());
    int errors;
    Kokkos::parallel_reduce(Kokkos::TeamPolicy<TEST_EXECSPACE>(1, team_size),
                            *this, errors);
    ASSERT_EQ(errors, 0);
  }

  KOKKOS_FUNCTION void operator()(const TeamMember& team_member,
                                  int& error_count) const {
    int* ptr = nullptr;
    Kokkos::single(
        Kokkos::PerTeam(team_member),
        [&](int*& my_ptr) { my_ptr = view.data(); }, ptr);
    if (ptr != view.data()) error_count++;
  }

  Kokkos::View<int, TEST_EXECSPACE> view;
};

TEST(TEST_CATEGORY, team_single_team_int_ptr) {
  TeamSingleTeamIntPtrFunctor team_single_team_functor;
  team_single_team_functor.run();
}

TEST(TEST_CATEGORY, team_handle_by_value) {
  { TestTeamPolicyHandleByValue<TEST_EXECSPACE>(); }
}

}  // namespace Test

#ifndef KOKKOS_ENABLE_OPENMPTARGET
#include <TestTeamVector.hpp>
#endif
#endif
