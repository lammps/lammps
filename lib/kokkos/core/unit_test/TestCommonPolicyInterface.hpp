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

#include <Kokkos_Core.hpp>

namespace {

template <typename Policy, typename ExpectedExecutionSpace,
          typename ExpectedIndexType, typename ExpectedScheduleType,
          typename ExpectedWorkTag>
constexpr bool compile_time_test() {
  using execution_space = typename Policy::execution_space;
  using index_type      = typename Policy::index_type;
  using schedule_type   = typename Policy::schedule_type;
  using work_tag        = typename Policy::work_tag;

  static_assert(std::is_same_v<execution_space, ExpectedExecutionSpace>);
  static_assert(std::is_same_v<index_type, ExpectedIndexType>);
  static_assert(std::is_same_v<schedule_type, ExpectedScheduleType>);
  static_assert(std::is_same_v<work_tag, ExpectedWorkTag>);

  return true;
}

// Separate class type from class template args so that different
// combinations of template args can be used, while still including
// any necessary templates args (stored in "Args...").
// Example: MDRangePolicy required an iteration pattern be included.
template <template <class...> class PolicyType, class... Args>
constexpr bool test_compile_time_parameters() {
  struct SomeTag {};

  using TestExecSpace    = TEST_EXECSPACE;
  using DefaultExecSpace = Kokkos::DefaultExecutionSpace;
  using TestIndex        = TestExecSpace::size_type;
  using DefaultIndex     = DefaultExecSpace::size_type;
  using LongIndex        = Kokkos::IndexType<long>;
  using StaticSchedule   = Kokkos::Schedule<Kokkos::Static>;
  using DynamicSchedule  = Kokkos::Schedule<Kokkos::Dynamic>;

  // clang-format off
  compile_time_test<PolicyType<                                                            Args...>, DefaultExecSpace, DefaultIndex, StaticSchedule,  void   >();
  compile_time_test<PolicyType<TestExecSpace,                                              Args...>, TestExecSpace,    TestIndex,    StaticSchedule,  void   >();
  compile_time_test<PolicyType<DynamicSchedule,                                            Args...>, DefaultExecSpace, DefaultIndex, DynamicSchedule, void   >();
  compile_time_test<PolicyType<TestExecSpace,   DynamicSchedule,                           Args...>, TestExecSpace,    TestIndex,    DynamicSchedule, void   >();
  compile_time_test<PolicyType<DynamicSchedule, LongIndex,                                 Args...>, DefaultExecSpace, long,         DynamicSchedule, void   >();
  compile_time_test<PolicyType<LongIndex,       DynamicSchedule,                           Args...>, DefaultExecSpace, long,         DynamicSchedule, void   >();
  compile_time_test<PolicyType<TestExecSpace,   DynamicSchedule, LongIndex,                Args...>, TestExecSpace,    long,         DynamicSchedule, void   >();
  compile_time_test<PolicyType<LongIndex,       TestExecSpace,   DynamicSchedule,          Args...>, TestExecSpace,    long,         DynamicSchedule, void   >();
  compile_time_test<PolicyType<DynamicSchedule, LongIndex,       SomeTag,                  Args...>, DefaultExecSpace, long,         DynamicSchedule, SomeTag>();
  compile_time_test<PolicyType<SomeTag,         DynamicSchedule, LongIndex,                Args...>, DefaultExecSpace, long,         DynamicSchedule, SomeTag>();
  compile_time_test<PolicyType<TestExecSpace,   DynamicSchedule, LongIndex, SomeTag,       Args...>, TestExecSpace,    long,         DynamicSchedule, SomeTag>();
  compile_time_test<PolicyType<DynamicSchedule, TestExecSpace,   LongIndex, SomeTag,       Args...>, TestExecSpace,    long,         DynamicSchedule, SomeTag>();
  compile_time_test<PolicyType<SomeTag,         DynamicSchedule, LongIndex, TestExecSpace, Args...>, TestExecSpace,    long,         DynamicSchedule, SomeTag>();
  // clang-format on

  return true;
}

static_assert(test_compile_time_parameters<Kokkos::RangePolicy>());
static_assert(test_compile_time_parameters<Kokkos::TeamPolicy>());
static_assert(
    test_compile_time_parameters<Kokkos::MDRangePolicy, Kokkos::Rank<2>>());

// Asserts that worktag conversion works properly.
template <class Policy>
constexpr bool test_worktag() {
  struct WorkTag1 {};
  struct WorkTag2 {};

  // Apply WorkTag1
  using PolicyWithWorkTag1 =
      Kokkos::Impl::WorkTagTrait::policy_with_trait<Policy, WorkTag1>;
  // Swap for WorkTag2
  using PolicyWithWorkTag2 =
      Kokkos::Impl::WorkTagTrait::policy_with_trait<PolicyWithWorkTag1,
                                                    WorkTag2>;

  static_assert(std::is_void_v<typename Policy::work_tag>);
  static_assert(
      std::is_same_v<typename PolicyWithWorkTag1::work_tag, WorkTag1>);
  static_assert(
      std::is_same_v<typename PolicyWithWorkTag2::work_tag, WorkTag2>);

  // Currently not possible to remove the work tag from a policy.
  // Uncomment the line below to see the compile error.
  // using PolicyRemoveWorkTag =
  // Kokkos::Impl::WorkTagTrait::policy_with_trait<PolicyWithWorkTag2, void>;
  // static_assert(std::is_void_v<PolicyRemoveWorkTag::work_tag>);

  return true;
}

static_assert(test_worktag<Kokkos::RangePolicy<>>());
static_assert(test_worktag<Kokkos::TeamPolicy<>>());
static_assert(test_worktag<Kokkos::MDRangePolicy<Kokkos::Rank<2>>>());

}  // namespace
