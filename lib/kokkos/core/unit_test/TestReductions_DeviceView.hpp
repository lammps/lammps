#include <Kokkos_Core.hpp>

namespace Test {
namespace {

struct TestIsAsynchFunctor {
  Kokkos::View<double, TEST_EXECSPACE> atomic_test;
  TestIsAsynchFunctor(Kokkos::View<double, TEST_EXECSPACE> atomic_test_)
      : atomic_test(atomic_test_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int) const { Kokkos::atomic_add(&atomic_test(), 1.0); }
};

template <class PolicyType, class ReduceFunctor>
void test_reduce_device_view(int64_t N, PolicyType policy,
                             ReduceFunctor functor) {
  using ExecSpace = TEST_EXECSPACE;

  Kokkos::View<int64_t, TEST_EXECSPACE> result("Result");
  Kokkos::View<double, TEST_EXECSPACE> atomic_test("Atomic");
  int64_t reducer_result, view_result, scalar_result;

  Kokkos::Timer timer;

  // Establish whether execspace is asynchronous
  Kokkos::parallel_for("Test::ReduceDeviceView::TestIsAsynch",
                       Kokkos::RangePolicy<TEST_EXECSPACE>(0, 1000000),
                       TestIsAsynchFunctor(atomic_test));
  double time0 = timer.seconds();
  timer.reset();
  typename ExecSpace::execution_space().fence();
  double time_fence0 = timer.seconds();
  Kokkos::deep_copy(result, 0);
  timer.reset();
  bool is_async = time0 < time_fence0;

  // Test Reducer

  Kokkos::parallel_reduce("Test::ReduceDeviceView::TestReducer", policy,
                          functor,
                          Kokkos::Sum<int64_t, TEST_EXECSPACE>(result));
  double time1 = timer.seconds();
  // Check whether it was asyncronous
  timer.reset();
  typename ExecSpace::execution_space().fence();
  double time_fence1 = timer.seconds();
  Kokkos::deep_copy(reducer_result, result);
  Kokkos::deep_copy(result, 0);
  ASSERT_EQ(N, reducer_result);
  timer.reset();

  // Test View
  Kokkos::parallel_reduce("Test::ReduceDeviceView::TestView", policy, functor,
                          result);
  double time2 = timer.seconds();
  // Check whether it was asyncronous
  timer.reset();
  typename ExecSpace::execution_space().fence();
  double time_fence2 = timer.seconds();
  Kokkos::deep_copy(view_result, result);
  Kokkos::deep_copy(result, 0);
  ASSERT_EQ(N, view_result);
  timer.reset();

  // Test Scalar
  Kokkos::parallel_reduce("Test::ReduceDeviceView::TestScalar", policy, functor,
                          scalar_result);
  double time3 = timer.seconds();

  // Check whether it was asyncronous
  timer.reset();
  typename ExecSpace::execution_space().fence();
  double time_fence3 = timer.seconds();

  ASSERT_EQ(N, scalar_result);
  if (is_async) {
    ASSERT_TRUE(time1 < time_fence1);
  }
  if (is_async) {
    ASSERT_TRUE(time2 < time_fence2);
    ASSERT_TRUE(time3 > time_fence3);
  }
}

struct RangePolicyFunctor {
  KOKKOS_INLINE_FUNCTION
  void operator()(const int, int64_t& lsum) const { lsum += 1; }
};

struct MDRangePolicyFunctor {
  KOKKOS_INLINE_FUNCTION
  void operator()(const int, const int, const int, int64_t& lsum) const {
    lsum += 1;
  }
};

struct TeamPolicyFunctor {
  int M;
  TeamPolicyFunctor(int M_) : M(M_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const Kokkos::TeamPolicy<TEST_EXECSPACE>::member_type& team,
                  int64_t& lsum) const {
    for (int i = team.team_rank(); i < M; i += team.team_size()) lsum += 1;
  }
};

}  // namespace

TEST(TEST_CATEGORY, reduce_device_view_range_policy) {
  // Avoid running out of memory
#ifdef KOKKOS_ENABLE_SYCL
  int N = 100 * 1024 * 1024;
#else
  int N = 1000 * 1024 * 1024;
#endif
  test_reduce_device_view(N, Kokkos::RangePolicy<TEST_EXECSPACE>(0, N),
                          RangePolicyFunctor());
}

TEST(TEST_CATEGORY, reduce_device_view_mdrange_policy) {
  int N = 1000 * 1024 * 1024;
  test_reduce_device_view(
      N,
      Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<3>>(
          {0, 0, 0}, {1000, 1024, 1024}),
      MDRangePolicyFunctor());
}

// FIXME_HIP
#ifndef KOKKOS_ENABLE_HIP
TEST(TEST_CATEGORY, reduce_device_view_team_policy) {
// FIXME_SYCL The number of workgroups on CUDA devices can not be larger than
// 65535
#ifdef KOKKOS_ENABLE_SYCL
  int N = 63 * 1024 * 1024;
  test_reduce_device_view(
      N, Kokkos::TeamPolicy<TEST_EXECSPACE>(63 * 1024, Kokkos::AUTO),
      TeamPolicyFunctor(1024));
#else
  int N = 1000 * 1024 * 1024;
  test_reduce_device_view(
      N, Kokkos::TeamPolicy<TEST_EXECSPACE>(1000 * 1024, Kokkos::AUTO),
      TeamPolicyFunctor(1024));
#endif
}
#endif
}  // namespace Test
