// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <sstream>
#include <iostream>
#include <limits>

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include <TestNonTrivialScalarTypes.hpp>

//--------------------------------------------------------------------------

namespace Test {
struct MyPair : Kokkos::pair<int, int> {};
}  // namespace Test

template <>
struct Kokkos::reduction_identity<Test::MyPair> {
  KOKKOS_FUNCTION static Test::MyPair min() {
    return Test::MyPair{{INT_MAX, INT_MAX}};
  }
};

namespace Test {

struct ReducerTag {};

template <class Scalar, class ExecSpace = Kokkos::DefaultExecutionSpace>
struct TestReducers {
  struct SumFunctor {
    Kokkos::View<const Scalar*, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(const int& i, Scalar& value) const { value += values(i); }
  };

  struct TeamSumFunctor {
    using member_type = typename Kokkos::TeamPolicy<ExecSpace>::member_type;

    KOKKOS_INLINE_FUNCTION
    void operator()(const member_type& m, Scalar& value) const {
      if (m.team_rank() == m.team_size() - 1) value += Scalar(1);
    }
  };

  struct TeamSumNestedFunctor {
    using member_type = typename Kokkos::TeamPolicy<ExecSpace>::member_type;

    SumFunctor f;
    int M, N;
    Kokkos::View<Scalar*, ExecSpace> result;

    TeamSumNestedFunctor(SumFunctor& f_, const int M_, const int N_,
                         Kokkos::View<Scalar*, ExecSpace> result_)
        : f(f_), M(M_), N(N_), result(result_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const member_type& m) const {
      const int i = m.league_rank();
      Scalar local_scalar;
      Kokkos::Sum<Scalar, typename ExecSpace::memory_space> reducer_scalar(
          local_scalar);
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(m, N), f, reducer_scalar);
      result(i) = local_scalar;
    }
  };

  struct ProdFunctor {
    Kokkos::View<const Scalar*, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(const int& i, Scalar& value) const { value *= values(i); }
  };

  struct MinFunctor {
    Kokkos::View<const Scalar*, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(const int& i, Scalar& value) const {
      if (values(i) < value) value = values(i);
    }
  };

  struct MaxFunctor {
    Kokkos::View<const Scalar*, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(const int& i, Scalar& value) const {
      if (values(i) > value) value = values(i);
    }
  };

  struct MinLocFunctor {
    Kokkos::View<const Scalar*, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(
        const int& i,
        typename Kokkos::MinLoc<Scalar, int>::value_type& value) const {
      if (values(i) < value.val) {
        value.val = values(i);
        value.loc = i;
      }
    }
  };

  struct MinLocFunctor2D {
    Kokkos::View<const Scalar**, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(
        const int& i, const int& j,
        typename Kokkos::MinLoc<Scalar, MyPair>::value_type& value) const {
      if (values(i, j) < value.val) {
        value.val = values(i, j);
        value.loc = {{i, j}};
      }
    }
  };

  struct MaxLocFunctor {
    Kokkos::View<const Scalar*, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(
        const int& i,
        typename Kokkos::MaxLoc<Scalar, int>::value_type& value) const {
      if (values(i) > value.val) {
        value.val = values(i);
        value.loc = i;
      }
    }
  };

  struct MaxLocFunctor2D {
    Kokkos::View<const Scalar**, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(
        const int& i, const int& j,
        typename Kokkos::MaxLoc<Scalar, MyPair>::value_type& value) const {
      if (values(i, j) > value.val) {
        value.val = values(i, j);
        value.loc = {{i, j}};
      }
    }
  };

  struct MinMaxLocFunctor {
    Kokkos::View<const Scalar*, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(
        const int& i,
        typename Kokkos::MinMaxLoc<Scalar, int>::value_type& value) const {
      if (values(i) > value.max_val) {
        value.max_val = values(i);
        value.max_loc = i;
      }

      if (values(i) < value.min_val) {
        value.min_val = values(i);
        value.min_loc = i;
      }
    }
  };

  struct MinMaxLocFunctor2D {
    Kokkos::View<const Scalar**, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(
        const int& i, const int& j,
        typename Kokkos::MinMaxLoc<Scalar, MyPair>::value_type& value) const {
      if (values(i, j) > value.max_val) {
        value.max_val = values(i, j);
        value.max_loc = {{i, j}};
      }

      if (values(i, j) < value.min_val) {
        value.min_val = values(i, j);
        value.min_loc = {{i, j}};
      }
    }
  };

  struct BAndFunctor {
    Kokkos::View<const Scalar*, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(const int& i, Scalar& value) const {
      value = value & values(i);
    }
  };

  struct BOrFunctor {
    Kokkos::View<const Scalar*, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(const int& i, Scalar& value) const {
      value = value | values(i);
    }
  };

  struct LAndFunctor {
    Kokkos::View<const Scalar*, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(const int& i, Scalar& value) const {
      value = value && values(i);
    }
  };

  struct LOrFunctor {
    Kokkos::View<const Scalar*, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(const int& i, Scalar& value) const {
      value = value || values(i);
    }
  };

  struct SumFunctorTag {
    Kokkos::View<const Scalar*, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(const ReducerTag, const int& i, Scalar& value) const {
      value += values(i);
    }
  };

  struct ProdFunctorTag {
    Kokkos::View<const Scalar*, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(const ReducerTag, const int& i, Scalar& value) const {
      value *= values(i);
    }
  };

  struct MinFunctorTag {
    Kokkos::View<const Scalar*, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(const ReducerTag, const int& i, Scalar& value) const {
      if (values(i) < value) value = values(i);
    }
  };

  struct MaxFunctorTag {
    Kokkos::View<const Scalar*, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(const ReducerTag, const int& i, Scalar& value) const {
      if (values(i) > value) value = values(i);
    }
  };

  struct MinLocFunctorTag {
    Kokkos::View<const Scalar*, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(
        const ReducerTag, const int& i,
        typename Kokkos::MinLoc<Scalar, int>::value_type& value) const {
      if (values(i) < value.val) {
        value.val = values(i);
        value.loc = i;
      }
    }
  };

  struct MaxLocFunctorTag {
    Kokkos::View<const Scalar*, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(
        const ReducerTag, const int& i,
        typename Kokkos::MaxLoc<Scalar, int>::value_type& value) const {
      if (values(i) > value.val) {
        value.val = values(i);
        value.loc = i;
      }
    }
  };

  struct MinMaxLocFunctorTag {
    Kokkos::View<const Scalar*, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(
        const ReducerTag, const int& i,
        typename Kokkos::MinMaxLoc<Scalar, int>::value_type& value) const {
      if (values(i) > value.max_val) {
        value.max_val = values(i);
        value.max_loc = i;
      }

      if (values(i) < value.min_val) {
        value.min_val = values(i);
        value.min_loc = i;
      }
    }
  };

  struct BAndFunctorTag {
    Kokkos::View<const Scalar*, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(const ReducerTag, const int& i, Scalar& value) const {
      value = value & values(i);
    }
  };

  struct BOrFunctorTag {
    Kokkos::View<const Scalar*, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(const ReducerTag, const int& i, Scalar& value) const {
      value = value | values(i);
    }
  };

  struct LAndFunctorTag {
    Kokkos::View<const Scalar*, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(const ReducerTag, const int& i, Scalar& value) const {
      value = value && values(i);
    }
  };

  struct LOrFunctorTag {
    Kokkos::View<const Scalar*, ExecSpace> values;

    KOKKOS_INLINE_FUNCTION
    void operator()(const ReducerTag, const int& i, Scalar& value) const {
      value = value || values(i);
    }
  };

  // get number of teams for TeamPolicy depending on the tested type
  constexpr static int get_num_teams() {
    if constexpr (sizeof(Scalar) == 1) {
      return 126;
    } else if constexpr (std::is_same_v<Scalar,
                                        Kokkos::Experimental::bhalf_t>) {
      return 256;
    }

    return 1024;
  }

  static void test_sum_team_policy(int N, SumFunctor f, Scalar reference_sum) {
    Scalar sum_scalar;
    Kokkos::View<Scalar, ExecSpace> sum_view("result");
    Kokkos::deep_copy(sum_view, Scalar(1));

    // Test team policy reduction
    {
      constexpr int num_teams = get_num_teams();
      TeamSumFunctor tf;
      // FIXME_OPENMPTARGET temporary restriction for team size to be at least
      // 32
#ifdef KOKKOS_ENABLE_OPENMPTARGET
      int team_size =
          std::is_same<ExecSpace, Kokkos::Experimental::OpenMPTarget>::value
              ? 32
              : 1;
#else
      int team_size         = 1;
#endif
      auto team_pol = Kokkos::TeamPolicy<ExecSpace>(num_teams, team_size);
      Kokkos::parallel_reduce(team_pol, tf, sum_view);
      Kokkos::deep_copy(sum_scalar, sum_view);
      ASSERT_EQ(sum_scalar, Scalar{num_teams}) << "num_teams: " << num_teams;
    }

    // Test TeamThreadRange level reduction with 0 work produces 0 result
    {
      const int league_size = 1;
      Kokkos::View<Scalar*, ExecSpace> result("result", league_size);
      TeamSumNestedFunctor tnf(f, league_size, 0, result);
      // FIXME_OPENMPTARGET temporary restriction for team size to be at least
      // 32
#ifdef KOKKOS_ENABLE_OPENMPTARGET
      int team_size =
          std::is_same<ExecSpace, Kokkos::Experimental::OpenMPTarget>::value
              ? 32
              : 1;
#else
      int team_size         = 1;
#endif
      auto team_pol = Kokkos::TeamPolicy<ExecSpace>(1, team_size);
      Kokkos::parallel_for(team_pol, tnf);
      auto result_h =
          Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), result);
      ASSERT_EQ(result_h(0), Scalar{0}) << "N: " << N;
    }

    // Same test as above, but with inner reduction over N, and league_size=10
    {
      const int league_size = 10;
      Kokkos::View<Scalar*, ExecSpace> result("result", league_size);
      TeamSumNestedFunctor tnf(f, league_size, N, result);
      // FIXME_OPENMPTARGET temporary restriction for team size to be at least
      // 32
#ifdef KOKKOS_ENABLE_OPENMPTARGET
      int initial_team_size =
          std::is_same_v<ExecSpace, Kokkos::Experimental::OpenMPTarget> ? 32
                                                                        : 1;
#else
      int initial_team_size = 1;
#endif
      auto team_size_max =
          Kokkos::TeamPolicy<ExecSpace>(league_size, initial_team_size)
              .team_size_max(tnf, Kokkos::ParallelForTag());
      auto team_size = std::min(team_size_max, TEST_EXECSPACE().concurrency());
      auto team_pol  = Kokkos::TeamPolicy<ExecSpace>(league_size, team_size);
      Kokkos::parallel_for(team_pol, tnf);
      auto result_h =
          Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), result);
      for (int i = 0; i < result_h.extent_int(0); ++i) {
        ASSERT_EQ(result_h(i), reference_sum) << "N: " << N;
      }
    }
  }

  // Test that reducers return correct results with LaunchBounds value smaller
  // than 32.
  template <int N>
  static void test_launch_bounds() {
    Kokkos::View<Scalar*> v("", N);
    Kokkos::deep_copy(v, Scalar(2));

    // Functor
    auto suml  = KOKKOS_LAMBDA(const int i, Scalar& ret) { ret += v(i); };
    auto prodl = KOKKOS_LAMBDA(const int i, Scalar& ret) { ret *= v(i); };

    // This first reduction is necessary to ensure we trigger the bug #8443
    Scalar ret;
    Kokkos::parallel_reduce(Kokkos::RangePolicy(0, 1), suml,
                            Kokkos::Sum<Scalar>(ret));
    EXPECT_EQ(ret, 2) << "N=" << N;

    // Test that LaunchBounds<N> works with Sum
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Kokkos::LaunchBounds<N>>(0, N),
                            suml, Kokkos::Sum<Scalar>(ret));
    EXPECT_EQ(ret, 2 * N) << "N=" << N;

    // This test can't be run on int32 with N>= 31 because of overflow
    if constexpr (!(std::is_integral_v<Scalar> && sizeof(Scalar) <= 32) ||
                  N < 31) {
      // Test that LaunchBounds<N> works with Prod
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<Kokkos::LaunchBounds<N>>(0, N), prodl,
          Kokkos::Prod<Scalar>(ret));

      Scalar expected = Scalar(1ull << N);
      EXPECT_EQ(ret, expected) << "N=" << N;
    }
  }

  static void test_launch_bounds() {
    test_launch_bounds<1>();
    test_launch_bounds<5>();
    test_launch_bounds<31>();
  }

  static void test_sum(int N) {
    Kokkos::View<Scalar*, ExecSpace> values("Values", N);
    auto h_values        = Kokkos::create_mirror_view(values);
    Scalar reference_sum = 0;

// FIXME spurious warnings like
// error: 'SR.14123' may be used uninitialized [-Werror=maybe-uninitialized]
#if defined(KOKKOS_COMPILER_GNU) && KOKKOS_COMPILER_GNU >= 1500
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

    for (int i = 0; i < N; i++) {
      int denom = sizeof(Scalar) <= 2 ? 10 : 100;
      // clang-format off
      // For bhalf, we start overflowing integer values at 2^8.
      //            after 2^8,  we lose representation of odd numbers;
      //            after 2^9,  we lose representation of odd and even numbers in position 1.
      //            after 2^10, we lose representation of odd and even numbers in position 1-3.
      //            after 2^11, we lose representation of odd and even numbers in position 1-7.
      //            ...
      // Generally, for IEEE 754 floating point numbers, we start this overflow pattern at: 2^(num_fraction_bits+1).
      // brain float has num_fraction_bits = 7.
      // This mask addresses #4719 for N <= 51.
      // The mask is not needed for N <= 25.
      // clang-format on
      int mask = std::is_same_v<Scalar, Kokkos::Experimental::bhalf_t> && N > 25
                     ? (int)0xfffffffe
                     : (int)0xffffffff;
      h_values(i) = (Scalar)((rand() % denom) & mask);
      reference_sum += h_values(i);
    }

// FIXME spurious warnings like
// error: 'SR.14123' may be used uninitialized [-Werror=maybe-uninitialized]
#if defined(KOKKOS_COMPILER_GNU) && KOKKOS_COMPILER_GNU >= 1500
#pragma GCC diagnostic pop
#endif

    Kokkos::deep_copy(values, h_values);

    SumFunctor f;
    f.values = values;
    SumFunctorTag f_tag;
    f_tag.values = values;
    Scalar init  = 0;

    {
      Scalar sum_scalar = Scalar(1);
      Kokkos::Sum<Scalar> reducer_scalar(sum_scalar);

      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, 0), f,
                              reducer_scalar);
      ASSERT_EQ(sum_scalar, init) << "N: " << N;

      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_scalar);
      ASSERT_EQ(sum_scalar, reference_sum) << "N: " << N;

      sum_scalar = init;
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace, ReducerTag>(0, N),
                              f_tag, reducer_scalar);
      ASSERT_EQ(sum_scalar, reference_sum) << "N: " << N;

      Scalar sum_scalar_view = reducer_scalar.reference();
      ASSERT_EQ(sum_scalar_view, reference_sum) << "N: " << N;
    }

// FIXME_OPENACC - custom reduction with TeamPolicy is not yet implemented.
#ifdef KOKKOS_ENABLE_OPENACC
    if constexpr (!std::is_same_v<ExecSpace, Kokkos::Experimental::OpenACC>)
      test_sum_team_policy(N, f, reference_sum);
#endif

    {
      Kokkos::View<Scalar, Kokkos::HostSpace> sum_view("View");
      sum_view() = Scalar(1);
      Kokkos::Sum<Scalar> reducer_view(sum_view);
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, 0), f,
                              reducer_view);
      Kokkos::fence();
      Scalar sum_view_scalar = sum_view();
      ASSERT_EQ(sum_view_scalar, init) << "N: " << N;

      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_view);
      Kokkos::fence();
      sum_view_scalar = sum_view();
      ASSERT_EQ(sum_view_scalar, reference_sum) << "N: " << N;

      Scalar sum_view_view = reducer_view.reference();
      ASSERT_EQ(sum_view_view, reference_sum) << "N: " << N;
    }

    {
      Kokkos::View<Scalar, typename ExecSpace::memory_space> sum_view("View");
      Kokkos::deep_copy(sum_view, Scalar(1));
      Kokkos::Sum<Scalar, typename ExecSpace::memory_space> reducer_view(
          sum_view);
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, 0), f,
                              reducer_view);
      Kokkos::fence();
      Scalar sum_view_scalar;
      Kokkos::deep_copy(sum_view_scalar, sum_view);
      ASSERT_EQ(sum_view_scalar, init) << "N: " << N;

      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_view);
      Kokkos::fence();
      Kokkos::deep_copy(sum_view_scalar, sum_view);
      ASSERT_EQ(sum_view_scalar, reference_sum) << "N: " << N;
    }
  }

  static void test_prod(int N) {
    Kokkos::View<Scalar*, ExecSpace> values("Values", N);
    auto h_values         = Kokkos::create_mirror_view(values);
    Scalar reference_prod = 1;

    for (int i = 0; i < N; i++) {
      h_values(i) = (Scalar)(rand() % 4 + 1);
      reference_prod *= h_values(i);
    }
    Kokkos::deep_copy(values, h_values);

    ProdFunctor f;
    f.values = values;
    ProdFunctorTag f_tag;
    f_tag.values = values;
    Scalar init  = 1;

    {
      Scalar prod_scalar = Scalar(0);
      Kokkos::Prod<Scalar> reducer_scalar(prod_scalar);
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, 0), f,
                              reducer_scalar);
      ASSERT_EQ(prod_scalar, init);

      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_scalar);
      ASSERT_EQ(prod_scalar, reference_prod);

      prod_scalar = init;
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace, ReducerTag>(0, N),
                              f_tag, reducer_scalar);
      ASSERT_EQ(prod_scalar, reference_prod);

      Scalar prod_scalar_view = reducer_scalar.reference();
      ASSERT_EQ(prod_scalar_view, reference_prod);
    }

    {
      Kokkos::View<Scalar, Kokkos::HostSpace> prod_view("View");
      prod_view() = Scalar(0);
      Kokkos::Prod<Scalar> reducer_view(prod_view);
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, 0), f,
                              reducer_view);
      Kokkos::fence();
      Scalar prod_view_scalar = prod_view();
      ASSERT_EQ(prod_view_scalar, init);

      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_view);
      Kokkos::fence();
      prod_view_scalar = prod_view();
      ASSERT_EQ(prod_view_scalar, reference_prod);

      Scalar prod_view_view = reducer_view.reference();
      ASSERT_EQ(prod_view_view, reference_prod);
    }

    {
      Kokkos::View<Scalar, typename ExecSpace::memory_space> prod_view("View");
      Kokkos::deep_copy(prod_view, Scalar(0));
      Kokkos::Prod<Scalar, typename ExecSpace::memory_space> reducer_view(
          prod_view);
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, 0), f,
                              reducer_view);
      Kokkos::fence();
      Scalar prod_view_scalar;
      Kokkos::deep_copy(prod_view_scalar, prod_view);
      ASSERT_EQ(prod_view_scalar, init);

      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_view);
      Kokkos::fence();
      Kokkos::deep_copy(prod_view_scalar, prod_view);
      ASSERT_EQ(prod_view_scalar, reference_prod);
    }
  }

  static void test_min(int N) {
    Kokkos::View<Scalar*, ExecSpace> values("Values", N);
    auto h_values        = Kokkos::create_mirror_view(values);
    Scalar reference_min = std::numeric_limits<Scalar>::max();

    for (int i = 0; i < N; i++) {
      h_values(i) = (Scalar)(rand() % 100000);

      if (h_values(i) < reference_min) reference_min = h_values(i);
    }
    Kokkos::deep_copy(values, h_values);

    MinFunctor f;
    f.values = values;
    MinFunctorTag f_tag;
    f_tag.values = values;
    Scalar init  = std::numeric_limits<Scalar>::max();

    {
      Scalar min_scalar = init;
      Kokkos::Min<Scalar> reducer_scalar(min_scalar);

      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_scalar);
      ASSERT_EQ(min_scalar, reference_min);

      min_scalar = init;
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace, ReducerTag>(0, N),
                              f_tag, reducer_scalar);
      ASSERT_EQ(min_scalar, reference_min);

      Scalar min_scalar_view = reducer_scalar.reference();
      ASSERT_EQ(min_scalar_view, reference_min);
    }

    {
      Kokkos::View<Scalar, Kokkos::HostSpace> min_view("View");
      min_view() = init;
      Kokkos::Min<Scalar> reducer_view(min_view);
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_view);
      Kokkos::fence();

      Scalar min_view_scalar = min_view();
      ASSERT_EQ(min_view_scalar, reference_min);

      Scalar min_view_view = reducer_view.reference();
      ASSERT_EQ(min_view_view, reference_min);
    }
  }

  static void test_max(int N) {
    Kokkos::View<Scalar*, ExecSpace> values("Values", N);
    auto h_values        = Kokkos::create_mirror_view(values);
    Scalar reference_max = std::numeric_limits<Scalar>::min();

    for (int i = 0; i < N; i++) {
      h_values(i) = (Scalar)(rand() % 100000 + 1);

      if (h_values(i) > reference_max) reference_max = h_values(i);
    }
    Kokkos::deep_copy(values, h_values);

    MaxFunctor f;
    f.values = values;
    MaxFunctorTag f_tag;
    f_tag.values = values;
    Scalar init  = std::numeric_limits<Scalar>::min();

    {
      Scalar max_scalar = init;
      Kokkos::Max<Scalar> reducer_scalar(max_scalar);

      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_scalar);
      ASSERT_EQ(max_scalar, reference_max);

      max_scalar = init;
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace, ReducerTag>(0, N),
                              f_tag, reducer_scalar);
      ASSERT_EQ(max_scalar, reference_max);

      Scalar max_scalar_view = reducer_scalar.reference();
      ASSERT_EQ(max_scalar_view, reference_max);
    }

    {
      Kokkos::View<Scalar, Kokkos::HostSpace> max_view("View");
      max_view() = init;
      Kokkos::Max<Scalar> reducer_view(max_view);
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_view);
      Kokkos::fence();

      Scalar max_view_scalar = max_view();
      ASSERT_EQ(max_view_scalar, reference_max);

      Scalar max_view_view = reducer_view.reference();
      ASSERT_EQ(max_view_view, reference_max);
    }
  }

  static void test_minloc(int N) {
    using value_type = typename Kokkos::MinLoc<Scalar, int>::value_type;

    Kokkos::View<Scalar*, ExecSpace> values("Values", N);
    auto h_values        = Kokkos::create_mirror_view(values);
    Scalar reference_min = std::numeric_limits<Scalar>::max();
    int reference_loc    = -1;

    for (int i = 0; i < N; i++) {
      h_values(i) = (Scalar)(rand() % 100000 + 2);

      if (h_values(i) < reference_min) {
        reference_min = h_values(i);
        reference_loc = i;
      } else if (h_values(i) == reference_min) {
        // Make min unique.
        h_values(i) += Scalar(1);
      }
    }
    Kokkos::deep_copy(values, h_values);

    MinLocFunctor f;
    f.values = values;
    MinLocFunctorTag f_tag;
    f_tag.values = values;

    {
      value_type min_scalar;
      Kokkos::MinLoc<Scalar, int> reducer_scalar(min_scalar);

      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_scalar);
      ASSERT_EQ(min_scalar.val, reference_min);
      ASSERT_EQ(min_scalar.loc, reference_loc);

      min_scalar = value_type();
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace, ReducerTag>(0, N),
                              f_tag, reducer_scalar);
      ASSERT_EQ(min_scalar.val, reference_min);
      ASSERT_EQ(min_scalar.loc, reference_loc);

      value_type min_scalar_view = reducer_scalar.reference();
      ASSERT_EQ(min_scalar_view.val, reference_min);
      ASSERT_EQ(min_scalar_view.loc, reference_loc);
    }

    {
      Kokkos::View<value_type, Kokkos::HostSpace> min_view("View");
      Kokkos::MinLoc<Scalar, int> reducer_view(min_view);
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_view);
      Kokkos::fence();

      value_type min_view_scalar = min_view();
      ASSERT_EQ(min_view_scalar.val, reference_min);
      ASSERT_EQ(min_view_scalar.loc, reference_loc);

      value_type min_view_view = reducer_view.reference();
      ASSERT_EQ(min_view_view.val, reference_min);
      ASSERT_EQ(min_view_view.loc, reference_loc);
    }
  }

  static void test_minloc_2d(int N) {
    using reducer_type = Kokkos::MinLoc<Scalar, MyPair>;
    using value_type   = typename reducer_type::value_type;

    Kokkos::View<Scalar**, ExecSpace> values("Values", N, N);
    auto h_values        = Kokkos::create_mirror_view(values);
    Scalar reference_min = std::numeric_limits<Scalar>::max();
    MyPair reference_loc = {{-1, -1}};

    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) {
        h_values(i, j) = (Scalar)(rand() % 100000 + 2);

        if (h_values(i, j) < reference_min) {
          reference_min = h_values(i, j);
          reference_loc = {{i, j}};
        } else if (h_values(i, j) == reference_min) {
          // Make min unique.
          h_values(i, j) += Scalar(1);
        }
      }
    Kokkos::deep_copy(values, h_values);

    MinLocFunctor2D f;
    f.values = values;

    {
      value_type min_scalar;
      reducer_type reducer_scalar(min_scalar);

      Kokkos::parallel_reduce(
          Kokkos::MDRangePolicy<Kokkos::Rank<2>, ExecSpace>({0, 0}, {N, N}), f,
          reducer_scalar);
      ASSERT_EQ(min_scalar.val, reference_min);
      ASSERT_EQ(min_scalar.loc, reference_loc);
    }
  }

  static void test_minloc_loc_init(int N) {
    using reducer_type       = Kokkos::MinLoc<Scalar, int>;
    using reducer_value_type = typename reducer_type::value_type;

    Kokkos::View<Scalar*, ExecSpace> values("Values", N);
    auto h_values = Kokkos::create_mirror_view(values);

    for (int i = 0; i < N; ++i) {
      h_values(i) = Kokkos::reduction_identity<Scalar>::min();
    }
    Kokkos::deep_copy(values, h_values);

    reducer_value_type value_loc{0, -1};

    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<ExecSpace>(0, N),
        KOKKOS_LAMBDA(const int i, reducer_value_type& update) {
          auto x = values(i);
          if (i % 2 == 0)
            return;
          else if (x <= update.val) {
            update.val = x;
            update.loc = i;
          }
        },
        reducer_type(value_loc));

    ASSERT_EQ(value_loc.val, h_values(0));
    ASSERT_GE(value_loc.loc, 0);
    ASSERT_LT(value_loc.loc, N);
  }

  static void test_maxloc(int N) {
    using value_type = typename Kokkos::MaxLoc<Scalar, int>::value_type;

    Kokkos::View<Scalar*, ExecSpace> values("Values", N);
    auto h_values        = Kokkos::create_mirror_view(values);
    Scalar reference_max = std::numeric_limits<Scalar>::min();
    int reference_loc    = -1;

    for (int i = 0; i < N; i++) {
      h_values(i) = (Scalar)(rand() % 100000 + 2);

      if (h_values(i) > reference_max) {
        reference_max = h_values(i);
        reference_loc = i;
      } else if (h_values(i) == reference_max) {
        // Make max unique.
        h_values(i) -= Scalar(1);
      }
    }
    Kokkos::deep_copy(values, h_values);

    MaxLocFunctor f;
    f.values = values;
    MaxLocFunctorTag f_tag;
    f_tag.values = values;

    {
      value_type max_scalar;
      Kokkos::MaxLoc<Scalar, int> reducer_scalar(max_scalar);

      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_scalar);
      ASSERT_EQ(max_scalar.val, reference_max);
      ASSERT_EQ(max_scalar.loc, reference_loc);

      max_scalar = value_type();
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace, ReducerTag>(0, N),
                              f_tag, reducer_scalar);
      ASSERT_EQ(max_scalar.val, reference_max);
      ASSERT_EQ(max_scalar.loc, reference_loc);

      value_type max_scalar_view = reducer_scalar.reference();
      ASSERT_EQ(max_scalar_view.val, reference_max);
      ASSERT_EQ(max_scalar_view.loc, reference_loc);
    }

    {
      Kokkos::View<value_type, Kokkos::HostSpace> max_view("View");
      Kokkos::MaxLoc<Scalar, int> reducer_view(max_view);
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_view);
      Kokkos::fence();

      value_type max_view_scalar = max_view();
      ASSERT_EQ(max_view_scalar.val, reference_max);
      ASSERT_EQ(max_view_scalar.loc, reference_loc);

      value_type max_view_view = reducer_view.reference();
      ASSERT_EQ(max_view_view.val, reference_max);
      ASSERT_EQ(max_view_view.loc, reference_loc);
    }
  }

  static void test_maxloc_2d(int N) {
    using reducer_type = Kokkos::MaxLoc<Scalar, MyPair>;
    using value_type   = typename reducer_type::value_type;

    Kokkos::View<Scalar**, ExecSpace> values("Values", N, N);
    auto h_values        = Kokkos::create_mirror_view(values);
    Scalar reference_max = std::numeric_limits<Scalar>::min();
    MyPair reference_loc = {{-1, -1}};

    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j) {
        h_values(i, j) = (Scalar)(rand() % 100000 + 2);

        if (h_values(i, j) > reference_max) {
          reference_max = h_values(i, j);
          reference_loc = {{i, j}};
        } else if (h_values(i, j) == reference_max) {
          // Make max unique.
          h_values(i, j) -= Scalar(1);
        }
      }
    Kokkos::deep_copy(values, h_values);

    MaxLocFunctor2D f;
    f.values = values;

    {
      value_type max_scalar;
      reducer_type reducer_scalar(max_scalar);

      Kokkos::parallel_reduce(
          Kokkos::MDRangePolicy<Kokkos::Rank<2>, ExecSpace>({0, 0}, {N, N}), f,
          reducer_scalar);
      ASSERT_EQ(max_scalar.val, reference_max);
      ASSERT_EQ(max_scalar.loc, reference_loc);
    }
  }

  static void test_maxloc_loc_init(int N) {
    using reducer_type       = Kokkos::MaxLoc<Scalar, int>;
    using reducer_value_type = typename reducer_type::value_type;

    Kokkos::View<Scalar*, ExecSpace> values("Values", N);
    auto h_values = Kokkos::create_mirror_view(values);

    for (int i = 0; i < N; ++i) {
      h_values(i) = Kokkos::reduction_identity<Scalar>::max();
    }
    Kokkos::deep_copy(values, h_values);

    reducer_value_type value_loc{0, -1};

    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<ExecSpace>(0, N),
        KOKKOS_LAMBDA(const int i, reducer_value_type& update) {
          auto x = values(i);
          if (i % 2 == 0)
            return;
          else if (x >= update.val) {
            update.val = x;
            update.loc = i;
          }
        },
        reducer_type(value_loc));

    ASSERT_EQ(value_loc.val, h_values(0));
    ASSERT_GE(value_loc.loc, 0);
    ASSERT_LT(value_loc.loc, N);
  }

  static void test_minmaxloc(int N) {
    using value_type = typename Kokkos::MinMaxLoc<Scalar, int>::value_type;

    Kokkos::View<Scalar*, ExecSpace> values("Values", N);
    auto h_values        = Kokkos::create_mirror_view(values);
    Scalar reference_max = std::numeric_limits<Scalar>::min();
    Scalar reference_min = std::numeric_limits<Scalar>::max();
    int reference_minloc = -1;
    int reference_maxloc = -1;

    for (int i = 0; i < N; i++) {
      h_values(i) = (Scalar)(rand() % 100000 + 2);
    }

    for (int i = 0; i < N; i++) {
      if (h_values(i) > reference_max) {
        reference_max    = h_values(i);
        reference_maxloc = i;
      } else if (h_values(i) == reference_max) {
        // Make max unique.
        h_values(i) -= Scalar(1);
      }
    }

    for (int i = 0; i < N; i++) {
      if (h_values(i) < reference_min) {
        reference_min    = h_values(i);
        reference_minloc = i;
      } else if (h_values(i) == reference_min) {
        // Make min unique.
        h_values(i) += Scalar(1);
      }
    }

    Kokkos::deep_copy(values, h_values);

    MinMaxLocFunctor f;
    f.values = values;
    MinMaxLocFunctorTag f_tag;
    f_tag.values = values;

    {
      value_type minmax_scalar;
      Kokkos::MinMaxLoc<Scalar, int> reducer_scalar(minmax_scalar);

      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_scalar);
      ASSERT_EQ(minmax_scalar.min_val, reference_min);

      for (int i = 0; i < N; i++) {
        if ((i == minmax_scalar.min_loc) && (h_values(i) == reference_min)) {
          reference_minloc = i;
        }
      }

      ASSERT_EQ(minmax_scalar.min_loc, reference_minloc);
      ASSERT_EQ(minmax_scalar.max_val, reference_max);

      for (int i = 0; i < N; i++) {
        if ((i == minmax_scalar.max_loc) && (h_values(i) == reference_max)) {
          reference_maxloc = i;
        }
      }

      ASSERT_EQ(minmax_scalar.max_loc, reference_maxloc);

      minmax_scalar = value_type();
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace, ReducerTag>(0, N),
                              f_tag, reducer_scalar);
      ASSERT_EQ(minmax_scalar.min_val, reference_min);

      for (int i = 0; i < N; i++) {
        if ((i == minmax_scalar.min_loc) && (h_values(i) == reference_min)) {
          reference_minloc = i;
        }
      }

      ASSERT_EQ(minmax_scalar.min_loc, reference_minloc);
      ASSERT_EQ(minmax_scalar.max_val, reference_max);

      for (int i = 0; i < N; i++) {
        if ((i == minmax_scalar.max_loc) && (h_values(i) == reference_max)) {
          reference_maxloc = i;
        }
      }

      ASSERT_EQ(minmax_scalar.max_loc, reference_maxloc);

      value_type minmax_scalar_view = reducer_scalar.reference();
      ASSERT_EQ(minmax_scalar_view.min_val, reference_min);
      ASSERT_EQ(minmax_scalar_view.min_loc, reference_minloc);
      ASSERT_EQ(minmax_scalar_view.max_val, reference_max);
      ASSERT_EQ(minmax_scalar_view.max_loc, reference_maxloc);
    }

    {
      Kokkos::View<value_type, Kokkos::HostSpace> minmax_view("View");
      Kokkos::MinMaxLoc<Scalar, int> reducer_view(minmax_view);
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_view);
      Kokkos::fence();

      value_type minmax_view_scalar = minmax_view();
      ASSERT_EQ(minmax_view_scalar.min_val, reference_min);
      ASSERT_EQ(minmax_view_scalar.min_loc, reference_minloc);
      ASSERT_EQ(minmax_view_scalar.max_val, reference_max);
      ASSERT_EQ(minmax_view_scalar.max_loc, reference_maxloc);

      value_type minmax_view_view = reducer_view.reference();
      ASSERT_EQ(minmax_view_view.min_val, reference_min);
      ASSERT_EQ(minmax_view_view.min_loc, reference_minloc);
      ASSERT_EQ(minmax_view_view.max_val, reference_max);
      ASSERT_EQ(minmax_view_view.max_loc, reference_maxloc);
    }
  }

  static void test_minmaxloc_2d(int N) {
    using reducer_type = Kokkos::MinMaxLoc<Scalar, MyPair>;
    using value_type   = typename reducer_type::value_type;

    Kokkos::View<Scalar**, ExecSpace> values("Values", N, N);
    auto h_values           = Kokkos::create_mirror_view(values);
    Scalar reference_max    = std::numeric_limits<Scalar>::min();
    Scalar reference_min    = std::numeric_limits<Scalar>::max();
    MyPair reference_minloc = {{-1, -1}};
    MyPair reference_maxloc = {{-1, -1}};

    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) {
        h_values(i, j) = (Scalar)(rand() % 100000 + 2);
      }

    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) {
        if (h_values(i, j) > reference_max) {
          reference_max    = h_values(i, j);
          reference_maxloc = {{i, j}};
        } else if (h_values(i, j) == reference_max) {
          // Make max unique.
          h_values(i, j) -= Scalar(1);
        }
      }

    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) {
        if (h_values(i, j) < reference_min) {
          reference_min    = h_values(i, j);
          reference_minloc = {{i, j}};
        } else if (h_values(i, j) == reference_min) {
          // Make min unique.
          h_values(i, j) += Scalar(1);
        }
      }

    Kokkos::deep_copy(values, h_values);

    MinMaxLocFunctor2D f;
    f.values = values;
    {
      value_type minmax_scalar;
      reducer_type reducer_scalar(minmax_scalar);

      Kokkos::parallel_reduce(
          Kokkos::MDRangePolicy<Kokkos::Rank<2>, ExecSpace>({0, 0}, {N, N}), f,
          reducer_scalar);

      ASSERT_EQ(minmax_scalar.min_val, reference_min);
      for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
          if ((minmax_scalar.min_loc == MyPair{{i, j}}) &&
              (h_values(i, j) == reference_min)) {
            reference_minloc = {{i, j}};
          }
        }
      ASSERT_EQ(minmax_scalar.min_loc, reference_minloc);

      ASSERT_EQ(minmax_scalar.max_val, reference_max);
      for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
          if ((minmax_scalar.max_loc == MyPair{{i, j}}) &&
              (h_values(i, j) == reference_max)) {
            reference_maxloc = {{i, j}};
          }
        }
      ASSERT_EQ(minmax_scalar.max_loc, reference_maxloc);
    }
  }

  static void test_minmaxloc_loc_init(int N) {
    using reducer_type       = Kokkos::MinMaxLoc<Scalar, int>;
    using reducer_value_type = typename reducer_type::value_type;

    Kokkos::View<Scalar*, ExecSpace> values("Values", N);
    auto h_values = Kokkos::create_mirror_view(values);

    auto functor = KOKKOS_LAMBDA(const int i, reducer_value_type& update) {
      auto x = values(i);
      if (i % 2 == 0) return;
      if (x <= update.min_val) {
        update.min_val = x;
        update.min_loc = i;
      }
      if (x >= update.max_val) {
        update.max_val = x;
        update.max_loc = i;
      }
    };

    {
      for (int i = 0; i < N; ++i) {
        h_values(i) = Kokkos::reduction_identity<Scalar>::min();
      }
      Kokkos::deep_copy(values, h_values);

      reducer_value_type value_loc{0, 0, -1, -1};

      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), functor,
                              reducer_type(value_loc));

      ASSERT_EQ(value_loc.min_val, h_values(0));
      ASSERT_EQ(value_loc.max_val, h_values(0));
      ASSERT_GE(value_loc.min_loc, 0);
      ASSERT_LT(value_loc.min_loc, N);
      ASSERT_GE(value_loc.max_loc, 0);
      ASSERT_LT(value_loc.max_loc, N);
    }

    {
      for (int i = 0; i < N; ++i) {
        h_values(i) = Kokkos::reduction_identity<Scalar>::max();
      }
      Kokkos::deep_copy(values, h_values);

      reducer_value_type value_loc{0, 0, -1, -1};

      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), functor,
                              reducer_type(value_loc));

      ASSERT_EQ(value_loc.min_val, h_values(0));
      ASSERT_EQ(value_loc.max_val, h_values(0));
      ASSERT_GE(value_loc.min_loc, 0);
      ASSERT_LT(value_loc.min_loc, N);
      ASSERT_GE(value_loc.max_loc, 0);
      ASSERT_LT(value_loc.max_loc, N);
    }
  }

  static void test_minmaxfirstlastloc_loc_init(int N) {
    using reducer_type       = Kokkos::MinMaxFirstLastLoc<Scalar, int>;
    using reducer_value_type = typename reducer_type::value_type;

    Kokkos::View<Scalar*, ExecSpace> values("Values", N);
    auto h_values = Kokkos::create_mirror_view(values);

    auto functor = KOKKOS_LAMBDA(const int i, reducer_value_type& update) {
      auto x = values(i);
      if (i % 2 == 0) return;
      if (x <= update.min_val) {
        update.min_val = x;
        update.min_loc = i;
      }
      if (x >= update.max_val) {
        update.max_val = x;
        update.max_loc = i;
      }
    };

    {
      for (int i = 0; i < N; ++i) {
        h_values(i) = Kokkos::reduction_identity<Scalar>::min();
      }
      Kokkos::deep_copy(values, h_values);

      reducer_value_type value_loc{0, 0, -1, -1};

      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), functor,
                              reducer_type(value_loc));

      ASSERT_EQ(value_loc.min_val, h_values(0));
      ASSERT_EQ(value_loc.max_val, h_values(0));
      ASSERT_GE(value_loc.min_loc, 0);
      ASSERT_LT(value_loc.min_loc, N);
      ASSERT_GE(value_loc.max_loc, 0);
      ASSERT_LT(value_loc.max_loc, N);
    }

    {
      for (int i = 0; i < N; ++i) {
        h_values(i) = Kokkos::reduction_identity<Scalar>::max();
      }
      Kokkos::deep_copy(values, h_values);

      reducer_value_type value_loc{0, 0, -1, -1};

      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), functor,
                              reducer_type(value_loc));

      ASSERT_EQ(value_loc.min_val, h_values(0));
      ASSERT_EQ(value_loc.max_val, h_values(0));
      ASSERT_GE(value_loc.min_loc, 0);
      ASSERT_LT(value_loc.min_loc, N);
      ASSERT_GE(value_loc.max_loc, 0);
      ASSERT_LT(value_loc.max_loc, N);
    }
  }

  static void test_minfirstloc_loc_init(int N) {
    using reducer_type       = Kokkos::MinFirstLoc<Scalar, int>;
    using reducer_value_type = typename reducer_type::value_type;

    Kokkos::View<Scalar*, ExecSpace> values("Values", N);
    auto h_values = Kokkos::create_mirror_view(values);

    for (int i = 0; i < N; ++i) {
      h_values(i) = Kokkos::reduction_identity<Scalar>::min();
    }
    Kokkos::deep_copy(values, h_values);

    reducer_value_type value_loc{0, -1};

    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<ExecSpace>(0, N),
        KOKKOS_LAMBDA(const int i, reducer_value_type& update) {
          auto x = values(i);
          if (i % 2 == 0)
            return;
          else if (x <= update.val) {
            update.val = x;
            update.loc = i;
          }
        },
        reducer_type(value_loc));

    ASSERT_EQ(value_loc.val, h_values(0));
    ASSERT_GE(value_loc.loc, 0);
    ASSERT_LT(value_loc.loc, N);
  }

  static void test_maxfirstloc_loc_init(int N) {
    using reducer_type       = Kokkos::MaxFirstLoc<Scalar, int>;
    using reducer_value_type = typename reducer_type::value_type;

    Kokkos::View<Scalar*, ExecSpace> values("Values", N);
    auto h_values = Kokkos::create_mirror_view(values);

    for (int i = 0; i < N; ++i) {
      h_values(i) = Kokkos::reduction_identity<Scalar>::max();
    }
    Kokkos::deep_copy(values, h_values);

    reducer_value_type value_loc{0, -1};

    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<ExecSpace>(0, N),
        KOKKOS_LAMBDA(const int i, reducer_value_type& update) {
          auto x = values(i);
          if (i % 2 == 0)
            return;
          else if (x >= update.val) {
            update.val = x;
            update.loc = i;
          }
        },
        reducer_type(value_loc));

    ASSERT_EQ(value_loc.val, h_values(0));
    ASSERT_GE(value_loc.loc, 0);
    ASSERT_LT(value_loc.loc, N);
  }

  static void test_BAnd(int N) {
    Kokkos::View<Scalar*, ExecSpace> values("Values", N);
    auto h_values         = Kokkos::create_mirror_view(values);
    Scalar reference_band = Scalar() | (~Scalar());

    for (int i = 0; i < N; i++) {
      h_values(i)    = (Scalar)(rand() % 100000 + 1);
      reference_band = reference_band & h_values(i);
    }
    Kokkos::deep_copy(values, h_values);

    BAndFunctor f;
    f.values = values;
    BAndFunctorTag f_tag;
    f_tag.values = values;
    Scalar init  = Scalar() | (~Scalar());

    {
      Scalar band_scalar = init;
      Kokkos::BAnd<Scalar> reducer_scalar(band_scalar);

      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_scalar);
      ASSERT_EQ(band_scalar, reference_band);

      band_scalar = init;
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace, ReducerTag>(0, N),
                              f_tag, reducer_scalar);
      ASSERT_EQ(band_scalar, reference_band);

      Scalar band_scalar_view = reducer_scalar.reference();

      ASSERT_EQ(band_scalar_view, reference_band);
    }

    {
      Kokkos::View<Scalar, Kokkos::HostSpace> band_view("View");
      band_view() = init;
      Kokkos::BAnd<Scalar> reducer_view(band_view);
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_view);
      Kokkos::fence();

      Scalar band_view_scalar = band_view();
      ASSERT_EQ(band_view_scalar, reference_band);

      Scalar band_view_view = reducer_view.reference();
      ASSERT_EQ(band_view_view, reference_band);
    }
  }

  static void test_BOr(int N) {
    Kokkos::View<Scalar*, ExecSpace> values("Values", N);
    auto h_values        = Kokkos::create_mirror_view(values);
    Scalar reference_bor = Scalar() & (~Scalar());

    for (int i = 0; i < N; i++) {
      h_values(i)   = (Scalar)((rand() % 100000 + 1) * 2);
      reference_bor = reference_bor | h_values(i);
    }
    Kokkos::deep_copy(values, h_values);

    BOrFunctor f;
    f.values = values;
    BOrFunctorTag f_tag;
    f_tag.values = values;
    Scalar init  = Scalar() & (~Scalar());

    {
      Scalar bor_scalar = init;
      Kokkos::BOr<Scalar> reducer_scalar(bor_scalar);

      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_scalar);
      ASSERT_EQ(bor_scalar, reference_bor);

      bor_scalar = init;
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace, ReducerTag>(0, N),
                              f_tag, reducer_scalar);
      ASSERT_EQ(bor_scalar, reference_bor);

      Scalar bor_scalar_view = reducer_scalar.reference();
      ASSERT_EQ(bor_scalar_view, reference_bor);
    }

    {
      Kokkos::View<Scalar, Kokkos::HostSpace> bor_view("View");
      bor_view() = init;
      Kokkos::BOr<Scalar> reducer_view(bor_view);
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_view);
      Kokkos::fence();

      Scalar bor_view_scalar = bor_view();
      ASSERT_EQ(bor_view_scalar, reference_bor);

      Scalar bor_view_view = reducer_view.reference();
      ASSERT_EQ(bor_view_view, reference_bor);
    }
  }

  static void test_LAnd(int N) {
    Kokkos::View<Scalar*, ExecSpace> values("Values", N);
    auto h_values         = Kokkos::create_mirror_view(values);
    Scalar reference_land = 1;

    for (int i = 0; i < N; i++) {
      h_values(i)    = (Scalar)(rand() % 2);
      reference_land = reference_land && h_values(i);
    }
    Kokkos::deep_copy(values, h_values);

    LAndFunctor f;
    f.values = values;
    LAndFunctorTag f_tag;
    f_tag.values = values;
    Scalar init  = 1;

    {
      Scalar land_scalar = init;
      Kokkos::LAnd<Scalar> reducer_scalar(land_scalar);

      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_scalar);
      ASSERT_EQ(land_scalar, reference_land);

      land_scalar = init;
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace, ReducerTag>(0, N),
                              f_tag, reducer_scalar);
      ASSERT_EQ(land_scalar, reference_land);

      Scalar land_scalar_view = reducer_scalar.reference();
      ASSERT_EQ(land_scalar_view, reference_land);
    }

    {
      Kokkos::View<Scalar, Kokkos::HostSpace> land_view("View");
      land_view() = init;
      Kokkos::LAnd<Scalar> reducer_view(land_view);
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_view);
      Kokkos::fence();

      Scalar land_view_scalar = land_view();
      ASSERT_EQ(land_view_scalar, reference_land);

      Scalar land_view_view = reducer_view.reference();
      ASSERT_EQ(land_view_view, reference_land);
    }
  }

  static void test_LOr(int N) {
    Kokkos::View<Scalar*, ExecSpace> values("Values", N);
    auto h_values        = Kokkos::create_mirror_view(values);
    Scalar reference_lor = 0;

    for (int i = 0; i < N; i++) {
      h_values(i)   = (Scalar)(rand() % 2);
      reference_lor = reference_lor || h_values(i);
    }
    Kokkos::deep_copy(values, h_values);

    LOrFunctor f;
    f.values = values;
    LOrFunctorTag f_tag;
    f_tag.values = values;
    Scalar init  = 0;

    {
      Scalar lor_scalar = init;
      Kokkos::LOr<Scalar> reducer_scalar(lor_scalar);

      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_scalar);
      ASSERT_EQ(lor_scalar, reference_lor);

      lor_scalar = init;
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace, ReducerTag>(0, N),
                              f_tag, reducer_scalar);
      ASSERT_EQ(lor_scalar, reference_lor);

      Scalar lor_scalar_view = reducer_scalar.reference();
      ASSERT_EQ(lor_scalar_view, reference_lor);
    }

    {
      Kokkos::View<Scalar, Kokkos::HostSpace> lor_view("View");
      lor_view() = init;
      Kokkos::LOr<Scalar> reducer_view(lor_view);
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0, N), f,
                              reducer_view);
      Kokkos::fence();

      Scalar lor_view_scalar = lor_view();
      ASSERT_EQ(lor_view_scalar, reference_lor);

      Scalar lor_view_view = reducer_view.reference();
      ASSERT_EQ(lor_view_view, reference_lor);
    }
  }

  static void execute_float() {
    test_launch_bounds();
    test_sum(10001);
    test_prod(35);
    test_min(10003);
    test_minloc(10003);
    test_minloc_loc_init(3);
// FIXME_OPENACC - custom reduction with MDRangePolicy is not yet implemented.
#if !defined(KOKKOS_ENABLE_OPENACC)
// FIXME_OPENMPTARGET requires custom reductions.
#if !defined(KOKKOS_ENABLE_OPENMPTARGET)
    test_minloc_2d(100);
#endif
#endif
    test_max(10007);
    test_maxloc(10007);
    test_maxloc_loc_init(3);
// FIXME_OPENACC - custom reduction with MDRangePolicy is not yet implemented.
#if !defined(KOKKOS_ENABLE_OPENACC)
// FIXME_OPENMPTARGET requires custom reductions.
#if !defined(KOKKOS_ENABLE_OPENMPTARGET)
    test_maxloc_2d(100);
#endif
#endif
#if defined(KOKKOS_ENABLE_OPENMPTARGET)  // FIXME_OPENMPTARGET custom reducers
    test_minmaxloc(10007);
#else
    test_minmaxloc(10007);
    test_minmaxloc_loc_init(3);
// FIXME_OPENACC - custom reduction with MDRangePolicy is not yet implemented.
#if !defined(KOKKOS_ENABLE_OPENACC)
    test_minmaxloc_2d(100);
#endif

    test_minmaxfirstlastloc_loc_init(3);
    test_minfirstloc_loc_init(3);
    test_maxfirstloc_loc_init(3);
#endif
  }

  // NOTE test_prod generates N random numbers between 1 and 4.
  // Although unlikely, the test below could still in principle overflow.
  // For reference log(numeric_limits<int>)/log(4) is 15.5
  static void execute_integer() {
    test_launch_bounds();
    test_sum(10001);
    test_prod(sizeof(Scalar) > 4 ? 35 : 19);  // avoid int overflow (see above)
    test_min(10003);
    test_minloc(10003);
    test_minloc_loc_init(3);
#if defined(KOKKOS_ENABLE_CUDA)
    if (!std::is_same_v<ExecSpace, Kokkos::Cuda>)
#endif
    // FIXME_OPENACC - custom reduction with MDRangePolicy is not yet
    // implemented.
#if !defined(KOKKOS_ENABLE_OPENACC)
    // FIXME_OPENMPTARGET requires custom reductions.
#if !defined(KOKKOS_ENABLE_OPENMPTARGET)
      test_minloc_2d(100);
#endif
#endif
    test_max(10007);
    test_maxloc(10007);
    test_maxloc_loc_init(3);
#if defined(KOKKOS_ENABLE_CUDA)
    if (!std::is_same_v<ExecSpace, Kokkos::Cuda>)
#endif
// FIXME_OPENACC - custom reduction with MDRangePolicy is not yet implemented.
#if !defined(KOKKOS_ENABLE_OPENACC)
// FIXME_OPENMPTARGET requires custom reductions.
#if !defined(KOKKOS_ENABLE_OPENMPTARGET)
      test_maxloc_2d(100);
#endif
#endif
#if defined(KOKKOS_ENABLE_OPENMPTARGET)  // FIXME_OPENMPTARGET custom reducers
    test_minmaxloc(10007);
#else
    test_minmaxloc(10007);
    test_minmaxloc_loc_init(3);
// FIXME_OPENACC - custom reduction with MDRangePolicy is not yet implemented.
#if !defined(KOKKOS_ENABLE_OPENACC)
    test_minmaxloc_2d(100);
#endif

    test_minmaxfirstlastloc_loc_init(3);
    test_minfirstloc_loc_init(3);
    test_maxfirstloc_loc_init(3);
#endif
    test_BAnd(35);
    test_BOr(35);
    test_LAnd(35);
    test_LOr(35);
  }

  static void execute_basic() {
    test_sum(10001);
    test_prod(35);
  }

  static void execute_bool() {
    test_LAnd(10001);
    test_LOr(35);
  }
};

}  // namespace Test
