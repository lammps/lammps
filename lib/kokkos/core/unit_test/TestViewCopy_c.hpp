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

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

namespace {
// Do not rely on deep_copy(0) as we want to test it!
template <class ViewType, class ExecSpace>
void reset_view(const ExecSpace& space, ViewType& a, int magic) {
  auto policy = Kokkos::RangePolicy<ExecSpace>(space, 0, a.span());

  assert(a.span_is_contiguous());

  Kokkos::parallel_for(
      "TestViewCopy::ResetView", policy,
      KOKKOS_LAMBDA(int i) { a.data()[i] = magic; });
}

template <class ViewType, class ExecSpace>
size_t compute_overall_sum(const ExecSpace& space, ViewType& a) {
  auto policy = Kokkos::RangePolicy<ExecSpace>(space, 0, a.span());

  assert(a.span_is_contiguous());

  typename ViewType::value_type sum = 0;
  Kokkos::parallel_reduce(
      "TestViewCopy::ComputeSum", policy,
      KOKKOS_LAMBDA(int i, int& lcl_sum) { lcl_sum += a.data()[i]; }, sum);

  return static_cast<size_t>(sum);
}

template <typename ExecSpace, typename DT, typename... DP>
bool check_magic_value(
    const ExecSpace& space, const Kokkos::View<DT, DP...>& a, int magic,
    std::enable_if_t<Kokkos::ViewTraits<DT, DP...>::rank == 0>* = nullptr) {
  auto policy = Kokkos::RangePolicy<ExecSpace>(space, 0, 1);

  bool all_elements_are_set;  // Uninitialized, set by parallel_reduce

  Kokkos::parallel_reduce(
      "TestViewCopy::CheckMagicValueRank0", policy,
      KOKKOS_LAMBDA(int, bool& local_check) { local_check &= (a() == magic); },
      Kokkos::LAnd<bool>(all_elements_are_set));

  return all_elements_are_set;
}

template <typename ExecSpace, typename DT, typename... DP>
bool check_magic_value(
    const ExecSpace& space, const Kokkos::View<DT, DP...>& a, int magic,
    std::enable_if_t<Kokkos::ViewTraits<DT, DP...>::rank == 1>* = nullptr) {
  auto policy = Kokkos::RangePolicy<ExecSpace>(space, 0, a.extent(0));

  bool all_elements_are_set;  // Uninitialized, set by parallel_reduce

  Kokkos::parallel_reduce(
      "TestViewCopy::CheckMagicValueRank1", policy,
      KOKKOS_LAMBDA(int i, bool& local_check) {
        local_check &= (a(i) == magic);
      },
      Kokkos::LAnd<bool>(all_elements_are_set));

  return all_elements_are_set;
}

template <typename ExecSpace, typename DT, typename... DP>
bool check_magic_value(
    const ExecSpace& space, const Kokkos::View<DT, DP...>& a, int magic,
    std::enable_if_t<Kokkos::ViewTraits<DT, DP...>::rank == 2>* = nullptr) {
  auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>, ExecSpace>(
      space, {0, 0}, {a.extent(0), a.extent(1)});

  bool all_elements_are_set;  // Uninitialized, set by parallel_reduce

  Kokkos::parallel_reduce(
      "TestViewCopy::CheckMagicValueRank2", policy,
      KOKKOS_LAMBDA(int i0, int i1, bool& local_check) {
        local_check &= (a(i0, i1) == magic);
      },
      Kokkos::LAnd<bool>(all_elements_are_set));

  return all_elements_are_set;
}

template <typename ExecSpace, typename DT, typename... DP>
bool check_magic_value(
    const ExecSpace& space, const Kokkos::View<DT, DP...>& a, int magic,
    std::enable_if_t<Kokkos::ViewTraits<DT, DP...>::rank == 3>* = nullptr) {
  auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecSpace>(
      space, {0, 0, 0}, {a.extent(0), a.extent(1), a.extent(2)});

  bool all_elements_are_set;  // Uninitialized, set by parallel_reduce

  Kokkos::parallel_reduce(
      "TestViewCopy::CheckMagicValueRank3", policy,
      KOKKOS_LAMBDA(int i0, int i1, int i2, bool& local_check) {
        local_check &= (a(i0, i1, i2) == magic);
      },
      Kokkos::LAnd<bool>(all_elements_are_set));

  return all_elements_are_set;
}

template <typename ExecSpace, typename DT, typename... DP>
bool check_magic_value(
    const ExecSpace& space, const Kokkos::View<DT, DP...>& a, int magic,
    std::enable_if_t<Kokkos::ViewTraits<DT, DP...>::rank == 4>* = nullptr) {
  auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<4>, ExecSpace>(
      space, {0, 0, 0, 0},
      {a.extent(0), a.extent(1), a.extent(2), a.extent(3)});

  bool all_elements_are_set;  // Uninitialized, set by parallel_reduce

  Kokkos::parallel_reduce(
      "TestViewCopy::CheckMagicValueRank4", policy,
      KOKKOS_LAMBDA(int i0, int i1, int i2, int i3, bool& local_check) {
        local_check &= (a(i0, i1, i2, i3) == magic);
      },
      Kokkos::LAnd<bool>(all_elements_are_set));

  return all_elements_are_set;
}

template <typename ExecSpace, typename DT, typename... DP>
bool check_magic_value(
    const ExecSpace& space, const Kokkos::View<DT, DP...>& a, int magic,
    std::enable_if_t<Kokkos::ViewTraits<DT, DP...>::rank == 5>* = nullptr) {
  auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<5>, ExecSpace>(
      space, {0, 0, 0, 0, 0},
      {a.extent(0), a.extent(1), a.extent(2), a.extent(3), a.extent(4)});

  bool all_elements_are_set;  // Uninitialized, set by parallel_reduce

  Kokkos::parallel_reduce(
      "TestViewCopy::CheckMagicValueRank5", policy,
      KOKKOS_LAMBDA(int i0, int i1, int i2, int i3, int i4, bool& local_check) {
        local_check &= (a(i0, i1, i2, i3, i4) == magic);
      },
      Kokkos::LAnd<bool>(all_elements_are_set));

  return all_elements_are_set;
}

template <typename ExecSpace, typename DT, typename... DP>
bool check_magic_value(
    const ExecSpace& space, const Kokkos::View<DT, DP...>& a, int magic,
    std::enable_if_t<Kokkos::ViewTraits<DT, DP...>::rank == 6>* = nullptr) {
  auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<6>, ExecSpace>(
      space, {0, 0, 0, 0, 0, 0},
      {a.extent(0), a.extent(1), a.extent(2), a.extent(3), a.extent(4),
       a.extent(5)});

  bool all_elements_are_set;  // Uninitialized, set by parallel_reduce

  Kokkos::parallel_reduce(
      "TestViewCopy::CheckMagicValueRank6", policy,
      KOKKOS_LAMBDA(int i0, int i1, int i2, int i3, int i4, int i5,
                    bool& local_check) {
        local_check &= (a(i0, i1, i2, i3, i4, i5) == magic);
      },
      Kokkos::LAnd<bool>(all_elements_are_set));

  return all_elements_are_set;
}

template <typename ExecSpace, typename DT, typename... DP>
bool check_magic_value(
    const ExecSpace& space, const Kokkos::View<DT, DP...>& a, int magic,
    std::enable_if_t<Kokkos::ViewTraits<DT, DP...>::rank == 7>* = nullptr) {
  auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<6>, ExecSpace>(
      space, {0, 0, 0, 0, 0, 0},
      {a.extent(0), a.extent(1), a.extent(2), a.extent(3), a.extent(4),
       a.extent(5)});

  bool all_elements_are_set = true;

  for (size_t outer = 0; outer < a.extent(6); ++outer) {
    bool all_local_elements_are_set;  // Uninitialized, set by parallel_reduce
    Kokkos::parallel_reduce(
        "TestViewCopy::CheckMagicValueRank7", policy,
        KOKKOS_LAMBDA(int i0, int i1, int i2, int i3, int i4, int i5,
                      bool& local_check) {
          local_check &= (a(i0, i1, i2, i3, i4, i5, outer) == magic);
        },
        Kokkos::LAnd<bool>(all_local_elements_are_set));

    all_elements_are_set = all_elements_are_set && all_local_elements_are_set;
  }
  return all_elements_are_set;
}

template <typename ExecSpace, typename DT, typename... DP>
bool check_magic_value(
    const ExecSpace& space, const Kokkos::View<DT, DP...>& a, int magic,
    std::enable_if_t<Kokkos::ViewTraits<DT, DP...>::rank == 8>* = nullptr) {
  auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<6>, ExecSpace>(
      space, {0, 0, 0, 0, 0, 0},
      {a.extent(0), a.extent(1), a.extent(2), a.extent(3), a.extent(4),
       a.extent(5)});

  bool all_elements_are_set = true;

  for (size_t outer = 0; outer < a.extent(7); ++outer) {
    for (size_t inner = 0; inner < a.extent(6); ++inner) {
      bool all_local_elements_are_set;  // Uninitialized, set by parallel_reduce
      Kokkos::parallel_reduce(
          "TestViewCopy::CheckMagicValueRank8", policy,
          KOKKOS_LAMBDA(int i0, int i1, int i2, int i3, int i4, int i5,
                        bool& local_check) {
            local_check &= (a(i0, i1, i2, i3, i4, i5, inner, outer) == magic);
          },
          Kokkos::LAnd<bool>(all_local_elements_are_set));

      all_elements_are_set = all_elements_are_set && all_local_elements_are_set;
    }
  }
  return all_elements_are_set;
}

template <class ExecSpace, class ViewType>
bool view_fill_test(const ExecSpace& space, ViewType& a, int magic) {
  Kokkos::deep_copy(space, a, magic);
#if defined(KOKKOS_ENABLE_OPENMPTARGET)
  // FIXME_OPENMPTARGET Does not work with Land reducer
  return true;
#else   // KOKKOS_ENABLE_OPENMPTARGET
  return check_magic_value(space, a, magic);
#endif  // KOKKOS_ENABLE_OPENMPTARGET
}

template <class Layout, class Space>
void run_test() {
  int magic = 19;

  using ViewType = Kokkos::View<int********, Layout, Space>;
  // Create views with different lengths for each dimension
  // We want to test if all loops are over the correct dimensions
  // We use prime numbers to make sure that the strides are different
  ViewType a_decreasing("a", 23, 19, 17, 13, 11, 7, 5, 3);
  // We also test with increasing strides to catch more "out-of-bounds" errors
  // within subviews.
  ViewType a_increasing("a", 3, 5, 7, 11, 13, 17, 19, 23);

  using exec_space = typename Space::execution_space;
  auto space       = exec_space();

  // Use subviews in the tests to have cases with different ranks and
  // non-contiguous memory
  // Tests have two parts:
  // 1. Fill the subview with a magic value and check that all elements are set
  // 2. Check if only the subview is set by summing all elements in the view and
  // comparing to the subview size times the magic value

  // Rank 0
  {
    auto sub_dec = Kokkos::subview(a_decreasing, 0, 0, 0, 0, 0, 0, 0, 0);
    EXPECT_TRUE(view_fill_test(space, sub_dec, magic));
    EXPECT_EQ(compute_overall_sum(space, a_decreasing),
              static_cast<size_t>(magic));

    auto sub_inc = Kokkos::subview(a_increasing, 0, 0, 0, 0, 0, 0, 0, 0);
    EXPECT_TRUE(view_fill_test(space, sub_inc, magic));
    EXPECT_EQ(compute_overall_sum(space, a_increasing),
              static_cast<size_t>(magic));
  }
  reset_view(space, a_decreasing, 0);
  reset_view(space, a_increasing, 0);

  // Rank 1
  {
    auto sub_dec =
        Kokkos::subview(a_decreasing, Kokkos::ALL, 0, 0, 0, 0, 0, 0, 0);
    EXPECT_TRUE(view_fill_test(space, sub_dec, magic));
    EXPECT_EQ(compute_overall_sum(space, a_decreasing), sub_dec.size() * magic);

    auto sub_inc =
        Kokkos::subview(a_increasing, Kokkos::ALL, 0, 0, 0, 0, 0, 0, 0);
    EXPECT_TRUE(view_fill_test(space, sub_inc, magic));
    EXPECT_EQ(compute_overall_sum(space, a_increasing), sub_inc.size() * magic);
  }
  reset_view(space, a_decreasing, 0);
  reset_view(space, a_increasing, 0);

  // Rank 2
  {
    auto sub_dec = Kokkos::subview(a_decreasing, Kokkos::ALL, Kokkos::ALL, 0, 0,
                                   0, 0, 0, 0);
    EXPECT_TRUE(view_fill_test(space, sub_dec, magic));
    EXPECT_EQ(compute_overall_sum(space, a_decreasing), sub_dec.size() * magic);

    auto sub_inc = Kokkos::subview(a_increasing, Kokkos::ALL, Kokkos::ALL, 0, 0,
                                   0, 0, 0, 0);
    EXPECT_TRUE(view_fill_test(space, sub_inc, magic));
    EXPECT_EQ(compute_overall_sum(space, a_increasing), sub_inc.size() * magic);
  }
  reset_view(space, a_decreasing, 0);
  reset_view(space, a_increasing, 0);
  space.fence();

  // Rank 3
  {
    auto sub_dec = Kokkos::subview(a_decreasing, Kokkos::ALL, Kokkos::ALL,
                                   Kokkos::ALL, 0, 0, 0, 0, 0);
    EXPECT_TRUE(view_fill_test(space, sub_dec, magic));
    EXPECT_EQ(
        compute_overall_sum(space, a_decreasing),
        sub_dec.extent(0) * sub_dec.extent(1) * sub_dec.extent(2) * magic);

    auto sub_inc = Kokkos::subview(a_increasing, Kokkos::ALL, Kokkos::ALL,
                                   Kokkos::ALL, 0, 0, 0, 0, 0);
    EXPECT_TRUE(view_fill_test(space, sub_inc, magic));
    EXPECT_EQ(compute_overall_sum(space, a_increasing), sub_inc.size() * magic);
  }
  reset_view(space, a_decreasing, 0);
  reset_view(space, a_increasing, 0);
  space.fence();

  // Rank 4
  {
    auto sub_dec = Kokkos::subview(a_decreasing, Kokkos::ALL, Kokkos::ALL,
                                   Kokkos::ALL, Kokkos::ALL, 0, 0, 0, 0);
    EXPECT_TRUE(view_fill_test(space, sub_dec, magic));
    EXPECT_EQ(compute_overall_sum(space, a_decreasing),
              sub_dec.extent(0) * sub_dec.extent(1) * sub_dec.extent(2) *
                  sub_dec.extent(3) * magic);

    auto sub_inc = Kokkos::subview(a_increasing, Kokkos::ALL, Kokkos::ALL,
                                   Kokkos::ALL, Kokkos::ALL, 0, 0, 0, 0);
    EXPECT_TRUE(view_fill_test(space, sub_inc, magic));
    EXPECT_EQ(compute_overall_sum(space, a_increasing), sub_inc.size() * magic);
  }
  reset_view(space, a_decreasing, 0);
  reset_view(space, a_increasing, 0);
  space.fence();

  // Rank 5
  {
    auto sub_dec =
        Kokkos::subview(a_decreasing, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                        Kokkos::ALL, Kokkos::ALL, 0, 0, 0);
    EXPECT_TRUE(view_fill_test(space, sub_dec, magic));
    EXPECT_EQ(compute_overall_sum(space, a_decreasing), sub_dec.size() * magic);

    auto sub_inc =
        Kokkos::subview(a_increasing, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                        Kokkos::ALL, Kokkos::ALL, 0, 0, 0);
    EXPECT_TRUE(view_fill_test(space, sub_inc, magic));
    EXPECT_EQ(compute_overall_sum(space, a_increasing), sub_inc.size() * magic);
  }
  reset_view(space, a_decreasing, 0);
  reset_view(space, a_increasing, 0);
  space.fence();

  // Rank 6
  {
    auto sub_dec =
        Kokkos::subview(a_decreasing, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                        Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, 0, 0);
    EXPECT_TRUE(view_fill_test(space, sub_dec, magic));
    EXPECT_EQ(compute_overall_sum(space, a_decreasing), sub_dec.size() * magic);

    auto sub_inc =
        Kokkos::subview(a_increasing, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                        Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, 0, 0);
    EXPECT_TRUE(view_fill_test(space, sub_inc, magic));
    EXPECT_EQ(compute_overall_sum(space, a_increasing), sub_inc.size() * magic);
  }
  reset_view(space, a_decreasing, 0);
  reset_view(space, a_increasing, 0);
  space.fence();

  // Rank 7
  {
    auto sub_dec =
        Kokkos::subview(a_decreasing, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                        Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, 0);
    EXPECT_TRUE(view_fill_test(space, sub_dec, magic));
    EXPECT_EQ(compute_overall_sum(space, a_decreasing), sub_dec.size() * magic);

    auto sub_inc =
        Kokkos::subview(a_increasing, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                        Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, 0);
    EXPECT_TRUE(view_fill_test(space, sub_inc, magic));
    EXPECT_EQ(compute_overall_sum(space, a_increasing), sub_inc.size() * magic);
  }
  reset_view(space, a_decreasing, 0);
  reset_view(space, a_increasing, 0);
  space.fence();

  // Rank 8
  {
    auto sub_dec = Kokkos::subview(
        a_decreasing, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
        Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, std::make_pair(0, 2));
    EXPECT_TRUE(view_fill_test(space, sub_dec, magic));
    EXPECT_EQ(compute_overall_sum(space, a_decreasing), sub_dec.size() * magic);

    auto sub_inc = Kokkos::subview(
        a_increasing, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
        Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, std::make_pair(0, 2));
    EXPECT_TRUE(view_fill_test(space, sub_inc, magic));
    EXPECT_EQ(compute_overall_sum(space, a_increasing), sub_inc.size() * magic);
  }
}

TEST(TEST_CATEGORY, view_fill_tests_layout_right) {
  using Space  = TEST_EXECSPACE;
  using Layout = Kokkos::LayoutRight;
  run_test<Layout, Space>();
}

TEST(TEST_CATEGORY, view_fill_tests_layout_left) {
  using Space  = TEST_EXECSPACE;
  using Layout = Kokkos::LayoutLeft;
  run_test<Layout, Space>();
}

}  // namespace
