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
#include <Kokkos_StdAlgorithms.hpp>

namespace Test {
namespace stdalgos {

TEST(std_algorithms, is_admissible_to_std_algorithms) {
  namespace KE     = Kokkos::Experimental;
  using value_type = double;

  static constexpr size_t extent0 = 13;
  static constexpr size_t extent1 = 18;
  static constexpr size_t extent2 = 18;

  //-------------
  // 1d views
  //-------------
  using static_view_1d_t = Kokkos::View<value_type[extent0]>;
  static_view_1d_t static_view_1d{"std-algo-test-1d-contiguous-view-static"};

  using dyn_view_1d_t = Kokkos::View<value_type*>;
  dyn_view_1d_t dynamic_view_1d{"std-algo-test-1d-contiguous-view-dynamic",
                                extent0};

  using strided_view_1d_t = Kokkos::View<value_type*, Kokkos::LayoutStride>;
  Kokkos::LayoutStride layout1d{extent0, 2};
  strided_view_1d_t strided_view_1d{"std-algo-test-1d-strided-view", layout1d};
  ASSERT_EQ(layout1d.dimension[0], 13u);
  ASSERT_EQ(layout1d.stride[0], 2u);
  // they are admissible
  KE::Impl::static_assert_is_admissible_to_kokkos_std_algorithms(
      static_view_1d);
  KE::Impl::static_assert_is_admissible_to_kokkos_std_algorithms(
      dynamic_view_1d);
  KE::Impl::static_assert_is_admissible_to_kokkos_std_algorithms(
      strided_view_1d);

  //-------------
  // 2d views
  //-------------
  using static_view_2d_t  = Kokkos::View<value_type[extent0][extent1]>;
  using dyn_view_2d_t     = Kokkos::View<value_type**>;
  using strided_view_2d_t = Kokkos::View<value_type**, Kokkos::LayoutStride>;
  // non admissible
  EXPECT_FALSE(KE::Impl::is_admissible_to_kokkos_std_algorithms<
               static_view_2d_t>::value);
  EXPECT_FALSE(
      KE::Impl::is_admissible_to_kokkos_std_algorithms<dyn_view_2d_t>::value);
  EXPECT_FALSE(KE::Impl::is_admissible_to_kokkos_std_algorithms<
               strided_view_2d_t>::value);

  //-------------
  // 3d views
  //-------------
  using static_view_3d_t  = Kokkos::View<value_type[extent0][extent1][extent2]>;
  using dyn_view_3d_t     = Kokkos::View<value_type***>;
  using strided_view_3d_t = Kokkos::View<value_type***, Kokkos::LayoutStride>;
  // non admissible
  EXPECT_FALSE(KE::Impl::is_admissible_to_kokkos_std_algorithms<
               static_view_3d_t>::value);
  EXPECT_FALSE(
      KE::Impl::is_admissible_to_kokkos_std_algorithms<dyn_view_3d_t>::value);
  EXPECT_FALSE(KE::Impl::is_admissible_to_kokkos_std_algorithms<
               strided_view_3d_t>::value);
}

TEST(std_algorithms, expect_no_overlap) {
  namespace KE     = Kokkos::Experimental;
  using value_type = double;

  static constexpr size_t extent0 = 13;

  //-------------
  // 1d views
  //-------------
  using static_view_1d_t = Kokkos::View<value_type[extent0]>;
  [[maybe_unused]] static_view_1d_t static_view_1d{
      "std-algo-test-1d-contiguous-view-static"};

  using dyn_view_1d_t = Kokkos::View<value_type*>;
  [[maybe_unused]] dyn_view_1d_t dynamic_view_1d{
      "std-algo-test-1d-contiguous-view-dynamic", extent0};

  using strided_view_1d_t = Kokkos::View<value_type*, Kokkos::LayoutStride>;
  Kokkos::LayoutStride layout1d{extent0, 2};
  strided_view_1d_t strided_view_1d{"std-algo-test-1d-strided-view", layout1d};

// Overlapping because iterators are identical
#if defined(KOKKOS_ENABLE_DEBUG)
  auto first_s = KE::begin(static_view_1d);
  auto last_s  = first_s + extent0;
  EXPECT_DEATH({ KE::Impl::expect_no_overlap(first_s, last_s, first_s); },
               "Kokkos contract violation:.*");

  auto first_d = KE::begin(dynamic_view_1d);
  auto last_d  = first_d + extent0;
  EXPECT_DEATH({ KE::Impl::expect_no_overlap(first_d, last_d, first_d); },
               "Kokkos contract violation:.*");

  auto first_st = KE::begin(strided_view_1d);
  auto last_st  = first_st + extent0;
  EXPECT_DEATH({ KE::Impl::expect_no_overlap(first_st, last_st, first_st); },
               "Kokkos contract violation:.*");
#endif

  // Ranges are overlapped
  static constexpr size_t sub_extent0 = 6, offset0 = 3;
  std::pair<size_t, size_t> range0(0, sub_extent0),
      range1(offset0, offset0 + sub_extent0);
#if defined(KOKKOS_ENABLE_DEBUG)
  auto static_view_1d_0 = Kokkos::subview(static_view_1d, range0);
  auto static_view_1d_1 = Kokkos::subview(static_view_1d, range1);
  auto first_s0         = KE::begin(static_view_1d_0);  // [0, 6)
  auto last_s0          = first_s0 + static_view_1d_0.extent(0);
  auto first_s1         = KE::begin(static_view_1d_1);  // [3, 9)
  EXPECT_DEATH({ KE::Impl::expect_no_overlap(first_s0, last_s0, first_s1); },
               "Kokkos contract violation:.*");

  auto dynamic_view_1d_0 = Kokkos::subview(dynamic_view_1d, range0);
  auto dynamic_view_1d_1 = Kokkos::subview(dynamic_view_1d, range1);
  auto first_d0          = KE::begin(dynamic_view_1d_0);  // [0, 6)
  auto last_d0           = first_d0 + dynamic_view_1d_0.extent(0);
  auto first_d1          = KE::begin(dynamic_view_1d_1);  // [3, 9)
  EXPECT_DEATH({ KE::Impl::expect_no_overlap(first_d0, last_d0, first_d1); },
               "Kokkos contract violation:.*");
#endif

  auto strided_view_1d_0 = Kokkos::subview(strided_view_1d, range0);
  auto strided_view_1d_1 = Kokkos::subview(strided_view_1d, range1);
  auto first_st0         = KE::begin(strided_view_1d_0);  // [0, 12)
  auto last_st0          = first_st0 + strided_view_1d_0.extent(0);
  auto first_st1         = KE::begin(strided_view_1d_1);  // [3, 15)
  // Does not overlap since offset (=3) is not divisible by stride (=2)
  EXPECT_NO_THROW(
      { KE::Impl::expect_no_overlap(first_st0, last_st0, first_st1); });

  // Iterating over the same range without overlapping
  Kokkos::View<value_type[2][extent0], Kokkos::LayoutLeft> static_view_2d{
      "std-algo-test-2d-contiguous-view-static"};
  auto sub_static_view_1d_0 = Kokkos::subview(static_view_2d, 0, Kokkos::ALL);
  auto sub_static_view_1d_1 = Kokkos::subview(static_view_2d, 1, Kokkos::ALL);
  auto sub_first_s0         = KE::begin(sub_static_view_1d_0);  // 0, 2, 4, ...
  auto sub_last_s0          = sub_first_s0 + sub_static_view_1d_0.extent(0);
  auto sub_first_s1         = KE::begin(sub_static_view_1d_1);  // 1, 3, 5, ...

  EXPECT_NO_THROW({
    KE::Impl::expect_no_overlap(sub_first_s0, sub_last_s0, sub_first_s1);
  });

  Kokkos::View<value_type**, Kokkos::LayoutLeft> dynamic_view_2d{
      "std-algo-test-2d-contiguous-view-dynamic", 2, extent0};
  auto sub_dynamic_view_1d_0 = Kokkos::subview(dynamic_view_2d, 0, Kokkos::ALL);
  auto sub_dynamic_view_1d_1 = Kokkos::subview(dynamic_view_2d, 1, Kokkos::ALL);
  auto sub_first_d0 = KE::begin(sub_dynamic_view_1d_0);  // 0, 2, 4, ...
  auto sub_last_d0  = sub_first_d0 + sub_dynamic_view_1d_0.extent(0);
  auto sub_first_d1 = KE::begin(sub_dynamic_view_1d_1);  // 1, 3, 5, ...

  EXPECT_NO_THROW({
    KE::Impl::expect_no_overlap(sub_first_d0, sub_last_d0, sub_first_d1);
  });

  Kokkos::LayoutStride layout2d{2, 3, extent0, 2 * 3};
  Kokkos::View<value_type**, Kokkos::LayoutStride> strided_view_2d{
      "std-algo-test-2d-contiguous-view-strided", layout2d};
  auto sub_strided_view_1d_0 = Kokkos::subview(strided_view_2d, 0, Kokkos::ALL);
  auto sub_strided_view_1d_1 = Kokkos::subview(strided_view_2d, 1, Kokkos::ALL);
  auto sub_first_st0 = KE::begin(sub_strided_view_1d_0);  // 0, 6, 12, ...
  auto sub_last_st0  = sub_first_st0 + sub_strided_view_1d_0.extent(0);
  auto sub_first_st1 = KE::begin(sub_strided_view_1d_1);  // 1, 7, 13, ...

  EXPECT_NO_THROW({
    KE::Impl::expect_no_overlap(sub_first_st0, sub_last_st0, sub_first_st1);
  });
}

}  // namespace stdalgos
}  // namespace Test
