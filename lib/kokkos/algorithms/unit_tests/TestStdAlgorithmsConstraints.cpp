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

}  // namespace stdalgos
}  // namespace Test
