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

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <std_algorithms/Kokkos_Constraints.hpp>

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
  EXPECT_EQ(layout1d.dimension[0], 13u);
  EXPECT_EQ(layout1d.stride[0], 2u);
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
