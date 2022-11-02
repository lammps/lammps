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

#include <TestStdAlgorithmsCommon.hpp>
#include <utility>

namespace Test {
namespace stdalgos {
namespace TransformUnaryOp {

namespace KE = Kokkos::Experimental;

std::string value_type_to_string(int) { return "int"; }
std::string value_type_to_string(double) { return "double"; }

template <class ViewType>
void fill_view(ViewType dest_view) {
  using value_type      = typename ViewType::value_type;
  using exe_space       = typename ViewType::execution_space;
  const std::size_t ext = dest_view.extent(0);
  using aux_view_t      = Kokkos::View<value_type*, exe_space>;
  aux_view_t aux_view("aux_view", ext);
  auto v_h = create_mirror_view(Kokkos::HostSpace(), aux_view);

  for (std::size_t i = 0; i < ext; ++i) {
    v_h(i) = static_cast<value_type>(i);
  }

  Kokkos::deep_copy(aux_view, v_h);
  CopyFunctor<aux_view_t, ViewType> F1(aux_view, dest_view);
  Kokkos::parallel_for("copy", dest_view.extent(0), F1);
}

template <class ViewTypeFrom, class ViewTypeTest>
void verify_data(ViewTypeFrom view_from, ViewTypeTest view_test) {
  using value_type = typename ViewTypeFrom::value_type;

  //! always careful because views might not be deep copyable
  auto view_test_dc = create_deep_copyable_compatible_clone(view_test);
  auto view_test_h =
      create_mirror_view_and_copy(Kokkos::HostSpace(), view_test_dc);

  auto view_from_dc = create_deep_copyable_compatible_clone(view_from);
  auto view_from_h =
      create_mirror_view_and_copy(Kokkos::HostSpace(), view_from_dc);

  for (std::size_t i = 0; i < view_test_h.extent(0); ++i) {
    EXPECT_EQ(view_test_h(i), view_from_h(i) + value_type(1));
  }
}

template <class ValueType>
struct TransformFunctor {
  KOKKOS_INLINE_FUNCTION
  ValueType operator()(const ValueType& val) const {
    return val + ValueType(1);
  }
};

template <class Tag, class ValueType, class InfoType>
void run_single_scenario(const InfoType& scenario_info) {
  const auto name            = std::get<0>(scenario_info);
  const std::size_t view_ext = std::get<1>(scenario_info);
  // std::cout << "transform_unary_op: " << name << ", "
  //           << view_tag_to_string(Tag{}) << ", "
  //           << value_type_to_string(ValueType()) << std::endl;

  auto view_from =
      create_view<ValueType>(Tag{}, view_ext, "transform_uop_from");
  fill_view(view_from);
  TransformFunctor<ValueType> unOp;

  {
    auto view_dest =
        create_view<ValueType>(Tag{}, view_ext, "transform_uop_dest");
    auto r1 = KE::transform(exespace(), KE::begin(view_from),
                            KE::end(view_from), KE::begin(view_dest), unOp);
    verify_data(view_from, view_dest);
    EXPECT_EQ(r1, KE::end(view_dest));
  }

  {
    auto view_dest =
        create_view<ValueType>(Tag{}, view_ext, "transform_uop_dest");
    auto r1 = KE::transform("label", exespace(), KE::begin(view_from),
                            KE::end(view_from), KE::begin(view_dest), unOp);
    verify_data(view_from, view_dest);
    EXPECT_EQ(r1, KE::end(view_dest));
  }

  {
    auto view_dest =
        create_view<ValueType>(Tag{}, view_ext, "transform_uop_dest");
    auto r1 = KE::transform(exespace(), view_from, view_dest, unOp);
    verify_data(view_from, view_dest);
    EXPECT_EQ(r1, KE::end(view_dest));
  }

  {
    auto view_dest =
        create_view<ValueType>(Tag{}, view_ext, "transform_uop_dest");
    auto r1 = KE::transform("label", exespace(), view_from, view_dest, unOp);
    verify_data(view_from, view_dest);
    EXPECT_EQ(r1, KE::end(view_dest));
  }

  Kokkos::fence();
}

template <class Tag, class ValueType>
void run_all_scenarios() {
  const std::map<std::string, std::size_t> scenarios = {
      {"empty", 0},          {"one-element-a", 1},  {"one-element-b", 1},
      {"two-elements-a", 2}, {"two-elements-b", 2}, {"small-a", 9},
      {"small-b", 13},       {"medium", 1103},      {"large", 101513}};

  for (const auto& it : scenarios) {
    run_single_scenario<Tag, ValueType>(it);
  }
}

TEST(std_algorithms_transform_ops_test, transform_unary_op) {
  run_all_scenarios<DynamicTag, int>();
  run_all_scenarios<StridedThreeTag, int>();
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedThreeTag, double>();
}

}  // namespace TransformUnaryOp
}  // namespace stdalgos
}  // namespace Test
