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
#include <std_algorithms/Kokkos_BeginEnd.hpp>
#include <std_algorithms/Kokkos_ModifyingSequenceOperations.hpp>
#include <utility>
#include <algorithm>

namespace Test {
namespace stdalgos {
namespace RotateCopy {

namespace KE = Kokkos::Experimental;

template <class ValueType>
struct UnifDist;

template <>
struct UnifDist<int> {
  using dist_type = std::uniform_int_distribution<int>;
  std::mt19937 m_gen;
  dist_type m_dist;

  UnifDist() : m_dist(-50, 50) { m_gen.seed(1034343); }
  int operator()() { return m_dist(m_gen); }
};

template <>
struct UnifDist<double> {
  using dist_type = std::uniform_real_distribution<double>;
  std::mt19937 m_gen;
  dist_type m_dist;

  UnifDist() : m_dist(-90., 100.) { m_gen.seed(1034343); }

  double operator()() { return m_dist(m_gen); }
};

template <class ViewType>
void fill_view(ViewType dest_view, const std::string& name) {
  using value_type      = typename ViewType::value_type;
  using exe_space       = typename ViewType::execution_space;
  const std::size_t ext = dest_view.extent(0);
  using aux_view_t      = Kokkos::View<value_type*, exe_space>;
  aux_view_t aux_view("aux_view", ext);
  auto v_h = create_mirror_view(Kokkos::HostSpace(), aux_view);

  if (name == "empty") {
    // no op
  }

  else if (name == "one-element-a") {
    v_h(0) = static_cast<value_type>(1);
  }

  else if (name == "one-element-b") {
    v_h(0) = static_cast<value_type>(2);
  }

  else if (name == "two-elements-a") {
    v_h(0) = static_cast<value_type>(1);
    v_h(1) = static_cast<value_type>(2);
  }

  else if (name == "two-elements-b") {
    v_h(0) = static_cast<value_type>(2);
    v_h(1) = static_cast<value_type>(-1);
  }

  else if (name == "small-a") {
    v_h(0)  = static_cast<value_type>(0);
    v_h(1)  = static_cast<value_type>(1);
    v_h(2)  = static_cast<value_type>(1);
    v_h(3)  = static_cast<value_type>(2);
    v_h(4)  = static_cast<value_type>(3);
    v_h(5)  = static_cast<value_type>(4);
    v_h(6)  = static_cast<value_type>(4);
    v_h(7)  = static_cast<value_type>(4);
    v_h(8)  = static_cast<value_type>(5);
    v_h(9)  = static_cast<value_type>(6);
    v_h(10) = static_cast<value_type>(6);
  }

  else if (name == "small-b") {
    v_h(0)  = static_cast<value_type>(1);
    v_h(1)  = static_cast<value_type>(1);
    v_h(2)  = static_cast<value_type>(1);
    v_h(3)  = static_cast<value_type>(2);
    v_h(4)  = static_cast<value_type>(3);
    v_h(5)  = static_cast<value_type>(4);
    v_h(6)  = static_cast<value_type>(4);
    v_h(7)  = static_cast<value_type>(4);
    v_h(8)  = static_cast<value_type>(5);
    v_h(9)  = static_cast<value_type>(6);
    v_h(10) = static_cast<value_type>(8);
    v_h(11) = static_cast<value_type>(9);
    v_h(12) = static_cast<value_type>(8);
  }

  else if (name == "medium" || name == "large") {
    UnifDist<value_type> randObj;
    for (std::size_t i = 0; i < ext; ++i) {
      v_h(i) = randObj();
    }
  }

  else {
    throw std::runtime_error("invalid choice");
  }

  Kokkos::deep_copy(aux_view, v_h);
  CopyFunctor<aux_view_t, ViewType> F1(aux_view, dest_view);
  Kokkos::parallel_for("copy", dest_view.extent(0), F1);
}

template <class ViewTypeFrom, class ViewTypeTest>
void verify_data(ViewTypeFrom view_from, ViewTypeTest view_test,
                 std::size_t rotation_point) {
  auto view_from_h      = create_host_space_copy(view_from);
  auto view_test_h      = create_host_space_copy(view_test);
  const std::size_t ext = view_test_h.extent(0);

  using value_type = typename ViewTypeTest::value_type;
  std::vector<value_type> std_gold_h(ext);
  auto first_n = KE::cbegin(view_from_h) + rotation_point;
  std::rotate_copy(KE::cbegin(view_from_h), first_n, KE::cend(view_from_h),
                   std_gold_h.begin());

  for (std::size_t i = 0; i < ext; ++i) {
    EXPECT_TRUE(view_test_h(i) == std_gold_h[i]);
    // std::cout << "i= " << i << " "
    // 	      << "from: " << view_from_h(i) << " "
    // 	      << "mine: " << view_test_h(i) << " "
    // 	      << "std: " << std_gold_h[i]
    // 	      << '\n';
  }
}

std::string value_type_to_string(int) { return "int"; }
std::string value_type_to_string(double) { return "double"; }

template <class Tag, class ValueType>
void print_scenario_details(const std::string& name,
                            std::size_t rotation_point) {
  std::cout << "rotate_copy: "
            << " at " << rotation_point << ", " << name << ", "
            << view_tag_to_string(Tag{}) << ", "
            << value_type_to_string(ValueType()) << std::endl;
}

template <class Tag, class ValueType, class InfoType>
void run_single_scenario(const InfoType& scenario_info,
                         std::size_t rotation_point) {
  const auto name            = std::get<0>(scenario_info);
  const std::size_t view_ext = std::get<1>(scenario_info);
  // print_scenario_details<Tag, ValueType>(name, rotation_point);

  auto view_from = create_view<ValueType>(Tag{}, view_ext, "rotate_copy_from");
  fill_view(view_from, name);

  {
    auto view_dest =
        create_view<ValueType>(Tag{}, view_ext, "rotate_copy_dest");
    auto n_it = KE::cbegin(view_from) + rotation_point;
    auto rit  = KE::rotate_copy(exespace(), KE::cbegin(view_from), n_it,
                               KE::cend(view_from), KE::begin(view_dest));
    verify_data(view_from, view_dest, rotation_point);
    EXPECT_TRUE(rit == (KE::begin(view_dest) + view_ext));
  }

  {
    auto view_dest =
        create_view<ValueType>(Tag{}, view_ext, "rotate_copy_dest");
    auto n_it = KE::cbegin(view_from) + rotation_point;
    auto rit = KE::rotate_copy("label", exespace(), KE::cbegin(view_from), n_it,
                               KE::cend(view_from), KE::begin(view_dest));
    verify_data(view_from, view_dest, rotation_point);
    EXPECT_TRUE(rit == (KE::begin(view_dest) + view_ext));
  }

  {
    auto view_dest =
        create_view<ValueType>(Tag{}, view_ext, "rotate_copy_dest");
    auto rit =
        KE::rotate_copy(exespace(), view_from, rotation_point, view_dest);
    verify_data(view_from, view_dest, rotation_point);
    EXPECT_TRUE(rit == (KE::begin(view_dest) + view_ext));
  }

  {
    auto view_dest =
        create_view<ValueType>(Tag{}, view_ext, "rotate_copy_dest");
    auto rit = KE::rotate_copy("label", exespace(), view_from, rotation_point,
                               view_dest);
    verify_data(view_from, view_dest, rotation_point);
    EXPECT_TRUE(rit == (KE::begin(view_dest) + view_ext));
  }

  Kokkos::fence();
}

template <class Tag, class ValueType>
void run_all_scenarios() {
  const std::map<std::string, std::size_t> scenarios = {
      {"empty", 0},          {"one-element-a", 1},  {"one-element-b", 1},
      {"two-elements-a", 2}, {"two-elements-b", 2}, {"small-a", 11},
      {"small-b", 13},       {"medium", 21103},     {"large", 101513}};

  std::vector<std::size_t> rotation_points = {0,  1,   2,    3,     8,
                                              56, 101, 1003, 101501};

  for (const auto& it : scenarios) {
    for (const auto& it2 : rotation_points) {
      // for each view scenario, we rotate at multiple points
      // but only if the view has an extent that is >= rotation point
      const auto view_ext = it.second;
      if (view_ext >= it2) {
        run_single_scenario<Tag, ValueType>(it, it2);
      }
    }
  }
}

TEST(std_algorithms_mod_seq_ops, rotate_copy) {
  run_all_scenarios<DynamicTag, int>();
  run_all_scenarios<StridedThreeTag, int>();
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedThreeTag, double>();
}

}  // namespace RotateCopy
}  // namespace stdalgos
}  // namespace Test
