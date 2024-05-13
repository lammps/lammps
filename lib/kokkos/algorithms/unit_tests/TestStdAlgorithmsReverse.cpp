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

#include <TestStdAlgorithmsCommon.hpp>
#include <utility>

namespace Test {
namespace stdalgos {
namespace Reverse {

namespace KE = Kokkos::Experimental;

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

  else if (name == "small-a" || name == "small-b" || name == "medium" ||
           name == "large") {
    for (std::size_t i = 0; i < ext; ++i) {
      v_h(i) = static_cast<value_type>(-11) + static_cast<value_type>(i);
    }
  }

  else {
    throw std::runtime_error("invalid choice");
  }

  Kokkos::deep_copy(aux_view, v_h);
  CopyFunctor<aux_view_t, ViewType> F1(aux_view, dest_view);
  Kokkos::parallel_for("copy", dest_view.extent(0), F1);
}

template <class ViewType1, class ViewType2>
void verify_data(ViewType1 test_view, ViewType2 orig_view) {
  auto tv_h = create_host_space_copy(test_view);
  auto ov_h = create_host_space_copy(orig_view);

  const std::size_t ext = test_view.extent(0);
  for (std::size_t i = 0; i < ext; ++i) {
    ASSERT_EQ(tv_h(i), ov_h(ext - i - 1));
  }
}

std::string value_type_to_string(int) { return "int"; }
std::string value_type_to_string(double) { return "double"; }

template <class Tag, class ValueType, class InfoType>
void run_single_scenario(const InfoType& scenario_info) {
  const auto name            = std::get<0>(scenario_info);
  const std::size_t view_ext = std::get<1>(scenario_info);
  // std::cout << "reverse: " << name << ", " << view_tag_to_string(Tag{}) << ",
  // "
  //           << value_type_to_string(ValueType()) << std::endl;

  auto test_view = create_view<ValueType>(Tag{}, view_ext, "reverse");
  auto orig_view = create_view<ValueType>(Tag{}, view_ext, "reverse");

  {
    fill_view(test_view, name);
    fill_view(orig_view, name);
    KE::reverse(exespace(), KE::begin(test_view), KE::end(test_view));
    verify_data(test_view, orig_view);
  }

  {
    fill_view(test_view, name);
    fill_view(orig_view, name);
    KE::reverse("label", exespace(), KE::begin(test_view), KE::end(test_view));
    verify_data(test_view, orig_view);
  }

  {
    fill_view(test_view, name);
    fill_view(orig_view, name);
    KE::reverse(exespace(), test_view);
    verify_data(test_view, orig_view);
  }

  {
    fill_view(test_view, name);
    fill_view(orig_view, name);
    KE::reverse("label", exespace(), test_view);
    verify_data(test_view, orig_view);
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

TEST(std_algorithms_modseq_test, reverse) {
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedThreeTag, double>();
  run_all_scenarios<DynamicTag, int>();
  run_all_scenarios<StridedThreeTag, int>();
}

}  // namespace Reverse
}  // namespace stdalgos
}  // namespace Test
