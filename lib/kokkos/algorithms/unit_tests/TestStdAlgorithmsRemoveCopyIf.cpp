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
#include <algorithm>

namespace Test {
namespace stdalgos {
namespace RemoveCopyIf {

namespace KE = Kokkos::Experimental;

template <class ValueType>
struct UnifDist;

template <>
struct UnifDist<int> {
  using dist_type = std::uniform_int_distribution<int>;
  std::mt19937 m_gen;
  dist_type m_dist;

  UnifDist() : m_dist(-100, 100) { m_gen.seed(1034343); }

  int operator()() { return m_dist(m_gen); }
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
    for (std::size_t i = 0; i < ext; ++i) {
      v_h(i) = static_cast<value_type>(i + 1);
    }
  }

  else if (name == "small-b") {
    for (std::size_t i = 0; i < ext; ++i) {
      if (i % 2 == 0) {
        v_h(i) = static_cast<value_type>(22);
      } else {
        v_h(i) = static_cast<value_type>(-12);
      }
    }
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

template <class ViewTypeFrom, class ViewTypeDest, class MyItResult,
          class PredicateType>
void verify_data(ViewTypeFrom view_from, ViewTypeDest view_dest,
                 MyItResult my_result, PredicateType pred) {
  // make a host copy of the view_from
  auto view_from_h      = create_host_space_copy(view_from);
  const std::size_t ext = view_from_h.extent(0);
  using value_type      = typename ViewTypeFrom::value_type;

  // run std::remove_copy_if
  std::vector<value_type> gold_dest_std(ext);
  auto std_result =
      std::remove_copy_if(KE::cbegin(view_from_h), KE::cend(view_from_h),
                          gold_dest_std.begin(), pred);

  // check that returned iterators are correct
  const std::size_t std_diff = std_result - gold_dest_std.begin();
  const std::size_t my_diff  = my_result - KE::begin(view_dest);
  ASSERT_EQ(std_diff, my_diff);

  // check the actual data after algo has been applied
  auto view_dest_h = create_host_space_copy(view_dest);
  for (std::size_t i = 0; i < my_diff; ++i) {
    ASSERT_EQ(view_dest_h(i), gold_dest_std[i]);
    // std::cout << "i= " << i << " "
    // 	      << "mine: " << view_dest_h(i) << " "
    // 	      << "std: " << gold_dest_std[i]
    // 	      << '\n';
  }
}

std::string value_type_to_string(int) { return "int"; }
std::string value_type_to_string(double) { return "double"; }

template <class Tag, class ValueType, class InfoType>
void run_single_scenario(const InfoType& scenario_info) {
  const auto name            = std::get<0>(scenario_info);
  const std::size_t view_ext = std::get<1>(scenario_info);
  // std::cout << "remove_copy_if: " << name << ", " <<
  // view_tag_to_string(Tag{})
  //           << ", " << value_type_to_string(ValueType()) << std::endl;

  using pred_type = IsEvenFunctor<ValueType>;
  pred_type remove_if_even;

  {
    auto view_from =
        create_view<ValueType>(Tag{}, view_ext, "remove_copy_if_view_from");
    auto view_dest =
        create_view<ValueType>(Tag{}, view_ext, "remove_copy_if_view_dest");
    fill_view(view_from, name);
    auto rit = KE::remove_copy_if(exespace(), KE::cbegin(view_from),
                                  KE::cend(view_from), KE::begin(view_dest),
                                  remove_if_even);
    verify_data(view_from, view_dest, rit, remove_if_even);
  }

  {
    auto view_from =
        create_view<ValueType>(Tag{}, view_ext, "remove_copy_if_view_from");
    auto view_dest =
        create_view<ValueType>(Tag{}, view_ext, "remove_copy_if_view_dest");
    fill_view(view_from, name);
    auto rit = KE::remove_copy_if("label", exespace(), KE::cbegin(view_from),
                                  KE::cend(view_from), KE::begin(view_dest),
                                  remove_if_even);
    verify_data(view_from, view_dest, rit, remove_if_even);
  }

  {
    auto view_from =
        create_view<ValueType>(Tag{}, view_ext, "remove_copy_if_view_from");
    auto view_dest =
        create_view<ValueType>(Tag{}, view_ext, "remove_copy_if_view_dest");
    fill_view(view_from, name);
    auto rit =
        KE::remove_copy_if(exespace(), view_from, view_dest, remove_if_even);
    verify_data(view_from, view_dest, rit, remove_if_even);
  }

  {
    auto view_from =
        create_view<ValueType>(Tag{}, view_ext, "remove_copy_if_view_from");
    auto view_dest =
        create_view<ValueType>(Tag{}, view_ext, "remove_copy_if_view_dest");
    fill_view(view_from, name);
    auto rit = KE::remove_copy_if("label", exespace(), view_from, view_dest,
                                  remove_if_even);
    verify_data(view_from, view_dest, rit, remove_if_even);
  }

  Kokkos::fence();
}

template <class Tag, class ValueType>
void run_all_scenarios() {
  const std::map<std::string, std::size_t> scenarios = {
      {"empty", 0},          {"one-element-a", 1},  {"one-element-b", 1},
      {"two-elements-a", 2}, {"two-elements-b", 2}, {"small-a", 9},
      {"small-b", 13},       {"medium", 23103},     {"large", 101513}};

  for (const auto& it : scenarios) {
    run_single_scenario<Tag, ValueType>(it);
  }
}

TEST(std_algorithms_mod_seq_ops, remove_copy_if) {
  run_all_scenarios<DynamicTag, int>();
  run_all_scenarios<StridedThreeTag, int>();
}

}  // namespace RemoveCopyIf
}  // namespace stdalgos
}  // namespace Test
