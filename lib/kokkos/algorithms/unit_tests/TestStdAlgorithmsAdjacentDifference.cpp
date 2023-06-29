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
#include <std_algorithms/Kokkos_BeginEnd.hpp>
#include <std_algorithms/Kokkos_AdjacentDifference.hpp>
#include <utility>
#include <numeric>

namespace Test {
namespace stdalgos {
namespace AdjacentDifference {

namespace KE = Kokkos::Experimental;

template <class DestViewType>
void fill_view(DestViewType dest_view, const std::string& name) {
  // we need to be careful because dest_view might not be deep copyable
  // for instance strided layout

  using value_type      = typename DestViewType::value_type;
  const std::size_t ext = dest_view.extent(0);
  auto aux_view =
      create_deep_copyable_compatible_view_with_same_extent(dest_view);
  auto aux_v_h = create_mirror_view(Kokkos::HostSpace(), aux_view);

  if (name == "empty") {
    // no op
  }

  else if (name == "one-element") {
    aux_v_h(0) = static_cast<value_type>(1);
  }

  else if (name == "two-elements-a") {
    aux_v_h(0) = static_cast<value_type>(1);
    aux_v_h(1) = static_cast<value_type>(2);
  }

  else if (name == "two-elements-b") {
    aux_v_h(0) = static_cast<value_type>(2);
    aux_v_h(1) = static_cast<value_type>(-1);
  }

  else if (name == "small-a") {
    for (std::size_t i = 0; i < ext; ++i) {
      aux_v_h(i) = static_cast<value_type>(i) * 2;
    }
  }

  else if (name == "small-b") {
    for (std::size_t i = 0; i < ext; ++i) {
      aux_v_h(i) = static_cast<value_type>(i) * 3;
    }
    aux_v_h(5) = static_cast<value_type>(-15);
  }

  else if (name == "medium-a") {
    for (std::size_t i = 0; i < ext; ++i) {
      aux_v_h(i) = static_cast<value_type>(i) * 2;
    }
  }

  else if (name == "medium-b") {
    for (std::size_t i = 0; i < ext; ++i) {
      aux_v_h(i) = static_cast<value_type>(i) * 2;
    }
    aux_v_h(4) = static_cast<value_type>(-1);
  }

  else if (name == "large-a") {
    for (std::size_t i = 0; i < ext; ++i) {
      aux_v_h(i) = static_cast<value_type>(-100) + static_cast<value_type>(i);
    }
  }

  else if (name == "large-b") {
    for (std::size_t i = 0; i < ext; ++i) {
      aux_v_h(i) = static_cast<value_type>(-100) + static_cast<value_type>(i);
    }
    aux_v_h(156) = static_cast<value_type>(-250);

  }

  else {
    throw std::runtime_error("invalid choice");
  }

  Kokkos::deep_copy(aux_view, aux_v_h);
  CopyFunctor<decltype(aux_view), DestViewType> F1(aux_view, dest_view);
  Kokkos::parallel_for("copy", dest_view.extent(0), F1);
}

template <class TestViewType, class... Args>
auto compute_gold(TestViewType test_view, const std::string& name,
                  Args... args /* copy on purpose */) {
  // we need to be careful because test_view might not be deep copyable
  // for instance strided layout

  const std::size_t ext = test_view.extent(0);

  // create a deep copyable clone of test_view
  auto test_view_dc = create_deep_copyable_compatible_clone(test_view);
  auto test_view_dc_h =
      create_mirror_view_and_copy(Kokkos::HostSpace(), test_view_dc);

  // create gold deep copyable view
  auto gold_view =
      create_deep_copyable_compatible_view_with_same_extent(test_view);
  auto gold_view_h = create_mirror_view(Kokkos::HostSpace(), gold_view);

  // compute gold solution on host and deep copy to device
  if (name == "empty") {
    return gold_view;
  } else {
    using value_type = typename TestViewType::value_type;
    std::vector<value_type> tmp(ext);
    for (std::size_t i = 0; i < ext; ++i) {
      tmp[i] = test_view_dc_h(i);
    }
    // run adj-diff on tmp directly
    std::adjacent_difference(tmp.begin(), tmp.end(), tmp.begin(),
                             std::forward<Args>(args)...);

    // copy from tmp to gold_h
    for (std::size_t i = 0; i < ext; ++i) {
      gold_view_h(i) = tmp[i];
    }
    // deep_copy to device
    Kokkos::deep_copy(gold_view, gold_view_h);
    return gold_view;
  }
}

template <class TestViewType, class GoldViewType>
void verify_data(TestViewType test_view, GoldViewType gold) {
  // we need to be careful because test_view might not be deep copyable
  // for instance strided layout

  auto test_view_dc = create_deep_copyable_compatible_clone(test_view);
  auto test_view_dc_h =
      create_mirror_view_and_copy(Kokkos::HostSpace(), test_view_dc);
  // gold is deep_copyable for sure
  const auto gold_h = create_mirror_view_and_copy(Kokkos::HostSpace(), gold);

  for (std::size_t i = 0; i < test_view.extent(0); ++i) {
    ASSERT_EQ(gold_h(i), test_view_dc_h(i));
  }
}

template <class ValueType1, class ValueType2 = ValueType1,
          class RetType = ValueType2>
struct CustomBinaryOpFunctor {
  KOKKOS_INLINE_FUNCTION
  RetType operator()(const ValueType1& a, const ValueType2& b) const {
    return a * b;
  }
};

template <class ValueType1, class ValueType2 = ValueType1,
          class RetType = ValueType2>
struct DefaultBinaryOpFunctor {
  KOKKOS_INLINE_FUNCTION
  RetType operator()(const ValueType1& a, const ValueType2& b) const {
    return a - b;
  }
};

template <class Tag, class ValueType, class InfoType, class... Args>
void run_single_scenario(const InfoType& scenario_info,
                         Args... args /* copy on purpose */) {
  const auto name            = std::get<0>(scenario_info);
  const std::size_t view_ext = std::get<1>(scenario_info);

  auto view_from =
      create_view<ValueType>(Tag{}, view_ext, "adj_diff_from_view");
  fill_view(view_from, name);

  const auto gold = compute_gold(view_from, name, args...);

  {
    auto view_dest =
        create_view<ValueType>(Tag{}, view_ext, "adj_diff_dest_view");
    auto res1 = KE::adjacent_difference(exespace(), KE::cbegin(view_from),
                                        KE::cend(view_from),
                                        KE::begin(view_dest), args...);
    ASSERT_EQ(res1, KE::end(view_dest));
    verify_data(view_dest, gold);
  }

  {
    auto view_dest =
        create_view<ValueType>(Tag{}, view_ext, "adj_diff_dest_view");
    auto res2 = KE::adjacent_difference(
        "label", exespace(), KE::cbegin(view_from), KE::cend(view_from),
        KE::begin(view_dest), args...);
    ASSERT_EQ(res2, KE::end(view_dest));
    verify_data(view_dest, gold);
  }

  {
    auto view_dest =
        create_view<ValueType>(Tag{}, view_ext, "adj_diff_dest_view");
    auto res3 =
        KE::adjacent_difference(exespace(), view_from, view_dest, args...);
    ASSERT_EQ(res3, KE::end(view_dest));
    verify_data(view_dest, gold);
  }

  {
    auto view_dest =
        create_view<ValueType>(Tag{}, view_ext, "adj_diff_dest_view");
    auto res4 = KE::adjacent_difference("label", exespace(), view_from,
                                        view_dest, args...);
    ASSERT_EQ(res4, KE::end(view_dest));
    verify_data(view_dest, gold);
  }

  Kokkos::fence();
}

template <class Tag, class ValueType, class... Args>
void run_all_scenarios(Args... args /* copy on purpose */) {
  // if (0 < sizeof...(args)) {
  //   std::cout << "adjacent_difference: " << view_tag_to_string(Tag{})
  //             << ", custom binary op, all overloads \n";
  // } else {
  //   std::cout << "adjacent_difference: " << view_tag_to_string(Tag{})
  //             << ", default binary op, all overloads \n";
  // }

  for (const auto& it : default_scenarios) {
    run_single_scenario<Tag, ValueType>(it, args...);
  }
}

TEST(std_algorithms_numerics_ops_test, adjecent_difference) {
  using value_type = double;

  run_all_scenarios<DynamicTag, value_type>();
  run_all_scenarios<StridedTwoTag, value_type>();
  run_all_scenarios<StridedThreeTag, value_type>();

  using custom_binary_op = CustomBinaryOpFunctor<value_type>;
  run_all_scenarios<DynamicTag, value_type>(custom_binary_op{});
  run_all_scenarios<StridedTwoTag, value_type>(custom_binary_op{});
  run_all_scenarios<StridedThreeTag, value_type>(custom_binary_op{});
}

}  // namespace AdjacentDifference
}  // namespace stdalgos
}  // namespace Test
