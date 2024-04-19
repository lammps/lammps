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
namespace IsSorted {

namespace KE = Kokkos::Experimental;

template <class ViewType>
void fill_view(ViewType dest_view, const std::string& name) {
  using value_type = typename ViewType::value_type;
  using exe_space  = typename ViewType::execution_space;

  const std::size_t ext = dest_view.extent(0);
  using aux_view_t      = Kokkos::View<value_type*, exe_space>;
  aux_view_t aux_view("aux_view", ext);
  auto v_h = create_mirror_view(Kokkos::HostSpace(), aux_view);

  if (name == "empty") {
    // no op
  }

  else if (name == "one-element") {
    v_h(0) = static_cast<value_type>(1);
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
      v_h(i) = static_cast<value_type>(i);
    }
  }

  else if (name == "small-b") {
    for (std::size_t i = 0; i < ext; ++i) {
      v_h(i) = static_cast<value_type>(i);
    }
    v_h(5) = static_cast<value_type>(-15);
  }

  else if (name == "medium-a") {
    for (std::size_t i = 0; i < ext; ++i) {
      v_h(i) = static_cast<value_type>(i);
    }
  }

  else if (name == "medium-b") {
    for (std::size_t i = 0; i < ext; ++i) {
      v_h(i) = static_cast<value_type>(i);
    }
    v_h(4) = static_cast<value_type>(-1);
  }

  else if (name == "large-a") {
    for (std::size_t i = 0; i < ext; ++i) {
      v_h(i) = static_cast<value_type>(-100) + static_cast<value_type>(i);
    }
  }

  else if (name == "large-b") {
    for (std::size_t i = 0; i < ext; ++i) {
      v_h(i) = static_cast<value_type>(-100) + static_cast<value_type>(i);
    }
    v_h(156) = static_cast<value_type>(-250);

  }

  else {
    throw std::runtime_error("invalid choice");
  }

  Kokkos::deep_copy(aux_view, v_h);
  CopyFunctor<aux_view_t, ViewType> F1(aux_view, dest_view);
  Kokkos::parallel_for("copy", dest_view.extent(0), F1);
}

bool compute_gold(const std::string& name) {
  if (name == "empty") {
    return true;
  } else if (name == "one-element") {
    return true;
  } else if (name == "two-elements-a") {
    return true;
  } else if (name == "two-elements-b") {
    return false;
  } else if (name == "small-a") {
    return true;
  } else if (name == "small-b") {
    return false;
  } else if (name == "medium-a") {
    return true;
  } else if (name == "medium-b") {
    return false;
  } else if (name == "large-a") {
    return true;
  } else if (name == "large-b") {
    return false;
  } else {
    throw std::runtime_error("invalid choice");
  }
}

template <class Tag, class ValueType, class InfoType>
void run_single_scenario(const InfoType& scenario_info) {
  const auto name            = std::get<0>(scenario_info);
  const std::size_t view_ext = std::get<1>(scenario_info);

  // std::cout << "is-sorted: " << name << ", " << view_tag_to_string(Tag{})
  //           << std::endl;

  auto view = create_view<ValueType>(Tag{}, view_ext, "is_sorted");
  fill_view(view, name);
  const auto gold = compute_gold(name);

  std::vector<bool> resultsA(4);
  resultsA[0] = KE::is_sorted(exespace(), KE::cbegin(view), KE::cend(view));
  resultsA[1] =
      KE::is_sorted("label", exespace(), KE::cbegin(view), KE::cend(view));
  resultsA[2]     = KE::is_sorted(exespace(), view);
  resultsA[3]     = KE::is_sorted("label", exespace(), view);
  const auto allA = std::all_of(resultsA.cbegin(), resultsA.cend(),
                                [=](bool v) { return v == gold; });
  EXPECT_TRUE(allA) << name << ", " << view_tag_to_string(Tag{});

#if !defined KOKKOS_ENABLE_OPENMPTARGET
  CustomLessThanComparator<ValueType, ValueType> comp;
  std::vector<bool> resultsB(4);
  resultsB[0] =
      KE::is_sorted(exespace(), KE::cbegin(view), KE::cend(view), comp);
  resultsB[1]     = KE::is_sorted("label", exespace(), KE::cbegin(view),
                              KE::cend(view), comp);
  resultsB[2]     = KE::is_sorted(exespace(), view, comp);
  resultsB[3]     = KE::is_sorted("label", exespace(), view, comp);
  const auto allB = std::all_of(resultsB.cbegin(), resultsB.cend(),
                                [=](bool v) { return v == gold; });
  EXPECT_TRUE(allB) << name << ", " << view_tag_to_string(Tag{});
#endif

  Kokkos::fence();
}

template <class Tag, class ValueType>
void run_is_sorted_all_scenarios() {
  const std::map<std::string, std::size_t> scenarios = {
      {"empty", 0},          {"one-element", 1}, {"two-elements-a", 2},
      {"two-elements-b", 2}, {"small-a", 9},     {"small-b", 13},
      {"medium-a", 1003},    {"medium-b", 1003}, {"large-a", 101513},
      {"large-b", 101513}};

  for (const auto& it : scenarios) {
    run_single_scenario<Tag, ValueType>(it);
  }
}

TEST(std_algorithms_sorting_ops_test, is_sorted) {
  run_is_sorted_all_scenarios<DynamicTag, double>();
  run_is_sorted_all_scenarios<StridedTwoTag, double>();
  run_is_sorted_all_scenarios<StridedThreeTag, double>();
}

}  // namespace IsSorted
}  // namespace stdalgos
}  // namespace Test
