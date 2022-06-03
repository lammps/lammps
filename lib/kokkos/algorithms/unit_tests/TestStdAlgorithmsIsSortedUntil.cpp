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
#include <std_algorithms/Kokkos_SortingOperations.hpp>
#include <utility>

namespace Test {
namespace stdalgos {
namespace IsSortedUntil {

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
    v_h(5) = static_cast<value_type>(15);
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

template <class ViewType>
auto compute_gold(ViewType view, const std::string& name) {
  if (name == "empty") {
    return KE::end(view);
  } else if (name == "one-element") {
    return KE::end(view);
  } else if (name == "two-elements-a") {
    return KE::end(view);
  } else if (name == "two-elements-b") {
    return KE::begin(view) + 1;
  } else if (name == "small-a") {
    return KE::end(view);
  } else if (name == "small-b") {
    return KE::begin(view) + 6;
  } else if (name == "medium-a") {
    return KE::end(view);
  } else if (name == "medium-b") {
    return KE::begin(view) + 4;
  } else if (name == "large-a") {
    return KE::end(view);
  } else if (name == "large-b") {
    return KE::begin(view) + 156;
  } else {
    throw std::runtime_error("invalid choice");
  }
}

template <class Tag, class ValueType, class InfoType>
void run_single_scenario(const InfoType& scenario_info) {
  const auto name            = std::get<0>(scenario_info);
  const std::size_t view_ext = std::get<1>(scenario_info);

  // std::cout << "is-sorted-until: " << name << ", " <<
  // view_tag_to_string(Tag{})
  //           << std::endl;

  auto view = create_view<ValueType>(Tag{}, view_ext, "is_sorted_until");
  fill_view(view, name);
  const auto gold = compute_gold(view, name);

  auto r1 = KE::is_sorted_until(exespace(), KE::begin(view), KE::end(view));
  auto r2 =
      KE::is_sorted_until("label", exespace(), KE::begin(view), KE::end(view));
  auto r3 = KE::is_sorted_until(exespace(), view);
  auto r4 = KE::is_sorted_until("label", exespace(), view);
  EXPECT_TRUE(r1 == gold);
  EXPECT_TRUE(r2 == gold);
  EXPECT_TRUE(r3 == gold);
  EXPECT_TRUE(r4 == gold);

#if not defined KOKKOS_ENABLE_OPENMPTARGET
  CustomLessThanComparator<ValueType, ValueType> comp;
  auto r5 =
      KE::is_sorted_until(exespace(), KE::cbegin(view), KE::cend(view), comp);
  auto r6 = KE::is_sorted_until("label", exespace(), KE::cbegin(view),
                                KE::cend(view), comp);
  auto r7 = KE::is_sorted_until(exespace(), view, comp);
  auto r8 = KE::is_sorted_until("label", exespace(), view, comp);
#endif

  EXPECT_TRUE(r1 == gold);
  EXPECT_TRUE(r2 == gold);
  EXPECT_TRUE(r3 == gold);
  EXPECT_TRUE(r4 == gold);

  Kokkos::fence();
}

template <class Tag, class ValueType>
void run_is_sorted_until_all_scenarios() {
  const std::map<std::string, std::size_t> scenarios = {
      {"empty", 0},          {"one-element", 1}, {"two-elements-a", 2},
      {"two-elements-b", 2}, {"small-a", 9},     {"small-b", 13},
      {"medium-a", 1003},    {"medium-b", 1003}, {"large-a", 101513},
      {"large-b", 101513}};

  std::cout << "is_sorted_until: " << view_tag_to_string(Tag{})
            << ", all overloads \n";

  for (const auto& it : scenarios) {
    run_single_scenario<Tag, ValueType>(it);
  }
}

TEST(std_algorithms_sorting_ops_test, is_sorted_until) {
  run_is_sorted_until_all_scenarios<DynamicTag, double>();
  run_is_sorted_until_all_scenarios<StridedTwoTag, double>();
  run_is_sorted_until_all_scenarios<StridedThreeTag, double>();
}

}  // namespace IsSortedUntil
}  // namespace stdalgos
}  // namespace Test
