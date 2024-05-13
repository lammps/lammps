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
namespace CopyIf {

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

template <class ViewType, class PredicateType>
std::size_t fill_view(ViewType dest_view, const std::string& name,
                      PredicateType pred) {
  using value_type      = typename ViewType::value_type;
  using exe_space       = typename ViewType::execution_space;
  const std::size_t ext = dest_view.extent(0);
  using aux_view_t      = Kokkos::View<value_type*, exe_space>;
  aux_view_t aux_view("aux_view", ext);
  auto v_h = create_mirror_view(Kokkos::HostSpace(), aux_view);

  std::size_t count = 0;
  if (name == "empty") {
    // no op
  }

  else if (name == "one-element-a") {
    v_h(0) = static_cast<value_type>(1);
    // 1 is not even, so count is not incremented
  }

  else if (name == "one-element-b") {
    v_h(0) = static_cast<value_type>(2);
    count++;
  }

  else if (name == "two-elements-a") {
    v_h(0) = static_cast<value_type>(1);
    v_h(1) = static_cast<value_type>(2);
    count++;
  }

  else if (name == "two-elements-b") {
    v_h(0) = static_cast<value_type>(2);
    v_h(1) = static_cast<value_type>(-1);
    count++;
  }

  else if (name == "small-a") {
    for (std::size_t i = 0; i < ext; ++i) {
      v_h(i) = value_type{-5} + static_cast<value_type>(i + 1);
      if (pred(v_h(i))) {
        count++;
      }
    }
  }

  else if (name == "small-b") {
    for (std::size_t i = 0; i < ext; ++i) {
      if (i % 2 == 0) {
        v_h(i) = static_cast<value_type>(22);
      } else {
        v_h(i) = static_cast<value_type>(-12);
      }
      if (pred(v_h(i))) {
        count++;
      }
    }
  }

  else if (name == "medium" || name == "large") {
    UnifDist<value_type> randObj;
    for (std::size_t i = 0; i < ext; ++i) {
      v_h(i) = randObj();
      if (pred(v_h(i))) {
        count++;
      }
    }
  }

  else {
    throw std::runtime_error("invalid choice");
  }

  Kokkos::deep_copy(aux_view, v_h);
  CopyFunctor<aux_view_t, ViewType> F1(aux_view, dest_view);
  Kokkos::parallel_for("copy", dest_view.extent(0), F1);
  return count;
}

template <class ViewTypeFrom, class ViewTypeTest, class PredType>
void verify_data(const std::string& name, ViewTypeFrom view_from,
                 ViewTypeTest view_test, PredType pred) {
  using value_type = typename ViewTypeTest::value_type;

  //! always careful because views might not be deep copyable
  auto view_test_dc = create_deep_copyable_compatible_clone(view_test);
  auto view_test_h =
      create_mirror_view_and_copy(Kokkos::HostSpace(), view_test_dc);

  auto view_from_dc = create_deep_copyable_compatible_clone(view_from);
  auto view_from_h =
      create_mirror_view_and_copy(Kokkos::HostSpace(), view_from_dc);

  if (name == "empty") {
    // no op
  }

  else if (name == "one-element-a") {
    ASSERT_EQ(view_test_h(0), static_cast<value_type>(0));
  }

  else if (name == "one-element-b") {
    ASSERT_EQ(view_test_h(0), static_cast<value_type>(2));
  }

  else if (name == "two-elements-a") {
    ASSERT_EQ(view_test_h(0), static_cast<value_type>(2));
    ASSERT_EQ(view_test_h(1), static_cast<value_type>(0));
  }

  else if (name == "two-elements-b") {
    ASSERT_EQ(view_test_h(0), static_cast<value_type>(2));
    ASSERT_EQ(view_test_h(1), static_cast<value_type>(0));
  }

  else if (name == "small-a") {
    ASSERT_EQ(view_test_h(0), static_cast<value_type>(-4));
    ASSERT_EQ(view_test_h(1), static_cast<value_type>(-2));
    ASSERT_EQ(view_test_h(2), static_cast<value_type>(0));
    ASSERT_EQ(view_test_h(3), static_cast<value_type>(2));
    ASSERT_EQ(view_test_h(4), static_cast<value_type>(4));
    ASSERT_EQ(view_test_h(5), static_cast<value_type>(0));
    ASSERT_EQ(view_test_h(6), static_cast<value_type>(0));
    ASSERT_EQ(view_test_h(7), static_cast<value_type>(0));
    ASSERT_EQ(view_test_h(8), static_cast<value_type>(0));
  }

  else if (name == "small-b") {
    ASSERT_EQ(view_test_h(0), static_cast<value_type>(22));
    ASSERT_EQ(view_test_h(1), static_cast<value_type>(-12));
    ASSERT_EQ(view_test_h(2), static_cast<value_type>(22));
    ASSERT_EQ(view_test_h(3), static_cast<value_type>(-12));
    ASSERT_EQ(view_test_h(4), static_cast<value_type>(22));
    ASSERT_EQ(view_test_h(5), static_cast<value_type>(-12));
    ASSERT_EQ(view_test_h(6), static_cast<value_type>(22));
    ASSERT_EQ(view_test_h(7), static_cast<value_type>(-12));
    ASSERT_EQ(view_test_h(8), static_cast<value_type>(22));
    ASSERT_EQ(view_test_h(9), static_cast<value_type>(-12));
    ASSERT_EQ(view_test_h(10), static_cast<value_type>(22));
    ASSERT_EQ(view_test_h(11), static_cast<value_type>(-12));
    ASSERT_EQ(view_test_h(12), static_cast<value_type>(22));
  }

  else if (name == "medium" || name == "large") {
    // for (std::size_t i = 0; i < view_from_h.extent(0); ++i){
    //   std::cout << "i= " << i << " "
    // 		<< " vf = " << view_from_h(i) << " "
    // 		<< " vt = " << view_test_h(i) << '\n';
    // }

    std::size_t count = 0;
    for (std::size_t i = 0; i < view_from_h.extent(0); ++i) {
      if (pred(view_from_h(i))) {
        ASSERT_EQ(view_test_h(count), view_from_h(i));
        count++;
      }
    }
    // all other entries of test view should be zero
    for (; count < view_test_h.extent(0); ++count) {
      // std::cout << count << '\n';
      ASSERT_EQ(view_test_h(count), value_type(0));
    }
  }

  else {
    throw std::runtime_error("invalid choice");
  }
}

std::string value_type_to_string(int) { return "int"; }
std::string value_type_to_string(double) { return "double"; }

template <class Tag, class ValueType, class InfoType>
void run_single_scenario(const InfoType& scenario_info) {
  const auto name            = std::get<0>(scenario_info);
  const std::size_t view_ext = std::get<1>(scenario_info);
  // std::cout << "copy_if: " << name << ", " << view_tag_to_string(Tag{}) << ",
  // "
  //           << value_type_to_string(ValueType()) << std::endl;

  auto view_from = create_view<ValueType>(Tag{}, view_ext, "copy_if_from");
  IsEvenFunctor<ValueType> pred;

  {
    auto n         = fill_view(view_from, name, pred);
    auto view_dest = create_view<ValueType>(Tag{}, view_ext, "copy_if_dest");
    auto rit       = KE::copy_if(exespace(), KE::cbegin(view_from),
                           KE::cend(view_from), KE::begin(view_dest), pred);
    verify_data(name, view_from, view_dest, pred);
    ASSERT_EQ(rit, (KE::begin(view_dest) + n));
  }

  {
    auto n         = fill_view(view_from, name, pred);
    auto view_dest = create_view<ValueType>(Tag{}, view_ext, "copy_if_dest");
    auto rit       = KE::copy_if("label", exespace(), KE::cbegin(view_from),
                           KE::cend(view_from), KE::begin(view_dest), pred);
    verify_data(name, view_from, view_dest, pred);
    ASSERT_EQ(rit, (KE::begin(view_dest) + n));
  }

  {
    auto n         = fill_view(view_from, name, pred);
    auto view_dest = create_view<ValueType>(Tag{}, view_ext, "copy_if_dest");
    auto rit       = KE::copy_if(exespace(), view_from, view_dest, pred);
    verify_data(name, view_from, view_dest, pred);
    ASSERT_EQ(rit, (KE::begin(view_dest) + n));
  }

  {
    auto n         = fill_view(view_from, name, pred);
    auto view_dest = create_view<ValueType>(Tag{}, view_ext, "copy_if_dest");
    auto rit = KE::copy_if("label", exespace(), view_from, view_dest, pred);
    verify_data(name, view_from, view_dest, pred);
    ASSERT_EQ(rit, (KE::begin(view_dest) + n));
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

TEST(std_algorithms_mod_seq_ops, copy_if) {
  run_all_scenarios<DynamicTag, int>();
  run_all_scenarios<StridedThreeTag, int>();
}

}  // namespace CopyIf
}  // namespace stdalgos
}  // namespace Test
