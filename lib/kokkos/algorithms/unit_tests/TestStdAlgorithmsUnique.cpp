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
namespace Unique {

namespace KE = Kokkos::Experimental;

template <class ValueType>
struct UnifDist;

template <>
struct UnifDist<int> {
  using dist_type = std::uniform_int_distribution<int>;
  std::mt19937 m_gen;
  dist_type m_dist;

  UnifDist() : m_dist(2, 7) { m_gen.seed(1034343); }
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
    v_h(0)  = static_cast<value_type>(0);
    v_h(1)  = static_cast<value_type>(1);
    v_h(2)  = static_cast<value_type>(3);
    v_h(3)  = static_cast<value_type>(2);
    v_h(4)  = static_cast<value_type>(2);
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

  else if (name == "medium-a") {
    // beginning just contains increasing values
    for (std::size_t i = 0; i < 1000; ++i) {
      v_h(i) = static_cast<value_type>(i);
    }

    // then use random
    UnifDist<value_type> randObj;
    for (std::size_t i = 1000; i < ext; ++i) {
      v_h(i) = randObj();
    }
  }

  else if (name == "medium-b") {
    for (std::size_t i = 0; i < 1000; ++i) {
      v_h(i) = static_cast<value_type>(22);
    }
    for (std::size_t i = 1000; i < ext; ++i) {
      v_h(i) = static_cast<value_type>(44);
    }
  }

  else if (name == "large-a") {
    // put equal elements at the end
    for (std::size_t i = 0; i < ext; ++i) {
      v_h(i) = static_cast<value_type>(i);
    }
    v_h(ext - 3) = static_cast<value_type>(44);
    v_h(ext - 2) = static_cast<value_type>(44);
    v_h(ext - 1) = static_cast<value_type>(44);
  }

  else if (name == "large-b") {
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

template <class ResultIt, class ViewType1, class ViewType2, class... Args>
void verify_data(const std::string& name, ResultIt my_result_it,
                 ViewType1 view_test, ViewType2 data_v_h, Args... args) {
  // run std unique on host data
  auto std_r = std::unique(KE::begin(data_v_h), KE::end(data_v_h), args...);

  //
  // check the returned iterator is correct
  //
  const auto std_diff = (std::size_t)(std_r - KE::begin(data_v_h));
  const auto my_diff  = (std::size_t)(my_result_it - KE::begin(view_test));
  EXPECT_TRUE(my_diff == std_diff);

  //
  // check the data in the view
  //
  // Note that we need to stop after std_diff because
  // after that values are unspecified, see std::unique
  auto view_test_h = create_host_space_copy(view_test);
  for (std::size_t i = 0; i < std_diff; ++i) {
    // std::cout << "i = " << i
    // 		<< " my  = " << view_test_h(i) << " "
    // 		<< " std = " << data_v_h(i)
    // 		<< '\n';
    EXPECT_TRUE(view_test_h(i) == data_v_h(i));
  }

  if (name == "medium-b") {
    using value_type = typename ViewType1::value_type;
    EXPECT_TRUE(my_diff == (std::size_t)2);
    EXPECT_TRUE(view_test_h(0) == (value_type)22);
    EXPECT_TRUE(view_test_h(1) == (value_type)44);
  }
}

std::string value_type_to_string(int) { return "int"; }
std::string value_type_to_string(double) { return "double"; }

template <class Tag, class ValueType>
void print_scenario_details(const std::string& name) {
  std::cout << "unique: default predicate: " << name << ", "
            << view_tag_to_string(Tag{}) << " "
            << value_type_to_string(ValueType()) << '\n';
}

template <class Tag, class ValueType, class Predicate>
void print_scenario_details(const std::string& name, Predicate pred) {
  (void)pred;
  std::cout << "unique: custom  predicate: " << name << ", "
            << view_tag_to_string(Tag{}) << " "
            << value_type_to_string(ValueType()) << '\n';
}

template <class Tag, class ValueType, class InfoType, class... Args>
void run_single_scenario(const InfoType& scenario_info, Args... args) {
  const auto name            = std::get<0>(scenario_info);
  const std::size_t view_ext = std::get<1>(scenario_info);
  // print_scenario_details<Tag, ValueType>(name, args...);

  auto test_view = create_view<ValueType>(Tag{}, view_ext, "unique_test_view");

  {
    fill_view(test_view, name);
    // make host copy BEFORE running unique on it since unique modifies it
    auto data_h = create_host_space_copy(test_view);

    // run unique
    auto rit = KE::unique(exespace(), KE::begin(test_view), KE::end(test_view),
                          args...);
    // verify
    verify_data(name, rit, test_view, data_h, args...);
  }

  {
    fill_view(test_view, name);
    // make host copy BEFORE running unique on it since unique modifies it
    auto data_h = create_host_space_copy(test_view);

    // run unique
    auto rit = KE::unique("label", exespace(), KE::begin(test_view),
                          KE::end(test_view), args...);
    // verify
    verify_data(name, rit, test_view, data_h, args...);
  }

  {
    fill_view(test_view, name);
    // make host copy BEFORE running unique on it since unique modifies it
    auto data_h = create_host_space_copy(test_view);

    // run unique
    auto rit = KE::unique(exespace(), test_view, args...);
    // verify
    verify_data(name, rit, test_view, data_h, args...);
  }

  {
    fill_view(test_view, name);
    // make host copy BEFORE running unique on it since unique modifies it
    auto data_h = create_host_space_copy(test_view);

    // run unique
    auto rit = KE::unique("label", exespace(), test_view, args...);
    // verify
    verify_data(name, rit, test_view, data_h, args...);
  }

  Kokkos::fence();
}

template <class Tag, class ValueType>
void run_all_scenarios() {
  const std::map<std::string, std::size_t> scenarios = {
      {"empty", 0},          {"one-element-a", 1},  {"one-element-b", 1},
      {"two-elements-a", 2}, {"two-elements-b", 2}, {"small-a", 11},
      {"small-b", 13},       {"medium-a", 11103},   {"medium-b", 21103},
      {"large-a", 101513},   {"large-b", 100111}};

  for (const auto& it : scenarios) {
    run_single_scenario<Tag, ValueType>(it);

    using func_t = IsEqualFunctor<ValueType>;
    run_single_scenario<Tag, ValueType>(it, func_t());
  }
}

TEST(std_algorithms_mod_seq_ops, unique) {
  run_all_scenarios<DynamicTag, int>();
  run_all_scenarios<StridedThreeTag, int>();
}

}  // namespace Unique
}  // namespace stdalgos
}  // namespace Test
