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

namespace Test {
namespace stdalgos {
namespace ReplaceIf {

namespace KE = Kokkos::Experimental;

template <class ValueType>
struct UnifDist;

template <>
struct UnifDist<double> {
  using dist_type = std::uniform_real_distribution<double>;
  std::mt19937 m_gen;
  dist_type m_dist;

  UnifDist() : m_dist(-20., 20.) { m_gen.seed(1034343); }

  double operator()() { return m_dist(m_gen); }
};

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

  UnifDist<value_type> randObj;
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
      v_h(i) = value_type{-5} + static_cast<value_type>(i + 1);
    }
  }

  else if (name == "small-b") {
    for (std::size_t i = 0; i < ext; ++i) {
      v_h(i) = randObj();
    }
    v_h(5) = static_cast<value_type>(-2);
  }

  else if (name == "medium" || name == "large") {
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

// my own because std::replace_if is ONLY found with std=c++20
template <class ForwardIt, class UnaryPredicate, class T>
void my_host_replace_if(ForwardIt first, ForwardIt last, UnaryPredicate p,
                        const T& new_value) {
  for (; first != last; ++first) {
    if (p(*first)) {
      *first = new_value;
    }
  }
}

template <class ViewType1, class ViewType2, class ValueType,
          class PredicateType>
void verify_data(ViewType1 data_view,  // contains data
                 ViewType2 test_view,  // the view to test
                 ValueType new_value, PredicateType pred) {
  //! always careful because views might not be deep copyable

  auto data_view_dc = create_deep_copyable_compatible_clone(data_view);
  auto data_view_h =
      create_mirror_view_and_copy(Kokkos::HostSpace(), data_view_dc);
  my_host_replace_if(KE::begin(data_view_h), KE::end(data_view_h), pred,
                     new_value);

  auto test_view_dc = create_deep_copyable_compatible_clone(test_view);
  auto test_view_h =
      create_mirror_view_and_copy(Kokkos::HostSpace(), test_view_dc);

  if (test_view_h.extent(0) > 0) {
    for (std::size_t i = 0; i < test_view_h.extent(0); ++i) {
      // std::cout << i << " " << std::setprecision(15)
      // 		<< data_view_dc(i) << " "
      // 		<< data_view_h(i) << " "
      // 		<< test_view_h(i) << std::endl;
      EXPECT_TRUE(data_view_h(i) == test_view_h(i));
    }
  }
}

std::string value_type_to_string(int) { return "int"; }
std::string value_type_to_string(double) { return "double"; }

template <class Tag, class ValueType, class InfoType, class PredicateType>
void run_single_scenario(const InfoType& scenario_info, PredicateType pred) {
  const auto name            = std::get<0>(scenario_info);
  const std::size_t view_ext = std::get<1>(scenario_info);
  // std::cout << "replace_if: " << name << ", " << view_tag_to_string(Tag{})
  //           << ", " << value_type_to_string(ValueType()) << std::endl;

  ValueType new_value{23};
  auto view_with_data =
      create_view<ValueType>(Tag{}, view_ext, "replace_if_v2");
  auto view_to_test = create_view<ValueType>(Tag{}, view_ext, "replace_if_v1");
  fill_view(view_with_data, name);

  {
    CopyFunctor<decltype(view_with_data), decltype(view_to_test)> F1(
        view_with_data, view_to_test);
    Kokkos::parallel_for("copy", view_to_test.extent(0), F1);

    KE::replace_if(exespace(), KE::begin(view_to_test), KE::end(view_to_test),
                   pred, new_value);
    verify_data(view_with_data, view_to_test, new_value, pred);
  }

  {
    CopyFunctor<decltype(view_with_data), decltype(view_to_test)> F1(
        view_with_data, view_to_test);
    Kokkos::parallel_for("copy", view_to_test.extent(0), F1);

    KE::replace_if("label", exespace(), KE::begin(view_to_test),
                   KE::end(view_to_test), pred, new_value);
    verify_data(view_with_data, view_to_test, new_value, pred);
  }

  {
    CopyFunctor<decltype(view_with_data), decltype(view_to_test)> F1(
        view_with_data, view_to_test);
    Kokkos::parallel_for("copy", view_to_test.extent(0), F1);

    KE::replace_if(exespace(), view_to_test, pred, new_value);
    verify_data(view_with_data, view_to_test, new_value, pred);
  }

  {
    CopyFunctor<decltype(view_with_data), decltype(view_to_test)> F1(
        view_with_data, view_to_test);
    Kokkos::parallel_for("copy", view_to_test.extent(0), F1);

    KE::replace_if("label", exespace(), view_to_test, pred, new_value);
    verify_data(view_with_data, view_to_test, new_value, pred);
  }

  Kokkos::fence();
}

template <class Tag, class ValueType>
void run_all_scenarios() {
  const std::map<std::string, std::size_t> scenarios = {
      {"empty", 0},          {"one-element", 1}, {"two-elements-a", 2},
      {"two-elements-b", 2}, {"small-a", 9},     {"small-b", 13},
      {"medium", 1103},      {"large", 101513}};

  for (const auto& it : scenarios) {
    using pred_p_t = IsPositiveFunctor<ValueType>;
    run_single_scenario<Tag, ValueType>(it, pred_p_t{});
    using pred_n_t = IsNegativeFunctor<ValueType>;
    run_single_scenario<Tag, ValueType>(it, pred_n_t{});
  }
}

TEST(std_algorithms_replace_ops_test, replace_if) {
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedThreeTag, double>();
  run_all_scenarios<DynamicTag, int>();
  run_all_scenarios<StridedThreeTag, int>();
}

}  // namespace ReplaceIf
}  // namespace stdalgos
}  // namespace Test
