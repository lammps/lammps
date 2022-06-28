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
#include <std_algorithms/Kokkos_Numeric.hpp>
#include <utility>

namespace Test {
namespace stdalgos {
namespace EScan {

namespace KE = Kokkos::Experimental;

template <class ValueType>
struct UnifDist;

template <>
struct UnifDist<double> {
  using dist_type = std::uniform_real_distribution<double>;
  std::mt19937 m_gen;
  dist_type m_dist;

  UnifDist() : m_dist(0.05, 1.2) { m_gen.seed(1034343); }

  double operator()() { return m_dist(m_gen); }
};

template <>
struct UnifDist<int> {
  using dist_type = std::uniform_int_distribution<int>;
  std::mt19937 m_gen;
  dist_type m_dist;

  UnifDist() : m_dist(1, 3) { m_gen.seed(1034343); }

  int operator()() { return m_dist(m_gen); }
};

template <class ViewType>
void fill_zero(ViewType view) {
  Kokkos::parallel_for(view.extent(0), FillZeroFunctor<ViewType>(view));
}

template <class ViewType>
void fill_view(ViewType dest_view, const std::string& name) {
  using value_type = typename ViewType::value_type;
  using exe_space  = typename ViewType::execution_space;

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
      v_h(i) = static_cast<value_type>(i + 1);
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

// I had to write my own because std::exclusive_scan is ONLY found with
// std=c++17
template <class it1, class it2, class ValType, class BopType>
void my_host_exclusive_scan(it1 first, it1 last, it2 dest, ValType init,
                            BopType bop) {
  const auto num_elements = last - first;
  if (num_elements > 0) {
    while (first < last - 1) {
      *(dest++) = init;
      init      = bop(*first++, init);
    }
    *dest = init;
  }
}

template <class ViewType1, class ViewType2, class ValueType, class BinaryOp>
void verify_data(ViewType1 data_view,  // contains data
                 ViewType2 test_view,  // the view to test
                 ValueType init_value, BinaryOp bop) {
  //! always careful because views might not be deep copyable

  auto data_view_dc = create_deep_copyable_compatible_clone(data_view);
  auto data_view_h =
      create_mirror_view_and_copy(Kokkos::HostSpace(), data_view_dc);

  using gold_view_value_type = typename ViewType2::value_type;
  Kokkos::View<gold_view_value_type*, Kokkos::HostSpace> gold_h(
      "goldh", data_view.extent(0));
  my_host_exclusive_scan(KE::cbegin(data_view_h), KE::cend(data_view_h),
                         KE::begin(gold_h), init_value, bop);

  auto test_view_dc = create_deep_copyable_compatible_clone(test_view);
  auto test_view_h =
      create_mirror_view_and_copy(Kokkos::HostSpace(), test_view_dc);
  if (test_view_h.extent(0) > 0) {
    for (std::size_t i = 0; i < test_view_h.extent(0); ++i) {
      // std::cout << i << " " << std::setprecision(15) << data_view_h(i) << " "
      //           << gold_h(i) << " " << test_view_h(i) << " "
      //           << std::abs(gold_h(i) - test_view_h(i)) << std::endl;
      if (std::is_same<gold_view_value_type, int>::value) {
        EXPECT_TRUE(gold_h(i) == test_view_h(i));
      } else {
        const auto error = std::abs(gold_h(i) - test_view_h(i));
        if (error > 1e-10) {
          std::cout << i << " " << std::setprecision(15) << data_view_h(i)
                    << " " << gold_h(i) << " " << test_view_h(i) << " "
                    << std::abs(gold_h(i) - test_view_h(i)) << std::endl;
        }
        EXPECT_TRUE(error < 1e-10);
      }
    }
  }
}

template <class ValueType>
struct MultiplyFunctor {
  KOKKOS_INLINE_FUNCTION
  ValueType operator()(const ValueType& a, const ValueType& b) const {
    return (a * b);
  }

  KOKKOS_INLINE_FUNCTION
  ValueType operator()(const volatile ValueType& a,
                       const volatile ValueType& b) const {
    return (a * b);
  }
};

template <class ValueType>
struct SumFunctor {
  KOKKOS_INLINE_FUNCTION
  ValueType operator()(const ValueType& a, const ValueType& b) const {
    return (a + b);
  }

  KOKKOS_INLINE_FUNCTION
  ValueType operator()(const volatile ValueType& a,
                       const volatile ValueType& b) const {
    return (a + b);
  }
};

std::string value_type_to_string(int) { return "int"; }

std::string value_type_to_string(double) { return "double"; }

template <class Tag, class ValueType, class InfoType>
void run_single_scenario_default_op(const InfoType& scenario_info,
                                    ValueType init_value) {
  using default_op           = SumFunctor<ValueType>;
  const auto name            = std::get<0>(scenario_info);
  const std::size_t view_ext = std::get<1>(scenario_info);
  // std::cout << "exclusive_scan default op: " << name << ", "
  //           << view_tag_to_string(Tag{}) << ", "
  //           << value_type_to_string(ValueType()) << ", "
  //           << "init = " << init_value << std::endl;

  auto view_dest = create_view<ValueType>(Tag{}, view_ext, "exclusive_scan");
  auto view_from = create_view<ValueType>(Tag{}, view_ext, "exclusive_scan");
  fill_view(view_from, name);

  {
    fill_zero(view_dest);
    auto r = KE::exclusive_scan(exespace(), KE::cbegin(view_from),
                                KE::cend(view_from), KE::begin(view_dest),
                                init_value);
    EXPECT_TRUE(r == KE::end(view_dest));
    verify_data(view_from, view_dest, init_value, default_op());
  }

  {
    fill_zero(view_dest);
    auto r = KE::exclusive_scan("label", exespace(), KE::cbegin(view_from),
                                KE::cend(view_from), KE::begin(view_dest),
                                init_value);
    EXPECT_TRUE(r == KE::end(view_dest));
    verify_data(view_from, view_dest, init_value, default_op());
  }

  {
    fill_zero(view_dest);
    auto r = KE::exclusive_scan(exespace(), view_from, view_dest, init_value);
    EXPECT_TRUE(r == KE::end(view_dest));
    verify_data(view_from, view_dest, init_value, default_op());
  }

  {
    fill_zero(view_dest);
    auto r = KE::exclusive_scan("label", exespace(), view_from, view_dest,
                                init_value);
    EXPECT_TRUE(r == KE::end(view_dest));
    verify_data(view_from, view_dest, init_value, default_op());
  }

  Kokkos::fence();
}

template <class Tag, class ValueType, class InfoType, class BinaryOp>
void run_single_scenario_custom_op(const InfoType& scenario_info,
                                   ValueType init_value, BinaryOp bop) {
  const auto name            = std::get<0>(scenario_info);
  const std::size_t view_ext = std::get<1>(scenario_info);
  // std::cout << "exclusive_scan custom op: " << name << ", "
  //           << view_tag_to_string(Tag{}) << ", "
  //           << value_type_to_string(ValueType()) << ", "
  //           << "init = " << init_value << std::endl;

  auto view_dest = create_view<ValueType>(Tag{}, view_ext, "exclusive_scan");
  auto view_from = create_view<ValueType>(Tag{}, view_ext, "exclusive_scan");
  fill_view(view_from, name);

  {
    fill_zero(view_dest);
    auto r = KE::exclusive_scan(exespace(), KE::cbegin(view_from),
                                KE::cend(view_from), KE::begin(view_dest),
                                init_value, bop);
    EXPECT_TRUE(r == KE::end(view_dest));
    verify_data(view_from, view_dest, init_value, bop);
  }

  {
    fill_zero(view_dest);
    auto r = KE::exclusive_scan("label", exespace(), KE::cbegin(view_from),
                                KE::cend(view_from), KE::begin(view_dest),
                                init_value, bop);
    EXPECT_TRUE(r == KE::end(view_dest));
    verify_data(view_from, view_dest, init_value, bop);
  }

  {
    fill_zero(view_dest);
    auto r =
        KE::exclusive_scan(exespace(), view_from, view_dest, init_value, bop);
    EXPECT_TRUE(r == KE::end(view_dest));
    verify_data(view_from, view_dest, init_value, bop);
  }

  {
    fill_zero(view_dest);
    auto r = KE::exclusive_scan("label", exespace(), view_from, view_dest,
                                init_value, bop);
    EXPECT_TRUE(r == KE::end(view_dest));
    verify_data(view_from, view_dest, init_value, bop);
  }

  Kokkos::fence();
}

template <class Tag, class ValueType>
void run_exclusive_scan_all_scenarios() {
  const std::map<std::string, std::size_t> scenarios = {
      {"empty", 0},          {"one-element", 1}, {"two-elements-a", 2},
      {"two-elements-b", 2}, {"small-a", 9},     {"small-b", 13},
      {"medium", 1103},      {"large", 10513}};

  for (const auto& it : scenarios) {
    run_single_scenario_default_op<Tag, ValueType>(it, ValueType{0});
    run_single_scenario_default_op<Tag, ValueType>(it, ValueType{1});
    run_single_scenario_default_op<Tag, ValueType>(it, ValueType{-2});
    run_single_scenario_default_op<Tag, ValueType>(it, ValueType{3});

#if not defined KOKKOS_ENABLE_OPENMPTARGET
    // custom multiply op is only run for small views otherwise it overflows
    if (it.first == "small-a" || it.first == "small-b") {
      using custom_bop_t = MultiplyFunctor<ValueType>;
      run_single_scenario_custom_op<Tag, ValueType>(it, ValueType{0},
                                                    custom_bop_t());
      run_single_scenario_custom_op<Tag, ValueType>(it, ValueType{1},
                                                    custom_bop_t());
      run_single_scenario_custom_op<Tag, ValueType>(it, ValueType{-2},
                                                    custom_bop_t());
      run_single_scenario_custom_op<Tag, ValueType>(it, ValueType{3},
                                                    custom_bop_t());
    }

    using custom_bop_t = SumFunctor<ValueType>;
    run_single_scenario_custom_op<Tag, ValueType>(it, ValueType{0},
                                                  custom_bop_t());
    run_single_scenario_custom_op<Tag, ValueType>(it, ValueType{1},
                                                  custom_bop_t());
    run_single_scenario_custom_op<Tag, ValueType>(it, ValueType{-2},
                                                  custom_bop_t());
    run_single_scenario_custom_op<Tag, ValueType>(it, ValueType{3},
                                                  custom_bop_t());
#endif
  }
}

TEST(std_algorithms_numeric_ops_test, exclusive_scan) {
  run_exclusive_scan_all_scenarios<DynamicTag, double>();
  run_exclusive_scan_all_scenarios<StridedThreeTag, double>();
  run_exclusive_scan_all_scenarios<DynamicTag, int>();
  run_exclusive_scan_all_scenarios<StridedThreeTag, int>();
}

}  // namespace EScan
}  // namespace stdalgos
}  // namespace Test
