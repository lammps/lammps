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
namespace TransformEScan {

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
    assert(v_h.extent(0) == 1);
    v_h(0) = static_cast<value_type>(1);
  }

  else if (name == "two-elements-a") {
    assert(v_h.extent(0) == 2);
    v_h(0) = static_cast<value_type>(1);
    v_h(1) = static_cast<value_type>(2);
  }

  else if (name == "two-elements-b") {
    assert(v_h.extent(0) == 2);
    v_h(0) = static_cast<value_type>(2);
    v_h(1) = static_cast<value_type>(-1);
  }

  else if (name == "small-a") {
    assert(v_h.extent(0) == 9);
    v_h(0) = static_cast<value_type>(3);
    v_h(1) = static_cast<value_type>(1);
    v_h(2) = static_cast<value_type>(4);
    v_h(3) = static_cast<value_type>(1);
    v_h(4) = static_cast<value_type>(5);
    v_h(5) = static_cast<value_type>(9);
    v_h(6) = static_cast<value_type>(2);
    v_h(7) = static_cast<value_type>(6);
    v_h(8) = static_cast<value_type>(2);
  }

  else if (name == "small-b") {
    assert(v_h.extent(0) >= 6);
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

// I had to write my own because std::transform_exclusive_scan is ONLY found
// with std=c++17
template <class it1, class it2, class ValType, class BopType, class UopType>
void my_host_transform_exclusive_scan(it1 first, it1 last, it2 dest,
                                      ValType init, BopType bop, UopType uop) {
  const auto num_elements = last - first;
  if (num_elements > 0) {
    while (first < last - 1) {
      *(dest++) = init;
      init      = bop(uop(*(first++)), init);
    }
    *dest = init;
  }
}

template <class ViewType1, class ViewType2, class ValueType, class BinaryOp,
          class UnaryOp>
void verify_data(ViewType1 data_view,  // contains data
                 ViewType2 test_view,  // the view to test
                 ValueType init_value, BinaryOp bop, UnaryOp uop) {
  //! always careful because views might not be deep copyable

  auto data_view_dc = create_deep_copyable_compatible_clone(data_view);
  auto data_view_h =
      create_mirror_view_and_copy(Kokkos::HostSpace(), data_view_dc);

  using gold_view_value_type = typename ViewType2::value_type;
  Kokkos::View<gold_view_value_type*, Kokkos::HostSpace> gold_h(
      "goldh", data_view.extent(0));
  my_host_transform_exclusive_scan(KE::cbegin(data_view_h),
                                   KE::cend(data_view_h), KE::begin(gold_h),
                                   init_value, bop, uop);

  auto test_view_dc = create_deep_copyable_compatible_clone(test_view);
  auto test_view_h =
      create_mirror_view_and_copy(Kokkos::HostSpace(), test_view_dc);
  if (test_view_h.extent(0) > 0) {
    for (std::size_t i = 0; i < test_view_h.extent(0); ++i) {
      // std::cout << i << " " << std::setprecision(15) << data_view_h(i) << " "
      //           << gold_h(i) << " " << test_view_h(i) << " "
      //           << std::abs(gold_h(i) - test_view_h(i)) << std::endl;

      if (std::is_same<gold_view_value_type, int>::value) {
        ASSERT_EQ(gold_h(i), test_view_h(i));
      } else {
        const auto error = std::abs(gold_h(i) - test_view_h(i));
        if (error > 1e-10) {
          std::cout << i << " " << std::setprecision(15) << data_view_h(i)
                    << " " << gold_h(i) << " " << test_view_h(i) << " "
                    << std::abs(gold_h(i) - test_view_h(i)) << std::endl;
        }
        EXPECT_LT(error, 1e-10);
      }
    }
    // std::cout << " last el: " << test_view_h(test_view_h.extent(0)-1) <<
    // std::endl;
  }
}

template <class ValueType>
struct TimesTwoUnaryFunctor {
  KOKKOS_INLINE_FUNCTION
  ValueType operator()(const ValueType& a) const { return (a * ValueType(2)); }
};

template <class ValueType>
struct SumBinaryFunctor {
  KOKKOS_INLINE_FUNCTION
  ValueType operator()(const ValueType& a, const ValueType& b) const {
    return (a + b);
  }
};

std::string value_type_to_string(int) { return "int"; }

std::string value_type_to_string(double) { return "double"; }

template <class Tag, class ValueType, class InfoType, class BinaryOp,
          class UnaryOp>
void run_single_scenario(const InfoType& scenario_info, ValueType init_value,
                         BinaryOp bop, UnaryOp uop) {
  const auto name            = std::get<0>(scenario_info);
  const std::size_t view_ext = std::get<1>(scenario_info);
  // std::cout << "transform_exclusive_scan custom op: " << name << ", "
  //           << view_tag_to_string(Tag{}) << ", "
  //           << value_type_to_string(ValueType()) << ", "
  //           << "init = " << init_value << std::endl;

  auto view_dest =
      create_view<ValueType>(Tag{}, view_ext, "transform_exclusive_scan");
  auto view_from =
      create_view<ValueType>(Tag{}, view_ext, "transform_exclusive_scan");
  fill_view(view_from, name);

  {
    fill_zero(view_dest);
    auto r = KE::transform_exclusive_scan(
        exespace(), KE::cbegin(view_from), KE::cend(view_from),
        KE::begin(view_dest), init_value, bop, uop);
    ASSERT_EQ(r, KE::end(view_dest));
    verify_data(view_from, view_dest, init_value, bop, uop);
  }

  {
    fill_zero(view_dest);
    auto r = KE::transform_exclusive_scan(
        "label", exespace(), KE::cbegin(view_from), KE::cend(view_from),
        KE::begin(view_dest), init_value, bop, uop);
    ASSERT_EQ(r, KE::end(view_dest));
    verify_data(view_from, view_dest, init_value, bop, uop);
  }

  {
    fill_zero(view_dest);
    auto r = KE::transform_exclusive_scan(exespace(), view_from, view_dest,
                                          init_value, bop, uop);
    ASSERT_EQ(r, KE::end(view_dest));
    verify_data(view_from, view_dest, init_value, bop, uop);
  }

  {
    fill_zero(view_dest);
    auto r = KE::transform_exclusive_scan("label", exespace(), view_from,
                                          view_dest, init_value, bop, uop);
    ASSERT_EQ(r, KE::end(view_dest));
    verify_data(view_from, view_dest, init_value, bop, uop);
  }

  Kokkos::fence();
}

template <class Tag, class ValueType>
void run_all_scenarios() {
  const std::map<std::string, std::size_t> scenarios = {
      {"empty", 0},          {"one-element", 1}, {"two-elements-a", 2},
      {"two-elements-b", 2}, {"small-a", 9},     {"small-b", 13},
      {"medium", 1103},      {"large", 10513}};

  for (const auto& it : scenarios) {
    using uop_t = TimesTwoUnaryFunctor<ValueType>;
    using bop_t = SumBinaryFunctor<ValueType>;
    run_single_scenario<Tag, ValueType>(it, ValueType{0}, bop_t(), uop_t());
    run_single_scenario<Tag, ValueType>(it, ValueType{1}, bop_t(), uop_t());
    run_single_scenario<Tag, ValueType>(it, ValueType{-2}, bop_t(), uop_t());
    run_single_scenario<Tag, ValueType>(it, ValueType{3}, bop_t(), uop_t());
  }
}

#if !defined KOKKOS_ENABLE_OPENMPTARGET
TEST(std_algorithms_numeric_ops_test, transform_exclusive_scan) {
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedThreeTag, double>();
  run_all_scenarios<DynamicTag, int>();
  run_all_scenarios<StridedThreeTag, int>();
}
#endif

template <class ValueType>
struct MultiplyFunctor {
  KOKKOS_INLINE_FUNCTION
  ValueType operator()(const ValueType& a, const ValueType& b) const {
    return (a * b);
  }
};

TEST(std_algorithms_numeric_ops_test, transform_exclusive_scan_functor) {
  int dummy       = 0;
  using view_type = Kokkos::View<int*, exespace>;
  view_type dummy_view("dummy_view", 0);
  using unary_op_type =
      Kokkos::Experimental::Impl::StdNumericScanIdentityReferenceUnaryFunctor<
          int>;
  using functor_type =
      Kokkos::Experimental::Impl::TransformExclusiveScanFunctor<
          exespace, int, int, view_type, view_type, MultiplyFunctor<int>,
          unary_op_type>;
  functor_type functor(dummy, dummy_view, dummy_view, {}, {});
  using value_type = functor_type::value_type;

  value_type value1;
  functor.init(value1);
  ASSERT_EQ(value1.val, 0);
  ASSERT_EQ(value1.is_initial, true);

  value_type value2;
  value2.val        = 1;
  value2.is_initial = false;
  functor.join(value1, value2);
  ASSERT_EQ(value1.val, 1);
  ASSERT_EQ(value1.is_initial, false);

  functor.init(value1);
  functor.join(value2, value1);
  ASSERT_EQ(value2.val, 1);
  ASSERT_EQ(value2.is_initial, false);

  functor.init(value2);
  functor.join(value2, value1);
  ASSERT_EQ(value2.val, 0);
  ASSERT_EQ(value2.is_initial, true);

  value1.val        = 3;
  value1.is_initial = false;
  value2.val        = 2;
  value2.is_initial = false;
  functor.join(value2, value1);
  ASSERT_EQ(value2.val, 6);
  ASSERT_EQ(value2.is_initial, false);
}

}  // namespace TransformEScan
}  // namespace stdalgos
}  // namespace Test
