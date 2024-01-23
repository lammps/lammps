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
namespace IncScan {

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

template <>
struct UnifDist<CustomValueType> {
  using dist_type = std::uniform_real_distribution<double>;
  std::mt19937 m_gen;
  dist_type m_dist;

  UnifDist() : m_dist(0.05, 1.2) { m_gen.seed(1034343); }

  CustomValueType operator()() { return m_dist(m_gen); }
};

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

  else if (name == "medium-a" || name == "medium-b" || name == "large") {
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

// my own because std::inclusive_scan is ONLY found with std=c++17
template <class it1, class it2, class BinOp>
void my_host_inclusive_scan(it1 first, it1 last, it2 dest, BinOp bop) {
  if (first != last) {
    auto init = *first;
    *dest     = init;
    while (++first < last) {
      init      = bop(*first, init);
      *(++dest) = init;
    }
  }
}

template <class it1, class it2, class BinOp, class ValType>
void my_host_inclusive_scan(it1 first, it1 last, it2 dest, BinOp bop,
                            ValType init) {
  if (first != last) {
    init  = bop(*first, init);
    *dest = init;
    while (++first < last) {
      init      = bop(*first, init);
      *(++dest) = init;
    }
  }
}

template <class ViewType1, class ViewType2, class BinaryOp, class... Args>
void verify_data(ViewType1 data_view,  // contains data
                 ViewType2 test_view,  // the view to test
                 BinaryOp bop, Args... args /* copy on purpose */) {
  //! always careful because views might not be deep copyable

  auto data_view_dc = create_deep_copyable_compatible_clone(data_view);
  auto data_view_h =
      create_mirror_view_and_copy(Kokkos::HostSpace(), data_view_dc);

  using gold_view_value_type = typename ViewType2::value_type;
  Kokkos::View<gold_view_value_type*, Kokkos::HostSpace> gold_h(
      "goldh", data_view.extent(0));
  my_host_inclusive_scan(KE::cbegin(data_view_h), KE::cend(data_view_h),
                         KE::begin(gold_h), bop, args...);

  auto test_view_dc = create_deep_copyable_compatible_clone(test_view);
  auto test_view_h =
      create_mirror_view_and_copy(Kokkos::HostSpace(), test_view_dc);

  const auto ext = test_view_h.extent(0);
  if (ext > 0) {
    for (std::size_t i = 0; i < ext; ++i) {
      // std::cout << i << " " << std::setprecision(15) << data_view_h(i) << " "
      //           << gold_h(i) << " " << test_view_h(i) << " "
      //           << std::abs(gold_h(i) - test_view_h(i)) << std::endl;

      if (std::is_same<gold_view_value_type, int>::value) {
        ASSERT_EQ(gold_h(i), test_view_h(i));
      } else {
        const auto error =
            std::abs(static_cast<double>(gold_h(i) - test_view_h(i)));
        if (error > 1e-10) {
          std::cout << i << " " << std::setprecision(15) << data_view_h(i)
                    << " " << gold_h(i) << " " << test_view_h(i) << " "
                    << std::abs(static_cast<double>(gold_h(i) - test_view_h(i)))
                    << std::endl;
        }
        EXPECT_LT(error, 1e-10);
      }
    }
    // std::cout << " last el: " << test_view_h(ext-1) << std::endl;
  }
}

template <class ValueType>
struct MultiplyFunctor {
  KOKKOS_INLINE_FUNCTION
  ValueType operator()(const ValueType& a, const ValueType& b) const {
    return (a * b);
  }
};

template <class ValueType>
struct SumFunctor {
  KOKKOS_INLINE_FUNCTION
  ValueType operator()(const ValueType& a, const ValueType& b) const {
    return (a + b);
  }
};

std::string value_type_to_string(int) { return "int"; }
std::string value_type_to_string(double) { return "double"; }

template <class Tag, class ValueType, class InfoType>
void run_single_scenario_default_op(const InfoType& scenario_info) {
  using default_op           = SumFunctor<ValueType>;
  const auto name            = std::get<0>(scenario_info);
  const std::size_t view_ext = std::get<1>(scenario_info);
  // std::cout << "inclusive_scan default op: " << name << ", "
  //           << view_tag_to_string(Tag{}) << ", "
  //           << value_type_to_string(ValueType()) << std::endl;

  auto view_dest = create_view<ValueType>(Tag{}, view_ext, "inclusive_scan");
  auto view_from = create_view<ValueType>(Tag{}, view_ext, "inclusive_scan");
  fill_view(view_from, name);

  {
    fill_zero(view_dest);
    auto r = KE::inclusive_scan(exespace(), KE::cbegin(view_from),
                                KE::cend(view_from), KE::begin(view_dest));
    ASSERT_EQ(r, KE::end(view_dest));
    verify_data(view_from, view_dest, default_op());
  }

  {
    fill_zero(view_dest);
    auto r = KE::inclusive_scan("label", exespace(), KE::cbegin(view_from),
                                KE::cend(view_from), KE::begin(view_dest));
    ASSERT_EQ(r, KE::end(view_dest));
    verify_data(view_from, view_dest, default_op());
  }

  {
    fill_zero(view_dest);
    auto r = KE::inclusive_scan(exespace(), view_from, view_dest);
    ASSERT_EQ(r, KE::end(view_dest));
    verify_data(view_from, view_dest, default_op());
  }

  {
    fill_zero(view_dest);
    auto r = KE::inclusive_scan("label", exespace(), view_from, view_dest);
    ASSERT_EQ(r, KE::end(view_dest));
    verify_data(view_from, view_dest, default_op());
  }

  Kokkos::fence();
}

template <class Tag, class ValueType, class InfoType, class BinaryOp,
          class... Args>
void run_single_scenario_custom_op(const InfoType& scenario_info, BinaryOp bop,
                                   Args... args /* copy on purpose */) {
  const auto name            = std::get<0>(scenario_info);
  const std::size_t view_ext = std::get<1>(scenario_info);

  // if (1 == sizeof...(Args)) {
  //   std::cout << "inclusive_scan custom op and init value: " << name << ", "
  //             << view_tag_to_string(Tag{}) << ", "
  //             << value_type_to_string(ValueType()) << ", " << std::endl;
  // } else {
  //   std::cout << "inclusive_scan custom op: " << name << ", "
  //             << view_tag_to_string(Tag{}) << ", "
  //             << value_type_to_string(ValueType()) << ", " << std::endl;
  // }

  auto view_dest = create_view<ValueType>(Tag{}, view_ext, "inclusive_scan");
  auto view_from = create_view<ValueType>(Tag{}, view_ext, "inclusive_scan");
  fill_view(view_from, name);

  {
    fill_zero(view_dest);
    auto r = KE::inclusive_scan(exespace(), KE::cbegin(view_from),
                                KE::cend(view_from), KE::begin(view_dest), bop,
                                args...);
    ASSERT_EQ(r, KE::end(view_dest));
    verify_data(view_from, view_dest, bop, args...);
  }

  {
    fill_zero(view_dest);
    auto r = KE::inclusive_scan("label", exespace(), KE::cbegin(view_from),
                                KE::cend(view_from), KE::begin(view_dest), bop,
                                args...);
    ASSERT_EQ(r, KE::end(view_dest));
    verify_data(view_from, view_dest, bop, args...);
  }

  {
    fill_zero(view_dest);
    auto r = KE::inclusive_scan(exespace(), view_from, view_dest, bop, args...);
    ASSERT_EQ(r, KE::end(view_dest));
    verify_data(view_from, view_dest, bop, args...);
  }

  {
    fill_zero(view_dest);
    auto r = KE::inclusive_scan("label", exespace(), view_from, view_dest, bop,
                                args...);
    ASSERT_EQ(r, KE::end(view_dest));
    verify_data(view_from, view_dest, bop, args...);
  }

  Kokkos::fence();
}

template <class Tag, class ValueType>
void run_inclusive_scan_all_scenarios() {
  const std::map<std::string, std::size_t> scenarios = {
      {"empty", 0},          {"one-element", 1}, {"two-elements-a", 2},
      {"two-elements-b", 2}, {"small-a", 9},     {"small-b", 13},
      {"medium-a", 313},     {"medium-b", 1103}, {"large", 10513}};

  for (const auto& it : scenarios) {
    run_single_scenario_default_op<Tag, ValueType>(it);

#if !defined KOKKOS_ENABLE_OPENMPTARGET
    // the sum custom op is always run
    using sum_binary_op = SumFunctor<ValueType>;
    sum_binary_op sbop;
    run_single_scenario_custom_op<Tag, ValueType>(it, sbop);
    run_single_scenario_custom_op<Tag, ValueType>(it, sbop, ValueType{0});
    run_single_scenario_custom_op<Tag, ValueType>(it, sbop, ValueType{1});
    run_single_scenario_custom_op<Tag, ValueType>(it, sbop, ValueType{-2});
    run_single_scenario_custom_op<Tag, ValueType>(it, sbop, ValueType{3});

    // custom multiply only for small views to avoid overflows
    if (it.first == "small-a" || it.first == "small-b") {
      using mult_binary_op = MultiplyFunctor<ValueType>;
      mult_binary_op mbop;
      run_single_scenario_custom_op<Tag, ValueType>(it, mbop);
      run_single_scenario_custom_op<Tag, ValueType>(it, mbop, ValueType{0});
      run_single_scenario_custom_op<Tag, ValueType>(it, mbop, ValueType{1});
      run_single_scenario_custom_op<Tag, ValueType>(it, mbop, ValueType{-2});
      run_single_scenario_custom_op<Tag, ValueType>(it, mbop, ValueType{3});
    }
#endif
  }
}

TEST(std_algorithms_numeric_ops_test, inclusive_scan) {
  run_inclusive_scan_all_scenarios<DynamicTag, double>();
  run_inclusive_scan_all_scenarios<StridedThreeTag, double>();
  run_inclusive_scan_all_scenarios<DynamicTag, int>();
  run_inclusive_scan_all_scenarios<StridedThreeTag, int>();
  run_inclusive_scan_all_scenarios<DynamicTag, CustomValueType>();
  run_inclusive_scan_all_scenarios<StridedThreeTag, CustomValueType>();
}

TEST(std_algorithms_numeric_ops_test, inclusive_scan_functor) {
  using view_type = Kokkos::View<int*, exespace>;
  view_type dummy_view("dummy_view", 0);
  using functor_type = Kokkos::Experimental::Impl::InclusiveScanDefaultFunctor<
      exespace, int, int, view_type, view_type>;
  functor_type functor(dummy_view, dummy_view);
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

  value1.val        = 1;
  value1.is_initial = false;
  value2.val        = 2;
  value2.is_initial = false;
  functor.join(value2, value1);
  ASSERT_EQ(value2.val, 3);
  ASSERT_EQ(value2.is_initial, false);
}

}  // namespace IncScan
}  // namespace stdalgos
}  // namespace Test
