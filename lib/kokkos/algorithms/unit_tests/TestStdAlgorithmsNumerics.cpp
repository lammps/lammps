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

namespace KE = Kokkos::Experimental;

namespace Test {
namespace stdalgos {

template <class ValueType>
struct TimesTwoUnaryTransformFunctor {
  KOKKOS_INLINE_FUNCTION
  ValueType operator()(const ValueType& a) const { return (a * 2.); }
};

template <class ValueType>
struct MultiplyAndHalveBinaryTransformFunctor {
  KOKKOS_INLINE_FUNCTION
  ValueType operator()(const ValueType& a, const ValueType& b) const {
    return (a * b) * 0.5;
  }
};

template <class ValueType>
struct SumJoinFunctor {
  KOKKOS_INLINE_FUNCTION
  ValueType operator()(const ValueType& a, const ValueType& b) const {
    return a + b;
  }
};

struct std_algorithms_numerics_test : public ::testing::Test {
  Kokkos::LayoutStride layout{20, 2};

  // value_type
  using static_view_t  = Kokkos::View<value_type[20]>;
  using dyn_view_t     = Kokkos::View<value_type*>;
  using strided_view_t = Kokkos::View<value_type*, Kokkos::LayoutStride>;

  static_view_t m_static_view{"std-algo-test-1D-contiguous-view-static"};
  dyn_view_t m_dynamic_view{"std-algo-test-1D-contiguous-view-dyn", 20};
  strided_view_t m_strided_view{"std-algo-test-1D-strided-view", layout};

  // custom scalar (cs)
  using static_view_cs_t = Kokkos::View<CustomValueType[20]>;
  using dyn_view_cs_t    = Kokkos::View<CustomValueType*>;
  using strided_view_cs_t =
      Kokkos::View<CustomValueType*, Kokkos::LayoutStride>;

  static_view_cs_t m_static_view_cs{
      "std-algo-test-1D-contiguous-view-static-custom-scalar"};
  dyn_view_cs_t m_dynamic_view_cs{
      "std-algo-test-1D-contiguous-view-dyn-custom_scalar", 20};
  strided_view_cs_t m_strided_view_cs{
      "std-algo-test-1D-strided-view-custom-scalar", layout};

  template <class ViewFromType, class ViewToType>
  void copyPodViewToCustom(ViewFromType v_from, ViewToType v_to) {
    for (std::size_t i = 0; i < v_from.extent(0); ++i) {
      v_to(i)() = v_from(i);
    }
  }

  void fillFixtureViews() {
    static_view_t tmpView("tmpView");
    static_view_cs_t tmpViewCs("tmpViewCs");
    auto tmp_view_h = Kokkos::create_mirror_view(Kokkos::HostSpace(), tmpView);
    auto tmp_view_cs_h =
        Kokkos::create_mirror_view(Kokkos::HostSpace(), tmpViewCs);
    tmp_view_h(0)  = 0.;
    tmp_view_h(1)  = 0.;
    tmp_view_h(2)  = 0.;
    tmp_view_h(3)  = 2.;
    tmp_view_h(4)  = 2.;
    tmp_view_h(5)  = 1.;
    tmp_view_h(6)  = 1.;
    tmp_view_h(7)  = 1.;
    tmp_view_h(8)  = 1.;
    tmp_view_h(9)  = 0.;
    tmp_view_h(10) = -2.;
    tmp_view_h(11) = -2.;
    tmp_view_h(12) = 0.;
    tmp_view_h(13) = 2.;
    tmp_view_h(14) = 2.;
    tmp_view_h(15) = 1.;
    tmp_view_h(16) = 1.;
    tmp_view_h(17) = 1.;
    tmp_view_h(18) = 1.;
    tmp_view_h(19) = 0.;

    copyPodViewToCustom(tmp_view_h, tmp_view_cs_h);

    Kokkos::deep_copy(tmpView, tmp_view_h);
    Kokkos::deep_copy(tmpViewCs, tmp_view_cs_h);

    CopyFunctor<static_view_t, static_view_t> F1(tmpView, m_static_view);
    Kokkos::parallel_for("_std_algo_copy1", 20, F1);

    CopyFunctor<static_view_t, dyn_view_t> F2(tmpView, m_dynamic_view);
    Kokkos::parallel_for("_std_algo_copy2", 20, F2);

    CopyFunctor<static_view_t, strided_view_t> F3(tmpView, m_strided_view);
    Kokkos::parallel_for("_std_algo_copy3", 20, F3);

    CopyFunctor<static_view_cs_t, static_view_cs_t> F4(tmpViewCs,
                                                       m_static_view_cs);
    Kokkos::parallel_for("_std_algo_copy4", 20, F4);

    CopyFunctor<static_view_cs_t, dyn_view_cs_t> F5(tmpViewCs,
                                                    m_dynamic_view_cs);
    Kokkos::parallel_for("_std_algo_copy5", 20, F5);

    CopyFunctor<static_view_cs_t, strided_view_cs_t> F6(tmpViewCs,
                                                        m_strided_view_cs);
    Kokkos::parallel_for("_std_algo_copy6", 20, F6);
  }
};

#if !defined KOKKOS_ENABLE_OPENMPTARGET

// -------------------------------------------------------------------
// test default case of transform_reduce
//
// test for both POD types and custom scalar types
// -------------------------------------------------------------------
template <class ExecutionSpace, class ViewType1, class ViewType2,
          class ValueType>
void run_and_check_transform_reduce_default(ViewType1 first_view,
                                            ViewType2 second_view,
                                            ValueType init_value,
                                            ValueType result_value) {
  // trivial cases
  const auto r1 = KE::transform_reduce(ExecutionSpace(), KE::cbegin(first_view),
                                       KE::cbegin(first_view),
                                       KE::cbegin(second_view), init_value);

  const auto r2 = KE::transform_reduce(
      "MYLABEL", ExecutionSpace(), KE::cbegin(first_view),
      KE::cbegin(first_view), KE::cbegin(second_view), init_value);
  ASSERT_EQ(r1, init_value);
  ASSERT_EQ(r2, init_value);

  // non-trivial cases
  const auto r3 = KE::transform_reduce(ExecutionSpace(), KE::cbegin(first_view),
                                       KE::cend(first_view),
                                       KE::cbegin(second_view), init_value);

  const auto r4 = KE::transform_reduce(
      "MYLABEL", ExecutionSpace(), KE::cbegin(first_view), KE::cend(first_view),
      KE::cbegin(second_view), init_value);

  const auto r5 = KE::transform_reduce(ExecutionSpace(), first_view,
                                       second_view, init_value);
  const auto r6 = KE::transform_reduce("MYLABEL", ExecutionSpace(), first_view,
                                       second_view, init_value);

  ASSERT_EQ(r3, result_value);
  ASSERT_EQ(r4, result_value);
  ASSERT_EQ(r5, result_value);
  ASSERT_EQ(r6, result_value);
}

TEST_F(std_algorithms_numerics_test,
       transform_reduce_default_functors_using_pod_value_type) {
  fillFixtureViews();
  const value_type init0 = 0.;
  const value_type init5 = 5.;
  const value_type gold0 = 32.;
  const value_type gold5 = 37.;

  run_and_check_transform_reduce_default<exespace>(
      m_static_view, m_dynamic_view, init0, gold0);
  run_and_check_transform_reduce_default<exespace>(
      m_static_view, m_dynamic_view, init5, gold5);

  run_and_check_transform_reduce_default<exespace>(
      m_static_view, m_strided_view, init0, gold0);
  run_and_check_transform_reduce_default<exespace>(
      m_static_view, m_strided_view, init5, gold5);

  run_and_check_transform_reduce_default<exespace>(
      m_dynamic_view, m_strided_view, init0, gold0);
  run_and_check_transform_reduce_default<exespace>(
      m_dynamic_view, m_strided_view, init5, gold5);
}

TEST_F(std_algorithms_numerics_test,
       transform_reduce_default_functors_using_custom_value_type) {
  fillFixtureViews();
  const CustomValueType init0{0.};
  const CustomValueType init5{5.};
  const CustomValueType gold0{32.};
  const CustomValueType gold5{37.};

  run_and_check_transform_reduce_default<exespace>(
      m_static_view_cs, m_dynamic_view_cs, init0, gold0);
  run_and_check_transform_reduce_default<exespace>(
      m_static_view_cs, m_dynamic_view_cs, init5, gold5);

  run_and_check_transform_reduce_default<exespace>(
      m_static_view_cs, m_strided_view_cs, init0, gold0);
  run_and_check_transform_reduce_default<exespace>(
      m_static_view_cs, m_strided_view_cs, init5, gold5);

  run_and_check_transform_reduce_default<exespace>(
      m_dynamic_view_cs, m_strided_view_cs, init0, gold0);
  run_and_check_transform_reduce_default<exespace>(
      m_dynamic_view_cs, m_strided_view_cs, init5, gold5);
}

// -------------------------------------------------------------------
// transform_reduce for custom joiner and custom transform op
// test for both POD types and custom scalar types
//
// test overload1 accepting two intervals
//
// Note that in the std, the reducer is called BinaryReductionOp
// but in the Kokkos naming convention, it corresponds to a "joiner"
// that knows how to join two values.
// the "joiner" is assumed to be commutative:
//
// https://en.cppreference.com/w/cpp/algorithm/transform_reduce
//
// -------------------------------------------------------------------

template <class ExecutionSpace, class ViewType1, class ViewType2,
          class ValueType, class... Args>
void run_and_check_transform_reduce_overloadA(ViewType1 first_view,
                                              ViewType2 second_view,
                                              ValueType init_value,
                                              ValueType result_value,
                                              Args&&... args) {
  // trivial cases
  const auto r1 = KE::transform_reduce(
      ExecutionSpace(), KE::cbegin(first_view), KE::cbegin(first_view),
      KE::cbegin(second_view), init_value, std::forward<Args>(args)...);

  const auto r2 =
      KE::transform_reduce("MYLABEL", ExecutionSpace(), KE::cbegin(first_view),
                           KE::cbegin(first_view), KE::cbegin(second_view),
                           init_value, std::forward<Args>(args)...);

  ASSERT_EQ(r1, init_value);
  ASSERT_EQ(r2, init_value);

  // non trivial cases
  const auto r3 = KE::transform_reduce(
      ExecutionSpace(), KE::cbegin(first_view), KE::cend(first_view),
      KE::cbegin(second_view), init_value, std::forward<Args>(args)...);

  const auto r4 = KE::transform_reduce(
      "MYLABEL", ExecutionSpace(), KE::cbegin(first_view), KE::cend(first_view),
      KE::cbegin(second_view), init_value, std::forward<Args>(args)...);

  const auto r5 =
      KE::transform_reduce(ExecutionSpace(), first_view, second_view,
                           init_value, std::forward<Args>(args)...);
  const auto r6 =
      KE::transform_reduce("MYLABEL", ExecutionSpace(), first_view, second_view,
                           init_value, std::forward<Args>(args)...);

  ASSERT_EQ(r3, result_value);
  ASSERT_EQ(r4, result_value);
  ASSERT_EQ(r5, result_value);
  ASSERT_EQ(r6, result_value);
}

TEST_F(std_algorithms_numerics_test,
       transform_reduce_custom_functors_overloadA_using_pod_value_type) {
  using joiner_type = SumJoinFunctor<value_type>;
  using transf_type = MultiplyAndHalveBinaryTransformFunctor<value_type>;

  const value_type init0 = 0.;
  const value_type init5 = 5.;
  const value_type gold0 = 16.;
  const value_type gold5 = 21.;

  fillFixtureViews();
  run_and_check_transform_reduce_overloadA<exespace>(
      m_static_view, m_dynamic_view, init0, gold0, joiner_type(),
      transf_type());
  run_and_check_transform_reduce_overloadA<exespace>(
      m_static_view, m_dynamic_view, init5, gold5, joiner_type(),
      transf_type());

  run_and_check_transform_reduce_overloadA<exespace>(
      m_static_view, m_strided_view, init0, gold0, joiner_type(),
      transf_type());
  run_and_check_transform_reduce_overloadA<exespace>(
      m_static_view, m_strided_view, init5, gold5, joiner_type(),
      transf_type());
  run_and_check_transform_reduce_overloadA<exespace>(
      m_dynamic_view, m_strided_view, init0, gold0, joiner_type(),
      transf_type());
  run_and_check_transform_reduce_overloadA<exespace>(
      m_dynamic_view, m_strided_view, init5, gold5, joiner_type(),
      transf_type());
}

TEST_F(std_algorithms_numerics_test,
       transform_reduce_custom_functors_overloadA_using_custom_value_type) {
  using joiner_type = SumJoinFunctor<CustomValueType>;
  using transf_type = MultiplyAndHalveBinaryTransformFunctor<CustomValueType>;

  const CustomValueType init0{0.};
  const CustomValueType init5{5.};
  const CustomValueType gold0{16.};
  const CustomValueType gold5{21.};

  fillFixtureViews();
  run_and_check_transform_reduce_overloadA<exespace>(
      m_static_view_cs, m_dynamic_view_cs, init0, gold0, joiner_type(),
      transf_type());
  run_and_check_transform_reduce_overloadA<exespace>(
      m_static_view_cs, m_dynamic_view_cs, init5, gold5, joiner_type(),
      transf_type());

  run_and_check_transform_reduce_overloadA<exespace>(
      m_static_view_cs, m_strided_view_cs, init0, gold0, joiner_type(),
      transf_type());
  run_and_check_transform_reduce_overloadA<exespace>(
      m_static_view_cs, m_strided_view_cs, init5, gold5, joiner_type(),
      transf_type());

  run_and_check_transform_reduce_overloadA<exespace>(
      m_dynamic_view_cs, m_strided_view_cs, init0, gold0, joiner_type(),
      transf_type());
  run_and_check_transform_reduce_overloadA<exespace>(
      m_dynamic_view_cs, m_strided_view_cs, init5, gold5, joiner_type(),
      transf_type());
}

// -------------------------------------------------------------------
// transform_reduce for custom joiner and custom transform op
// test for both POD types and custom scalar types
//
// test overload1 accepting single interval/view
//
// Note that in the std, the reducer is called BinaryReductionOp
// but in the Kokkos naming convention, it corresponds to a "joiner"
// that knows how to join two values.
// the "joiner" is assumed to be commutative:
//
// https://en.cppreference.com/w/cpp/algorithm/transform_reduce
//
// -------------------------------------------------------------------

template <class ExecutionSpace, class ViewType, class ValueType, class... Args>
void run_and_check_transform_reduce_overloadB(ViewType view,
                                              ValueType init_value,
                                              ValueType result_value,
                                              Args&&... args) {
  // trivial
  const auto r1 =
      KE::transform_reduce(ExecutionSpace(), KE::cbegin(view), KE::cbegin(view),
                           init_value, std::forward<Args>(args)...);

  const auto r2 = KE::transform_reduce("MYLABEL", ExecutionSpace(),
                                       KE::cbegin(view), KE::cbegin(view),
                                       init_value, std::forward<Args>(args)...);

  ASSERT_EQ(r1, init_value);
  ASSERT_EQ(r2, init_value);

  // non trivial
  const auto r3 =
      KE::transform_reduce(ExecutionSpace(), KE::cbegin(view), KE::cend(view),
                           init_value, std::forward<Args>(args)...);

  const auto r4 = KE::transform_reduce("MYLABEL", ExecutionSpace(),
                                       KE::cbegin(view), KE::cend(view),
                                       init_value, std::forward<Args>(args)...);
  const auto r5 = KE::transform_reduce(ExecutionSpace(), view, init_value,
                                       std::forward<Args>(args)...);

  const auto r6 = KE::transform_reduce("MYLABEL", ExecutionSpace(), view,
                                       init_value, std::forward<Args>(args)...);

  ASSERT_EQ(r3, result_value);
  ASSERT_EQ(r4, result_value);
  ASSERT_EQ(r5, result_value);
  ASSERT_EQ(r6, result_value);
}

TEST_F(std_algorithms_numerics_test,
       transform_reduce_custom_functors_overloadB_using_pod_value_type) {
  using joiner_type = SumJoinFunctor<value_type>;
  using transf_type = TimesTwoUnaryTransformFunctor<value_type>;

  const value_type init0 = 0.;
  const value_type init5 = 5.;
  const value_type gold0 = 24.;
  const value_type gold5 = 29.;

  fillFixtureViews();
  run_and_check_transform_reduce_overloadB<exespace>(
      m_static_view, init0, gold0, joiner_type(), transf_type());
  run_and_check_transform_reduce_overloadB<exespace>(
      m_dynamic_view, init5, gold5, joiner_type(), transf_type());
  run_and_check_transform_reduce_overloadB<exespace>(
      m_strided_view, init0, gold0, joiner_type(), transf_type());
}

TEST_F(std_algorithms_numerics_test,
       transform_reduce_custom_functors_overloadB_using_custom_value_type) {
  using joiner_type = SumJoinFunctor<CustomValueType>;
  using transf_type = TimesTwoUnaryTransformFunctor<CustomValueType>;

  const CustomValueType init0{0.};
  const CustomValueType init5{5.};
  const CustomValueType gold0{24.};
  const CustomValueType gold5{29.};

  fillFixtureViews();
  run_and_check_transform_reduce_overloadB<exespace>(
      m_static_view_cs, init0, gold0, joiner_type(), transf_type());
  run_and_check_transform_reduce_overloadB<exespace>(
      m_dynamic_view_cs, init5, gold5, joiner_type(), transf_type());
  run_and_check_transform_reduce_overloadB<exespace>(
      m_strided_view_cs, init0, gold0, joiner_type(), transf_type());
}

// -------------------------------------------------------------------
// test reduce overload1
//
// test for both POD types and custom scalar types
// -------------------------------------------------------------------
template <class ExecutionSpace, class ViewType, class ValueType>
void run_and_check_reduce_overloadA(ViewType view, ValueType non_trivial_result,
                                    ValueType trivial_result) {
  // trivial cases
  const auto r1 =
      KE::reduce(ExecutionSpace(), KE::cbegin(view), KE::cbegin(view));
  const auto r2 = KE::reduce("MYLABEL", ExecutionSpace(), KE::cbegin(view),
                             KE::cbegin(view));
  ASSERT_EQ(r1, trivial_result);
  ASSERT_EQ(r2, trivial_result);

  // non trivial cases
  const auto r3 =
      KE::reduce(ExecutionSpace(), KE::cbegin(view), KE::cend(view));
  const auto r4 =
      KE::reduce("MYLABEL", ExecutionSpace(), KE::cbegin(view), KE::cend(view));
  const auto r5 = KE::reduce(ExecutionSpace(), view);
  const auto r6 = KE::reduce("MYLABEL", ExecutionSpace(), view);

  ASSERT_EQ(r3, non_trivial_result);
  ASSERT_EQ(r4, non_trivial_result);
  ASSERT_EQ(r5, non_trivial_result);
  ASSERT_EQ(r6, non_trivial_result);
}

TEST_F(std_algorithms_numerics_test,
       reduce_default_functors_overloadA_using_pod_value_type) {
  fillFixtureViews();
  const value_type trivial_gold     = 0.;
  const value_type non_trivial_gold = 12.;
  run_and_check_reduce_overloadA<exespace>(m_static_view, non_trivial_gold,
                                           trivial_gold);
  run_and_check_reduce_overloadA<exespace>(m_dynamic_view, non_trivial_gold,
                                           trivial_gold);
  run_and_check_reduce_overloadA<exespace>(m_strided_view, non_trivial_gold,
                                           trivial_gold);
}

TEST_F(std_algorithms_numerics_test,
       reduce_default_functors_overloadA_using_custom_value_type) {
  fillFixtureViews();
  const CustomValueType trivial_gold{0.};
  const CustomValueType non_trivial_gold{12.};
  run_and_check_reduce_overloadA<exespace>(m_static_view_cs, non_trivial_gold,
                                           trivial_gold);
  run_and_check_reduce_overloadA<exespace>(m_dynamic_view_cs, non_trivial_gold,
                                           trivial_gold);
  run_and_check_reduce_overloadA<exespace>(m_strided_view_cs, non_trivial_gold,
                                           trivial_gold);
}

// -------------------------------------------------------------------
// test reduce overload2 with init value
//
// test for both POD types and custom scalar types
// -------------------------------------------------------------------
template <class ExecutionSpace, class ViewType, class ValueType>
void run_and_check_reduce_overloadB(ViewType view, ValueType result_value,
                                    ValueType init_value) {
  // trivial cases
  const auto r1 = KE::reduce(ExecutionSpace(), KE::cbegin(view),
                             KE::cbegin(view), init_value);
  const auto r2 = KE::reduce("MYLABEL", ExecutionSpace(), KE::cbegin(view),
                             KE::cbegin(view), init_value);
  ASSERT_EQ(r1, init_value);
  ASSERT_EQ(r2, init_value);

  // non trivial cases
  const auto r3 = KE::reduce(ExecutionSpace(), KE::cbegin(view), KE::cend(view),
                             init_value);
  const auto r4 = KE::reduce("MYLABEL", ExecutionSpace(), KE::cbegin(view),
                             KE::cend(view), init_value);
  const auto r5 = KE::reduce(ExecutionSpace(), view, init_value);
  const auto r6 = KE::reduce("MYLABEL", ExecutionSpace(), view, init_value);

  ASSERT_EQ(r3, result_value);
  ASSERT_EQ(r4, result_value);
  ASSERT_EQ(r5, result_value);
  ASSERT_EQ(r6, result_value);
}

TEST_F(std_algorithms_numerics_test,
       reduce_default_functors_overloadB_using_pod_value_type) {
  fillFixtureViews();
  const value_type init = 5.;
  const value_type gold = 17.;
  run_and_check_reduce_overloadB<exespace>(m_static_view, gold, init);
  run_and_check_reduce_overloadB<exespace>(m_dynamic_view, gold, init);
  run_and_check_reduce_overloadB<exespace>(m_strided_view, gold, init);
}

TEST_F(std_algorithms_numerics_test,
       reduce_default_functors_overloadB_using_custom_value_type) {
  fillFixtureViews();
  const CustomValueType init{5.};
  const CustomValueType gold{17.};
  run_and_check_reduce_overloadB<exespace>(m_static_view_cs, gold, init);
  run_and_check_reduce_overloadB<exespace>(m_dynamic_view_cs, gold, init);
  run_and_check_reduce_overloadB<exespace>(m_strided_view_cs, gold, init);
}

// -------------------------------------------------------------------
// test reduce overload3 with init value
//
// test for both POD types and custom scalar types
// -------------------------------------------------------------------
template <class ExecutionSpace, class ViewType, class ValueType, class BinaryOp>
void run_and_check_reduce_overloadC(ViewType view, ValueType result_value,
                                    ValueType init_value, BinaryOp joiner) {
  // trivial cases
  const auto r1 = KE::reduce(ExecutionSpace(), KE::cbegin(view),
                             KE::cbegin(view), init_value, joiner);
  const auto r2 = KE::reduce("MYLABEL", ExecutionSpace(), KE::cbegin(view),
                             KE::cbegin(view), init_value, joiner);
  ASSERT_EQ(r1, init_value);
  ASSERT_EQ(r2, init_value);

  // non trivial cases
  const auto r3 = KE::reduce(ExecutionSpace(), KE::cbegin(view), KE::cend(view),
                             init_value, joiner);
  const auto r4 = KE::reduce("MYLABEL", ExecutionSpace(), KE::cbegin(view),
                             KE::cend(view), init_value, joiner);
  const auto r5 = KE::reduce(ExecutionSpace(), view, init_value, joiner);
  const auto r6 =
      KE::reduce("MYLABEL", ExecutionSpace(), view, init_value, joiner);

  ASSERT_EQ(r3, result_value);
  ASSERT_EQ(r4, result_value);
  ASSERT_EQ(r5, result_value);
  ASSERT_EQ(r6, result_value);
}

TEST_F(std_algorithms_numerics_test,
       reduce_custom_functors_using_pod_value_type) {
  using joiner_type = SumJoinFunctor<value_type>;

  fillFixtureViews();
  const value_type init = 5.;
  const value_type gold = 17.;
  run_and_check_reduce_overloadC<exespace>(m_static_view, gold, init,
                                           joiner_type());
  run_and_check_reduce_overloadC<exespace>(m_dynamic_view, gold, init,
                                           joiner_type());
  run_and_check_reduce_overloadC<exespace>(m_strided_view, gold, init,
                                           joiner_type());
}

TEST_F(std_algorithms_numerics_test,
       reduce_custom_functors_using_custom_value_type) {
  using joiner_type = SumJoinFunctor<CustomValueType>;

  fillFixtureViews();
  const CustomValueType init{5.};
  const CustomValueType gold{17.};
  run_and_check_reduce_overloadC<exespace>(m_static_view_cs, gold, init,
                                           joiner_type());
  run_and_check_reduce_overloadC<exespace>(m_dynamic_view_cs, gold, init,
                                           joiner_type());
  run_and_check_reduce_overloadC<exespace>(m_strided_view_cs, gold, init,
                                           joiner_type());
}

#endif  // not defined KOKKOS_ENABLE_OPENMPTARGET

}  // namespace stdalgos
}  // namespace Test
