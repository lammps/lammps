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

#ifndef KOKKOS_ALGORITHMS_UNITTESTS_TEST_STD_ALGOS_COMMON_HPP
#define KOKKOS_ALGORITHMS_UNITTESTS_TEST_STD_ALGOS_COMMON_HPP

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_StdAlgorithms.hpp>
#include <Kokkos_Random.hpp>
#include <TestStdAlgorithmsHelperFunctors.hpp>
#include <utility>
#include <numeric>
#include <random>

namespace Test {
namespace stdalgos {

using exespace = Kokkos::DefaultExecutionSpace;

//
// tags
//
struct DynamicTag {};
struct DynamicLayoutLeftTag {};
struct DynamicLayoutRightTag {};

// these are for rank-1
struct StridedTwoTag {};
struct StridedThreeTag {};

// these are for rank-2
struct StridedTwoRowsTag {};
struct StridedThreeRowsTag {};

#ifndef _WIN32
const std::vector<int> teamSizesToTest = {1, 2, 23, 77, 123};
#else
// avoid timeouts in AppVeyor CI
const std::vector<int> teamSizesToTest = {1, 2, 23};
#endif

// map of scenarios where the key is a description
// and the value is the extent
const std::map<std::string, std::size_t> default_scenarios = {
    {"empty", 0},          {"one-element", 1}, {"two-elements-a", 2},
    {"two-elements-b", 2}, {"small-a", 9},     {"small-b", 13},
    {"medium-a", 1003},    {"medium-b", 1003}, {"large-a", 101513},
    {"large-b", 101513}};

// see cpp file for these functions
std::string view_tag_to_string(DynamicTag);
std::string view_tag_to_string(DynamicLayoutLeftTag);
std::string view_tag_to_string(DynamicLayoutRightTag);
std::string view_tag_to_string(StridedTwoTag);
std::string view_tag_to_string(StridedThreeTag);
std::string view_tag_to_string(StridedTwoRowsTag);
std::string view_tag_to_string(StridedThreeRowsTag);

//
// overload set for create_view for rank1
//

// dynamic
template <class ValueType,
          class MemSpace = typename Kokkos::DefaultExecutionSpace::memory_space>
auto create_view(DynamicTag, std::size_t ext, const std::string label) {
  using view_t = Kokkos::View<ValueType*, MemSpace>;
  view_t view{label + "_" + view_tag_to_string(DynamicTag{}), ext};
  return view;
}

// dynamic layout left
template <class ValueType,
          class MemSpace = typename Kokkos::DefaultExecutionSpace::memory_space>
auto create_view(DynamicLayoutLeftTag, std::size_t ext,
                 const std::string label) {
  using view_t = Kokkos::View<ValueType*, Kokkos::LayoutLeft, MemSpace>;
  view_t view{label + "_" + view_tag_to_string(DynamicLayoutLeftTag{}), ext};
  return view;
}

// dynamic layout right
template <class ValueType,
          class MemSpace = typename Kokkos::DefaultExecutionSpace::memory_space>
auto create_view(DynamicLayoutRightTag, std::size_t ext,
                 const std::string label) {
  using view_t = Kokkos::View<ValueType*, Kokkos::LayoutRight, MemSpace>;
  view_t view{label + "_" + view_tag_to_string(DynamicLayoutRightTag{}), ext};
  return view;
}

// stride2
template <class ValueType,
          class MemSpace = typename Kokkos::DefaultExecutionSpace::memory_space>
auto create_view(StridedTwoTag, std::size_t ext, const std::string label) {
  using view_t = Kokkos::View<ValueType*, Kokkos::LayoutStride, MemSpace>;
  Kokkos::LayoutStride layout{ext, 2};
  view_t view{label + "_" + view_tag_to_string(StridedTwoTag{}), layout};
  return view;
}

// stride3
template <class ValueType,
          class MemSpace = typename Kokkos::DefaultExecutionSpace::memory_space>
auto create_view(StridedThreeTag, std::size_t ext, const std::string label) {
  using view_t = Kokkos::View<ValueType*, Kokkos::LayoutStride, MemSpace>;
  Kokkos::LayoutStride layout{ext, 3};
  view_t view{label + "_" + view_tag_to_string(StridedThreeTag{}), layout};
  return view;
}

//
// overload set for create_view for rank2
//

// dynamic
template <class ValueType,
          class MemSpace = typename Kokkos::DefaultExecutionSpace::memory_space>
auto create_view(DynamicTag, std::size_t ext0, std::size_t ext1,
                 const std::string label) {
  using view_t = Kokkos::View<ValueType**, MemSpace>;
  view_t view{label + "_" + view_tag_to_string(DynamicTag{}), ext0, ext1};
  return view;
}

// dynamic layout left
template <class ValueType,
          class MemSpace = typename Kokkos::DefaultExecutionSpace::memory_space>
auto create_view(DynamicLayoutLeftTag, std::size_t ext0, std::size_t ext1,
                 const std::string label) {
  using view_t = Kokkos::View<ValueType**, Kokkos::LayoutLeft, MemSpace>;
  view_t view{label + "_" + view_tag_to_string(DynamicLayoutLeftTag{}), ext0,
              ext1};
  return view;
}

// dynamic layout right
template <class ValueType,
          class MemSpace = typename Kokkos::DefaultExecutionSpace::memory_space>
auto create_view(DynamicLayoutRightTag, std::size_t ext0, std::size_t ext1,
                 const std::string label) {
  using view_t = Kokkos::View<ValueType**, Kokkos::LayoutRight, MemSpace>;
  view_t view{label + "_" + view_tag_to_string(DynamicLayoutRightTag{}), ext0,
              ext1};
  return view;
}

// stride2rows
template <class ValueType,
          class MemSpace = typename Kokkos::DefaultExecutionSpace::memory_space>
auto create_view(StridedTwoRowsTag, std::size_t ext0, std::size_t ext1,
                 const std::string label) {
  using view_t = Kokkos::View<ValueType**, Kokkos::LayoutStride, MemSpace>;
  Kokkos::LayoutStride layout{ext0, 2, ext1, ext0 * 2};
  view_t view{label + "_" + view_tag_to_string(StridedTwoRowsTag{}), layout};
  return view;
}

// stride3rows
template <class ValueType,
          class MemSpace = typename Kokkos::DefaultExecutionSpace::memory_space>
auto create_view(StridedThreeRowsTag, std::size_t ext0, std::size_t ext1,
                 const std::string label) {
  using view_t = Kokkos::View<ValueType**, Kokkos::LayoutStride, MemSpace>;
  Kokkos::LayoutStride layout{ext0, 3, ext1, ext0 * 3};
  view_t view{label + "_" + view_tag_to_string(StridedThreeRowsTag{}), layout};
  return view;
}

template <class ViewType>
auto create_deep_copyable_compatible_view_with_same_extent(ViewType view) {
  using view_value_type  = typename ViewType::value_type;
  using view_exespace    = typename ViewType::execution_space;
  const std::size_t ext0 = view.extent(0);
  if constexpr (ViewType::rank == 1) {
    using view_deep_copyable_t = Kokkos::View<view_value_type*, view_exespace>;
    return view_deep_copyable_t{"view_dc", ext0};
  } else {
    static_assert(ViewType::rank == 2, "Only rank 1 or 2 supported.");
    using view_deep_copyable_t = Kokkos::View<view_value_type**, view_exespace>;
    const std::size_t ext1     = view.extent(1);
    return view_deep_copyable_t{"view_dc", ext0, ext1};
  }

  // this is needed for intel to avoid
  // error #1011: missing return statement at end of non-void function
#if defined KOKKOS_COMPILER_INTEL ||                                  \
    (defined(KOKKOS_COMPILER_NVCC) && KOKKOS_COMPILER_NVCC >= 1130 && \
     !defined(KOKKOS_COMPILER_MSVC))
  __builtin_unreachable();
#endif
}

template <class ViewType>
auto create_deep_copyable_compatible_clone(ViewType view) {
  auto view_dc    = create_deep_copyable_compatible_view_with_same_extent(view);
  using view_dc_t = decltype(view_dc);
  using exe_space = typename view_dc_t::execution_space;
  if constexpr (ViewType::rank == 1) {
    CopyFunctor<ViewType, view_dc_t> F1(view, view_dc);
    Kokkos::RangePolicy<exe_space> policy(0, view.extent(0));
    Kokkos::parallel_for("copy", policy, F1);

  } else {
    static_assert(ViewType::rank == 2, "Only rank 1 or 2 supported.");
    CopyFunctorRank2<ViewType, view_dc_t> F1(view, view_dc);
    Kokkos::RangePolicy<exe_space> policy(0, view.extent(0) * view.extent(1));
    Kokkos::parallel_for("copy", policy, F1);
  }
  return view_dc;
}

//
// others
//

template <class TeamHandleType, class ValueType1, class ValueType2>
KOKKOS_FUNCTION bool team_members_have_matching_result(
    const TeamHandleType& teamHandle, const ValueType1 memberValueIn,
    const ValueType2 targetIn) {
  using T             = std::common_type_t<ValueType1, ValueType2>;
  const T memberValue = memberValueIn;
  const T target      = targetIn;

  // set accum to 1 if a mismach is found
  const bool mismatch = memberValue != target;
  int accum           = static_cast<int>(mismatch);
  // FIXME_OPENMPTARGET: team API does not meet the TeamHandle concept and
  // ignores the reducer passed
#if defined KOKKOS_ENABLE_OPENMPTARGET
  Kokkos::Sum<int> dummyReducer(accum);
  const auto result = teamHandle.team_reduce(accum, dummyReducer);
  return (result == 0);
#else
  teamHandle.team_reduce(Kokkos::Sum<int>(accum));
  return (accum == 0);
#endif
}

template <class ValueType1, class ValueType2>
auto make_bounds(const ValueType1& lower, const ValueType2 upper) {
  return Kokkos::pair<ValueType1, ValueType2>{lower, upper};
}

// libstdc++ as provided by GCC 8 does not have reduce, transform_reduce,
// exclusive_scan, inclusive_scan, transform_exclusive_scan,
// transform_inclusive_scan and for GCC 9.1, 9.2 fails to compile them for
// missing overload not accepting policy so use here simplified versions of
// them, only for testing purpose
#if defined(_GLIBCXX_RELEASE) && (_GLIBCXX_RELEASE <= 9)

template <class InputIterator, class ValueType, class BinaryOp>
ValueType testing_reduce(InputIterator first, InputIterator last,
                         ValueType initIn, BinaryOp binOp) {
  using value_type = std::remove_const_t<ValueType>;
  value_type init  = initIn;

  while (last - first >= 4) {
    ValueType v1 = binOp(first[0], first[1]);
    ValueType v2 = binOp(first[2], first[3]);
    ValueType v3 = binOp(v1, v2);
    init         = binOp(init, v3);
    first += 4;
  }

  for (; first != last; ++first) {
    init = binOp(init, *first);
  }

  return init;
}

template <class InputIterator, class ValueType>
ValueType testing_reduce(InputIterator first, InputIterator last,
                         ValueType init) {
  return testing_reduce(
      first, last, init,
      [](const ValueType& lhs, const ValueType& rhs) { return lhs + rhs; });
}

template <class InputIterator>
auto testing_reduce(InputIterator first, InputIterator last) {
  using ValueType = typename InputIterator::value_type;
  return testing_reduce(
      first, last, ValueType{},
      [](const ValueType& lhs, const ValueType& rhs) { return lhs + rhs; });
}

template <class InputIterator1, class InputIterator2, class ValueType,
          class BinaryJoiner, class BinaryTransform>
ValueType testing_transform_reduce(InputIterator1 first1, InputIterator1 last1,
                                   InputIterator2 first2, ValueType initIn,
                                   BinaryJoiner binJoiner,
                                   BinaryTransform binTransform) {
  using value_type = std::remove_const_t<ValueType>;
  value_type init  = initIn;

  while (last1 - first1 >= 4) {
    ValueType v1 = binJoiner(binTransform(first1[0], first2[0]),
                             binTransform(first1[1], first2[1]));

    ValueType v2 = binJoiner(binTransform(first1[2], first2[2]),
                             binTransform(first1[3], first2[3]));

    ValueType v3 = binJoiner(v1, v2);
    init         = binJoiner(init, v3);

    first1 += 4;
    first2 += 4;
  }

  for (; first1 != last1; ++first1, ++first2) {
    init = binJoiner(init, binTransform(*first1, *first2));
  }

  return init;
}

template <class InputIterator1, class InputIterator2, class ValueType>
ValueType testing_transform_reduce(InputIterator1 first1, InputIterator1 last1,
                                   InputIterator2 first2, ValueType init) {
  return testing_transform_reduce(
      first1, last1, first2, init,
      [](const ValueType& lhs, const ValueType& rhs) { return lhs + rhs; },
      [](const ValueType& lhs, const ValueType& rhs) { return lhs * rhs; });
}

template <class InputIterator, class ValueType, class BinaryJoiner,
          class UnaryTransform>
ValueType testing_transform_reduce(InputIterator first, InputIterator last,
                                   ValueType initIn, BinaryJoiner binJoiner,
                                   UnaryTransform unaryTransform) {
  using value_type = std::remove_const_t<ValueType>;
  value_type init  = initIn;

  while (last - first >= 4) {
    ValueType v1 =
        binJoiner(unaryTransform(first[0]), unaryTransform(first[1]));
    ValueType v2 =
        binJoiner(unaryTransform(first[2]), unaryTransform(first[3]));
    ValueType v3 = binJoiner(v1, v2);
    init         = binJoiner(init, v3);
    first += 4;
  }

  for (; first != last; ++first) {
    init = binJoiner(init, unaryTransform(*first));
  }

  return init;
}

/*
   EXCLUSIVE_SCAN
 */
template <class InputIterator, class OutputIterator, class ValueType,
          class BinaryOp>
OutputIterator testing_exclusive_scan(InputIterator first, InputIterator last,
                                      OutputIterator result, ValueType initIn,
                                      BinaryOp binOp) {
  using value_type = std::remove_const_t<ValueType>;
  value_type init  = initIn;

  while (first != last) {
    auto v = init;
    init   = binOp(init, *first);
    ++first;
    *result++ = v;
  }

  return result;
}

template <class InputIterator, class OutputIterator, class ValueType>
OutputIterator testing_exclusive_scan(InputIterator first, InputIterator last,
                                      OutputIterator result, ValueType init) {
  return testing_exclusive_scan(
      first, last, result, init,
      [](const ValueType& lhs, const ValueType& rhs) { return lhs + rhs; });
}

/*
   INCLUSIVE_SCAN
 */
template <class InputIterator, class OutputIterator, class BinaryOp,
          class ValueType>
OutputIterator testing_inclusive_scan(InputIterator first, InputIterator last,
                                      OutputIterator result, BinaryOp binOp,
                                      ValueType initIn) {
  using value_type = std::remove_const_t<ValueType>;
  value_type init  = initIn;
  for (; first != last; ++first) {
    init      = binOp(init, *first);
    *result++ = init;
  }

  return result;
}

template <class InputIterator, class OutputIterator, class BinaryOp>
OutputIterator testing_inclusive_scan(InputIterator first, InputIterator last,
                                      OutputIterator result, BinaryOp bop) {
  if (first != last) {
    auto init = *first;
    *result++ = init;
    ++first;
    if (first != last) {
      result = testing_inclusive_scan(first, last, result, bop, init);
    }
  }
  return result;
}

template <class InputIterator, class OutputIterator>
OutputIterator testing_inclusive_scan(InputIterator first, InputIterator last,
                                      OutputIterator result) {
  using ValueType = typename InputIterator::value_type;
  return testing_inclusive_scan(
      first, last, result,
      [](const ValueType& lhs, const ValueType& rhs) { return lhs + rhs; });
}

/*
   TRANSFORM_EXCLUSIVE_SCAN
 */
template <class InputIterator, class OutputIterator, class ValueType,
          class BinaryOp, class UnaryOp>
OutputIterator testing_transform_exclusive_scan(
    InputIterator first, InputIterator last, OutputIterator result,
    ValueType initIn, BinaryOp binOp, UnaryOp unaryOp) {
  using value_type = std::remove_const_t<ValueType>;
  value_type init  = initIn;

  while (first != last) {
    auto v = init;
    init   = binOp(init, unaryOp(*first));
    ++first;
    *result++ = v;
  }

  return result;
}

template <class InputIterator, class OutputIterator, class BinaryOp,
          class UnaryOp, class ValueType>
OutputIterator testing_transform_inclusive_scan(InputIterator first,
                                                InputIterator last,
                                                OutputIterator result,
                                                BinaryOp binOp, UnaryOp unaryOp,
                                                ValueType initIn) {
  using value_type = std::remove_const_t<ValueType>;
  value_type init  = initIn;

  for (; first != last; ++first) {
    init      = binOp(init, unaryOp(*first));
    *result++ = init;
  }

  return result;
}

template <class InputIterator, class OutputIterator, class BinaryOp,
          class UnaryOp>
OutputIterator testing_transform_inclusive_scan(InputIterator first,
                                                InputIterator last,
                                                OutputIterator result,
                                                BinaryOp binOp,
                                                UnaryOp unaryOp) {
  if (first != last) {
    auto init = unaryOp(*first);
    *result++ = init;
    ++first;
    if (first != last) {
      result = testing_transform_inclusive_scan(first, last, result, binOp,
                                                unaryOp, init);
    }
  }

  return result;
}

#endif

template <class LayoutTagType, class ValueType>
auto create_random_view_and_host_clone(
    LayoutTagType LayoutTag, std::size_t numRows, std::size_t numCols,
    Kokkos::pair<ValueType, ValueType> bounds, const std::string& label,
    std::size_t seedIn = 12371) {
  // construct in memory space associated with default exespace
  auto dataView = create_view<ValueType>(LayoutTag, numRows, numCols, label);

  // dataView might not deep copyable (e.g. strided layout) so to
  // randomize it, we make a new view that is for sure deep copyable,
  // modify it on the host, deep copy to device and then launch
  // a kernel to copy to dataView
  auto dataView_dc =
      create_deep_copyable_compatible_view_with_same_extent(dataView);
  auto dataView_dc_h = create_mirror_view(Kokkos::HostSpace(), dataView_dc);

  // randomly fill the view
  Kokkos::Random_XorShift64_Pool<Kokkos::DefaultHostExecutionSpace> pool(
      seedIn);
  Kokkos::fill_random(dataView_dc_h, pool, bounds.first, bounds.second);

  // copy to dataView_dc and then to dataView
  Kokkos::deep_copy(dataView_dc, dataView_dc_h);
  // use CTAD
  CopyFunctorRank2 F1(dataView_dc, dataView);
  Kokkos::parallel_for("copy", dataView.extent(0) * dataView.extent(1), F1);

  return std::make_pair(dataView, dataView_dc_h);
}

template <class ViewType>
auto create_host_space_copy(ViewType view) {
  auto view_dc = create_deep_copyable_compatible_clone(view);
  return create_mirror_view_and_copy(Kokkos::HostSpace(), view_dc);
}

// fill the views with sequentially increasing values
template <class ViewType, class ViewHostType>
void fill_views_inc(ViewType view, ViewHostType host_view) {
  namespace KE = Kokkos::Experimental;

  Kokkos::parallel_for(view.extent(0), AssignIndexFunctor<ViewType>(view));
  std::iota(KE::begin(host_view), KE::end(host_view), 0);
  // compare_views(expected, view);
}

template <class ValueType, class ViewType>
std::enable_if_t<!std::is_same<typename ViewType::traits::array_layout,
                               Kokkos::LayoutStride>::value>
verify_values(ValueType expected, const ViewType view) {
  static_assert(std::is_same<ValueType, typename ViewType::value_type>::value,
                "Non-matching value types of view and reference value");
  auto view_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), view);
  for (std::size_t i = 0; i < view_h.extent(0); i++) {
    ASSERT_EQ(expected, view_h(i));
  }
}

template <class ValueType, class ViewType>
std::enable_if_t<std::is_same<typename ViewType::traits::array_layout,
                              Kokkos::LayoutStride>::value>
verify_values(ValueType expected, const ViewType view) {
  static_assert(std::is_same<ValueType, typename ViewType::value_type>::value,
                "Non-matching value types of view and reference value");

  using non_strided_view_t = Kokkos::View<typename ViewType::value_type*>;
  non_strided_view_t tmpView("tmpView", view.extent(0));

  Kokkos::parallel_for(
      "_std_algo_copy", view.extent(0),
      CopyFunctor<ViewType, non_strided_view_t>(view, tmpView));
  auto view_h =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), tmpView);
  for (std::size_t i = 0; i < view_h.extent(0); i++) {
    ASSERT_EQ(expected, view_h(i));
  }
}

template <class ViewType1, class ViewType2>
std::enable_if_t<!std::is_same<typename ViewType2::traits::array_layout,
                               Kokkos::LayoutStride>::value>
compare_views(ViewType1 expected, const ViewType2 actual) {
  static_assert(std::is_same<typename ViewType1::value_type,
                             typename ViewType2::value_type>::value,
                "Non-matching value types of expected and actual view");
  auto expected_h =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), expected);
  auto actual_h =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), actual);

  for (std::size_t i = 0; i < expected_h.extent(0); i++) {
    ASSERT_EQ(expected_h(i), actual_h(i));
  }
}

template <class ViewType1, class ViewType2>
std::enable_if_t<std::is_same<typename ViewType2::traits::array_layout,
                              Kokkos::LayoutStride>::value>
compare_views(ViewType1 expected, const ViewType2 actual) {
  static_assert(std::is_same<typename ViewType1::value_type,
                             typename ViewType2::value_type>::value,
                "Non-matching value types of expected and actual view");

  using non_strided_view_t = Kokkos::View<typename ViewType2::value_type*>;
  non_strided_view_t tmp_view("tmp_view", actual.extent(0));
  Kokkos::parallel_for(
      "_std_algo_copy", actual.extent(0),
      CopyFunctor<ViewType2, non_strided_view_t>(actual, tmp_view));

  auto actual_h =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), tmp_view);
  auto expected_h =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), expected);

  for (std::size_t i = 0; i < expected_h.extent(0); i++) {
    ASSERT_EQ(expected_h(i), actual_h(i));
  }
}

template <class ViewType1, class ViewType2>
void expect_equal_host_views(ViewType1 A, const ViewType2 B) {
  static_assert(
      ViewType1::rank == 2 && ViewType2::rank == 2 &&
          std::is_same_v<typename ViewType1::memory_space, Kokkos::HostSpace> &&
          std::is_same_v<typename ViewType2::memory_space, Kokkos::HostSpace>,
      "Expected 2-dimensional host view.");
  ASSERT_EQ(A.extent(0), B.extent(0));
  ASSERT_EQ(A.extent(1), B.extent(1));

  constexpr bool values_are_floast =
      std::is_floating_point_v<typename ViewType1::value_type> ||
      std::is_floating_point_v<typename ViewType2::value_type>;

  for (std::size_t i = 0; i < A.extent(0); i++) {
    for (std::size_t j = 0; j < A.extent(1); j++) {
      if constexpr (values_are_floast) {
        EXPECT_FLOAT_EQ(A(i, j), B(i, j));
      } else {
        ASSERT_EQ(A(i, j), B(i, j));
      }
    }
  }
}

template <class ViewType>
void fill_zero(ViewType a) {
  const auto functor = FillZeroFunctor<ViewType>(a);
  ::Kokkos::parallel_for(a.extent(0), std::move(functor));
}

template <class ViewType1, class ViewType2>
void fill_zero(ViewType1 a, ViewType2 b) {
  fill_zero(a);
  fill_zero(b);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// helpers for testing small views (extent = 10)
// prefer `default_scenarios` map for creating new tests
using value_type = double;

struct std_algorithms_test : public ::testing::Test {
  static constexpr size_t extent = 10;

  using static_view_t = Kokkos::View<value_type[extent]>;
  static_view_t m_static_view{"std-algo-test-1D-contiguous-view-static"};

  using dyn_view_t = Kokkos::View<value_type*>;
  dyn_view_t m_dynamic_view{"std-algo-test-1D-contiguous-view-dynamic", extent};

  using strided_view_t = Kokkos::View<value_type*, Kokkos::LayoutStride>;
  Kokkos::LayoutStride layout{extent, 2};
  strided_view_t m_strided_view{"std-algo-test-1D-strided-view", layout};

  using view_host_space_t = Kokkos::View<value_type[10], Kokkos::HostSpace>;

  template <class ViewFromType>
  void copyInputViewToFixtureViews(ViewFromType view) {
    CopyFunctor<ViewFromType, static_view_t> F1(view, m_static_view);
    Kokkos::parallel_for("_std_algo_copy1", view.extent(0), F1);

    CopyFunctor<ViewFromType, dyn_view_t> F2(view, m_dynamic_view);
    Kokkos::parallel_for("_std_algo_copy2", view.extent(0), F2);

    CopyFunctor<ViewFromType, strided_view_t> F3(view, m_strided_view);
    Kokkos::parallel_for("_std_algo_copy3", view.extent(0), F3);
  }
};

struct CustomValueType {
  KOKKOS_INLINE_FUNCTION
  CustomValueType(){};

  KOKKOS_INLINE_FUNCTION
  CustomValueType(value_type val) : value(val){};

  KOKKOS_INLINE_FUNCTION
  CustomValueType(const CustomValueType& other) { this->value = other.value; }

  KOKKOS_INLINE_FUNCTION
  explicit operator value_type() const { return value; }

  KOKKOS_INLINE_FUNCTION
  value_type& operator()() { return value; }

  KOKKOS_INLINE_FUNCTION
  const value_type& operator()() const { return value; }

  KOKKOS_INLINE_FUNCTION
  CustomValueType& operator+=(const CustomValueType& other) {
    this->value += other.value;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  CustomValueType& operator=(const CustomValueType& other) {
    this->value = other.value;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  CustomValueType operator+(const CustomValueType& other) const {
    CustomValueType result;
    result.value = this->value + other.value;
    return result;
  }

  KOKKOS_INLINE_FUNCTION
  CustomValueType operator-(const CustomValueType& other) const {
    CustomValueType result;
    result.value = this->value - other.value;
    return result;
  }

  KOKKOS_INLINE_FUNCTION
  CustomValueType operator*(const CustomValueType& other) const {
    CustomValueType result;
    result.value = this->value * other.value;
    return result;
  }

  KOKKOS_INLINE_FUNCTION
  bool operator==(const CustomValueType& other) const {
    return this->value == other.value;
  }

 private:
  friend std::ostream& operator<<(std::ostream& os,
                                  const CustomValueType& custom_value_type) {
    return os << custom_value_type.value;
  }
  value_type value = {};
};

}  // namespace stdalgos
}  // namespace Test

#endif
