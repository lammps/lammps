// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_STD_ALGORITHMS_CONSTRAINTS_HPP_
#define KOKKOS_STD_ALGORITHMS_CONSTRAINTS_HPP_

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
import kokkos.core_impl;
#else
#include <Kokkos_Core.hpp>
#include <Kokkos_Iterator.hpp>
#endif
#include <Kokkos_Assert.hpp>

#include <Kokkos_DetectionIdiom.hpp>

#include <iterator>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class T>
class RandomAccessIterator;

template <typename T, typename enable = void>
struct is_admissible_to_kokkos_std_algorithms : std::false_type {};

template <typename T>
struct is_admissible_to_kokkos_std_algorithms<
    T, std::enable_if_t<::Kokkos::is_view<T>::value && T::rank() == 1 &&
                        (std::is_same_v<typename T::traits::array_layout,
                                        Kokkos::LayoutLeft> ||
                         std::is_same_v<typename T::traits::array_layout,
                                        Kokkos::LayoutRight> ||
                         std::is_same_v<typename T::traits::array_layout,
                                        Kokkos::LayoutStride>)>>
    : std::true_type {};

template <class ViewType>
KOKKOS_INLINE_FUNCTION constexpr void
static_assert_is_admissible_to_kokkos_std_algorithms(
    const ViewType& /* view */) {
  static_assert(is_admissible_to_kokkos_std_algorithms<ViewType>::value,
                "Currently, Kokkos standard algorithms only accept 1D Views.");
}

//
// iterators_are_accessible_from
//
template <class... Args>
struct iterators_are_accessible_from;

template <class ExeSpace, class IteratorType>
struct iterators_are_accessible_from<ExeSpace, IteratorType> {
  using view_type = typename IteratorType::view_type;
  static constexpr bool value =
      SpaceAccessibility<ExeSpace,
                         typename view_type::memory_space>::accessible;
};

template <class ExeSpace, class Head, class... Tail>
struct iterators_are_accessible_from<ExeSpace, Head, Tail...> {
  static constexpr bool value =
      iterators_are_accessible_from<ExeSpace, Head>::value &&
      iterators_are_accessible_from<ExeSpace, Tail...>::value;
};

template <class ExecutionSpaceOrTeamHandleType, class... IteratorTypes>
KOKKOS_INLINE_FUNCTION constexpr void
static_assert_random_access_and_accessible(
    const ExecutionSpaceOrTeamHandleType& /* ex_or_th*/,
    IteratorTypes... /* iterators */) {
  static_assert(
      are_random_access_iterators<IteratorTypes...>::value,
      "Currently, Kokkos standard algorithms require random access iterators.");
  static_assert(iterators_are_accessible_from<
                    typename ExecutionSpaceOrTeamHandleType::execution_space,
                    IteratorTypes...>::value,
                "Incompatible view/iterator and execution space");
}

//
// have matching difference_type
//
template <class... Args>
struct iterators_have_matching_difference_type;

template <class T>
struct iterators_have_matching_difference_type<T> {
  static constexpr bool value = true;
};

template <class T1, class T2>
struct iterators_have_matching_difference_type<T1, T2> {
  static constexpr bool value = std::is_same_v<typename T1::difference_type,
                                               typename T2::difference_type>;
};

template <class T1, class T2, class... Tail>
struct iterators_have_matching_difference_type<T1, T2, Tail...> {
  static constexpr bool value =
      iterators_have_matching_difference_type<T1, T2>::value &&
      iterators_have_matching_difference_type<T2, Tail...>::value;
};

template <class IteratorType1, class IteratorType2>
KOKKOS_INLINE_FUNCTION constexpr void
static_assert_iterators_have_matching_difference_type(IteratorType1 /* it1 */,
                                                      IteratorType2 /* it2 */) {
  static_assert(iterators_have_matching_difference_type<IteratorType1,
                                                        IteratorType2>::value,
                "Iterators do not have matching difference_type");
}

template <class IteratorType1, class IteratorType2, class IteratorType3>
KOKKOS_INLINE_FUNCTION constexpr void
static_assert_iterators_have_matching_difference_type(IteratorType1 it1,
                                                      IteratorType2 it2,
                                                      IteratorType3 it3) {
  static_assert_iterators_have_matching_difference_type(it1, it2);
  static_assert_iterators_have_matching_difference_type(it2, it3);
}

//
// valid range
//
template <class IteratorType>
KOKKOS_INLINE_FUNCTION void expect_valid_range(IteratorType first,
                                               IteratorType last) {
  // this is a no-op for release
  KOKKOS_EXPECTS(last >= first);
  // avoid compiler complaining when KOKKOS_EXPECTS is no-op
  (void)first;
  (void)last;
}

//
// Check if kokkos iterators are overlapping
//
template <typename IteratorType1, typename IteratorType2>
KOKKOS_INLINE_FUNCTION void expect_no_overlap(
    [[maybe_unused]] IteratorType1 first, [[maybe_unused]] IteratorType1 last,
    [[maybe_unused]] IteratorType2 s_first) {
  if constexpr (is_iterator_v<IteratorType1> && is_iterator_v<IteratorType2>) {
    std::size_t stride1  = first.stride();
    std::size_t stride2  = s_first.stride();
    ptrdiff_t first_diff = first.data() - s_first.data();

    // FIXME If strides are not identical, checks may not be made
    // with the cost of O(1)
    // Currently, checks are made only if strides are identical
    // If first_diff == 0, there is already an overlap
    if (stride1 == stride2 || first_diff == 0) {
      [[maybe_unused]] bool is_no_overlap  = (first_diff % stride1);
      auto* first_pointer1                 = first.data();
      auto* first_pointer2                 = s_first.data();
      [[maybe_unused]] auto* last_pointer1 = first_pointer1 + (last - first);
      [[maybe_unused]] auto* last_pointer2 = first_pointer2 + (last - first);
      KOKKOS_EXPECTS(first_pointer1 >= last_pointer2 ||
                     last_pointer1 <= first_pointer2 || is_no_overlap);
    }
  }
}

template <typename DataType1, typename... Properties1, typename DataType2,
          typename... Properties2>
KOKKOS_INLINE_FUNCTION void expect_equal_extents(
    [[maybe_unused]] const ::Kokkos::View<DataType1, Properties1...>& a,
    [[maybe_unused]] const ::Kokkos::View<DataType2, Properties2...>& b) {
  KOKKOS_EXPECTS(a.extent(0) == b.extent(0));
}

// Returns whether two views have the same extent (rank-1).
// Use this when the caller needs a boolean (e.g. to early-return),
// rather than asserting via KOKKOS_EXPECTS.
template <typename DataType1, typename... Properties1, typename DataType2,
          typename... Properties2>
KOKKOS_INLINE_FUNCTION bool have_equal_extents(
    const ::Kokkos::View<DataType1, Properties1...>& a,
    const ::Kokkos::View<DataType2, Properties2...>& b) {
  return a.extent(0) == b.extent(0);
}

//
// Check if the destination view is large enough to hold the data from the
// source (i.e. a.extent(0) <= b.extent(0)).
//
template <typename DataType1, typename... Properties1, typename DataType2,
          typename... Properties2>
KOKKOS_INLINE_FUNCTION void expect_less_or_equal_extents(
    [[maybe_unused]] const ::Kokkos::View<DataType1, Properties1...>& a,
    [[maybe_unused]] const ::Kokkos::View<DataType2, Properties2...>& b) {
  KOKKOS_EXPECTS(a.extent(0) <= b.extent(0));
}

//
// Check that both views have at least `count` elements in their first
// extent (used by count-based algorithms such as copy_n).
//
template <typename DataType1, typename... Properties1, typename DataType2,
          typename... Properties2, typename Size>
KOKKOS_INLINE_FUNCTION void expect_extents_of_at_least(
    [[maybe_unused]] const ::Kokkos::View<DataType1, Properties1...>& a,
    [[maybe_unused]] const ::Kokkos::View<DataType2, Properties2...>& b,
    [[maybe_unused]] Size count) {
  KOKKOS_EXPECTS(a.extent(0) >= count);
  KOKKOS_EXPECTS(b.extent(0) >= count);
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
