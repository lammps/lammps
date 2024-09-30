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

#ifndef KOKKOS_STD_ALGORITHMS_CONSTRAINTS_HPP_
#define KOKKOS_STD_ALGORITHMS_CONSTRAINTS_HPP_

#include <Kokkos_DetectionIdiom.hpp>
#include <Kokkos_View.hpp>

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
                        (std::is_same<typename T::traits::array_layout,
                                      Kokkos::LayoutLeft>::value ||
                         std::is_same<typename T::traits::array_layout,
                                      Kokkos::LayoutRight>::value ||
                         std::is_same<typename T::traits::array_layout,
                                      Kokkos::LayoutStride>::value)>>
    : std::true_type {};

template <class ViewType>
KOKKOS_INLINE_FUNCTION constexpr void
static_assert_is_admissible_to_kokkos_std_algorithms(
    const ViewType& /* view */) {
  static_assert(is_admissible_to_kokkos_std_algorithms<ViewType>::value,
                "Currently, Kokkos standard algorithms only accept 1D Views.");
}

//
// is_iterator
//
template <class T>
using iterator_category_t = typename T::iterator_category;

template <class T>
using is_iterator = Kokkos::is_detected<iterator_category_t, T>;

template <class T>
inline constexpr bool is_iterator_v = is_iterator<T>::value;

template <typename ViewType>
struct is_kokkos_iterator : std::false_type {};

template <typename ViewType>
struct is_kokkos_iterator<RandomAccessIterator<ViewType>> {
  static constexpr bool value =
      is_admissible_to_kokkos_std_algorithms<ViewType>::value;
};

template <class T>
inline constexpr bool is_kokkos_iterator_v = is_kokkos_iterator<T>::value;

//
// are_iterators
//
template <class... Args>
struct are_iterators;

template <class T>
struct are_iterators<T> {
  static constexpr bool value = is_iterator_v<T>;
};

template <class Head, class... Tail>
struct are_iterators<Head, Tail...> {
  static constexpr bool value =
      are_iterators<Head>::value && (are_iterators<Tail>::value && ... && true);
};

template <class... Ts>
inline constexpr bool are_iterators_v = are_iterators<Ts...>::value;

//
// are_random_access_iterators
//
template <class... Args>
struct are_random_access_iterators;

template <class T>
struct are_random_access_iterators<T> {
  static constexpr bool value =
      is_iterator_v<T> && std::is_base_of<std::random_access_iterator_tag,
                                          typename T::iterator_category>::value;
};

template <class Head, class... Tail>
struct are_random_access_iterators<Head, Tail...> {
  static constexpr bool value =
      are_random_access_iterators<Head>::value &&
      (are_random_access_iterators<Tail>::value && ... && true);
};

template <class... Ts>
inline constexpr bool are_random_access_iterators_v =
    are_random_access_iterators<Ts...>::value;

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
  static constexpr bool value =
      std::is_same<typename T1::difference_type,
                   typename T2::difference_type>::value;
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
// not_openmptarget
//
template <class ExeSpace>
struct not_openmptarget {
#ifndef KOKKOS_ENABLE_OPENMPTARGET
  static constexpr bool value = true;
#else
  static constexpr bool value =
      !std::is_same<std::decay_t<ExeSpace>,
                    ::Kokkos::Experimental::OpenMPTarget>::value;
#endif
};

template <class ExecutionSpaceOrTeamHandleType>
KOKKOS_INLINE_FUNCTION constexpr void static_assert_is_not_openmptarget(
    const ExecutionSpaceOrTeamHandleType& /*ex_or_th*/) {
  static_assert(not_openmptarget<ExecutionSpaceOrTeamHandleType>::value,
                "Currently, Kokkos standard algorithms do not support custom "
                "comparators in OpenMPTarget");
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
  if constexpr (is_kokkos_iterator_v<IteratorType1> &&
                is_kokkos_iterator_v<IteratorType2>) {
    auto const view1 = first.view();
    auto const view2 = s_first.view();

    std::size_t stride1  = view1.stride(0);
    std::size_t stride2  = view2.stride(0);
    ptrdiff_t first_diff = view1.data() - view2.data();

    // FIXME If strides are not identical, checks may not be made
    // with the cost of O(1)
    // Currently, checks are made only if strides are identical
    // If first_diff == 0, there is already an overlap
    if (stride1 == stride2 || first_diff == 0) {
      [[maybe_unused]] bool is_no_overlap  = (first_diff % stride1);
      auto* first_pointer1                 = view1.data();
      auto* first_pointer2                 = view2.data();
      [[maybe_unused]] auto* last_pointer1 = first_pointer1 + (last - first);
      [[maybe_unused]] auto* last_pointer2 = first_pointer2 + (last - first);
      KOKKOS_EXPECTS(first_pointer1 >= last_pointer2 ||
                     last_pointer1 <= first_pointer2 || is_no_overlap);
    }
  }
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
