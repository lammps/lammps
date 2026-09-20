// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_ITERATOR_HPP
#define KOKKOS_ITERATOR_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif

#include <Kokkos_View.hpp>
#include <Kokkos_Macros.hpp>
#include "impl/Kokkos_RandomAccessIterator.hpp"

namespace Kokkos {
namespace Experimental {

namespace Impl {

template <typename T, typename enable = void>
struct is_iterable_view : std::false_type {};

template <typename T>
struct is_iterable_view<
    T, std::enable_if_t<::Kokkos::is_view<T>::value && T::rank() == 1 &&
                        (std::is_same_v<typename T::traits::array_layout,
                                        Kokkos::LayoutLeft> ||
                         std::is_same_v<typename T::traits::array_layout,
                                        Kokkos::LayoutRight> ||
                         std::is_same_v<typename T::traits::array_layout,
                                        Kokkos::LayoutStride>)>>
    : std::true_type {};

template <class ViewType>
KOKKOS_INLINE_FUNCTION constexpr void static_assert_is_iterable_view(
    const ViewType& /* view */) {
  static_assert(is_iterable_view<ViewType>::value,
                "Currently, Kokkos::(c)begin, (c)end only accept 1D Views with "
                "layout right, left or stride.");
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
      is_iterator_v<T> && std::is_base_of_v<std::random_access_iterator_tag,
                                            typename T::iterator_category>;
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
}  // namespace Impl

//
// begin, end
//
template <class DataType, class... Properties>
KOKKOS_INLINE_FUNCTION auto begin(
    const Kokkos::View<DataType, Properties...>& v) {
  Impl::static_assert_is_iterable_view(v);

  using it_t =
      Impl::RandomAccessIterator<Kokkos::View<DataType, Properties...>>;
  return it_t(v);
}

template <class DataType, class... Properties>
KOKKOS_INLINE_FUNCTION auto end(
    const Kokkos::View<DataType, Properties...>& v) {
  Impl::static_assert_is_iterable_view(v);

  using it_t =
      Impl::RandomAccessIterator<Kokkos::View<DataType, Properties...>>;
  return it_t(v, v.extent(0));
}

template <class DataType, class... Properties>
KOKKOS_INLINE_FUNCTION auto cbegin(
    const Kokkos::View<DataType, Properties...>& v) {
  Impl::static_assert_is_iterable_view(v);

  using ViewConstType =
      typename Kokkos::View<DataType, Properties...>::const_type;
  const ViewConstType cv = v;
  using it_t             = Impl::RandomAccessIterator<ViewConstType>;
  return it_t(cv);
}

template <class DataType, class... Properties>
KOKKOS_INLINE_FUNCTION auto cend(
    const Kokkos::View<DataType, Properties...>& v) {
  Impl::static_assert_is_iterable_view(v);

  using ViewConstType =
      typename Kokkos::View<DataType, Properties...>::const_type;
  const ViewConstType cv = v;
  using it_t             = Impl::RandomAccessIterator<ViewConstType>;
  return it_t(cv, cv.extent(0));
}

//
// distance
//
template <class IteratorType>
KOKKOS_INLINE_FUNCTION constexpr typename IteratorType::difference_type
distance(IteratorType first, IteratorType last) {
  static_assert(
      ::Kokkos::Experimental::Impl::are_random_access_iterators<
          IteratorType>::value,
      "Kokkos::Experimental::distance: only implemented for random access "
      "iterators.");

  return last - first;
}

}  // namespace Experimental
}  // namespace Kokkos

#endif /* #ifndef KOKKOS_ITERATOR_HPP */
