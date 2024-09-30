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

#ifndef KOKKOS_RANDOM_ACCESS_ITERATOR_IMPL_HPP
#define KOKKOS_RANDOM_ACCESS_ITERATOR_IMPL_HPP

#include <iterator>
#include <Kokkos_Macros.hpp>
#include <Kokkos_View.hpp>
#include "Kokkos_Constraints.hpp"

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class T>
class RandomAccessIterator;

template <class DataType, class... Args>
class RandomAccessIterator< ::Kokkos::View<DataType, Args...> > {
 public:
  using view_type     = ::Kokkos::View<DataType, Args...>;
  using iterator_type = RandomAccessIterator<view_type>;

  using iterator_category = std::random_access_iterator_tag;
  using value_type        = typename view_type::value_type;
  using difference_type   = ptrdiff_t;
  using pointer           = typename view_type::pointer_type;
  using reference         = typename view_type::reference_type;

  static_assert(view_type::rank == 1 &&
                    (std::is_same<typename view_type::traits::array_layout,
                                  Kokkos::LayoutLeft>::value ||
                     std::is_same<typename view_type::traits::array_layout,
                                  Kokkos::LayoutRight>::value ||
                     std::is_same<typename view_type::traits::array_layout,
                                  Kokkos::LayoutStride>::value),
                "RandomAccessIterator only supports 1D Views with LayoutLeft, "
                "LayoutRight, LayoutStride.");

  KOKKOS_DEFAULTED_FUNCTION RandomAccessIterator() = default;

  explicit KOKKOS_FUNCTION RandomAccessIterator(const view_type view)
      : m_view(view) {}
  explicit KOKKOS_FUNCTION RandomAccessIterator(const view_type view,
                                                ptrdiff_t current_index)
      : m_view(view), m_current_index(current_index) {}

#ifndef KOKKOS_ENABLE_CXX17  // C++20 and beyond
  template <class OtherViewType>
  requires(std::is_constructible_v<view_type, OtherViewType>) KOKKOS_FUNCTION
      explicit(!std::is_convertible_v<OtherViewType, view_type>)
          RandomAccessIterator(const RandomAccessIterator<OtherViewType>& other)
      : m_view(other.m_view), m_current_index(other.m_current_index) {}
#else
  template <
      class OtherViewType,
      std::enable_if_t<std::is_constructible_v<view_type, OtherViewType> &&
                           !std::is_convertible_v<OtherViewType, view_type>,
                       int> = 0>
  KOKKOS_FUNCTION explicit RandomAccessIterator(
      const RandomAccessIterator<OtherViewType>& other)
      : m_view(other.m_view), m_current_index(other.m_current_index) {}

  template <class OtherViewType,
            std::enable_if_t<std::is_convertible_v<OtherViewType, view_type>,
                             int> = 0>
  KOKKOS_FUNCTION RandomAccessIterator(
      const RandomAccessIterator<OtherViewType>& other)
      : m_view(other.m_view), m_current_index(other.m_current_index) {}
#endif

  KOKKOS_FUNCTION
  iterator_type& operator++() {
    ++m_current_index;
    return *this;
  }

  KOKKOS_FUNCTION
  iterator_type operator++(int) {
    auto tmp = *this;
    ++*this;
    return tmp;
  }

  KOKKOS_FUNCTION
  iterator_type& operator--() {
    --m_current_index;
    return *this;
  }

  KOKKOS_FUNCTION
  iterator_type operator--(int) {
    auto tmp = *this;
    --*this;
    return tmp;
  }

  KOKKOS_FUNCTION
  reference operator[](difference_type n) const {
    return m_view(m_current_index + n);
  }

  KOKKOS_FUNCTION
  iterator_type& operator+=(difference_type n) {
    m_current_index += n;
    return *this;
  }

  KOKKOS_FUNCTION
  iterator_type& operator-=(difference_type n) {
    m_current_index -= n;
    return *this;
  }

  KOKKOS_FUNCTION
  iterator_type operator+(difference_type n) const {
    return iterator_type(m_view, m_current_index + n);
  }

  KOKKOS_FUNCTION
  iterator_type operator-(difference_type n) const {
    return iterator_type(m_view, m_current_index - n);
  }

  KOKKOS_FUNCTION
  difference_type operator-(iterator_type it) const {
    return m_current_index - it.m_current_index;
  }

  KOKKOS_FUNCTION
  bool operator==(iterator_type other) const {
    return m_current_index == other.m_current_index &&
           m_view.data() == other.m_view.data();
  }

  KOKKOS_FUNCTION
  bool operator!=(iterator_type other) const {
    return m_current_index != other.m_current_index ||
           m_view.data() != other.m_view.data();
  }

  KOKKOS_FUNCTION
  bool operator<(iterator_type other) const {
    return m_current_index < other.m_current_index;
  }

  KOKKOS_FUNCTION
  bool operator<=(iterator_type other) const {
    return m_current_index <= other.m_current_index;
  }

  KOKKOS_FUNCTION
  bool operator>(iterator_type other) const {
    return m_current_index > other.m_current_index;
  }

  KOKKOS_FUNCTION
  bool operator>=(iterator_type other) const {
    return m_current_index >= other.m_current_index;
  }

  KOKKOS_FUNCTION
  reference operator*() const { return m_view(m_current_index); }

  KOKKOS_FUNCTION
  view_type view() const { return m_view; }

 private:
  view_type m_view;
  ptrdiff_t m_current_index = 0;

  // Needed for the converting constructor accepting another iterator
  template <class>
  friend class RandomAccessIterator;
};

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
