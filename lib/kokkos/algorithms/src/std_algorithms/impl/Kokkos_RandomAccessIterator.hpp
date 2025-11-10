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
#include <utility>  // declval
#include <Kokkos_Macros.hpp>
#include "Kokkos_Constraints.hpp"

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class T>
class RandomAccessIterator;

namespace {

template <typename ViewType>
struct is_always_strided {
  static_assert(is_view_v<ViewType>);

  constexpr static bool value =
#ifdef KOKKOS_ENABLE_IMPL_MDSPAN
      decltype(std::declval<ViewType>().to_mdspan())::is_always_strided();
#else
      (std::is_same_v<typename ViewType::traits::array_layout,
                      Kokkos::LayoutLeft> ||
       std::is_same_v<typename ViewType::traits::array_layout,
                      Kokkos::LayoutRight> ||
       std::is_same_v<typename ViewType::traits::array_layout,
                      Kokkos::LayoutStride>);
#endif
};

}  // namespace

template <class DataType, class... Args>
class RandomAccessIterator<::Kokkos::View<DataType, Args...>> {
 public:
  using view_type     = ::Kokkos::View<DataType, Args...>;
  using iterator_type = RandomAccessIterator<view_type>;

  using iterator_category = std::random_access_iterator_tag;
  using value_type        = typename view_type::non_const_value_type;
  using difference_type   = ptrdiff_t;
  using pointer           = typename view_type::pointer_type;
  using reference         = typename view_type::reference_type;

// oneDPL needs this alias in order not to assume the data is on the host but on
// the device, see
// https://github.com/uxlfoundation/oneDPL/blob/a045eac689f9107f50ba7b42235e9e927118e483/include/oneapi/dpl/pstl/hetero/dpcpp/utils_ranges_sycl.h#L210-L214
#ifdef KOKKOS_ENABLE_ONEDPL
  using is_passed_directly = std::true_type;
#endif

  static_assert(view_type::rank == 1 &&
                is_always_strided<::Kokkos::View<DataType, Args...>>::value);

  KOKKOS_DEFAULTED_FUNCTION RandomAccessIterator() = default;

  explicit KOKKOS_FUNCTION RandomAccessIterator(const view_type view)
      : m_data(view.data()), m_stride(view.stride_0()) {}
  explicit KOKKOS_FUNCTION RandomAccessIterator(const view_type view,
                                                ptrdiff_t current_index)
      : m_data(view.data() + current_index * view.stride_0()),
        m_stride(view.stride_0()) {}

#ifndef KOKKOS_ENABLE_CXX17  // C++20 and beyond
  template <class OtherViewType>
    requires(std::is_constructible_v<view_type, OtherViewType>)
  KOKKOS_FUNCTION explicit(!std::is_convertible_v<OtherViewType, view_type>)
      RandomAccessIterator(const RandomAccessIterator<OtherViewType>& other)
      : m_data(other.m_data), m_stride(other.m_stride) {}
#else
  template <
      class OtherViewType,
      std::enable_if_t<std::is_constructible_v<view_type, OtherViewType> &&
                           !std::is_convertible_v<OtherViewType, view_type>,
                       int> = 0>
  KOKKOS_FUNCTION explicit RandomAccessIterator(
      const RandomAccessIterator<OtherViewType>& other)
      : m_data(other.m_data), m_stride(other.m_stride) {}

  template <class OtherViewType,
            std::enable_if_t<std::is_convertible_v<OtherViewType, view_type>,
                             int> = 0>
  KOKKOS_FUNCTION RandomAccessIterator(
      const RandomAccessIterator<OtherViewType>& other)
      : m_data(other.m_data), m_stride(other.m_stride) {}
#endif

  KOKKOS_FUNCTION
  iterator_type& operator++() {
    if constexpr (is_always_contiguous)
      m_data++;
    else
      m_data += m_stride;
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
    if constexpr (is_always_contiguous)
      m_data--;
    else
      m_data -= m_stride;
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
    if constexpr (is_always_contiguous)
      return *(m_data + n);
    else
      return *(m_data + n * m_stride);
  }

  KOKKOS_FUNCTION
  iterator_type& operator+=(difference_type n) {
    if constexpr (is_always_contiguous)
      m_data += n;
    else
      m_data += n * m_stride;
    return *this;
  }

  KOKKOS_FUNCTION
  iterator_type& operator-=(difference_type n) {
    if constexpr (is_always_contiguous)
      m_data -= n;
    else
      m_data -= n * m_stride;
    return *this;
  }

  KOKKOS_FUNCTION
  iterator_type operator+(difference_type n) const {
    auto it = *this;
    it += n;
    return it;
  }

  friend iterator_type operator+(difference_type n, iterator_type other) {
    return other + n;
  }

  KOKKOS_FUNCTION
  iterator_type operator-(difference_type n) const {
    auto it = *this;
    it -= n;
    return it;
  }

  KOKKOS_FUNCTION
  difference_type operator-(iterator_type it) const {
    if constexpr (is_always_contiguous)
      return m_data - it.m_data;
    else
      return (m_data - it.m_data) / m_stride;
  }

  KOKKOS_FUNCTION
  bool operator==(iterator_type other) const {
    return m_data == other.m_data && m_stride == other.m_stride;
  }

  KOKKOS_FUNCTION
  bool operator!=(iterator_type other) const {
    return m_data != other.m_data || m_stride != other.m_stride;
  }

  KOKKOS_FUNCTION
  bool operator<(iterator_type other) const { return m_data < other.m_data; }

  KOKKOS_FUNCTION
  bool operator<=(iterator_type other) const { return m_data <= other.m_data; }

  KOKKOS_FUNCTION
  bool operator>(iterator_type other) const { return m_data > other.m_data; }

  KOKKOS_FUNCTION
  bool operator>=(iterator_type other) const { return m_data >= other.m_data; }

  KOKKOS_FUNCTION
  reference operator*() const { return *m_data; }

  KOKKOS_FUNCTION
  pointer data() const { return m_data; }

  KOKKOS_FUNCTION
  int stride() const { return m_stride; }

 private:
  pointer m_data;
  int m_stride;
  static constexpr bool is_always_contiguous =
      (std::is_same_v<typename view_type::traits::array_layout,
                      Kokkos::LayoutLeft> ||
       std::is_same_v<typename view_type::traits::array_layout,
                      Kokkos::LayoutRight>);

  // Needed for the converting constructor accepting another iterator
  template <class>
  friend class RandomAccessIterator;
};

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#ifdef KOKKOS_ENABLE_SYCL
template <class T>
struct sycl::is_device_copyable<
    Kokkos::Experimental::Impl::RandomAccessIterator<T>> : std::true_type {};
#endif

#endif
