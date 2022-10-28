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

 private:
  view_type m_view;
  ptrdiff_t m_current_index = 0;
};

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
