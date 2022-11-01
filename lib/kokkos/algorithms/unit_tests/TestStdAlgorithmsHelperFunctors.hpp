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

#ifndef KOKKOS_ALGORITHMS_UNITTESTS_TEST_STD_ALGOS_HELPERS_FUNCTORS_HPP
#define KOKKOS_ALGORITHMS_UNITTESTS_TEST_STD_ALGOS_HELPERS_FUNCTORS_HPP

#include <Kokkos_Core.hpp>
#include <type_traits>

namespace Test {
namespace stdalgos {

template <class ViewTypeFrom, class ViewTypeTo>
struct CopyFunctor {
  ViewTypeFrom m_view_from;
  ViewTypeTo m_view_to;

  CopyFunctor() = delete;

  CopyFunctor(const ViewTypeFrom view_from, const ViewTypeTo view_to)
      : m_view_from(view_from), m_view_to(view_to) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const { m_view_to(i) = m_view_from(i); }
};

template <class ItTypeFrom, class ViewTypeTo>
struct CopyFromIteratorFunctor {
  ItTypeFrom m_it_from;
  ViewTypeTo m_view_to;

  CopyFromIteratorFunctor(const ItTypeFrom it_from, const ViewTypeTo view_to)
      : m_it_from(it_from), m_view_to(view_to) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int) const { m_view_to() = *m_it_from; }
};

template <class ValueType>
struct IncrementElementWiseFunctor {
  KOKKOS_INLINE_FUNCTION
  void operator()(ValueType& val) const { ++val; }
};

template <class ViewType>
struct FillZeroFunctor {
  ViewType m_view;

  KOKKOS_INLINE_FUNCTION
  void operator()(int index) const {
    m_view(index) = static_cast<typename ViewType::value_type>(0);
  }

  KOKKOS_INLINE_FUNCTION
  FillZeroFunctor(ViewType viewIn) : m_view(viewIn) {}
};

template <class ValueType>
struct NoOpNonMutableFunctor {
  KOKKOS_INLINE_FUNCTION
  void operator()(const ValueType& val) const { (void)val; }
};

template <class ViewType>
struct AssignIndexFunctor {
  ViewType m_view;

  AssignIndexFunctor(ViewType view) : m_view(view) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const { m_view(i) = typename ViewType::value_type(i); }
};

template <class ValueType>
struct IsEvenFunctor {
  static_assert(std::is_integral<ValueType>::value,
                "IsEvenFunctor uses operator%, so ValueType must be int");

  KOKKOS_INLINE_FUNCTION
  bool operator()(const ValueType val) const { return (val % 2 == 0); }
};

template <class ValueType>
struct IsPositiveFunctor {
  KOKKOS_INLINE_FUNCTION
  bool operator()(const ValueType val) const { return (val > 0); }
};

template <class ValueType>
struct IsNegativeFunctor {
  KOKKOS_INLINE_FUNCTION
  bool operator()(const ValueType val) const { return (val < 0); }
};

template <class ValueType>
struct NotEqualsZeroFunctor {
  KOKKOS_INLINE_FUNCTION
  bool operator()(const ValueType val) const { return val != 0; }
};

template <class ValueType>
struct EqualsValFunctor {
  const ValueType m_value;

  EqualsValFunctor(ValueType value) : m_value(value) {}

  KOKKOS_INLINE_FUNCTION
  bool operator()(const ValueType val) const { return val == m_value; }
};

template <class ValueType1, class ValueType2>
struct CustomLessThanComparator {
  KOKKOS_INLINE_FUNCTION
  bool operator()(const ValueType1& a, const ValueType2& b) const {
    return a < b;
  }

  KOKKOS_INLINE_FUNCTION
  CustomLessThanComparator() {}
};

template <class ValueType>
struct CustomEqualityComparator {
  KOKKOS_INLINE_FUNCTION
  bool operator()(const ValueType& a, const ValueType& b) const {
    return a == b;
  }
};

template <class ValueType1, class ValueType2 = ValueType1>
struct IsEqualFunctor {
  KOKKOS_INLINE_FUNCTION
  bool operator()(const ValueType1& a, const ValueType2& b) const {
    return (a == b);
  }
};

}  // namespace stdalgos
}  // namespace Test

#endif
