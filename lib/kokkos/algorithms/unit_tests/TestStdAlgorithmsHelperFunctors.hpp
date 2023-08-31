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
