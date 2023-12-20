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

#ifndef KOKKOS_STD_ALGORITHMS_HELPER_PREDICATES_HPP
#define KOKKOS_STD_ALGORITHMS_HELPER_PREDICATES_HPP

#include <Kokkos_Macros.hpp>

// naming convetion:
// StdAlgoSomeExpressiveNameUnaryPredicate
// StdAlgoSomeExpressiveNameBinaryPredicate

namespace Kokkos {
namespace Experimental {
namespace Impl {

// ------------------
// UNARY PREDICATES
// ------------------
template <class T>
struct StdAlgoEqualsValUnaryPredicate {
  T m_value;

  KOKKOS_FUNCTION
  constexpr bool operator()(const T& val) const { return val == m_value; }

  KOKKOS_FUNCTION
  constexpr explicit StdAlgoEqualsValUnaryPredicate(const T& _value)
      : m_value(_value) {}
};

template <class T>
struct StdAlgoNotEqualsValUnaryPredicate {
  T m_value;

  KOKKOS_FUNCTION
  constexpr bool operator()(const T& val) const { return !(val == m_value); }

  KOKKOS_FUNCTION
  constexpr explicit StdAlgoNotEqualsValUnaryPredicate(const T& _value)
      : m_value(_value) {}
};

template <class ValueType, class PredicateType>
struct StdAlgoNegateUnaryPredicateWrapper {
  PredicateType m_pred;

  KOKKOS_FUNCTION
  constexpr bool operator()(const ValueType& val) const { return !m_pred(val); }

  KOKKOS_FUNCTION
  constexpr explicit StdAlgoNegateUnaryPredicateWrapper(
      const PredicateType& pred)
      : m_pred(pred) {}
};

// ------------------
// BINARY PREDICATES
// ------------------
template <class ValueType1, class ValueType2 = ValueType1>
struct StdAlgoEqualBinaryPredicate {
  KOKKOS_FUNCTION
  constexpr bool operator()(const ValueType1& a, const ValueType2& b) const {
    return a == b;
  }
};

template <class ValueType1, class ValueType2 = ValueType1>
struct StdAlgoLessThanBinaryPredicate {
  KOKKOS_FUNCTION
  constexpr bool operator()(const ValueType1& a, const ValueType2& b) const {
    return a < b;
  }
};

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos
#endif
