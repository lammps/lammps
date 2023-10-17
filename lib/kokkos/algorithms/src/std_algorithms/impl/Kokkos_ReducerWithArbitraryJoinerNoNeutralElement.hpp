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

#ifndef KOKKOS_STD_ALGORITHMS_REDUCER_WITH_ARBITRARY_JOINER_NONEUTRAL_ELEMENT_HPP
#define KOKKOS_STD_ALGORITHMS_REDUCER_WITH_ARBITRARY_JOINER_NONEUTRAL_ELEMENT_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_ValueWrapperForNoNeutralElement.hpp"

namespace Kokkos {
namespace Experimental {
namespace Impl {

// This reducer is here and not where all other reducers are
// because it is inside Impl and also because it would not work
// for OpenMPTarget backend. We can move this later.

template <class Scalar, class JoinerType, class Space = HostSpace>
struct ReducerWithArbitraryJoinerNoNeutralElement {
  using scalar_type = std::remove_cv_t<Scalar>;

 public:
  // Required
  using reducer =
      ReducerWithArbitraryJoinerNoNeutralElement<Scalar, JoinerType, Space>;
  using value_type = ValueWrapperForNoNeutralElement<scalar_type>;

  using result_view_type = Kokkos::View<value_type, Space>;

 private:
  JoinerType m_joiner;
  result_view_type m_value;
  bool m_references_scalar_v;

 public:
  KOKKOS_FUNCTION
  ReducerWithArbitraryJoinerNoNeutralElement(value_type& value_,
                                             JoinerType joiner_)
      : m_joiner(joiner_), m_value(&value_), m_references_scalar_v(true) {}

  KOKKOS_FUNCTION
  ReducerWithArbitraryJoinerNoNeutralElement(const result_view_type& value_,
                                             JoinerType joiner_)
      : m_joiner(joiner_), m_value(value_), m_references_scalar_v(false) {}

  // Required
  KOKKOS_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    dest.val = m_joiner(dest.val, src.val);
  }

  KOKKOS_FUNCTION
  void init(value_type& val) const {
    // I cannot call reduction_identity, so need to default this
    val = {};
  }

  KOKKOS_FUNCTION
  value_type& reference() const { return *m_value.data(); }

  KOKKOS_FUNCTION
  result_view_type view() const { return m_value; }

  KOKKOS_FUNCTION
  bool references_scalar() const { return m_references_scalar_v; }
};

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
