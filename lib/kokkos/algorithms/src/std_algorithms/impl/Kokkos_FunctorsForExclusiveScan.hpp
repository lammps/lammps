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

#ifndef KOKKOS_STD_ALGORITHMS_FUNCTORS_FOR_EXCLUSIVE_SCAN_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_FUNCTORS_FOR_EXCLUSIVE_SCAN_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_ValueWrapperForNoNeutralElement.hpp"

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <typename ValueType>
using ex_scan_has_reduction_identity_sum_t =
    decltype(Kokkos::reduction_identity<ValueType>::sum());

template <class ExeSpace, class IndexType, class ValueType, class FirstFrom,
          class FirstDest>
struct ExclusiveScanDefaultFunctorForKnownNeutralElement {
  using execution_space = ExeSpace;
  ValueType m_init_value;
  FirstFrom m_first_from;
  FirstDest m_first_dest;

  KOKKOS_FUNCTION
  ExclusiveScanDefaultFunctorForKnownNeutralElement(ValueType init,
                                                    FirstFrom first_from,
                                                    FirstDest first_dest)
      : m_init_value(std::move(init)),
        m_first_from(std::move(first_from)),
        m_first_dest(std::move(first_dest)) {}

  KOKKOS_FUNCTION
  void operator()(const IndexType i, ValueType& update,
                  const bool final_pass) const {
    const auto tmp = m_first_from[i];
    if (final_pass) m_first_dest[i] = update + m_init_value;
    update += tmp;
  }
};

template <class ExeSpace, class IndexType, class ValueType, class FirstFrom,
          class FirstDest>
struct ExclusiveScanDefaultFunctorWithValueWrapper {
  using execution_space = ExeSpace;
  using value_type =
      ::Kokkos::Experimental::Impl::ValueWrapperForNoNeutralElement<ValueType>;
  ValueType m_init_value;
  FirstFrom m_first_from;
  FirstDest m_first_dest;

  KOKKOS_FUNCTION
  ExclusiveScanDefaultFunctorWithValueWrapper(ValueType init,
                                              FirstFrom first_from,
                                              FirstDest first_dest)
      : m_init_value(std::move(init)),
        m_first_from(std::move(first_from)),
        m_first_dest(std::move(first_dest)) {}

  KOKKOS_FUNCTION
  void operator()(const IndexType i, value_type& update,
                  const bool final_pass) const {
    const auto tmp = value_type{m_first_from[i], false};
    if (final_pass) {
      if (i == 0) {
        m_first_dest[i] = m_init_value;
      } else {
        m_first_dest[i] = update.val + m_init_value;
      }
    }

    this->join(update, tmp);
  }

  KOKKOS_FUNCTION
  void init(value_type& update) const {
    update.val        = {};
    update.is_initial = true;
  }

  KOKKOS_FUNCTION
  void join(value_type& update, const value_type& input) const {
    if (input.is_initial) return;

    if (update.is_initial) {
      update.val        = input.val;
      update.is_initial = false;
    } else {
      update.val = update.val + input.val;
    }
  }
};

template <class ExeSpace, class IndexType, class ValueType, class FirstFrom,
          class FirstDest, class BinaryOpType, class UnaryOpType>
struct TransformExclusiveScanFunctorWithValueWrapper {
  using execution_space = ExeSpace;
  using value_type =
      ::Kokkos::Experimental::Impl::ValueWrapperForNoNeutralElement<ValueType>;

  ValueType m_init_value;
  FirstFrom m_first_from;
  FirstDest m_first_dest;
  BinaryOpType m_binary_op;
  UnaryOpType m_unary_op;

  KOKKOS_FUNCTION
  TransformExclusiveScanFunctorWithValueWrapper(ValueType init,
                                                FirstFrom first_from,
                                                FirstDest first_dest,
                                                BinaryOpType bop,
                                                UnaryOpType uop)
      : m_init_value(std::move(init)),
        m_first_from(std::move(first_from)),
        m_first_dest(std::move(first_dest)),
        m_binary_op(std::move(bop)),
        m_unary_op(std::move(uop)) {}

  KOKKOS_FUNCTION
  void operator()(const IndexType i, value_type& update,
                  const bool final_pass) const {
    const auto tmp = value_type{m_unary_op(m_first_from[i]), false};
    if (final_pass) {
      if (i == 0) {
        // for both ExclusiveScan and TransformExclusiveScan,
        // init is unmodified
        m_first_dest[i] = m_init_value;
      } else {
        m_first_dest[i] = m_binary_op(update.val, m_init_value);
      }
    }

    this->join(update, tmp);
  }

  KOKKOS_FUNCTION void init(value_type& value) const {
    value.val        = {};
    value.is_initial = true;
  }

  KOKKOS_FUNCTION
  void join(value_type& update, const value_type& input) const {
    if (input.is_initial) return;

    if (update.is_initial) {
      update.val = input.val;
    } else {
      update.val = m_binary_op(update.val, input.val);
    }
    update.is_initial = false;
  }
};

template <class ExeSpace, class IndexType, class ValueType, class FirstFrom,
          class FirstDest, class BinaryOpType, class UnaryOpType>
struct TransformExclusiveScanFunctorWithoutValueWrapper {
  using execution_space = ExeSpace;

  ValueType m_init_value;
  FirstFrom m_first_from;
  FirstDest m_first_dest;
  BinaryOpType m_binary_op;
  UnaryOpType m_unary_op;

  KOKKOS_FUNCTION
  TransformExclusiveScanFunctorWithoutValueWrapper(ValueType init,
                                                   FirstFrom first_from,
                                                   FirstDest first_dest,
                                                   BinaryOpType bop,
                                                   UnaryOpType uop)
      : m_init_value(std::move(init)),
        m_first_from(std::move(first_from)),
        m_first_dest(std::move(first_dest)),
        m_binary_op(std::move(bop)),
        m_unary_op(std::move(uop)) {}

  KOKKOS_FUNCTION
  void operator()(const IndexType i, ValueType& update,
                  const bool final_pass) const {
    const auto tmp = ValueType{m_unary_op(m_first_from[i])};
    if (final_pass) {
      if (i == 0) {
        // for both ExclusiveScan and TransformExclusiveScan,
        // init is unmodified
        m_first_dest[i] = m_init_value;
      } else {
        m_first_dest[i] = m_binary_op(update, m_init_value);
      }
    }

    this->join(update, tmp);
  }

  KOKKOS_FUNCTION
  void init(ValueType& update) const { update = {}; }

  KOKKOS_FUNCTION
  void join(ValueType& update, const ValueType& input) const {
    update = m_binary_op(update, input);
  }
};

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
