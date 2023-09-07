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

#ifndef KOKKOS_STD_ALGORITHMS_REPLACE_IF_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_REPLACE_IF_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class InputIterator, class PredicateType, class NewValueType>
struct StdReplaceIfFunctor {
  using index_type = typename InputIterator::difference_type;

  InputIterator m_first;
  PredicateType m_predicate;
  NewValueType m_new_value;

  KOKKOS_FUNCTION
  void operator()(index_type i) const {
    if (m_predicate(m_first[i])) {
      m_first[i] = m_new_value;
    }
  }

  KOKKOS_FUNCTION
  StdReplaceIfFunctor(InputIterator first, PredicateType pred,
                      NewValueType new_value)
      : m_first(std::move(first)),
        m_predicate(std::move(pred)),
        m_new_value(std::move(new_value)) {}
};

template <class ExecutionSpace, class IteratorType, class PredicateType,
          class ValueType>
void replace_if_impl(const std::string& label, const ExecutionSpace& ex,
                     IteratorType first, IteratorType last, PredicateType pred,
                     const ValueType& new_value) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  // aliases
  using func_t = StdReplaceIfFunctor<IteratorType, PredicateType, ValueType>;

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                         func_t(first, std::move(pred), new_value));
  ex.fence("Kokkos::replace_if: fence after operation");
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
