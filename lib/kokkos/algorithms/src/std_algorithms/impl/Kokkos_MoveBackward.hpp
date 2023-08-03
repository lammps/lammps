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

#ifndef KOKKOS_STD_ALGORITHMS_MOVE_BACKWARD_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_MOVE_BACKWARD_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IndexType, class IteratorType1, class IteratorType2>
struct StdMoveBackwardFunctor {
  static_assert(std::is_signed<IndexType>::value,
                "Kokkos: StdMoveBackwardFunctor requires signed index type");

  IteratorType1 m_last;
  IteratorType2 m_dest_last;

  KOKKOS_FUNCTION
  void operator()(IndexType i) const {
    m_dest_last[-i - 1] = std::move(m_last[-i - 1]);
  }

  StdMoveBackwardFunctor(IteratorType1 _last, IteratorType2 _dest_last)
      : m_last(std::move(_last)), m_dest_last(std::move(_dest_last)) {}
};

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
IteratorType2 move_backward_impl(const std::string& label,
                                 const ExecutionSpace& ex, IteratorType1 first,
                                 IteratorType1 last, IteratorType2 d_last) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first, d_last);
  Impl::static_assert_iterators_have_matching_difference_type(first, d_last);
  Impl::expect_valid_range(first, last);

  // aliases
  using index_type = typename IteratorType1::difference_type;
  using func_t =
      StdMoveBackwardFunctor<index_type, IteratorType1, IteratorType2>;

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                         func_t(last, d_last));
  ex.fence("Kokkos::move_backward: fence after operation");

  // return
  return d_last - num_elements;
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
