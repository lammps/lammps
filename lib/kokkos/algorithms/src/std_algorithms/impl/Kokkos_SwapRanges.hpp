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

#ifndef KOKKOS_STD_ALGORITHMS_SWAP_RANGES_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_SWAP_RANGES_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <std_algorithms/Kokkos_Swap.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IndexType, class IteratorType1, class IteratorType2>
struct StdSwapRangesFunctor {
  IteratorType1 m_first1;
  IteratorType2 m_first2;

  KOKKOS_FUNCTION
  void operator()(IndexType i) const {
    ::Kokkos::Experimental::swap(m_first1[i], m_first2[i]);
  }

  KOKKOS_FUNCTION
  StdSwapRangesFunctor(IteratorType1 _first1, IteratorType2 _first2)
      : m_first1(std::move(_first1)), m_first2(std::move(_first2)) {}
};

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
IteratorType2 swap_ranges_impl(const std::string& label,
                               const ExecutionSpace& ex, IteratorType1 first1,
                               IteratorType1 last1, IteratorType2 first2) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first1, first2);
  Impl::static_assert_iterators_have_matching_difference_type(first1, first2);
  Impl::expect_valid_range(first1, last1);

  // aliases
  using index_type = typename IteratorType1::difference_type;
  using func_t = StdSwapRangesFunctor<index_type, IteratorType1, IteratorType2>;

  // run
  const auto num_elements_to_swap =
      Kokkos::Experimental::distance(first1, last1);
  ::Kokkos::parallel_for(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements_to_swap),
      func_t(first1, first2));
  ex.fence("Kokkos::swap_ranges: fence after operation");

  // return
  return first2 + num_elements_to_swap;
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
