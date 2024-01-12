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

#ifndef KOKKOS_STD_ALGORITHMS_EQUAL_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_EQUAL_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IteratorType1, class IteratorType2, class BinaryPredicateType>
struct StdEqualFunctor {
  using index_type = typename IteratorType1::difference_type;

  IteratorType1 m_first1;
  IteratorType2 m_first2;
  BinaryPredicateType m_predicate;

  KOKKOS_FUNCTION
  void operator()(index_type i, std::size_t& lsum) const {
    if (!m_predicate(m_first1[i], m_first2[i])) {
      lsum = 1;
    }
  }

  KOKKOS_FUNCTION
  StdEqualFunctor(IteratorType1 _first1, IteratorType2 _first2,
                  BinaryPredicateType _predicate)
      : m_first1(std::move(_first1)),
        m_first2(std::move(_first2)),
        m_predicate(std::move(_predicate)) {}
};

//
// exespace impl
//
template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
bool equal_exespace_impl(const std::string& label, const ExecutionSpace& ex,
                         IteratorType1 first1, IteratorType1 last1,
                         IteratorType2 first2, BinaryPredicateType predicate) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first1, first2);
  Impl::static_assert_iterators_have_matching_difference_type(first1, first2);
  Impl::expect_valid_range(first1, last1);

  // run
  const auto num_elements = Kokkos::Experimental::distance(first1, last1);
  std::size_t different   = 0;
  ::Kokkos::parallel_reduce(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements),
      StdEqualFunctor(first1, first2, predicate), different);
  ex.fence("Kokkos::equal: fence after operation");

  return !different;
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
bool equal_exespace_impl(const std::string& label, const ExecutionSpace& ex,
                         IteratorType1 first1, IteratorType1 last1,
                         IteratorType2 first2) {
  using value_type1 = typename IteratorType1::value_type;
  using value_type2 = typename IteratorType2::value_type;
  using pred_t      = StdAlgoEqualBinaryPredicate<value_type1, value_type2>;
  return equal_exespace_impl(label, ex, first1, last1, first2, pred_t());
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
bool equal_exespace_impl(const std::string& label, const ExecutionSpace& ex,
                         IteratorType1 first1, IteratorType1 last1,
                         IteratorType2 first2, IteratorType2 last2,
                         BinaryPredicateType predicate) {
  const auto d1 = ::Kokkos::Experimental::distance(first1, last1);
  const auto d2 = ::Kokkos::Experimental::distance(first2, last2);
  if (d1 != d2) {
    return false;
  }

  return equal_exespace_impl(label, ex, first1, last1, first2, predicate);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
bool equal_exespace_impl(const std::string& label, const ExecutionSpace& ex,
                         IteratorType1 first1, IteratorType1 last1,
                         IteratorType2 first2, IteratorType2 last2) {
  Impl::expect_valid_range(first1, last1);
  Impl::expect_valid_range(first2, last2);

  using value_type1 = typename IteratorType1::value_type;
  using value_type2 = typename IteratorType2::value_type;
  using pred_t      = StdAlgoEqualBinaryPredicate<value_type1, value_type2>;
  return equal_exespace_impl(label, ex, first1, last1, first2, last2, pred_t());
}

//
// team impl
//
template <class TeamHandleType, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
KOKKOS_FUNCTION bool equal_team_impl(const TeamHandleType& teamHandle,
                                     IteratorType1 first1, IteratorType1 last1,
                                     IteratorType2 first2,
                                     BinaryPredicateType predicate) {
  // checks
  Impl::static_assert_random_access_and_accessible(teamHandle, first1, first2);
  Impl::static_assert_iterators_have_matching_difference_type(first1, first2);
  Impl::expect_valid_range(first1, last1);

  // run
  const auto num_elements = Kokkos::Experimental::distance(first1, last1);
  std::size_t different   = 0;
  ::Kokkos::parallel_reduce(TeamThreadRange(teamHandle, 0, num_elements),
                            StdEqualFunctor(first1, first2, predicate),
                            different);
  teamHandle.team_barrier();

  return !different;
}

template <class TeamHandleType, class IteratorType1, class IteratorType2>
KOKKOS_FUNCTION bool equal_team_impl(const TeamHandleType& teamHandle,
                                     IteratorType1 first1, IteratorType1 last1,
                                     IteratorType2 first2) {
  using value_type1 = typename IteratorType1::value_type;
  using value_type2 = typename IteratorType2::value_type;
  using pred_t      = StdAlgoEqualBinaryPredicate<value_type1, value_type2>;
  return equal_team_impl(teamHandle, first1, last1, first2, pred_t());
}

template <class TeamHandleType, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
KOKKOS_FUNCTION bool equal_team_impl(const TeamHandleType& teamHandle,
                                     IteratorType1 first1, IteratorType1 last1,
                                     IteratorType2 first2, IteratorType2 last2,
                                     BinaryPredicateType predicate) {
  const auto d1 = ::Kokkos::Experimental::distance(first1, last1);
  const auto d2 = ::Kokkos::Experimental::distance(first2, last2);
  if (d1 != d2) {
    return false;
  }

  return equal_team_impl(teamHandle, first1, last1, first2, predicate);
}

template <class TeamHandleType, class IteratorType1, class IteratorType2>
KOKKOS_FUNCTION bool equal_team_impl(const TeamHandleType& teamHandle,
                                     IteratorType1 first1, IteratorType1 last1,
                                     IteratorType2 first2,
                                     IteratorType2 last2) {
  Impl::expect_valid_range(first1, last1);
  Impl::expect_valid_range(first2, last2);

  using value_type1 = typename IteratorType1::value_type;
  using value_type2 = typename IteratorType2::value_type;
  using pred_t      = StdAlgoEqualBinaryPredicate<value_type1, value_type2>;
  return equal_team_impl(teamHandle, first1, last1, first2, last2, pred_t());
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
