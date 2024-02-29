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

#ifndef KOKKOS_STD_ALGORITHMS_UNIQUE_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_UNIQUE_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Move.hpp>
#include <std_algorithms/Kokkos_Distance.hpp>
#include <std_algorithms/Kokkos_AdjacentFind.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class InputIt, class OutputIt, class BinaryPredicateType>
struct StdUniqueFunctor {
  using index_type = typename InputIt::difference_type;

  InputIt m_first_from;
  InputIt m_last_from;
  OutputIt m_first_dest;
  BinaryPredicateType m_pred;

  KOKKOS_FUNCTION
  StdUniqueFunctor(InputIt first_from, InputIt last_from, OutputIt first_dest,
                   BinaryPredicateType pred)
      : m_first_from(std::move(first_from)),
        m_last_from(std::move(last_from)),
        m_first_dest(std::move(first_dest)),
        m_pred(std::move(pred)) {}

  KOKKOS_FUNCTION
  void operator()(const index_type i, index_type& update,
                  const bool final_pass) const {
    auto& val_i         = m_first_from[i];
    const auto& val_ip1 = m_first_from[i + 1];

    if (final_pass) {
      if (!m_pred(val_i, val_ip1)) {
        m_first_dest[update] = std::move(val_i);
      }
    }

    if (!m_pred(val_i, val_ip1)) {
      update += 1;
    }
  }
};

template <class ExecutionSpace, class IteratorType, class PredicateType>
IteratorType unique_exespace_impl(const std::string& label,
                                  const ExecutionSpace& ex, IteratorType first,
                                  IteratorType last, PredicateType pred) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  const auto num_elements = Kokkos::Experimental::distance(first, last);
  if (num_elements == 0) {
    return first;
  } else if (num_elements == 1) {
    return last;
  } else {
    // ----------
    // step 1:
    // find first location of adjacent equal elements
    // ----------
    auto it_found =
        ::Kokkos::Experimental::adjacent_find(ex, first, last, pred);

    // if none, all elements are unique, so nothing to do
    if (it_found == last) {
      return last;
    } else {
      // if here, we found some equal adjacent elements,
      // so count all preceeding unique elements
      const auto num_unique_found_in_step_one = it_found - first;

      // ----------
      // step 2:
      // ----------
      // since we found some unique elements, we don't need to explore
      // the full range [first, last), but only need to focus on the
      // remaining range [it_found, last)
      const auto num_elements_to_explore = last - it_found;

      // create a tmp view to use to *move* all unique elements
      // using the same algorithm used for unique_copy but we now move things
      using value_type    = typename IteratorType::value_type;
      using tmp_view_type = Kokkos::View<value_type*, ExecutionSpace>;
      tmp_view_type tmp_view("std_unique_tmp_view", num_elements_to_explore);

      // scan extent is: num_elements_to_explore - 1
      // for same reason as the one explained in unique_copy
      const auto scan_size = num_elements_to_explore - 1;
      auto tmp_first       = ::Kokkos::Experimental::begin(tmp_view);

      using index_type = typename IteratorType::difference_type;
      index_type count = 0;
      ::Kokkos::parallel_scan(
          label, RangePolicy<ExecutionSpace>(ex, 0, scan_size),
          StdUniqueFunctor(it_found, last, tmp_first, pred), count);

      // move last element too, for the same reason as the unique_copy
      [[maybe_unused]] auto unused_r = Impl::move_exespace_impl(
          "Kokkos::move_from_unique", ex, it_found + scan_size, last,
          tmp_first + count);

      // ----------
      // step 3
      // ----------
      // move back from tmp to original range,
      // ensuring we start overwriting after the original unique found
      using tmp_readwrite_iterator_type = decltype(begin(tmp_view));
      using step3_func_t =
          StdMoveFunctor<index_type, tmp_readwrite_iterator_type, IteratorType>;

      ::Kokkos::parallel_for(
          "unique_step3_parfor",
          RangePolicy<ExecutionSpace>(ex, 0, tmp_view.extent(0)),
          step3_func_t(begin(tmp_view),
                       (first + num_unique_found_in_step_one)));

      ex.fence("Kokkos::unique: fence after operation");

      // return iterator to one passed the last written
      // (the +1 is needed to account for the last element, see above)
      return (first + num_unique_found_in_step_one + count + 1);
    }
  }
}

template <class ExecutionSpace, class IteratorType>
IteratorType unique_exespace_impl(const std::string& label,
                                  const ExecutionSpace& ex, IteratorType first,
                                  IteratorType last) {
  using value_type    = typename IteratorType::value_type;
  using binary_pred_t = StdAlgoEqualBinaryPredicate<value_type>;
  return unique_exespace_impl(label, ex, first, last, binary_pred_t());
}

//
// team level
//
template <class TeamHandleType, class IteratorType, class PredicateType>
KOKKOS_FUNCTION IteratorType unique_team_impl(const TeamHandleType& teamHandle,
                                              IteratorType first,
                                              IteratorType last,
                                              PredicateType pred) {
  // checks
  Impl::static_assert_random_access_and_accessible(teamHandle, first);
  Impl::expect_valid_range(first, last);

  // branch for trivial vs non trivial case
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  if (num_elements == 0) {
    return first;
  } else if (num_elements == 1) {
    return last;
  } else {
    // FIXME: for the execution-space-based impl we used an auxiliary
    // allocation, but for the team level we cannot do the same, so do this
    // serially for now and later figure out if this can be done in parallel

    std::size_t count = 0;
    Kokkos::single(
        Kokkos::PerTeam(teamHandle),
        [=](std::size_t& lcount) {
          IteratorType result = first;
          IteratorType lfirst = first;
          while (++lfirst != last) {
            if (!pred(*result, *lfirst) && ++result != lfirst) {
              *result = std::move(*lfirst);
            }
          }
          lcount = Kokkos::Experimental::distance(first, result);
        },
        count);
    // no barrier needed since single above broadcasts to all members

    // +1 is needed because we want one element past the end
    return first + count + 1;
  }
}

template <class TeamHandleType, class IteratorType>
KOKKOS_FUNCTION IteratorType unique_team_impl(const TeamHandleType& teamHandle,
                                              IteratorType first,
                                              IteratorType last) {
  // checks
  Impl::static_assert_random_access_and_accessible(teamHandle, first);
  Impl::expect_valid_range(first, last);

  using binary_pred_t =
      StdAlgoEqualBinaryPredicate<typename IteratorType::value_type>;
  return unique_team_impl(teamHandle, first, last, binary_pred_t());
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
