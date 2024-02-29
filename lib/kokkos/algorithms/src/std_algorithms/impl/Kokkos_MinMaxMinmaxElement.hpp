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

#ifndef KOKKOS_STD_ALGORITHMS_MIN_MAX_MINMAX_ELEMENT_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_MIN_MAX_MINMAX_ELEMENT_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IteratorType, class ReducerType>
struct StdMinOrMaxElemFunctor {
  using index_type     = typename IteratorType::difference_type;
  using red_value_type = typename ReducerType::value_type;

  IteratorType m_first;
  ReducerType m_reducer;

  KOKKOS_FUNCTION
  void operator()(const index_type i, red_value_type& red_value) const {
    m_reducer.join(red_value, red_value_type{m_first[i], i});
  }

  KOKKOS_FUNCTION
  StdMinOrMaxElemFunctor(IteratorType first, ReducerType reducer)
      : m_first(std::move(first)), m_reducer(std::move(reducer)) {}
};

template <class IteratorType, class ReducerType>
struct StdMinMaxElemFunctor {
  using index_type     = typename IteratorType::difference_type;
  using red_value_type = typename ReducerType::value_type;
  IteratorType m_first;
  ReducerType m_reducer;

  KOKKOS_FUNCTION
  void operator()(const index_type i, red_value_type& red_value) const {
    const auto& my_value = m_first[i];
    m_reducer.join(red_value, red_value_type{my_value, my_value, i, i});
  }

  KOKKOS_FUNCTION
  StdMinMaxElemFunctor(IteratorType first, ReducerType reducer)
      : m_first(std::move(first)), m_reducer(std::move(reducer)) {}
};

//
// exespace impl
//
template <template <class... Args> class ReducerType, class ExecutionSpace,
          class IteratorType, class... Args>
IteratorType min_or_max_element_exespace_impl(const std::string& label,
                                              const ExecutionSpace& ex,
                                              IteratorType first,
                                              IteratorType last,
                                              Args&&... args) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  if (first == last) {
    return last;
  }

  // aliases
  using index_type           = typename IteratorType::difference_type;
  using value_type           = typename IteratorType::value_type;
  using reducer_type         = ReducerType<value_type, index_type, Args...>;
  using reduction_value_type = typename reducer_type::value_type;
  using func_t = StdMinOrMaxElemFunctor<IteratorType, reducer_type>;

  // run
  reduction_value_type red_result;
  reducer_type reducer(red_result, std::forward<Args>(args)...);
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_reduce(label,
                            RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                            func_t(first, reducer), reducer);

  // fence not needed because reducing into scalar

  // return
  return first + red_result.loc;
}

template <template <class... Args> class ReducerType, class ExecutionSpace,
          class IteratorType, class... Args>
::Kokkos::pair<IteratorType, IteratorType> minmax_element_exespace_impl(
    const std::string& label, const ExecutionSpace& ex, IteratorType first,
    IteratorType last, Args&&... args) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  if (first == last) {
    return {first, first};
  }

  // aliases
  using index_type           = typename IteratorType::difference_type;
  using value_type           = typename IteratorType::value_type;
  using reducer_type         = ReducerType<value_type, index_type, Args...>;
  using reduction_value_type = typename reducer_type::value_type;
  using func_t               = StdMinMaxElemFunctor<IteratorType, reducer_type>;

  // run
  reduction_value_type red_result;
  reducer_type reducer(red_result, std::forward<Args>(args)...);
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_reduce(label,
                            RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                            func_t(first, reducer), reducer);

  // fence not needed because reducing into scalar

  // return
  return {first + red_result.min_loc, first + red_result.max_loc};
}

//
// team level impl
//
template <template <class... Args> class ReducerType, class TeamHandleType,
          class IteratorType, class... Args>
KOKKOS_FUNCTION IteratorType min_or_max_element_team_impl(
    const TeamHandleType& teamHandle, IteratorType first, IteratorType last,
    Args&&... args) {
  // checks
  Impl::static_assert_random_access_and_accessible(teamHandle, first);
  Impl::expect_valid_range(first, last);

  if (first == last) {
    return last;
  }

  // aliases
  using index_type           = typename IteratorType::difference_type;
  using value_type           = typename IteratorType::value_type;
  using reducer_type         = ReducerType<value_type, index_type, Args...>;
  using reduction_value_type = typename reducer_type::value_type;
  using func_t = StdMinOrMaxElemFunctor<IteratorType, reducer_type>;

  // run
  reduction_value_type red_result;
  reducer_type reducer(red_result, std::forward<Args>(args)...);
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_reduce(TeamThreadRange(teamHandle, 0, num_elements),
                            func_t(first, reducer), reducer);
  teamHandle.team_barrier();
  // maybe the barrier is not needed since reducing into scalar?

  // return
  return first + red_result.loc;
}

template <template <class... Args> class ReducerType, class TeamHandleType,
          class IteratorType, class... Args>
KOKKOS_FUNCTION ::Kokkos::pair<IteratorType, IteratorType>
minmax_element_team_impl(const TeamHandleType& teamHandle, IteratorType first,
                         IteratorType last, Args&&... args) {
  // checks
  Impl::static_assert_random_access_and_accessible(teamHandle, first);
  Impl::expect_valid_range(first, last);

  if (first == last) {
    return {first, first};
  }

  // aliases
  using index_type           = typename IteratorType::difference_type;
  using value_type           = typename IteratorType::value_type;
  using reducer_type         = ReducerType<value_type, index_type, Args...>;
  using reduction_value_type = typename reducer_type::value_type;
  using func_t               = StdMinMaxElemFunctor<IteratorType, reducer_type>;

  // run
  reduction_value_type red_result;
  reducer_type reducer(red_result, std::forward<Args>(args)...);
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_reduce(TeamThreadRange(teamHandle, 0, num_elements),
                            func_t(first, reducer), reducer);
  teamHandle.team_barrier();
  // maybe the barrier is not needed since reducing into scalar?

  // return
  return {first + red_result.min_loc, first + red_result.max_loc};
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
