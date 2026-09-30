// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_STD_ALGORITHMS_PARTITION_COPY_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_PARTITION_COPY_IMPL_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include "Kokkos_MustUseKokkosSingleInTeam.hpp"
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class ValueType>
struct StdPartitionCopyScalar {
  ValueType true_count_  = {};
  ValueType false_count_ = {};

  KOKKOS_DEFAULTED_FUNCTION
  StdPartitionCopyScalar() = default;

  KOKKOS_FUNCTION StdPartitionCopyScalar& operator=(ValueType o) {
    true_count_  = o;
    false_count_ = o;
    return *this;
  }

  KOKKOS_FUNCTION StdPartitionCopyScalar(ValueType true_count,
                                         ValueType false_count)
      : true_count_(true_count), false_count_(false_count) {}

  // Non-explicit: team_scan scratch uses `type accum = 0` for generic ArgType.
  KOKKOS_FUNCTION StdPartitionCopyScalar(ValueType o)
      : StdPartitionCopyScalar(o, o) {}

  // Threads team_scan returns through volatile scratch; copy from volatile.
  KOKKOS_FUNCTION StdPartitionCopyScalar(
      volatile StdPartitionCopyScalar const& o)
      : true_count_(o.true_count_), false_count_(o.false_count_) {}

  KOKKOS_DEFAULTED_FUNCTION
  StdPartitionCopyScalar(StdPartitionCopyScalar const&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  StdPartitionCopyScalar& operator=(StdPartitionCopyScalar const&) = default;

  // Threads team_scan writes through volatile scratch. Return void so GCC
  // -Wvolatile does not warn on *volatile_ptr = rhs (discarded volatile ref).
  KOKKOS_FUNCTION
  void operator=(StdPartitionCopyScalar const& o) volatile {
    true_count_  = o.true_count_;
    false_count_ = o.false_count_;
  }

  KOKKOS_DEFAULTED_FUNCTION
  StdPartitionCopyScalar(StdPartitionCopyScalar&&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  StdPartitionCopyScalar& operator=(StdPartitionCopyScalar&&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  ~StdPartitionCopyScalar() = default;

  KOKKOS_FUNCTION
  StdPartitionCopyScalar& operator+=(StdPartitionCopyScalar const& o) {
    true_count_ += o.true_count_;
    false_count_ += o.false_count_;
    return *this;
  }

  KOKKOS_FUNCTION
  StdPartitionCopyScalar operator+(StdPartitionCopyScalar const& o) const {
    StdPartitionCopyScalar res;
    res.true_count_  = true_count_ + o.true_count_;
    res.false_count_ = false_count_ + o.false_count_;
    return res;
  }
};

template <class FirstFrom, class FirstDestTrue, class FirstDestFalse,
          class PredType>
struct StdPartitionCopyFunctor {
  using index_type = typename FirstFrom::difference_type;
  using value_type = StdPartitionCopyScalar<index_type>;

  FirstFrom m_first_from;
  FirstDestTrue m_first_dest_true;
  FirstDestFalse m_first_dest_false;
  PredType m_pred;

  KOKKOS_FUNCTION
  StdPartitionCopyFunctor(FirstFrom first_from, FirstDestTrue first_dest_true,
                          FirstDestFalse first_dest_false, PredType pred)
      : m_first_from(std::move(first_from)),
        m_first_dest_true(std::move(first_dest_true)),
        m_first_dest_false(std::move(first_dest_false)),
        m_pred(std::move(pred)) {}

  KOKKOS_FUNCTION
  void operator()(const index_type i, value_type& update,
                  const bool final_pass) const {
    const auto& myval = m_first_from[i];
    if (final_pass) {
      if (m_pred(myval)) {
        m_first_dest_true[update.true_count_] = myval;
      } else {
        m_first_dest_false[update.false_count_] = myval;
      }
    }

    if (m_pred(myval)) {
      update.true_count_ += 1;
    } else {
      update.false_count_ += 1;
    }
  }

  KOKKOS_FUNCTION
  void init(value_type& update) const { update = value_type{}; }

  KOKKOS_FUNCTION
  void join(value_type& update, const value_type& input) const {
    update.true_count_ += input.true_count_;
    update.false_count_ += input.false_count_;
  }
};

template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorTrueType, class OutputIteratorFalseType,
          class PredicateType>
::Kokkos::pair<OutputIteratorTrueType, OutputIteratorFalseType>
partition_copy_exespace_impl(const std::string& label, const ExecutionSpace& ex,
                             InputIteratorType from_first,
                             InputIteratorType from_last,
                             OutputIteratorTrueType to_first_true,
                             OutputIteratorFalseType to_first_false,
                             PredicateType pred) {
  // impl uses a scan, this is similar how we implemented copy_if

  // checks
  Impl::static_assert_random_access_and_accessible(
      ex, from_first, to_first_true, to_first_false);
  Impl::static_assert_iterators_have_matching_difference_type(
      from_first, to_first_true, to_first_false);
  Impl::expect_valid_range(from_first, from_last);

  if (from_first == from_last) {
    return {to_first_true, to_first_false};
  }

  using func_type =
      StdPartitionCopyFunctor<InputIteratorType, OutputIteratorTrueType,
                              OutputIteratorFalseType, PredicateType>;

  // run
  const auto num_elements =
      Kokkos::Experimental::distance(from_first, from_last);
  typename func_type::value_type counts;
  ::Kokkos::parallel_scan(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements),
      func_type(from_first, to_first_true, to_first_false, pred), counts);

  // fence not needed here because of the scan into counts

  return {to_first_true + counts.true_count_,
          to_first_false + counts.false_count_};
}

template <class TeamHandleType, class InputIteratorType,
          class OutputIteratorTrueType, class OutputIteratorFalseType,
          class PredicateType>
KOKKOS_FUNCTION ::Kokkos::pair<OutputIteratorTrueType, OutputIteratorFalseType>
partition_copy_team_impl(const TeamHandleType& teamHandle,
                         InputIteratorType from_first,
                         InputIteratorType from_last,
                         OutputIteratorTrueType to_first_true,
                         OutputIteratorFalseType to_first_false,
                         PredicateType pred) {
  // impl uses a scan, this is similar how we implemented copy_if

  // checks
  Impl::static_assert_random_access_and_accessible(
      teamHandle, from_first, to_first_true, to_first_false);
  Impl::static_assert_iterators_have_matching_difference_type(
      from_first, to_first_true, to_first_false);
  Impl::expect_valid_range(from_first, from_last);

  if (from_first == from_last) {
    return {to_first_true, to_first_false};
  }

  const std::size_t num_elements =
      Kokkos::Experimental::distance(from_first, from_last);

  if constexpr (stdalgo_must_use_kokkos_single_for_team_scan_v<
                    typename TeamHandleType::execution_space>) {
    using counts_t  = ::Kokkos::pair<std::size_t, std::size_t>;
    counts_t counts = {};
    Kokkos::single(
        Kokkos::PerTeam(teamHandle),
        [=](counts_t& lcounts) {
          lcounts = {};
          for (std::size_t i = 0; i < num_elements; ++i) {
            const auto& myval = from_first[i];
            if (pred(myval)) {
              to_first_true[lcounts.first++] = myval;
            } else {
              to_first_false[lcounts.second++] = myval;
            }
          }
        },
        counts);
    // no barrier needed since single above broadcasts to all members

    return {to_first_true + counts.first, to_first_false + counts.second};

  } else {
    using func_type =
        StdPartitionCopyFunctor<InputIteratorType, OutputIteratorTrueType,
                                OutputIteratorFalseType, PredicateType>;

    typename func_type::value_type counts;
    ::Kokkos::parallel_scan(
        TeamThreadRange(teamHandle, 0, num_elements),
        func_type(from_first, to_first_true, to_first_false, pred), counts);
    // no barrier needed since reducing into counts

    return {to_first_true + counts.true_count_,
            to_first_false + counts.false_count_};
  }
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
