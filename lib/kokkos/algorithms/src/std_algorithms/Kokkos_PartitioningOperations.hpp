/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_STD_PARTITIONING_OPERATIONS_HPP
#define KOKKOS_STD_PARTITIONING_OPERATIONS_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_BeginEnd.hpp"
#include "Kokkos_Constraints.hpp"
#include "Kokkos_ModifyingOperations.hpp"
#include "Kokkos_NonModifyingSequenceOperations.hpp"

namespace Kokkos {
namespace Experimental {
namespace Impl {

// -------------------------
//
// functors
//
// -------------------------

template <class IteratorType, class ReducerType, class PredicateType>
struct StdIsPartitionedFunctor {
  using red_value_type = typename ReducerType::value_type;
  using index_type     = typename IteratorType::difference_type;

  IteratorType m_first;
  ReducerType m_reducer;
  PredicateType m_p;

  KOKKOS_FUNCTION
  void operator()(const index_type i, red_value_type& redValue) const {
    const auto predicate_value = m_p(m_first[i]);
    constexpr index_type m_red_id_min =
        ::Kokkos::reduction_identity<index_type>::min();
    constexpr index_type m_red_id_max =
        ::Kokkos::reduction_identity<index_type>::max();
    auto rv = predicate_value ? red_value_type{i, m_red_id_min}
                              : red_value_type{m_red_id_max, i};

    m_reducer.join(redValue, rv);
  }

  KOKKOS_FUNCTION
  StdIsPartitionedFunctor(IteratorType first, ReducerType reducer,
                          PredicateType p)
      : m_first(std::move(first)),
        m_reducer(std::move(reducer)),
        m_p(std::move(p)) {}
};

template <class IteratorType, class ReducerType, class PredicateType>
struct StdPartitionPointFunctor {
  using red_value_type = typename ReducerType::value_type;
  using index_type     = typename IteratorType::difference_type;

  IteratorType m_first;
  ReducerType m_reducer;
  PredicateType m_p;

  KOKKOS_FUNCTION
  void operator()(const index_type i, red_value_type& redValue) const {
    const auto predicate_value = m_p(m_first[i]);
    auto rv =
        predicate_value
            ? red_value_type{::Kokkos::reduction_identity<index_type>::min()}
            : red_value_type{i};
    m_reducer.join(redValue, rv);
  }

  KOKKOS_FUNCTION
  StdPartitionPointFunctor(IteratorType first, ReducerType reducer,
                           PredicateType p)
      : m_first(std::move(first)),
        m_reducer(std::move(reducer)),
        m_p(std::move(p)) {}
};

template <class ValueType>
struct StdPartitionCopyScalar {
  ValueType true_count_;
  ValueType false_count_;

  // Here we implement the copy assignment operators explicitly for consistency
  // with how the Scalar structs are implemented inside
  // Kokkos_Parallel_Reduce.hpp.
  KOKKOS_FUNCTION
  void operator=(const StdPartitionCopyScalar& other) {
    true_count_  = other.true_count_;
    false_count_ = other.false_count_;
  }

  KOKKOS_FUNCTION
  void operator=(const volatile StdPartitionCopyScalar& other) volatile {
    true_count_  = other.true_count_;
    false_count_ = other.false_count_;
  }

  // this is needed for
  // OpenMPTarget/Kokkos_OpenMPTarget_Parallel.hpp:699:21: error: no viable
  // overloaded '=' m_returnvalue = 0;
  //
  KOKKOS_FUNCTION
  void operator=(const ValueType value) {
    true_count_  = value;
    false_count_ = value;
  }
};

template <class IndexType, class FirstFrom, class FirstDestTrue,
          class FirstDestFalse, class PredType>
struct StdPartitionCopyFunctor {
  using value_type = StdPartitionCopyScalar<IndexType>;

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
  void operator()(const IndexType i, value_type& update,
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
  void init(value_type& update) const {
    update.true_count_  = 0;
    update.false_count_ = 0;
  }

  KOKKOS_FUNCTION
  void join(volatile value_type& update,
            volatile const value_type& input) const {
    update.true_count_ += input.true_count_;
    update.false_count_ += input.false_count_;
  }
};

// ------------------------------------------
// is_partitioned_impl
// ------------------------------------------

template <class ExecutionSpace, class IteratorType, class PredicateType>
bool is_partitioned_impl(const std::string& label, const ExecutionSpace& ex,
                         IteratorType first, IteratorType last,
                         PredicateType pred) {
  // true if all elements in the range [first, last) that satisfy
  // the predicate "pred" appear before all elements that don't.
  // Also returns true if [first, last) is empty.
  // also true if all elements satisfy the predicate.

  // we implement it by finding:
  // - the max location where predicate is true  (max_loc_true)
  // - the min location where predicate is false (min_loc_false)
  // so the range is partitioned if max_loc_true < (min_loc_false)

  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  // trivial case
  if (first == last) {
    return true;
  }

  // aliases
  using index_type           = typename IteratorType::difference_type;
  using reducer_type         = StdIsPartitioned<index_type>;
  using reduction_value_type = typename reducer_type::value_type;
  using func_t =
      StdIsPartitionedFunctor<IteratorType, reducer_type, PredicateType>;

  // run
  reduction_value_type red_result;
  reducer_type reducer(red_result);
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_reduce(label,
                            RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                            func_t(first, reducer, pred), reducer);

  // fence not needed because reducing into scalar

  // decide and return
  constexpr index_type red_id_min =
      ::Kokkos::reduction_identity<index_type>::min();
  constexpr index_type red_id_max =
      ::Kokkos::reduction_identity<index_type>::max();

  if (red_result.max_loc_true != red_id_max &&
      red_result.min_loc_false != red_id_min) {
    return red_result.max_loc_true < red_result.min_loc_false;
  } else if (first + red_result.max_loc_true == --last) {
    return true;
  } else {
    return false;
  }
}

// ------------------------------------------
// partition_point_impl
// ------------------------------------------
template <class ExecutionSpace, class IteratorType, class PredicateType>
IteratorType partition_point_impl(const std::string& label,
                                  const ExecutionSpace& ex, IteratorType first,
                                  IteratorType last, PredicateType pred) {
  // locates the end of the first partition, that is, the first
  // element that does not satisfy p or last if all elements satisfy p.
  // Implementation below finds the first location where p is false.

  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  if (first == last) {
    return first;
  }

  // aliases
  using index_type           = typename IteratorType::difference_type;
  using reducer_type         = StdPartitionPoint<index_type>;
  using reduction_value_type = typename reducer_type::value_type;
  using func_t =
      StdPartitionPointFunctor<IteratorType, reducer_type, PredicateType>;

  // run
  reduction_value_type red_result;
  reducer_type reducer(red_result);
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_reduce(label,
                            RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                            func_t(first, reducer, pred), reducer);

  // fence not needed because reducing into scalar

  // decide and return
  if (red_result.min_loc_false ==
      ::Kokkos::reduction_identity<index_type>::min()) {
    // if all elements are true, return last
    return last;
  } else {
    return first + red_result.min_loc_false;
  }
}

// ------------------------------------------
// partition_copy_impl
// ------------------------------------------
template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorTrueType, class OutputIteratorFalseType,
          class PredicateType>
::Kokkos::pair<OutputIteratorTrueType, OutputIteratorFalseType>
partition_copy_impl(const std::string& label, const ExecutionSpace& ex,
                    InputIteratorType from_first, InputIteratorType from_last,
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

  // aliases
  using index_type = typename InputIteratorType::difference_type;
  using func_type =
      StdPartitionCopyFunctor<index_type, InputIteratorType,
                              OutputIteratorTrueType, OutputIteratorFalseType,
                              PredicateType>;

  // run
  const auto num_elements =
      Kokkos::Experimental::distance(from_first, from_last);
  typename func_type::value_type counts{0, 0};
  ::Kokkos::parallel_scan(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements),
      func_type(from_first, to_first_true, to_first_false, pred), counts);

  // fence not needed here because of the scan into counts

  return {to_first_true + counts.true_count_,
          to_first_false + counts.false_count_};
}

}  // end namespace Impl

// ----------------------
// is_partitioned public API
// ----------------------
template <class ExecutionSpace, class IteratorType, class PredicateType>
bool is_partitioned(const ExecutionSpace& ex, IteratorType first,
                    IteratorType last, PredicateType p) {
  return Impl::is_partitioned_impl(
      "Kokkos::is_partitioned_iterator_api_default", ex, first, last,
      std::move(p));
}

template <class ExecutionSpace, class IteratorType, class PredicateType>
bool is_partitioned(const std::string& label, const ExecutionSpace& ex,
                    IteratorType first, IteratorType last, PredicateType p) {
  return Impl::is_partitioned_impl(label, ex, first, last, std::move(p));
}

template <class ExecutionSpace, class PredicateType, class DataType,
          class... Properties>
bool is_partitioned(const ExecutionSpace& ex,
                    const ::Kokkos::View<DataType, Properties...>& v,
                    PredicateType p) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  return Impl::is_partitioned_impl("Kokkos::is_partitioned_view_api_default",
                                   ex, cbegin(v), cend(v), std::move(p));
}

template <class ExecutionSpace, class PredicateType, class DataType,
          class... Properties>
bool is_partitioned(const std::string& label, const ExecutionSpace& ex,
                    const ::Kokkos::View<DataType, Properties...>& v,
                    PredicateType p) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  return Impl::is_partitioned_impl(label, ex, cbegin(v), cend(v), std::move(p));
}

// ----------------------
// partition_copy
// ----------------------
template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorTrueType, class OutputIteratorFalseType,
          class PredicateType>
::Kokkos::pair<OutputIteratorTrueType, OutputIteratorFalseType> partition_copy(
    const ExecutionSpace& ex, InputIteratorType from_first,
    InputIteratorType from_last, OutputIteratorTrueType to_first_true,
    OutputIteratorFalseType to_first_false, PredicateType p) {
  return Impl::partition_copy_impl(
      "Kokkos::partition_copy_iterator_api_default", ex, from_first, from_last,
      to_first_true, to_first_false, std::move(p));
}

template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorTrueType, class OutputIteratorFalseType,
          class PredicateType>
::Kokkos::pair<OutputIteratorTrueType, OutputIteratorFalseType> partition_copy(
    const std::string& label, const ExecutionSpace& ex,
    InputIteratorType from_first, InputIteratorType from_last,
    OutputIteratorTrueType to_first_true,
    OutputIteratorFalseType to_first_false, PredicateType p) {
  return Impl::partition_copy_impl(label, ex, from_first, from_last,
                                   to_first_true, to_first_false, std::move(p));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class DataType3,
          class... Properties3, class PredicateType>
auto partition_copy(
    const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest_true,
    const ::Kokkos::View<DataType3, Properties3...>& view_dest_false,
    PredicateType p) {
  return Impl::partition_copy_impl("Kokkos::partition_copy_view_api_default",
                                   ex, cbegin(view_from), cend(view_from),
                                   begin(view_dest_true),
                                   begin(view_dest_false), std::move(p));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class DataType3,
          class... Properties3, class PredicateType>
auto partition_copy(
    const std::string& label, const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest_true,
    const ::Kokkos::View<DataType3, Properties3...>& view_dest_false,
    PredicateType p) {
  return Impl::partition_copy_impl(label, ex, cbegin(view_from),
                                   cend(view_from), begin(view_dest_true),
                                   begin(view_dest_false), std::move(p));
}

// ----------------------
// partition_point
// ----------------------
template <class ExecutionSpace, class IteratorType, class UnaryPredicate>
IteratorType partition_point(const ExecutionSpace& ex, IteratorType first,
                             IteratorType last, UnaryPredicate p) {
  return Impl::partition_point_impl(
      "Kokkos::partitioned_point_iterator_api_default", ex, first, last,
      std::move(p));
}

template <class ExecutionSpace, class IteratorType, class UnaryPredicate>
IteratorType partition_point(const std::string& label, const ExecutionSpace& ex,
                             IteratorType first, IteratorType last,
                             UnaryPredicate p) {
  return Impl::partition_point_impl(label, ex, first, last, std::move(p));
}

template <class ExecutionSpace, class UnaryPredicate, class DataType,
          class... Properties>
auto partition_point(const std::string& label, const ExecutionSpace& ex,
                     const ::Kokkos::View<DataType, Properties...>& v,
                     UnaryPredicate p) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);
  return Impl::partition_point_impl(label, ex, begin(v), end(v), std::move(p));
}

template <class ExecutionSpace, class UnaryPredicate, class DataType,
          class... Properties>
auto partition_point(const ExecutionSpace& ex,
                     const ::Kokkos::View<DataType, Properties...>& v,
                     UnaryPredicate p) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);
  return Impl::partition_point_impl("Kokkos::partition_point_view_api_default",
                                    ex, begin(v), end(v), std::move(p));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
