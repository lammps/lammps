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

#ifndef KOKKOS_STD_ALGORITHMS_REMOVE_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_REMOVE_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <std_algorithms/Kokkos_CountIf.hpp>
#include <std_algorithms/Kokkos_CopyIf.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IndexType, class FirstFrom, class FirstDest, class PredType>
struct StdRemoveIfStage1Functor {
  FirstFrom m_first_from;
  FirstDest m_first_dest;
  PredType m_must_remove;

  KOKKOS_FUNCTION
  StdRemoveIfStage1Functor(FirstFrom first_from, FirstDest first_dest,
                           PredType pred)
      : m_first_from(std::move(first_from)),
        m_first_dest(std::move(first_dest)),
        m_must_remove(std::move(pred)) {}

  KOKKOS_FUNCTION
  void operator()(const IndexType i, IndexType& update,
                  const bool final_pass) const {
    auto& myval = m_first_from[i];
    if (final_pass) {
      if (!m_must_remove(myval)) {
        // calling move here is ok because we are inside final pass
        // we are calling move assign as specified by the std
        m_first_dest[update] = std::move(myval);
      }
    }

    if (!m_must_remove(myval)) {
      update += 1;
    }
  }
};

template <class IndexType, class InputIteratorType, class OutputIteratorType>
struct StdRemoveIfStage2Functor {
  InputIteratorType m_first_from;
  OutputIteratorType m_first_to;

  KOKKOS_FUNCTION
  StdRemoveIfStage2Functor(InputIteratorType first_from,
                           OutputIteratorType first_to)
      : m_first_from(std::move(first_from)), m_first_to(std::move(first_to)) {}

  KOKKOS_FUNCTION
  void operator()(const IndexType i) const {
    m_first_to[i] = std::move(m_first_from[i]);
  }
};

template <class ExecutionSpace, class IteratorType, class UnaryPredicateType>
IteratorType remove_if_impl(const std::string& label, const ExecutionSpace& ex,
                            IteratorType first, IteratorType last,
                            UnaryPredicateType pred) {
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  if (first == last) {
    return last;
  } else {
    // create tmp buffer to use to *move* all elements that we need to keep.
    // note that the tmp buffer is just large enought to store
    // all elements to keep, because ideally we do not need/want one
    // as large as the original range.
    // To allocate the right tmp view, we need a call to count_if.
    // We could just do a "safe" allocation of a buffer as
    // large as (last-first), but I think a call to count_if is more afforable.

    // count how many elements we need to keep
    // note that the elements to remove are those that meet the predicate
    const auto remove_count =
        ::Kokkos::Experimental::count_if(ex, first, last, pred);
    const auto keep_count =
        Kokkos::Experimental::distance(first, last) - remove_count;

    // create helper tmp view
    using value_type    = typename IteratorType::value_type;
    using tmp_view_type = Kokkos::View<value_type*, ExecutionSpace>;
    tmp_view_type tmp_view("std_remove_if_tmp_view", keep_count);
    using tmp_readwrite_iterator_type = decltype(begin(tmp_view));

    // in stage 1, *move* all elements to keep from original range to tmp
    // we use similar impl as copy_if except that we *move* rather than copy
    using index_type = typename IteratorType::difference_type;
    using func1_type = StdRemoveIfStage1Functor<index_type, IteratorType,
                                                tmp_readwrite_iterator_type,
                                                UnaryPredicateType>;

    const auto scan_num_elements = Kokkos::Experimental::distance(first, last);
    index_type scan_count        = 0;
    ::Kokkos::parallel_scan(
        label, RangePolicy<ExecutionSpace>(ex, 0, scan_num_elements),
        func1_type(first, begin(tmp_view), pred), scan_count);

    // scan_count should be equal to keep_count
    assert(scan_count == keep_count);
    (void)scan_count;  // to avoid unused complaints

    // stage 2, we do parfor to move from tmp to original range
    using func2_type =
        StdRemoveIfStage2Functor<index_type, tmp_readwrite_iterator_type,
                                 IteratorType>;
    ::Kokkos::parallel_for(
        "remove_if_stage2_parfor",
        RangePolicy<ExecutionSpace>(ex, 0, tmp_view.extent(0)),
        func2_type(begin(tmp_view), first));
    ex.fence("Kokkos::remove_if: fence after stage2");

    // return
    return first + keep_count;
  }
}

template <class ExecutionSpace, class IteratorType, class ValueType>
auto remove_impl(const std::string& label, const ExecutionSpace& ex,
                 IteratorType first, IteratorType last,
                 const ValueType& value) {
  using predicate_type = StdAlgoEqualsValUnaryPredicate<ValueType>;
  return remove_if_impl(label, ex, first, last, predicate_type(value));
}

template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class ValueType>
auto remove_copy_impl(const std::string& label, const ExecutionSpace& ex,
                      InputIteratorType first_from, InputIteratorType last_from,
                      OutputIteratorType first_dest, const ValueType& value) {
  // this is like copy_if except that we need to *ignore* the elements
  // that match the value, so we can solve this as follows:

  using predicate_type = StdAlgoNotEqualsValUnaryPredicate<ValueType>;
  return ::Kokkos::Experimental::copy_if(label, ex, first_from, last_from,
                                         first_dest, predicate_type(value));
}

template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class UnaryPredicate>
auto remove_copy_if_impl(const std::string& label, const ExecutionSpace& ex,
                         InputIteratorType first_from,
                         InputIteratorType last_from,
                         OutputIteratorType first_dest,
                         const UnaryPredicate& pred) {
  // this is like copy_if except that we need to *ignore* the elements
  // satisfying the pred, so we can solve this as follows:

  using value_type = typename InputIteratorType::value_type;
  using pred_wrapper_type =
      StdAlgoNegateUnaryPredicateWrapper<value_type, UnaryPredicate>;
  return ::Kokkos::Experimental::copy_if(label, ex, first_from, last_from,
                                         first_dest, pred_wrapper_type(pred));
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
