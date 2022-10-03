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

#ifndef KOKKOS_NON_MODIFYING_SEQUENCE_OPERATIONS_HPP
#define KOKKOS_NON_MODIFYING_SEQUENCE_OPERATIONS_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_BeginEnd.hpp"
#include "Kokkos_Constraints.hpp"
#include "Kokkos_ModifyingOperations.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include "Kokkos_RandomAccessIterator.hpp"
#include "Kokkos_Distance.hpp"
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

// ------------------------------------------
//
// functors
//
// ------------------------------------------

template <bool is_find_if, class IndexType, class IteratorType,
          class ReducerType, class PredicateType>
struct StdFindIfOrNotFunctor {
  using red_value_type = typename ReducerType::value_type;

  IteratorType m_first;
  ReducerType m_reducer;
  PredicateType m_p;

  KOKKOS_FUNCTION
  void operator()(const IndexType i, red_value_type& red_value) const {
    const auto& my_value = m_first[i];

    // if doing find_if, look for when predicate is true
    // if doing find_if_not, look for when predicate is false
    const bool found_condition = is_find_if ? m_p(my_value) : !m_p(my_value);

    auto rv =
        found_condition
            ? red_value_type{i}
            : red_value_type{::Kokkos::reduction_identity<IndexType>::min()};

    m_reducer.join(red_value, rv);
  }

  KOKKOS_FUNCTION
  StdFindIfOrNotFunctor(IteratorType first, ReducerType reducer,
                        PredicateType p)
      : m_first(std::move(first)),
        m_reducer(std::move(reducer)),
        m_p(std::move(p)) {}
};

template <class IteratorType, class UnaryFunctorType>
struct StdForEachFunctor {
  using index_type = typename IteratorType::difference_type;
  IteratorType m_first;
  UnaryFunctorType m_functor;

  KOKKOS_FUNCTION
  void operator()(index_type i) const { m_functor(m_first[i]); }

  KOKKOS_FUNCTION
  StdForEachFunctor(IteratorType _first, UnaryFunctorType _functor)
      : m_first(std::move(_first)), m_functor(std::move(_functor)) {}
};

template <class IteratorType, class Predicate>
struct StdCountIfFunctor {
  using index_type = typename IteratorType::difference_type;
  IteratorType m_first;
  Predicate m_predicate;

  KOKKOS_FUNCTION
  void operator()(index_type i, index_type& lsum) const {
    if (m_predicate(m_first[i])) {
      lsum++;
    }
  }

  KOKKOS_FUNCTION
  StdCountIfFunctor(IteratorType _first, Predicate _predicate)
      : m_first(std::move(_first)), m_predicate(std::move(_predicate)) {}
};

template <class IndexType, class IteratorType1, class IteratorType2,
          class ReducerType, class BinaryPredicateType>
struct StdMismatchRedFunctor {
  using red_value_type = typename ReducerType::value_type;

  IteratorType1 m_first1;
  IteratorType2 m_first2;
  ReducerType m_reducer;
  BinaryPredicateType m_predicate;

  KOKKOS_FUNCTION
  void operator()(const IndexType i, red_value_type& red_value) const {
    const auto& my_value1 = m_first1[i];
    const auto& my_value2 = m_first2[i];

    auto rv =
        !m_predicate(my_value1, my_value2)
            ? red_value_type{i}
            : red_value_type{::Kokkos::reduction_identity<IndexType>::min()};

    m_reducer.join(red_value, rv);
  }

  KOKKOS_FUNCTION
  StdMismatchRedFunctor(IteratorType1 first1, IteratorType2 first2,
                        ReducerType reducer, BinaryPredicateType predicate)
      : m_first1(std::move(first1)),
        m_first2(std::move(first2)),
        m_reducer(std::move(reducer)),
        m_predicate(std::move(predicate)) {}
};

template <class IndexType, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
struct StdEqualFunctor {
  IteratorType1 m_first1;
  IteratorType2 m_first2;
  BinaryPredicateType m_predicate;

  KOKKOS_FUNCTION
  void operator()(IndexType i, std::size_t& lsum) const {
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

template <class IndexType, class IteratorType1, class IteratorType2,
          class ReducerType, class ComparatorType>
struct StdLexicographicalCompareFunctor {
  using red_value_type = typename ReducerType::value_type;
  IteratorType1 m_first1;
  IteratorType2 m_first2;
  ReducerType m_reducer;
  ComparatorType m_comparator;

  KOKKOS_FUNCTION
  void operator()(const IndexType i, red_value_type& red_value) const {
    const auto& my_value1 = m_first1[i];
    const auto& my_value2 = m_first2[i];

    bool different = m_comparator(my_value1, my_value2) ||
                     m_comparator(my_value2, my_value1);
    auto rv =
        different
            ? red_value_type{i}
            : red_value_type{::Kokkos::reduction_identity<IndexType>::min()};

    m_reducer.join(red_value, rv);
  }

  KOKKOS_FUNCTION
  StdLexicographicalCompareFunctor(IteratorType1 _first1, IteratorType2 _first2,
                                   ReducerType _reducer, ComparatorType _comp)
      : m_first1(std::move(_first1)),
        m_first2(std::move(_first2)),
        m_reducer(std::move(_reducer)),
        m_comparator(std::move(_comp)) {}
};

template <class IndexType, class IteratorType1, class IteratorType2,
          class ComparatorType>
struct StdCompareFunctor {
  IteratorType1 m_it1;
  IteratorType2 m_it2;
  ComparatorType m_predicate;

  KOKKOS_FUNCTION
  void operator()(IndexType /* i is unused */, int& lsum) const {
    if (m_predicate(*m_it1, *m_it2)) {
      lsum = 1;
    }
  }

  KOKKOS_FUNCTION
  StdCompareFunctor(IteratorType1 _it1, IteratorType2 _it2,
                    ComparatorType _predicate)
      : m_it1(std::move(_it1)),
        m_it2(std::move(_it2)),
        m_predicate(std::move(_predicate)) {}
};

template <class IndexType, class IteratorType, class ReducerType,
          class PredicateType>
struct StdAdjacentFindFunctor {
  using red_value_type = typename ReducerType::value_type;

  IteratorType m_first;
  ReducerType m_reducer;
  PredicateType m_p;

  KOKKOS_FUNCTION
  void operator()(const IndexType i, red_value_type& red_value) const {
    const auto& my_value   = m_first[i];
    const auto& next_value = m_first[i + 1];
    const bool are_equal   = m_p(my_value, next_value);

    auto rv =
        are_equal
            ? red_value_type{i}
            : red_value_type{::Kokkos::reduction_identity<IndexType>::min()};

    m_reducer.join(red_value, rv);
  }

  KOKKOS_FUNCTION
  StdAdjacentFindFunctor(IteratorType first, ReducerType reducer,
                         PredicateType p)
      : m_first(std::move(first)),
        m_reducer(std::move(reducer)),
        m_p(std::move(p)) {}
};

template <class IndexType, class IteratorType1, class IteratorType2,
          class ReducerType, class PredicateType>
struct StdSearchFunctor {
  using red_value_type = typename ReducerType::value_type;

  IteratorType1 m_first;
  IteratorType1 m_last;
  IteratorType2 m_s_first;
  IteratorType2 m_s_last;
  ReducerType m_reducer;
  PredicateType m_p;

  KOKKOS_FUNCTION
  void operator()(const IndexType i, red_value_type& red_value) const {
    namespace KE = ::Kokkos::Experimental;
    auto myit    = m_first + i;
    bool found   = true;

    const auto search_count = KE::distance(m_s_first, m_s_last);
    for (IndexType k = 0; k < search_count; ++k) {
      // note that we add this EXPECT to check if we are in a valid range
      // but I think we can remove this beceause the guarantee we don't go
      // out of bounds is taken care of at the calling site
      // where we launch the par-reduce.
      KOKKOS_EXPECTS((myit + k) < m_last);

      if (!m_p(myit[k], m_s_first[k])) {
        found = false;
        break;
      }
    }

    const auto rv =
        found ? red_value_type{i}
              : red_value_type{::Kokkos::reduction_identity<IndexType>::min()};

    m_reducer.join(red_value, rv);
  }

  KOKKOS_FUNCTION
  StdSearchFunctor(IteratorType1 first, IteratorType1 last,
                   IteratorType2 s_first, IteratorType2 s_last,
                   ReducerType reducer, PredicateType p)
      : m_first(std::move(first)),
        m_last(std::move(last)),
        m_s_first(std::move(s_first)),
        m_s_last(std::move(s_last)),
        m_reducer(std::move(reducer)),
        m_p(std::move(p)) {}
};

template <class IndexType, class IteratorType, class SizeType, class ValueType,
          class ReducerType, class PredicateType>
struct StdSearchNFunctor {
  using red_value_type = typename ReducerType::value_type;

  IteratorType m_first;
  IteratorType m_last;
  SizeType m_count;
  ValueType m_value;
  ReducerType m_reducer;
  PredicateType m_p;

  KOKKOS_FUNCTION
  void operator()(const IndexType i, red_value_type& red_value) const {
    namespace KE = ::Kokkos::Experimental;
    auto myit    = m_first + i;
    bool found   = true;

    for (SizeType k = 0; k < m_count; ++k) {
      // note that we add this EXPECT to check if we are in a valid range
      // but I think we can remove this beceause the guarantee we don't go
      // out of bounds is taken care of at the calling site
      // where we launch the par-reduce.
      KOKKOS_EXPECTS((myit + k) < m_last);

      if (!m_p(myit[k], m_value)) {
        found = false;
        break;
      }
    }

    const auto rv =
        found ? red_value_type{i}
              : red_value_type{::Kokkos::reduction_identity<IndexType>::min()};

    m_reducer.join(red_value, rv);
  }

  KOKKOS_FUNCTION
  StdSearchNFunctor(IteratorType first, IteratorType last, SizeType count,
                    ValueType value, ReducerType reducer, PredicateType p)
      : m_first(std::move(first)),
        m_last(std::move(last)),
        m_count(std::move(count)),
        m_value(std::move(value)),
        m_reducer(std::move(reducer)),
        m_p(std::move(p)) {}
};

template <class IndexType, class IteratorType1, class IteratorType2,
          class ReducerType, class PredicateType>
struct StdFindFirstOfFunctor {
  using red_value_type = typename ReducerType::value_type;

  IteratorType1 m_first;
  IteratorType2 m_s_first;
  IteratorType2 m_s_last;
  ReducerType m_reducer;
  PredicateType m_p;

  KOKKOS_FUNCTION
  void operator()(const IndexType i, red_value_type& red_value) const {
    namespace KE        = ::Kokkos::Experimental;
    const auto& myvalue = m_first[i];
    bool found          = false;

    const auto search_count = KE::distance(m_s_first, m_s_last);
    for (IndexType k = 0; k < search_count; ++k) {
      if (m_p(myvalue, m_s_first[k])) {
        found = true;
        break;
      }
    }

    const auto rv =
        found ? red_value_type{i}
              : red_value_type{::Kokkos::reduction_identity<IndexType>::min()};

    m_reducer.join(red_value, rv);
  }

  KOKKOS_FUNCTION
  StdFindFirstOfFunctor(IteratorType1 first, IteratorType2 s_first,
                        IteratorType2 s_last, ReducerType reducer,
                        PredicateType p)
      : m_first(std::move(first)),
        m_s_first(std::move(s_first)),
        m_s_last(std::move(s_last)),
        m_reducer(std::move(reducer)),
        m_p(std::move(p)) {}
};

template <class IndexType, class IteratorType1, class IteratorType2,
          class ReducerType, class PredicateType>
struct StdFindEndFunctor {
  using red_value_type = typename ReducerType::value_type;

  IteratorType1 m_first;
  IteratorType1 m_last;
  IteratorType2 m_s_first;
  IteratorType2 m_s_last;
  ReducerType m_reducer;
  PredicateType m_p;

  KOKKOS_FUNCTION
  void operator()(const IndexType i, red_value_type& red_value) const {
    namespace KE = ::Kokkos::Experimental;
    auto myit    = m_first + i;
    bool found   = true;

    const auto search_count = KE::distance(m_s_first, m_s_last);
    for (IndexType k = 0; k < search_count; ++k) {
      // note that we add this EXPECT to check if we are in a valid range
      // but I think we can remvoe this beceause the guarantee we don't go
      // out of bounds is taken care of at the calling site
      // where we launch the par-reduce.
      KOKKOS_EXPECTS((myit + k) < m_last);

      if (!m_p(myit[k], m_s_first[k])) {
        found = false;
        break;
      }
    }

    const auto rv =
        found ? red_value_type{i}
              : red_value_type{::Kokkos::reduction_identity<IndexType>::max()};

    m_reducer.join(red_value, rv);
  }

  KOKKOS_FUNCTION
  StdFindEndFunctor(IteratorType1 first, IteratorType1 last,
                    IteratorType2 s_first, IteratorType2 s_last,
                    ReducerType reducer, PredicateType p)
      : m_first(std::move(first)),
        m_last(std::move(last)),
        m_s_first(std::move(s_first)),
        m_s_last(std::move(s_last)),
        m_reducer(std::move(reducer)),
        m_p(std::move(p)) {}
};

// ------------------------------------------
// find_if_or_not_impl
// ------------------------------------------
template <bool is_find_if, class ExecutionSpace, class IteratorType,
          class PredicateType>
IteratorType find_if_or_not_impl(const std::string& label,
                                 const ExecutionSpace& ex, IteratorType first,
                                 IteratorType last, PredicateType pred) {
  // checks
  Impl::static_assert_random_access_and_accessible(
      ex, first);  // only need one It per type
  Impl::expect_valid_range(first, last);

  if (first == last) {
    return last;
  }

  // aliases
  using index_type           = typename IteratorType::difference_type;
  using reducer_type         = FirstLoc<index_type>;
  using reduction_value_type = typename reducer_type::value_type;
  using func_t = StdFindIfOrNotFunctor<is_find_if, index_type, IteratorType,
                                       reducer_type, PredicateType>;

  // run
  reduction_value_type red_result;
  reducer_type reducer(red_result);
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_reduce(label,
                            RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                            func_t(first, reducer, pred), reducer);

  // fence not needed because reducing into scalar

  // decide and return
  if (red_result.min_loc_true ==
      ::Kokkos::reduction_identity<index_type>::min()) {
    // here, it means a valid loc has not been found,
    return last;
  } else {
    // a location has been found
    return first + red_result.min_loc_true;
  }
}

// ------------------------------------------
// find_impl
// ------------------------------------------
template <class ExecutionSpace, class InputIterator, class T>
InputIterator find_impl(const std::string& label, ExecutionSpace ex,
                        InputIterator first, InputIterator last,
                        const T& value) {
  return find_if_or_not_impl<true>(
      label, ex, first, last,
      ::Kokkos::Experimental::Impl::StdAlgoEqualsValUnaryPredicate<T>(value));
}

// ------------------------------------------
// for_each_impl
// ------------------------------------------
template <class ExecutionSpace, class IteratorType, class UnaryFunctorType>
UnaryFunctorType for_each_impl(const std::string& label,
                               const ExecutionSpace& ex, IteratorType first,
                               IteratorType last, UnaryFunctorType functor) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements),
      StdForEachFunctor<IteratorType, UnaryFunctorType>(first, functor));
  ex.fence("Kokkos::for_each: fence after operation");

  return functor;
}

// ------------------------------------------
// for_each_n_impl
// ------------------------------------------
template <class ExecutionSpace, class IteratorType, class SizeType,
          class UnaryFunctorType>
IteratorType for_each_n_impl(const std::string& label, const ExecutionSpace& ex,
                             IteratorType first, SizeType n,
                             UnaryFunctorType functor) {
  auto last = first + n;
  Impl::static_assert_random_access_and_accessible(ex, first, last);
  Impl::expect_valid_range(first, last);

  if (n == 0) {
    return first;
  }

  for_each_impl(label, ex, first, last, std::move(functor));
  // no neeed to fence since for_each_impl fences already

  return last;
}

// ------------------------------------------
// count_if_impl
// ------------------------------------------
template <class ExecutionSpace, class IteratorType, class Predicate>
typename IteratorType::difference_type count_if_impl(const std::string& label,
                                                     const ExecutionSpace& ex,
                                                     IteratorType first,
                                                     IteratorType last,
                                                     Predicate predicate) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  // aliases
  using func_t = StdCountIfFunctor<IteratorType, Predicate>;

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  typename IteratorType::difference_type count = 0;
  ::Kokkos::parallel_reduce(label,
                            RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                            func_t(first, predicate), count);
  ex.fence("Kokkos::count_if: fence after operation");

  return count;
}

// ------------------------------------------
// count_impl
// ------------------------------------------
template <class ExecutionSpace, class IteratorType, class T>
auto count_impl(const std::string& label, const ExecutionSpace& ex,
                IteratorType first, IteratorType last, const T& value) {
  return count_if_impl(
      label, ex, first, last,
      ::Kokkos::Experimental::Impl::StdAlgoEqualsValUnaryPredicate<T>(value));
}

// ------------------------------------------
// mismatch_impl
// ------------------------------------------
template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
::Kokkos::pair<IteratorType1, IteratorType2> mismatch_impl(
    const std::string& label, const ExecutionSpace& ex, IteratorType1 first1,
    IteratorType1 last1, IteratorType2 first2, IteratorType2 last2,
    BinaryPredicateType predicate) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first1, first2);
  Impl::static_assert_iterators_have_matching_difference_type(first1, first2);
  Impl::expect_valid_range(first1, last1);
  Impl::expect_valid_range(first2, last2);

  // aliases
  using return_type          = ::Kokkos::pair<IteratorType1, IteratorType2>;
  using index_type           = typename IteratorType1::difference_type;
  using reducer_type         = FirstLoc<index_type>;
  using reduction_value_type = typename reducer_type::value_type;
  using functor_type =
      StdMismatchRedFunctor<index_type, IteratorType1, IteratorType2,
                            reducer_type, BinaryPredicateType>;

  // trivial case: note that this is important,
  // for OpenMPTarget, omitting special handling of
  // the trivial case was giving all sorts of strange stuff.
  const auto num_e1 = last1 - first1;
  const auto num_e2 = last2 - first2;
  if (num_e1 == 0 || num_e2 == 0) {
    return return_type(first1, first2);
  }

  // run
  const auto num_elemen_par_reduce = (num_e1 <= num_e2) ? num_e1 : num_e2;
  reduction_value_type red_result;
  reducer_type reducer(red_result);
  ::Kokkos::parallel_reduce(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elemen_par_reduce),
      functor_type(first1, first2, reducer, std::move(predicate)), reducer);

  // fence not needed because reducing into scalar

  // decide and return
  constexpr auto red_min = ::Kokkos::reduction_identity<index_type>::min();
  if (red_result.min_loc_true == red_min) {
    // in here means mismatch has not been found
    if (num_e1 == num_e2) {
      return return_type(last1, last2);
    } else if (num_e1 < num_e2) {
      return return_type(last1, first2 + num_e1);
    } else {
      return return_type(first1 + num_e2, last2);
    }
  } else {
    // in here means mismatch has been found
    return return_type(first1 + red_result.min_loc_true,
                       first2 + red_result.min_loc_true);
  }
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
::Kokkos::pair<IteratorType1, IteratorType2> mismatch_impl(
    const std::string& label, const ExecutionSpace& ex, IteratorType1 first1,
    IteratorType1 last1, IteratorType2 first2, IteratorType2 last2) {
  using value_type1 = typename IteratorType1::value_type;
  using value_type2 = typename IteratorType2::value_type;
  using pred_t      = StdAlgoEqualBinaryPredicate<value_type1, value_type2>;
  return mismatch_impl(label, ex, first1, last1, first2, last2, pred_t());
}

// ------------------------------------------
// all_of_impl, any_of_impl, none_of_impl
// ------------------------------------------
template <class ExecutionSpace, class InputIterator, class Predicate>
bool all_of_impl(const std::string& label, const ExecutionSpace& ex,
                 InputIterator first, InputIterator last, Predicate predicate) {
  return (find_if_or_not_impl<false>(label, ex, first, last, predicate) ==
          last);
}

template <class ExecutionSpace, class InputIterator, class Predicate>
bool any_of_impl(const std::string& label, const ExecutionSpace& ex,
                 InputIterator first, InputIterator last, Predicate predicate) {
  return (find_if_or_not_impl<true>(label, ex, first, last, predicate) != last);
}

template <class ExecutionSpace, class IteratorType, class Predicate>
bool none_of_impl(const std::string& label, const ExecutionSpace& ex,
                  IteratorType first, IteratorType last, Predicate predicate) {
  return (find_if_or_not_impl<true>(label, ex, first, last, predicate) == last);
}

// ------------------------------------------
// equal_impl
// ------------------------------------------
template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
bool equal_impl(const std::string& label, const ExecutionSpace& ex,
                IteratorType1 first1, IteratorType1 last1, IteratorType2 first2,
                BinaryPredicateType predicate) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first1, first2);
  Impl::static_assert_iterators_have_matching_difference_type(first1, first2);
  Impl::expect_valid_range(first1, last1);

  // aliases
  using index_type = typename IteratorType1::difference_type;
  using func_t     = StdEqualFunctor<index_type, IteratorType1, IteratorType2,
                                 BinaryPredicateType>;

  // run
  const auto num_elements = Kokkos::Experimental::distance(first1, last1);
  std::size_t different   = 0;
  ::Kokkos::parallel_reduce(label,
                            RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                            func_t(first1, first2, predicate), different);
  ex.fence("Kokkos::equal: fence after operation");

  return !different;
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
bool equal_impl(const std::string& label, const ExecutionSpace& ex,
                IteratorType1 first1, IteratorType1 last1,
                IteratorType2 first2) {
  using value_type1 = typename IteratorType1::value_type;
  using value_type2 = typename IteratorType2::value_type;
  using pred_t      = StdAlgoEqualBinaryPredicate<value_type1, value_type2>;
  return equal_impl(label, ex, first1, last1, first2, pred_t());
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
bool equal_impl(const std::string& label, const ExecutionSpace& ex,
                IteratorType1 first1, IteratorType1 last1, IteratorType2 first2,
                IteratorType2 last2, BinaryPredicateType predicate) {
  const auto d1 = ::Kokkos::Experimental::distance(first1, last1);
  const auto d2 = ::Kokkos::Experimental::distance(first2, last2);
  if (d1 != d2) {
    return false;
  }

  return equal_impl(label, ex, first1, last1, first2, predicate);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
bool equal_impl(const std::string& label, const ExecutionSpace& ex,
                IteratorType1 first1, IteratorType1 last1, IteratorType2 first2,
                IteratorType2 last2) {
  Impl::expect_valid_range(first1, last1);
  Impl::expect_valid_range(first2, last2);

  using value_type1 = typename IteratorType1::value_type;
  using value_type2 = typename IteratorType2::value_type;
  using pred_t      = StdAlgoEqualBinaryPredicate<value_type1, value_type2>;
  return equal_impl(label, ex, first1, last1, first2, last2, pred_t());
}

// ------------------------------------------
// lexicographical_compare_impl
// ------------------------------------------
template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class ComparatorType>
bool lexicographical_compare_impl(const std::string& label,
                                  const ExecutionSpace& ex,
                                  IteratorType1 first1, IteratorType1 last1,
                                  IteratorType2 first2, IteratorType2 last2,
                                  ComparatorType comp) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first1, first2);
  Impl::static_assert_iterators_have_matching_difference_type(first1, first2);
  Impl::expect_valid_range(first1, last1);
  Impl::expect_valid_range(first2, last2);

  // aliases
  using index_type           = typename IteratorType1::difference_type;
  using reducer_type         = FirstLoc<index_type>;
  using reduction_value_type = typename reducer_type::value_type;

  // run
  const auto d1    = Kokkos::Experimental::distance(first1, last1);
  const auto d2    = Kokkos::Experimental::distance(first2, last2);
  const auto range = Kokkos::Experimental::min(d1, d2);
  reduction_value_type red_result;
  reducer_type reducer(red_result);
  using func1_t =
      StdLexicographicalCompareFunctor<index_type, IteratorType1, IteratorType2,
                                       reducer_type, ComparatorType>;

  ::Kokkos::parallel_reduce(label, RangePolicy<ExecutionSpace>(ex, 0, range),
                            func1_t(first1, first2, reducer, comp), reducer);

  // fence not needed because reducing into scalar
  // no mismatch
  if (red_result.min_loc_true ==
      ::Kokkos::reduction_identity<index_type>::min()) {
    auto new_last1 = first1 + range;
    auto new_last2 = first2 + range;
    bool is_prefix = (new_last1 == last1) && (new_last2 != last2);
    return is_prefix;
  }

  // check mismatched
  int less      = 0;
  auto it1      = first1 + red_result.min_loc_true;
  auto it2      = first2 + red_result.min_loc_true;
  using func2_t = StdCompareFunctor<index_type, IteratorType1, IteratorType2,
                                    ComparatorType>;
  ::Kokkos::parallel_reduce(label, RangePolicy<ExecutionSpace>(ex, 0, 1),
                            func2_t(it1, it2, comp), less);

  // fence not needed because reducing into scalar
  return static_cast<bool>(less);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
bool lexicographical_compare_impl(const std::string& label,
                                  const ExecutionSpace& ex,
                                  IteratorType1 first1, IteratorType1 last1,
                                  IteratorType2 first2, IteratorType2 last2) {
  using value_type_1 = typename IteratorType1::value_type;
  using value_type_2 = typename IteratorType2::value_type;
  using predicate_t =
      Impl::StdAlgoLessThanBinaryPredicate<value_type_1, value_type_2>;
  return lexicographical_compare_impl(label, ex, first1, last1, first2, last2,
                                      predicate_t());
}

// ------------------------------------------
// adjacent_find_impl
// ------------------------------------------
template <class ExecutionSpace, class IteratorType, class PredicateType>
IteratorType adjacent_find_impl(const std::string& label,
                                const ExecutionSpace& ex, IteratorType first,
                                IteratorType last, PredicateType pred) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  const auto num_elements = Kokkos::Experimental::distance(first, last);

  if (num_elements <= 1) {
    return last;
  }

  using index_type           = typename IteratorType::difference_type;
  using reducer_type         = FirstLoc<index_type>;
  using reduction_value_type = typename reducer_type::value_type;
  using func_t = StdAdjacentFindFunctor<index_type, IteratorType, reducer_type,
                                        PredicateType>;

  reduction_value_type red_result;
  reducer_type reducer(red_result);

  // note that we use below num_elements-1 because
  // each index i in the reduction checks i and (i+1).
  ::Kokkos::parallel_reduce(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements - 1),
      func_t(first, reducer, pred), reducer);

  // fence not needed because reducing into scalar
  if (red_result.min_loc_true ==
      ::Kokkos::reduction_identity<index_type>::min()) {
    return last;
  } else {
    return first + red_result.min_loc_true;
  }
}

template <class ExecutionSpace, class IteratorType>
IteratorType adjacent_find_impl(const std::string& label,
                                const ExecutionSpace& ex, IteratorType first,
                                IteratorType last) {
  using value_type     = typename IteratorType::value_type;
  using default_pred_t = StdAlgoEqualBinaryPredicate<value_type>;
  return adjacent_find_impl(label, ex, first, last, default_pred_t());
}

// ------------------------------------------
// search_impl
// ------------------------------------------
template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
IteratorType1 search_impl(const std::string& label, const ExecutionSpace& ex,
                          IteratorType1 first, IteratorType1 last,
                          IteratorType2 s_first, IteratorType2 s_last,
                          const BinaryPredicateType& pred) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first, s_first);
  Impl::static_assert_iterators_have_matching_difference_type(first, s_first);
  Impl::expect_valid_range(first, last);
  Impl::expect_valid_range(s_first, s_last);

  // the target sequence should not be larger than the range [first, last)
  namespace KE            = ::Kokkos::Experimental;
  const auto num_elements = KE::distance(first, last);
  const auto s_count      = KE::distance(s_first, s_last);
  KOKKOS_EXPECTS(num_elements >= s_count);
  (void)s_count;  // needed when macro above is a no-op

  if (s_first == s_last) {
    return first;
  }

  if (first == last) {
    return last;
  }

  // special case where the two ranges have equal size
  if (num_elements == s_count) {
    const auto equal_result = equal_impl(label, ex, first, last, s_first, pred);
    return (equal_result) ? first : last;
  } else {
    using index_type           = typename IteratorType1::difference_type;
    using reducer_type         = FirstLoc<index_type>;
    using reduction_value_type = typename reducer_type::value_type;
    using func_t = StdSearchFunctor<index_type, IteratorType1, IteratorType2,
                                    reducer_type, BinaryPredicateType>;

    // run
    reduction_value_type red_result;
    reducer_type reducer(red_result);

    // decide the size of the range policy of the par_red:
    // note that the last feasible index to start looking is the index
    // whose distance from the "last" is equal to the sequence count.
    // the +1 is because we need to include that location too.
    const auto range_size = num_elements - s_count + 1;

    // run par reduce
    ::Kokkos::parallel_reduce(
        label, RangePolicy<ExecutionSpace>(ex, 0, range_size),
        func_t(first, last, s_first, s_last, reducer, pred), reducer);

    // fence not needed because reducing into scalar

    // decide and return
    if (red_result.min_loc_true ==
        ::Kokkos::reduction_identity<index_type>::min()) {
      // location has not been found
      return last;
    } else {
      // location has been found
      return first + red_result.min_loc_true;
    }
  }
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
IteratorType1 search_impl(const std::string& label, const ExecutionSpace& ex,
                          IteratorType1 first, IteratorType1 last,
                          IteratorType2 s_first, IteratorType2 s_last) {
  using value_type1    = typename IteratorType1::value_type;
  using value_type2    = typename IteratorType2::value_type;
  using predicate_type = StdAlgoEqualBinaryPredicate<value_type1, value_type2>;
  return search_impl(label, ex, first, last, s_first, s_last, predicate_type());
}

// ------------------------------------------
// search_n_impl
// ------------------------------------------
template <class ExecutionSpace, class IteratorType, class SizeType,
          class ValueType, class BinaryPredicateType>
IteratorType search_n_impl(const std::string& label, const ExecutionSpace& ex,
                           IteratorType first, IteratorType last,
                           SizeType count, const ValueType& value,
                           const BinaryPredicateType& pred) {
  // checks
  static_assert_random_access_and_accessible(ex, first);
  expect_valid_range(first, last);
  KOKKOS_EXPECTS((std::ptrdiff_t)count >= 0);

  // count should not be larger than the range [first, last)
  namespace KE            = ::Kokkos::Experimental;
  const auto num_elements = KE::distance(first, last);
  // cast things to avoid compiler warning
  KOKKOS_EXPECTS((std::size_t)num_elements >= (std::size_t)count);

  if (first == last) {
    return first;
  }

  // special case where num elements in [first, last) == count
  if ((std::size_t)num_elements == (std::size_t)count) {
    using equal_to_value = StdAlgoEqualsValUnaryPredicate<ValueType>;
    const auto satisfies =
        all_of_impl(label, ex, first, last, equal_to_value(value));
    return (satisfies) ? first : last;
  } else {
    // aliases
    using index_type           = typename IteratorType::difference_type;
    using reducer_type         = FirstLoc<index_type>;
    using reduction_value_type = typename reducer_type::value_type;
    using func_t =
        StdSearchNFunctor<index_type, IteratorType, SizeType, ValueType,
                          reducer_type, BinaryPredicateType>;

    // run
    reduction_value_type red_result;
    reducer_type reducer(red_result);

    // decide the size of the range policy of the par_red:
    // the last feasible index to start looking is the index
    // whose distance from the "last" is equal to count.
    // the +1 is because we need to include that location too.
    const auto range_size = num_elements - count + 1;

    // run par reduce
    ::Kokkos::parallel_reduce(
        label, RangePolicy<ExecutionSpace>(ex, 0, range_size),
        func_t(first, last, count, value, reducer, pred), reducer);

    // fence not needed because reducing into scalar

    // decide and return
    if (red_result.min_loc_true ==
        ::Kokkos::reduction_identity<index_type>::min()) {
      // location has not been found
      return last;
    } else {
      // location has been found
      return first + red_result.min_loc_true;
    }
  }
}

template <class ExecutionSpace, class IteratorType, class SizeType,
          class ValueType>
IteratorType search_n_impl(const std::string& label, const ExecutionSpace& ex,
                           IteratorType first, IteratorType last,
                           SizeType count, const ValueType& value) {
  using iter_value_type = typename IteratorType::value_type;
  using predicate_type =
      StdAlgoEqualBinaryPredicate<iter_value_type, ValueType>;

  /* above we use <iter_value_type, ValueType> for the predicate_type
     to be consistent with the standard, which says:

     "
     The signature of the predicate function should be equivalent to:

        bool pred(const Type1 &a, const Type2 &b);

     The type Type1 must be such that an object of type ForwardIt can be
     dereferenced and then implicitly converted to Type1. The type Type2 must be
     such that an object of type T can be implicitly converted to Type2.
     "

     In our case, IteratorType = ForwardIt, and ValueType = T.
   */

  return search_n_impl(label, ex, first, last, count, value, predicate_type());
}

// ------------------------------------------
// find_first_of_impl
// ------------------------------------------
template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
IteratorType1 find_first_of_impl(const std::string& label,
                                 const ExecutionSpace& ex, IteratorType1 first,
                                 IteratorType1 last, IteratorType2 s_first,
                                 IteratorType2 s_last,
                                 const BinaryPredicateType& pred) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first, s_first);
  Impl::static_assert_iterators_have_matching_difference_type(first, s_first);
  Impl::expect_valid_range(first, last);
  Impl::expect_valid_range(s_first, s_last);

  if ((s_first == s_last) || (first == last)) {
    return last;
  }

  using index_type           = typename IteratorType1::difference_type;
  using reducer_type         = FirstLoc<index_type>;
  using reduction_value_type = typename reducer_type::value_type;
  using func_t = StdFindFirstOfFunctor<index_type, IteratorType1, IteratorType2,
                                       reducer_type, BinaryPredicateType>;

  // run
  reduction_value_type red_result;
  reducer_type reducer(red_result);
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_reduce(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements),
      func_t(first, s_first, s_last, reducer, pred), reducer);

  // fence not needed because reducing into scalar

  // decide and return
  if (red_result.min_loc_true ==
      ::Kokkos::reduction_identity<index_type>::min()) {
    // if here, nothing found
    return last;
  } else {
    // a location has been found
    return first + red_result.min_loc_true;
  }
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
IteratorType1 find_first_of_impl(const std::string& label,
                                 const ExecutionSpace& ex, IteratorType1 first,
                                 IteratorType1 last, IteratorType2 s_first,
                                 IteratorType2 s_last) {
  using value_type1    = typename IteratorType1::value_type;
  using value_type2    = typename IteratorType2::value_type;
  using predicate_type = StdAlgoEqualBinaryPredicate<value_type1, value_type2>;
  return find_first_of_impl(label, ex, first, last, s_first, s_last,
                            predicate_type());
}

// ------------------------------------------
// find_end_impl
// ------------------------------------------
template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
IteratorType1 find_end_impl(const std::string& label, const ExecutionSpace& ex,
                            IteratorType1 first, IteratorType1 last,
                            IteratorType2 s_first, IteratorType2 s_last,
                            const BinaryPredicateType& pred) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first, s_first);
  Impl::static_assert_iterators_have_matching_difference_type(first, s_first);
  Impl::expect_valid_range(first, last);
  Impl::expect_valid_range(s_first, s_last);

  // the target sequence should not be larger than the range [first, last)
  namespace KE            = ::Kokkos::Experimental;
  const auto num_elements = KE::distance(first, last);
  const auto s_count      = KE::distance(s_first, s_last);
  KOKKOS_EXPECTS(num_elements >= s_count);
  (void)s_count;  // needed when macro above is a no-op

  if (s_first == s_last) {
    return last;
  }

  if (first == last) {
    return last;
  }

  // special case where the two ranges have equal size
  if (num_elements == s_count) {
    const auto equal_result = equal_impl(label, ex, first, last, s_first, pred);
    return (equal_result) ? first : last;
  } else {
    using index_type           = typename IteratorType1::difference_type;
    using reducer_type         = LastLoc<index_type>;
    using reduction_value_type = typename reducer_type::value_type;
    using func_t = StdFindEndFunctor<index_type, IteratorType1, IteratorType2,
                                     reducer_type, BinaryPredicateType>;

    // run
    reduction_value_type red_result;
    reducer_type reducer(red_result);

    // decide the size of the range policy of the par_red:
    // note that the last feasible index to start looking is the index
    // whose distance from the "last" is equal to the sequence count.
    // the +1 is because we need to include that location too.
    const auto range_size = num_elements - s_count + 1;

    // run par reduce
    ::Kokkos::parallel_reduce(
        label, RangePolicy<ExecutionSpace>(ex, 0, range_size),
        func_t(first, last, s_first, s_last, reducer, pred), reducer);

    // fence not needed because reducing into scalar

    // decide and return
    if (red_result.max_loc_true ==
        ::Kokkos::reduction_identity<index_type>::max()) {
      // if here, a subrange has not been found
      return last;
    } else {
      // a location has been found
      return first + red_result.max_loc_true;
    }
  }
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
IteratorType1 find_end_impl(const std::string& label, const ExecutionSpace& ex,
                            IteratorType1 first, IteratorType1 last,
                            IteratorType2 s_first, IteratorType2 s_last) {
  using value_type1    = typename IteratorType1::value_type;
  using value_type2    = typename IteratorType2::value_type;
  using predicate_type = StdAlgoEqualBinaryPredicate<value_type1, value_type2>;
  return find_end_impl(label, ex, first, last, s_first, s_last,
                       predicate_type());
}

}  // namespace Impl

// ----------------------------------
// find public API
// ----------------------------------
template <class ExecutionSpace, class InputIterator, class T>
InputIterator find(const ExecutionSpace& ex, InputIterator first,
                   InputIterator last, const T& value) {
  return Impl::find_impl("Kokkos::find_iterator_api_default", ex, first, last,
                         value);
}

template <class ExecutionSpace, class InputIterator, class T>
InputIterator find(const std::string& label, const ExecutionSpace& ex,
                   InputIterator first, InputIterator last, const T& value) {
  return Impl::find_impl(label, ex, first, last, value);
}

template <class ExecutionSpace, class DataType, class... Properties, class T>
auto find(const ExecutionSpace& ex,
          const ::Kokkos::View<DataType, Properties...>& view, const T& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::find_impl("Kokkos::find_view_api_default", ex, KE::begin(view),
                         KE::end(view), value);
}

template <class ExecutionSpace, class DataType, class... Properties, class T>
auto find(const std::string& label, const ExecutionSpace& ex,
          const ::Kokkos::View<DataType, Properties...>& view, const T& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::find_impl(label, ex, KE::begin(view), KE::end(view), value);
}

// -------------------
// find_if public API
// -------------------
template <class ExecutionSpace, class IteratorType, class PredicateType>
IteratorType find_if(const ExecutionSpace& ex, IteratorType first,
                     IteratorType last, PredicateType predicate) {
  return Impl::find_if_or_not_impl<true>("Kokkos::find_if_iterator_api_default",
                                         ex, first, last, std::move(predicate));
}

template <class ExecutionSpace, class IteratorType, class PredicateType>
IteratorType find_if(const std::string& label, const ExecutionSpace& ex,
                     IteratorType first, IteratorType last,
                     PredicateType predicate) {
  return Impl::find_if_or_not_impl<true>(label, ex, first, last,
                                         std::move(predicate));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class Predicate>
auto find_if(const ExecutionSpace& ex,
             const ::Kokkos::View<DataType, Properties...>& v,
             Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);
  namespace KE = ::Kokkos::Experimental;
  return Impl::find_if_or_not_impl<true>("Kokkos::find_if_view_api_default", ex,
                                         KE::begin(v), KE::end(v),
                                         std::move(predicate));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class Predicate>
auto find_if(const std::string& label, const ExecutionSpace& ex,
             const ::Kokkos::View<DataType, Properties...>& v,
             Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);
  namespace KE = ::Kokkos::Experimental;
  return Impl::find_if_or_not_impl<true>(label, ex, KE::begin(v), KE::end(v),
                                         std::move(predicate));
}

// ----------------------------------
// find_if_not public API
// ----------------------------------
template <class ExecutionSpace, class IteratorType, class Predicate>
IteratorType find_if_not(const ExecutionSpace& ex, IteratorType first,
                         IteratorType last, Predicate predicate) {
  return Impl::find_if_or_not_impl<false>(
      "Kokkos::find_if_not_iterator_api_default", ex, first, last,
      std::move(predicate));
}

template <class ExecutionSpace, class IteratorType, class Predicate>
IteratorType find_if_not(const std::string& label, const ExecutionSpace& ex,
                         IteratorType first, IteratorType last,
                         Predicate predicate) {
  return Impl::find_if_or_not_impl<false>(label, ex, first, last,
                                          std::move(predicate));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class Predicate>
auto find_if_not(const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& v,
                 Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::find_if_or_not_impl<false>(
      "Kokkos::find_if_not_view_api_default", ex, KE::begin(v), KE::end(v),
      std::move(predicate));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class Predicate>
auto find_if_not(const std::string& label, const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& v,
                 Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::find_if_or_not_impl<false>(label, ex, KE::begin(v), KE::end(v),
                                          std::move(predicate));
}

// ----------------------------------
// for_each public API
// ----------------------------------
template <class ExecutionSpace, class IteratorType, class UnaryFunctorType>
UnaryFunctorType for_each(const std::string& label, const ExecutionSpace& ex,
                          IteratorType first, IteratorType last,
                          UnaryFunctorType functor) {
  return Impl::for_each_impl(label, ex, first, last, std::move(functor));
}

template <class ExecutionSpace, class IteratorType, class UnaryFunctorType>
UnaryFunctorType for_each(const ExecutionSpace& ex, IteratorType first,
                          IteratorType last, UnaryFunctorType functor) {
  return Impl::for_each_impl("Kokkos::for_each_iterator_api_default", ex, first,
                             last, std::move(functor));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class UnaryFunctorType>
UnaryFunctorType for_each(const std::string& label, const ExecutionSpace& ex,
                          const ::Kokkos::View<DataType, Properties...>& v,
                          UnaryFunctorType functor) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::for_each_impl(label, ex, KE::begin(v), KE::end(v),
                             std::move(functor));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class UnaryFunctorType>
UnaryFunctorType for_each(const ExecutionSpace& ex,
                          const ::Kokkos::View<DataType, Properties...>& v,
                          UnaryFunctorType functor) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::for_each_impl("Kokkos::for_each_view_api_default", ex,
                             KE::begin(v), KE::end(v), std::move(functor));
}

// ----------------------------------
// for_each_n public API
// ----------------------------------
template <class ExecutionSpace, class IteratorType, class SizeType,
          class UnaryFunctorType>
IteratorType for_each_n(const std::string& label, const ExecutionSpace& ex,
                        IteratorType first, SizeType n,
                        UnaryFunctorType functor) {
  return Impl::for_each_n_impl(label, ex, first, n, std::move(functor));
}

template <class ExecutionSpace, class IteratorType, class SizeType,
          class UnaryFunctorType>
IteratorType for_each_n(const ExecutionSpace& ex, IteratorType first,
                        SizeType n, UnaryFunctorType functor) {
  return Impl::for_each_n_impl("Kokkos::for_each_n_iterator_api_default", ex,
                               first, n, std::move(functor));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class SizeType, class UnaryFunctorType>
auto for_each_n(const std::string& label, const ExecutionSpace& ex,
                const ::Kokkos::View<DataType, Properties...>& v, SizeType n,
                UnaryFunctorType functor) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::for_each_n_impl(label, ex, KE::begin(v), n, std::move(functor));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class SizeType, class UnaryFunctorType>
auto for_each_n(const ExecutionSpace& ex,
                const ::Kokkos::View<DataType, Properties...>& v, SizeType n,
                UnaryFunctorType functor) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::for_each_n_impl("Kokkos::for_each_n_view_api_default", ex,
                               KE::begin(v), n, std::move(functor));
}

// ----------------------------------
// count_if public API
// ----------------------------------
template <class ExecutionSpace, class IteratorType, class Predicate>
typename IteratorType::difference_type count_if(const ExecutionSpace& ex,
                                                IteratorType first,
                                                IteratorType last,
                                                Predicate predicate) {
  return Impl::count_if_impl("Kokkos::count_if_iterator_api_default", ex, first,
                             last, std::move(predicate));
}

template <class ExecutionSpace, class IteratorType, class Predicate>
typename IteratorType::difference_type count_if(const std::string& label,
                                                const ExecutionSpace& ex,
                                                IteratorType first,
                                                IteratorType last,
                                                Predicate predicate) {
  return Impl::count_if_impl(label, ex, first, last, std::move(predicate));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class Predicate>
auto count_if(const ExecutionSpace& ex,
              const ::Kokkos::View<DataType, Properties...>& v,
              Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::count_if_impl("Kokkos::count_if_view_api_default", ex,
                             KE::cbegin(v), KE::cend(v), std::move(predicate));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class Predicate>
auto count_if(const std::string& label, const ExecutionSpace& ex,
              const ::Kokkos::View<DataType, Properties...>& v,
              Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::count_if_impl(label, ex, KE::cbegin(v), KE::cend(v),
                             std::move(predicate));
}

// ----------------------------------
// count public API
// ----------------------------------
template <class ExecutionSpace, class IteratorType, class T>
typename IteratorType::difference_type count(const ExecutionSpace& ex,
                                             IteratorType first,
                                             IteratorType last,
                                             const T& value) {
  return Impl::count_impl("Kokkos::count_iterator_api_default", ex, first, last,
                          value);
}

template <class ExecutionSpace, class IteratorType, class T>
typename IteratorType::difference_type count(const std::string& label,
                                             const ExecutionSpace& ex,
                                             IteratorType first,
                                             IteratorType last,
                                             const T& value) {
  return Impl::count_impl(label, ex, first, last, value);
}

template <class ExecutionSpace, class DataType, class... Properties, class T>
auto count(const ExecutionSpace& ex,
           const ::Kokkos::View<DataType, Properties...>& v, const T& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::count_impl("Kokkos::count_view_api_default", ex, KE::cbegin(v),
                          KE::cend(v), value);
}

template <class ExecutionSpace, class DataType, class... Properties, class T>
auto count(const std::string& label, const ExecutionSpace& ex,
           const ::Kokkos::View<DataType, Properties...>& v, const T& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::count_impl(label, ex, KE::cbegin(v), KE::cend(v), value);
}

// ----------------------------------
// mismatch public API
// ----------------------------------
// FIXME: add mismatch overloads accepting 3 iterators.
// An overload consistent with other algorithms:
//
// auto mismatch(const ExecSpace& ex, It1 first1, It1 last1, It2 first2) {...}
//
// makes API ambiguous (with the overload accepting views).

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
::Kokkos::pair<IteratorType1, IteratorType2> mismatch(const ExecutionSpace& ex,
                                                      IteratorType1 first1,
                                                      IteratorType1 last1,
                                                      IteratorType2 first2,
                                                      IteratorType2 last2) {
  return Impl::mismatch_impl("Kokkos::mismatch_iterator_api_default", ex,
                             first1, last1, first2, last2);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
::Kokkos::pair<IteratorType1, IteratorType2> mismatch(
    const ExecutionSpace& ex, IteratorType1 first1, IteratorType1 last1,
    IteratorType2 first2, IteratorType2 last2,
    BinaryPredicateType&& predicate) {
  return Impl::mismatch_impl("Kokkos::mismatch_iterator_api_default", ex,
                             first1, last1, first2, last2,
                             std::forward<BinaryPredicateType>(predicate));
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
::Kokkos::pair<IteratorType1, IteratorType2> mismatch(
    const std::string& label, const ExecutionSpace& ex, IteratorType1 first1,
    IteratorType1 last1, IteratorType2 first2, IteratorType2 last2) {
  return Impl::mismatch_impl(label, ex, first1, last1, first2, last2);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
::Kokkos::pair<IteratorType1, IteratorType2> mismatch(
    const std::string& label, const ExecutionSpace& ex, IteratorType1 first1,
    IteratorType1 last1, IteratorType2 first2, IteratorType2 last2,
    BinaryPredicateType&& predicate) {
  return Impl::mismatch_impl(label, ex, first1, last1, first2, last2,
                             std::forward<BinaryPredicateType>(predicate));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto mismatch(const ExecutionSpace& ex,
              const ::Kokkos::View<DataType1, Properties1...>& view1,
              const ::Kokkos::View<DataType2, Properties2...>& view2) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::mismatch_impl("Kokkos::mismatch_view_api_default", ex,
                             KE::begin(view1), KE::end(view1), KE::begin(view2),
                             KE::end(view2));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class BinaryPredicateType>
auto mismatch(const ExecutionSpace& ex,
              const ::Kokkos::View<DataType1, Properties1...>& view1,
              const ::Kokkos::View<DataType2, Properties2...>& view2,
              BinaryPredicateType&& predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::mismatch_impl("Kokkos::mismatch_view_api_default", ex,
                             KE::begin(view1), KE::end(view1), KE::begin(view2),
                             KE::end(view2),
                             std::forward<BinaryPredicateType>(predicate));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto mismatch(const std::string& label, const ExecutionSpace& ex,
              const ::Kokkos::View<DataType1, Properties1...>& view1,
              const ::Kokkos::View<DataType2, Properties2...>& view2) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::mismatch_impl(label, ex, KE::begin(view1), KE::end(view1),
                             KE::begin(view2), KE::end(view2));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class BinaryPredicateType>
auto mismatch(const std::string& label, const ExecutionSpace& ex,
              const ::Kokkos::View<DataType1, Properties1...>& view1,
              const ::Kokkos::View<DataType2, Properties2...>& view2,
              BinaryPredicateType&& predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::mismatch_impl(label, ex, KE::begin(view1), KE::end(view1),
                             KE::begin(view2), KE::end(view2),
                             std::forward<BinaryPredicateType>(predicate));
}

// ----------------------------------
// all_of public API
// ----------------------------------
template <class ExecutionSpace, class InputIterator, class Predicate>
bool all_of(const ExecutionSpace& ex, InputIterator first, InputIterator last,
            Predicate predicate) {
  return Impl::all_of_impl("Kokkos::all_of_iterator_api_default", ex, first,
                           last, predicate);
}

template <class ExecutionSpace, class InputIterator, class Predicate>
bool all_of(const std::string& label, const ExecutionSpace& ex,
            InputIterator first, InputIterator last, Predicate predicate) {
  return Impl::all_of_impl(label, ex, first, last, predicate);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class Predicate>
bool all_of(const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& v,
            Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::all_of_impl("Kokkos::all_of_view_api_default", ex, KE::cbegin(v),
                           KE::cend(v), std::move(predicate));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class Predicate>
bool all_of(const std::string& label, const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& v,
            Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::all_of_impl(label, ex, KE::cbegin(v), KE::cend(v),
                           std::move(predicate));
}

// ----------------------------------
// any_of public API
// ----------------------------------
template <class ExecutionSpace, class InputIterator, class Predicate>
bool any_of(const ExecutionSpace& ex, InputIterator first, InputIterator last,
            Predicate predicate) {
  return Impl::any_of_impl("Kokkos::any_of_view_api_default", ex, first, last,
                           predicate);
}

template <class ExecutionSpace, class InputIterator, class Predicate>
bool any_of(const std::string& label, const ExecutionSpace& ex,
            InputIterator first, InputIterator last, Predicate predicate) {
  return Impl::any_of_impl(label, ex, first, last, predicate);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class Predicate>
bool any_of(const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& v,
            Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::any_of_impl("Kokkos::any_of_view_api_default", ex, KE::cbegin(v),
                           KE::cend(v), std::move(predicate));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class Predicate>
bool any_of(const std::string& label, const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& v,
            Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::any_of_impl(label, ex, KE::cbegin(v), KE::cend(v),
                           std::move(predicate));
}

// ----------------------------------
// none_of public API
// ----------------------------------
template <class ExecutionSpace, class IteratorType, class Predicate>
bool none_of(const ExecutionSpace& ex, IteratorType first, IteratorType last,
             Predicate predicate) {
  return Impl::none_of_impl("Kokkos::none_of_iterator_api_default", ex, first,
                            last, predicate);
}

template <class ExecutionSpace, class IteratorType, class Predicate>
bool none_of(const std::string& label, const ExecutionSpace& ex,
             IteratorType first, IteratorType last, Predicate predicate) {
  return Impl::none_of_impl(label, ex, first, last, predicate);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class Predicate>
bool none_of(const ExecutionSpace& ex,
             const ::Kokkos::View<DataType, Properties...>& v,
             Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::none_of_impl("Kokkos::none_of_view_api_default", ex,
                            KE::cbegin(v), KE::cend(v), std::move(predicate));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class Predicate>
bool none_of(const std::string& label, const ExecutionSpace& ex,
             const ::Kokkos::View<DataType, Properties...>& v,
             Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::none_of_impl(label, ex, KE::cbegin(v), KE::cend(v),
                            std::move(predicate));
}

// ----------------------------------
// equal public API
// ----------------------------------
template <class ExecutionSpace, class IteratorType1, class IteratorType2>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      IteratorType1, IteratorType2>::value,
                  bool>
equal(const ExecutionSpace& ex, IteratorType1 first1, IteratorType1 last1,
      IteratorType2 first2) {
  return Impl::equal_impl("Kokkos::equal_iterator_api_default", ex, first1,
                          last1, first2);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      IteratorType1, IteratorType2>::value,
                  bool>
equal(const std::string& label, const ExecutionSpace& ex, IteratorType1 first1,
      IteratorType1 last1, IteratorType2 first2) {
  return Impl::equal_impl(label, ex, first1, last1, first2);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      IteratorType1, IteratorType2>::value,
                  bool>
equal(const ExecutionSpace& ex, IteratorType1 first1, IteratorType1 last1,
      IteratorType2 first2, BinaryPredicateType predicate) {
  return Impl::equal_impl("Kokkos::equal_iterator_api_default", ex, first1,
                          last1, first2, std::move(predicate));
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      IteratorType1, IteratorType2>::value,
                  bool>
equal(const std::string& label, const ExecutionSpace& ex, IteratorType1 first1,
      IteratorType1 last1, IteratorType2 first2,
      BinaryPredicateType predicate) {
  return Impl::equal_impl(label, ex, first1, last1, first2,
                          std::move(predicate));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
bool equal(const ExecutionSpace& ex,
           const ::Kokkos::View<DataType1, Properties1...>& view1,
           ::Kokkos::View<DataType2, Properties2...>& view2) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::equal_impl("Kokkos::equal_view_api_default", ex,
                          KE::cbegin(view1), KE::cend(view1),
                          KE::cbegin(view2));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
bool equal(const std::string& label, const ExecutionSpace& ex,
           const ::Kokkos::View<DataType1, Properties1...>& view1,
           ::Kokkos::View<DataType2, Properties2...>& view2) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::equal_impl(label, ex, KE::cbegin(view1), KE::cend(view1),
                          KE::cbegin(view2));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class BinaryPredicateType>
bool equal(const ExecutionSpace& ex,
           const ::Kokkos::View<DataType1, Properties1...>& view1,
           ::Kokkos::View<DataType2, Properties2...>& view2,
           BinaryPredicateType predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::equal_impl("Kokkos::equal_view_api_default", ex,
                          KE::cbegin(view1), KE::cend(view1), KE::cbegin(view2),
                          std::move(predicate));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class BinaryPredicateType>
bool equal(const std::string& label, const ExecutionSpace& ex,
           const ::Kokkos::View<DataType1, Properties1...>& view1,
           ::Kokkos::View<DataType2, Properties2...>& view2,
           BinaryPredicateType predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::equal_impl(label, ex, KE::cbegin(view1), KE::cend(view1),
                          KE::cbegin(view2), std::move(predicate));
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      IteratorType1, IteratorType2>::value,
                  bool>
equal(const ExecutionSpace& ex, IteratorType1 first1, IteratorType1 last1,
      IteratorType2 first2, IteratorType2 last2) {
  return Impl::equal_impl("Kokkos::equal_iterator_api_default", ex, first1,
                          last1, first2, last2);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      IteratorType1, IteratorType2>::value,
                  bool>
equal(const std::string& label, const ExecutionSpace& ex, IteratorType1 first1,
      IteratorType1 last1, IteratorType2 first2, IteratorType2 last2) {
  return Impl::equal_impl(label, ex, first1, last1, first2, last2);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      IteratorType1, IteratorType2>::value,
                  bool>
equal(const ExecutionSpace& ex, IteratorType1 first1, IteratorType1 last1,
      IteratorType2 first2, IteratorType2 last2,
      BinaryPredicateType predicate) {
  return Impl::equal_impl("Kokkos::equal_iterator_api_default", ex, first1,
                          last1, first2, last2, std::move(predicate));
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      IteratorType1, IteratorType2>::value,
                  bool>
equal(const std::string& label, const ExecutionSpace& ex, IteratorType1 first1,
      IteratorType1 last1, IteratorType2 first2, IteratorType2 last2,
      BinaryPredicateType predicate) {
  return Impl::equal_impl(label, ex, first1, last1, first2, last2,
                          std::move(predicate));
}

// ----------------------------------
// lexicographical_compare public API
// ----------------------------------
template <class ExecutionSpace, class IteratorType1, class IteratorType2>
bool lexicographical_compare(const ExecutionSpace& ex, IteratorType1 first1,
                             IteratorType1 last1, IteratorType2 first2,
                             IteratorType2 last2) {
  return Impl::lexicographical_compare_impl(
      "Kokkos::lexicographical_compare_iterator_api_default", ex, first1, last1,
      first2, last2);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
bool lexicographical_compare(const std::string& label, const ExecutionSpace& ex,
                             IteratorType1 first1, IteratorType1 last1,
                             IteratorType2 first2, IteratorType2 last2) {
  return Impl::lexicographical_compare_impl(label, ex, first1, last1, first2,
                                            last2);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
bool lexicographical_compare(
    const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view1,
    ::Kokkos::View<DataType2, Properties2...>& view2) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::lexicographical_compare_impl(
      "Kokkos::lexicographical_compare_view_api_default", ex, KE::cbegin(view1),
      KE::cend(view1), KE::cbegin(view2), KE::cend(view2));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
bool lexicographical_compare(
    const std::string& label, const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view1,
    ::Kokkos::View<DataType2, Properties2...>& view2) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::lexicographical_compare_impl(label, ex, KE::cbegin(view1),
                                            KE::cend(view1), KE::cbegin(view2),
                                            KE::cend(view2));
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class ComparatorType>
bool lexicographical_compare(const ExecutionSpace& ex, IteratorType1 first1,
                             IteratorType1 last1, IteratorType2 first2,
                             IteratorType2 last2, ComparatorType comp) {
  return Impl::lexicographical_compare_impl(
      "Kokkos::lexicographical_compare_iterator_api_default", ex, first1, last1,
      first2, last2, comp);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class ComparatorType>
bool lexicographical_compare(const std::string& label, const ExecutionSpace& ex,
                             IteratorType1 first1, IteratorType1 last1,
                             IteratorType2 first2, IteratorType2 last2,
                             ComparatorType comp) {
  return Impl::lexicographical_compare_impl(label, ex, first1, last1, first2,
                                            last2, comp);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ComparatorType>
bool lexicographical_compare(
    const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view1,
    ::Kokkos::View<DataType2, Properties2...>& view2, ComparatorType comp) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::lexicographical_compare_impl(
      "Kokkos::lexicographical_compare_view_api_default", ex, KE::cbegin(view1),
      KE::cend(view1), KE::cbegin(view2), KE::cend(view2), comp);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ComparatorType>
bool lexicographical_compare(
    const std::string& label, const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view1,
    ::Kokkos::View<DataType2, Properties2...>& view2, ComparatorType comp) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::lexicographical_compare_impl(label, ex, KE::cbegin(view1),
                                            KE::cend(view1), KE::cbegin(view2),
                                            KE::cend(view2), comp);
}

// ----------------------------------
// adjacent_find
// ----------------------------------
// overload set1
template <class ExecutionSpace, class IteratorType>
IteratorType adjacent_find(const ExecutionSpace& ex, IteratorType first,
                           IteratorType last) {
  return Impl::adjacent_find_impl("Kokkos::adjacent_find_iterator_api_default",
                                  ex, first, last);
}

template <class ExecutionSpace, class IteratorType>
IteratorType adjacent_find(const std::string& label, const ExecutionSpace& ex,
                           IteratorType first, IteratorType last) {
  return Impl::adjacent_find_impl(label, ex, first, last);
}

template <class ExecutionSpace, class DataType, class... Properties>
auto adjacent_find(const ExecutionSpace& ex,
                   const ::Kokkos::View<DataType, Properties...>& v) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);
  namespace KE = ::Kokkos::Experimental;
  return Impl::adjacent_find_impl("Kokkos::adjacent_find_view_api_default", ex,
                                  KE::begin(v), KE::end(v));
}

template <class ExecutionSpace, class DataType, class... Properties>
auto adjacent_find(const std::string& label, const ExecutionSpace& ex,
                   const ::Kokkos::View<DataType, Properties...>& v) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);
  namespace KE = ::Kokkos::Experimental;
  return Impl::adjacent_find_impl(label, ex, KE::begin(v), KE::end(v));
}

// overload set2
template <class ExecutionSpace, class IteratorType, class BinaryPredicateType>
IteratorType adjacent_find(const ExecutionSpace& ex, IteratorType first,
                           IteratorType last, BinaryPredicateType pred) {
  return Impl::adjacent_find_impl("Kokkos::adjacent_find_iterator_api_default",
                                  ex, first, last, pred);
}

template <class ExecutionSpace, class IteratorType, class BinaryPredicateType>
IteratorType adjacent_find(const std::string& label, const ExecutionSpace& ex,
                           IteratorType first, IteratorType last,
                           BinaryPredicateType pred) {
  return Impl::adjacent_find_impl(label, ex, first, last, pred);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class BinaryPredicateType>
auto adjacent_find(const ExecutionSpace& ex,
                   const ::Kokkos::View<DataType, Properties...>& v,
                   BinaryPredicateType pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);
  namespace KE = ::Kokkos::Experimental;
  return Impl::adjacent_find_impl("Kokkos::adjacent_find_view_api_default", ex,
                                  KE::begin(v), KE::end(v), pred);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class BinaryPredicateType>
auto adjacent_find(const std::string& label, const ExecutionSpace& ex,
                   const ::Kokkos::View<DataType, Properties...>& v,
                   BinaryPredicateType pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);
  namespace KE = ::Kokkos::Experimental;
  return Impl::adjacent_find_impl(label, ex, KE::begin(v), KE::end(v), pred);
}

// ----------------------------------
// search
// ----------------------------------
// overload set 1: no binary predicate passed
template <class ExecutionSpace, class IteratorType1, class IteratorType2>
IteratorType1 search(const ExecutionSpace& ex, IteratorType1 first,
                     IteratorType1 last, IteratorType2 s_first,
                     IteratorType2 s_last) {
  return Impl::search_impl("Kokkos::search_iterator_api_default", ex, first,
                           last, s_first, s_last);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
IteratorType1 search(const std::string& label, const ExecutionSpace& ex,
                     IteratorType1 first, IteratorType1 last,
                     IteratorType2 s_first, IteratorType2 s_last) {
  return Impl::search_impl(label, ex, first, last, s_first, s_last);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto search(const ExecutionSpace& ex,
            const ::Kokkos::View<DataType1, Properties1...>& view,
            const ::Kokkos::View<DataType2, Properties2...>& s_view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(s_view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::search_impl("Kokkos::search_view_api_default", ex,
                           KE::begin(view), KE::end(view), KE::begin(s_view),
                           KE::end(s_view));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto search(const std::string& label, const ExecutionSpace& ex,
            const ::Kokkos::View<DataType1, Properties1...>& view,
            const ::Kokkos::View<DataType2, Properties2...>& s_view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(s_view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::search_impl(label, ex, KE::begin(view), KE::end(view),
                           KE::begin(s_view), KE::end(s_view));
}

// overload set 2: binary predicate passed
template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
IteratorType1 search(const ExecutionSpace& ex, IteratorType1 first,
                     IteratorType1 last, IteratorType2 s_first,
                     IteratorType2 s_last, const BinaryPredicateType& pred) {
  return Impl::search_impl("Kokkos::search_iterator_api_default", ex, first,
                           last, s_first, s_last, pred);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
IteratorType1 search(const std::string& label, const ExecutionSpace& ex,
                     IteratorType1 first, IteratorType1 last,
                     IteratorType2 s_first, IteratorType2 s_last,
                     const BinaryPredicateType& pred) {
  return Impl::search_impl(label, ex, first, last, s_first, s_last, pred);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class BinaryPredicateType>
auto search(const ExecutionSpace& ex,
            const ::Kokkos::View<DataType1, Properties1...>& view,
            const ::Kokkos::View<DataType2, Properties2...>& s_view,
            const BinaryPredicateType& pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(s_view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::search_impl("Kokkos::search_view_api_default", ex,
                           KE::begin(view), KE::end(view), KE::begin(s_view),
                           KE::end(s_view), pred);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class BinaryPredicateType>
auto search(const std::string& label, const ExecutionSpace& ex,
            const ::Kokkos::View<DataType1, Properties1...>& view,
            const ::Kokkos::View<DataType2, Properties2...>& s_view,
            const BinaryPredicateType& pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(s_view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::search_impl(label, ex, KE::begin(view), KE::end(view),
                           KE::begin(s_view), KE::end(s_view), pred);
}

// ----------------------------------
// find_first_of
// ----------------------------------
// overload set 1: no binary predicate passed
template <class ExecutionSpace, class IteratorType1, class IteratorType2>
IteratorType1 find_first_of(const ExecutionSpace& ex, IteratorType1 first,
                            IteratorType1 last, IteratorType2 s_first,
                            IteratorType2 s_last) {
  return Impl::find_first_of_impl("Kokkos::find_first_of_iterator_api_default",
                                  ex, first, last, s_first, s_last);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
IteratorType1 find_first_of(const std::string& label, const ExecutionSpace& ex,
                            IteratorType1 first, IteratorType1 last,
                            IteratorType2 s_first, IteratorType2 s_last) {
  return Impl::find_first_of_impl(label, ex, first, last, s_first, s_last);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto find_first_of(const ExecutionSpace& ex,
                   const ::Kokkos::View<DataType1, Properties1...>& view,
                   const ::Kokkos::View<DataType2, Properties2...>& s_view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(s_view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::find_first_of_impl("Kokkos::find_first_of_view_api_default", ex,
                                  KE::begin(view), KE::end(view),
                                  KE::begin(s_view), KE::end(s_view));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto find_first_of(const std::string& label, const ExecutionSpace& ex,
                   const ::Kokkos::View<DataType1, Properties1...>& view,
                   const ::Kokkos::View<DataType2, Properties2...>& s_view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(s_view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::find_first_of_impl(label, ex, KE::begin(view), KE::end(view),
                                  KE::begin(s_view), KE::end(s_view));
}

// overload set 2: binary predicate passed
template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
IteratorType1 find_first_of(const ExecutionSpace& ex, IteratorType1 first,
                            IteratorType1 last, IteratorType2 s_first,
                            IteratorType2 s_last,
                            const BinaryPredicateType& pred) {
  return Impl::find_first_of_impl("Kokkos::find_first_of_iterator_api_default",
                                  ex, first, last, s_first, s_last, pred);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
IteratorType1 find_first_of(const std::string& label, const ExecutionSpace& ex,
                            IteratorType1 first, IteratorType1 last,
                            IteratorType2 s_first, IteratorType2 s_last,
                            const BinaryPredicateType& pred) {
  return Impl::find_first_of_impl(label, ex, first, last, s_first, s_last,
                                  pred);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class BinaryPredicateType>
auto find_first_of(const ExecutionSpace& ex,
                   const ::Kokkos::View<DataType1, Properties1...>& view,
                   const ::Kokkos::View<DataType2, Properties2...>& s_view,
                   const BinaryPredicateType& pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(s_view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::find_first_of_impl("Kokkos::find_first_of_view_api_default", ex,
                                  KE::begin(view), KE::end(view),
                                  KE::begin(s_view), KE::end(s_view), pred);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class BinaryPredicateType>
auto find_first_of(const std::string& label, const ExecutionSpace& ex,
                   const ::Kokkos::View<DataType1, Properties1...>& view,
                   const ::Kokkos::View<DataType2, Properties2...>& s_view,
                   const BinaryPredicateType& pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(s_view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::find_first_of_impl(label, ex, KE::begin(view), KE::end(view),
                                  KE::begin(s_view), KE::end(s_view), pred);
}

// ----------------------------------
// search_n
// ----------------------------------
// overload set 1: no binary predicate passed
template <class ExecutionSpace, class IteratorType, class SizeType,
          class ValueType>
IteratorType search_n(const ExecutionSpace& ex, IteratorType first,
                      IteratorType last, SizeType count,
                      const ValueType& value) {
  return Impl::search_n_impl("Kokkos::search_n_iterator_api_default", ex, first,
                             last, count, value);
}

template <class ExecutionSpace, class IteratorType, class SizeType,
          class ValueType>
IteratorType search_n(const std::string& label, const ExecutionSpace& ex,
                      IteratorType first, IteratorType last, SizeType count,
                      const ValueType& value) {
  return Impl::search_n_impl(label, ex, first, last, count, value);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class SizeType, class ValueType>
auto search_n(const ExecutionSpace& ex,
              const ::Kokkos::View<DataType, Properties...>& view,
              SizeType count, const ValueType& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::search_n_impl("Kokkos::search_n_view_api_default", ex,
                             KE::begin(view), KE::end(view), count, value);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class SizeType, class ValueType>
auto search_n(const std::string& label, const ExecutionSpace& ex,
              const ::Kokkos::View<DataType, Properties...>& view,
              SizeType count, const ValueType& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::search_n_impl(label, ex, KE::begin(view), KE::end(view), count,
                             value);
}

// overload set 2: binary predicate passed
template <class ExecutionSpace, class IteratorType, class SizeType,
          class ValueType, class BinaryPredicateType>
IteratorType search_n(const ExecutionSpace& ex, IteratorType first,
                      IteratorType last, SizeType count, const ValueType& value,
                      const BinaryPredicateType& pred) {
  return Impl::search_n_impl("Kokkos::search_n_iterator_api_default", ex, first,
                             last, count, value, pred);
}

template <class ExecutionSpace, class IteratorType, class SizeType,
          class ValueType, class BinaryPredicateType>
IteratorType search_n(const std::string& label, const ExecutionSpace& ex,
                      IteratorType first, IteratorType last, SizeType count,
                      const ValueType& value, const BinaryPredicateType& pred) {
  return Impl::search_n_impl(label, ex, first, last, count, value, pred);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class SizeType, class ValueType, class BinaryPredicateType>
auto search_n(const ExecutionSpace& ex,
              const ::Kokkos::View<DataType, Properties...>& view,
              SizeType count, const ValueType& value,
              const BinaryPredicateType& pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::search_n_impl("Kokkos::search_n_view_api_default", ex,
                             KE::begin(view), KE::end(view), count, value,
                             pred);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class SizeType, class ValueType, class BinaryPredicateType>
auto search_n(const std::string& label, const ExecutionSpace& ex,
              const ::Kokkos::View<DataType, Properties...>& view,
              SizeType count, const ValueType& value,
              const BinaryPredicateType& pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::search_n_impl(label, ex, KE::begin(view), KE::end(view), count,
                             value, pred);
}

// ----------------------------------
// find_end
// ----------------------------------
// overload set 1: no binary predicate passed
template <class ExecutionSpace, class IteratorType1, class IteratorType2>
IteratorType1 find_end(const ExecutionSpace& ex, IteratorType1 first,
                       IteratorType1 last, IteratorType2 s_first,
                       IteratorType2 s_last) {
  return Impl::find_end_impl("Kokkos::find_end_iterator_api_default", ex, first,
                             last, s_first, s_last);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
IteratorType1 find_end(const std::string& label, const ExecutionSpace& ex,
                       IteratorType1 first, IteratorType1 last,
                       IteratorType2 s_first, IteratorType2 s_last) {
  return Impl::find_end_impl(label, ex, first, last, s_first, s_last);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto find_end(const ExecutionSpace& ex,
              const ::Kokkos::View<DataType1, Properties1...>& view,
              const ::Kokkos::View<DataType2, Properties2...>& s_view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(s_view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::find_end_impl("Kokkos::find_end_view_api_default", ex,
                             KE::begin(view), KE::end(view), KE::begin(s_view),
                             KE::end(s_view));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto find_end(const std::string& label, const ExecutionSpace& ex,
              const ::Kokkos::View<DataType1, Properties1...>& view,
              const ::Kokkos::View<DataType2, Properties2...>& s_view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(s_view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::find_end_impl(label, ex, KE::begin(view), KE::end(view),
                             KE::begin(s_view), KE::end(s_view));
}

// overload set 2: binary predicate passed
template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
IteratorType1 find_end(const ExecutionSpace& ex, IteratorType1 first,
                       IteratorType1 last, IteratorType2 s_first,
                       IteratorType2 s_last, const BinaryPredicateType& pred) {
  return Impl::find_end_impl("Kokkos::find_end_iterator_api_default", ex, first,
                             last, s_first, s_last, pred);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
IteratorType1 find_end(const std::string& label, const ExecutionSpace& ex,
                       IteratorType1 first, IteratorType1 last,
                       IteratorType2 s_first, IteratorType2 s_last,
                       const BinaryPredicateType& pred) {
  return Impl::find_end_impl(label, ex, first, last, s_first, s_last, pred);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class BinaryPredicateType>
auto find_end(const ExecutionSpace& ex,
              const ::Kokkos::View<DataType1, Properties1...>& view,
              const ::Kokkos::View<DataType2, Properties2...>& s_view,
              const BinaryPredicateType& pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(s_view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::find_end_impl("Kokkos::find_end_view_api_default", ex,
                             KE::begin(view), KE::end(view), KE::begin(s_view),
                             KE::end(s_view), pred);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class BinaryPredicateType>
auto find_end(const std::string& label, const ExecutionSpace& ex,
              const ::Kokkos::View<DataType1, Properties1...>& view,
              const ::Kokkos::View<DataType2, Properties2...>& s_view,
              const BinaryPredicateType& pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(s_view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::find_end_impl(label, ex, KE::begin(view), KE::end(view),
                             KE::begin(s_view), KE::end(s_view), pred);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
