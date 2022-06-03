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

#ifndef KOKKOS_MODIFYING_SEQUENCE_OPERATIONS_SET2_HPP
#define KOKKOS_MODIFYING_SEQUENCE_OPERATIONS_SET2_HPP

#include <Kokkos_Core.hpp>
#include "../Kokkos_BeginEnd.hpp"
#include "../Kokkos_Constraints.hpp"
#include "../Kokkos_ModifyingOperations.hpp"
#include "../Kokkos_NonModifyingSequenceOperations.hpp"

namespace Kokkos {
namespace Experimental {
namespace Impl {

//-------------------------
//
// functors
//
//-------------------------

template <class IndexType, class InputIt, class OutputIt,
          class BinaryPredicateType>
struct StdUniqueCopyFunctor {
  InputIt m_first_from;
  InputIt m_last_from;
  OutputIt m_first_dest;
  BinaryPredicateType m_pred;

  KOKKOS_FUNCTION
  StdUniqueCopyFunctor(InputIt first_from, InputIt last_from,
                       OutputIt first_dest, BinaryPredicateType pred)
      : m_first_from(std::move(first_from)),
        m_last_from(std::move(last_from)),
        m_first_dest(std::move(first_dest)),
        m_pred(std::move(pred)) {}

  KOKKOS_FUNCTION
  void operator()(const IndexType i, IndexType& update,
                  const bool final_pass) const {
    const auto& val_i   = m_first_from[i];
    const auto& val_ip1 = m_first_from[i + 1];

    if (final_pass) {
      if (!m_pred(val_i, val_ip1)) {
        m_first_dest[update] = val_i;
      }
    }

    if (!m_pred(val_i, val_ip1)) {
      update += 1;
    }
  }
};

template <class InputIterator>
struct StdReverseFunctor {
  using index_type = typename InputIterator::difference_type;
  static_assert(std::is_signed<index_type>::value,
                "Kokkos: StdReverseFunctor requires signed index type");

  InputIterator m_first;
  InputIterator m_last;

  KOKKOS_FUNCTION
  void operator()(index_type i) const {
    // the swap below is doing the same thing, but
    // for Intel 18.0.5 does not work.
    // But putting the impl directly here, it works.
#ifdef KOKKOS_COMPILER_INTEL
    typename InputIterator::value_type tmp = std::move(m_first[i]);
    m_first[i]                             = std::move(m_last[-i - 1]);
    m_last[-i - 1]                         = std::move(tmp);
#else
    ::Kokkos::Experimental::swap(m_first[i], m_last[-i - 1]);
#endif
  }

  StdReverseFunctor(InputIterator first, InputIterator last)
      : m_first(std::move(first)), m_last(std::move(last)) {}
};

template <class IndexType, class InputIterator, class OutputIterator>
struct StdReverseCopyFunctor {
  static_assert(std::is_signed<IndexType>::value,
                "Kokkos: StdReverseCopyFunctor requires signed index type");

  InputIterator m_last;
  OutputIterator m_dest_first;

  KOKKOS_FUNCTION
  void operator()(IndexType i) const { m_dest_first[i] = m_last[-1 - i]; }

  StdReverseCopyFunctor(InputIterator _last, OutputIterator _dest_first)
      : m_last(std::move(_last)), m_dest_first(std::move(_dest_first)) {}
};

template <class IndexType, class InputIterator, class OutputIterator>
struct StdMoveFunctor {
  InputIterator m_first;
  OutputIterator m_dest_first;

  KOKKOS_FUNCTION
  void operator()(IndexType i) const {
    m_dest_first[i] = std::move(m_first[i]);
  }

  StdMoveFunctor(InputIterator _first, OutputIterator _dest_first)
      : m_first(std::move(_first)), m_dest_first(std::move(_dest_first)) {}
};

template <class IndexType, class IteratorType1, class IteratorType2>
struct StdMoveBackwardFunctor {
  static_assert(std::is_signed<IndexType>::value,
                "Kokkos: StdMoveBackwardFunctor requires signed index type");

  IteratorType1 m_last;
  IteratorType2 m_dest_last;

  KOKKOS_FUNCTION
  void operator()(IndexType i) const {
    m_dest_last[-i] = std::move(m_last[-i]);
  }

  StdMoveBackwardFunctor(IteratorType1 _last, IteratorType2 _dest_last)
      : m_last(std::move(_last)), m_dest_last(std::move(_dest_last)) {}
};

template <class IndexType, class IteratorType1, class IteratorType2>
struct StdSwapRangesFunctor {
  IteratorType1 m_first1;
  IteratorType2 m_first2;

  KOKKOS_FUNCTION
  void operator()(IndexType i) const {
    // the swap below is doing the same thing, but
    // for Intel 18.0.5 does not work.
    // But putting the impl directly here, it works.
#ifdef KOKKOS_COMPILER_INTEL
    typename IteratorType1::value_type tmp = std::move(m_first1[i]);
    m_first1[i]                            = std::move(m_first2[i]);
    m_first2[i]                            = std::move(tmp);
#else
    ::Kokkos::Experimental::swap(m_first1[i], m_first2[i]);
#endif
  }

  KOKKOS_FUNCTION
  StdSwapRangesFunctor(IteratorType1 _first1, IteratorType2 _first2)
      : m_first1(std::move(_first1)), m_first2(std::move(_first2)) {}
};

template <class IndexType, class InputIt, class OutputIt,
          class BinaryPredicateType>
struct StdUniqueFunctor {
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
  void operator()(const IndexType i, IndexType& update,
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

template <class IndexType, class InputIterator, class OutputIterator>
struct StdRotateCopyFunctor {
  InputIterator m_first;
  InputIterator m_last;
  InputIterator m_first_n;
  OutputIterator m_dest_first;

  KOKKOS_FUNCTION
  void operator()(IndexType i) const {
    const IndexType shift = m_last - m_first_n;

    if (i < shift) {
      m_dest_first[i] = m_first_n[i];
    } else {
      m_dest_first[i] = m_first[i - shift];
    }
  }

  StdRotateCopyFunctor(InputIterator first, InputIterator last,
                       InputIterator first_n, OutputIterator dest_first)
      : m_first(std::move(first)),
        m_last(std::move(last)),
        m_first_n(std::move(first_n)),
        m_dest_first(std::move(dest_first)) {}
};

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

// ------------------------------------------
// unique_copy_impl
// ------------------------------------------
template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class PredicateType>
OutputIterator unique_copy_impl(const std::string& label,
                                const ExecutionSpace& ex, InputIterator first,
                                InputIterator last, OutputIterator d_first,
                                PredicateType pred) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first, d_first);
  Impl::static_assert_iterators_have_matching_difference_type(first, d_first);
  Impl::expect_valid_range(first, last);

  // branch for trivial vs non trivial case
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  if (num_elements == 0) {
    return d_first;
  } else if (num_elements == 1) {
    return Impl::copy_impl("Kokkos::copy_from_unique_copy", ex, first, last,
                           d_first);
  } else {
    // aliases
    using index_type = typename InputIterator::difference_type;
    using func_type  = StdUniqueCopyFunctor<index_type, InputIterator,
                                           OutputIterator, PredicateType>;

    // note here that we run scan for num_elements - 1
    // because of the way we implement this, the last element is always needed.
    // We avoid performing checks inside functor that we are within limits
    // and run a "safe" scan and then copy the last element.
    const auto scan_size = num_elements - 1;
    index_type count     = 0;
    ::Kokkos::parallel_scan(label,
                            RangePolicy<ExecutionSpace>(ex, 0, scan_size),
                            func_type(first, last, d_first, pred), count);

    return Impl::copy_impl("Kokkos::copy_from_unique_copy", ex,
                           first + scan_size, last, d_first + count);
  }
}

template <class ExecutionSpace, class InputIterator, class OutputIterator>
OutputIterator unique_copy_impl(const std::string& label,
                                const ExecutionSpace& ex, InputIterator first,
                                InputIterator last, OutputIterator d_first) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first, d_first);
  Impl::static_assert_iterators_have_matching_difference_type(first, d_first);
  Impl::expect_valid_range(first, last);

  // aliases
  using value_type1 = typename InputIterator::value_type;
  using value_type2 = typename OutputIterator::value_type;

  // default binary predicate uses ==
  using binary_pred_t = StdAlgoEqualBinaryPredicate<value_type1, value_type2>;

  // run
  return unique_copy_impl(label, ex, first, last, d_first, binary_pred_t());
}

// ------------------------------------------
// reverse_impl
// ------------------------------------------
template <class ExecutionSpace, class InputIterator>
void reverse_impl(const std::string& label, const ExecutionSpace& ex,
                  InputIterator first, InputIterator last) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  // aliases
  using func_t = StdReverseFunctor<InputIterator>;

  // run
  if (last >= first + 2) {
    // only need half
    const auto num_elements = Kokkos::Experimental::distance(first, last) / 2;
    ::Kokkos::parallel_for(label,
                           RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                           func_t(first, last));
    ex.fence("Kokkos::reverse: fence after operation");
  }
}

// ------------------------------------------
// reverse_copy_impl
// ------------------------------------------
template <class ExecutionSpace, class InputIterator, class OutputIterator>
OutputIterator reverse_copy_impl(const std::string& label,
                                 const ExecutionSpace& ex, InputIterator first,
                                 InputIterator last, OutputIterator d_first) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first, d_first);
  Impl::static_assert_iterators_have_matching_difference_type(first, d_first);
  Impl::expect_valid_range(first, last);

  // aliases
  using index_type = typename InputIterator::difference_type;
  using func_t =
      StdReverseCopyFunctor<index_type, InputIterator, OutputIterator>;

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                         func_t(last, d_first));
  ex.fence("Kokkos::reverse_copy: fence after operation");

  // return
  return d_first + num_elements;
}

// ------------------------------------------
// move_impl
// ------------------------------------------
template <class ExecutionSpace, class InputIterator, class OutputIterator>
OutputIterator move_impl(const std::string& label, const ExecutionSpace& ex,
                         InputIterator first, InputIterator last,
                         OutputIterator d_first) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first, d_first);
  Impl::static_assert_iterators_have_matching_difference_type(first, d_first);
  Impl::expect_valid_range(first, last);

  // aliases
  using index_type = typename InputIterator::difference_type;
  using func_t     = StdMoveFunctor<index_type, InputIterator, OutputIterator>;

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                         func_t(first, d_first));
  ex.fence("Kokkos::move: fence after operation");

  // return
  return d_first + num_elements;
}

// ------------------------------------------
// move_backward_impl
// ------------------------------------------
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

// ------------------------------------------
// swap_ranges_impl
// ------------------------------------------
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

// ------------------------------------------
// unique_impl
// ------------------------------------------
template <class ExecutionSpace, class IteratorType, class PredicateType>
IteratorType unique_impl(const std::string& label, const ExecutionSpace& ex,
                         IteratorType first, IteratorType last,
                         PredicateType pred) {
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
      using output_it      = decltype(tmp_first);

      using index_type = typename IteratorType::difference_type;
      using func_type =
          StdUniqueFunctor<index_type, IteratorType, output_it, PredicateType>;
      index_type count = 0;
      ::Kokkos::parallel_scan(
          label, RangePolicy<ExecutionSpace>(ex, 0, scan_size),
          func_type(it_found, last, tmp_first, pred), count);

      // move last element too, for the same reason as the unique_copy
      auto unused_r =
          Impl::move_impl("Kokkos::move_from_unique", ex, it_found + scan_size,
                          last, tmp_first + count);
      (void)unused_r;  // r1 not used

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
IteratorType unique_impl(const std::string& label, const ExecutionSpace& ex,
                         IteratorType first, IteratorType last) {
  using value_type    = typename IteratorType::value_type;
  using binary_pred_t = StdAlgoEqualBinaryPredicate<value_type>;
  return unique_impl(label, ex, first, last, binary_pred_t());
}

// ------------------------------------------
// rotate_copy_impl
// ------------------------------------------
template <class ExecutionSpace, class InputIterator, class OutputIterator>
OutputIterator rotate_copy_impl(const std::string& label,
                                const ExecutionSpace& ex, InputIterator first,
                                InputIterator n_first, InputIterator last,
                                OutputIterator d_first) {
  /*
    algorithm is implemented as follows:

    first 	   n_first		last
    |		      |                  |
    o  o  o  o  o  o  o  o  o  o  o  o

    dest+0 -> first_n
    dest+1 -> first_n+1
    dest+2 -> first_n+2
    dest+3 -> first
    dest+4 -> first+1
    dest+5 -> first+2
    dest+6 -> first+3
    dest+7 -> first+4
    dest+8 -> first+5
    ...
    let shift = last - first_n;

    then we have:
    if (i < shift){
      *(dest_first + i) = *(first_n + i);
    }
    else{
      *(dest_first + i) = *(from + i - shift);
    }
  */

  // checks
  Impl::static_assert_random_access_and_accessible(ex, first, d_first);
  Impl::static_assert_iterators_have_matching_difference_type(first, d_first);
  Impl::expect_valid_range(first, last);
  Impl::expect_valid_range(first, n_first);
  Impl::expect_valid_range(n_first, last);

  if (first == last) {
    return d_first;
  }

  // aliases
  using index_type = typename InputIterator::difference_type;
  using func_type =
      StdRotateCopyFunctor<index_type, InputIterator, OutputIterator>;

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                         func_type(first, last, n_first, d_first));

  ex.fence("Kokkos::rotate_copy: fence after operation");

  // return
  return d_first + num_elements;
}

// ------------------------------------------
// rotate_impl
// ------------------------------------------
template <class ExecutionSpace, class IteratorType>
IteratorType rotate_with_pivot_in_left_half(const std::string& label,
                                            const ExecutionSpace& ex,
                                            IteratorType first,
                                            IteratorType n_first,
                                            IteratorType last) {
  /*
    This impl is specific for when the n_first iterator points to
    an element that is before or equal to the middle of the range.

    If we have:

    | 0 | 1 | 2 | 1 | 4 | 5 | 2 | 2 | 10 | -3 | 1 | -6 | -5 | 8 | 9 | 11 | *
      ^           ^              mid					   ^
    first       n_first							  last

    In step 1, we create a temporary view with extent = distance(n_first, last)
    and *move* the elements from [n_first, last) to tmp view, such that
    tmp view becomes:

    | 1 | 4 | 5 | 2 | 2 | 10 | -3 | 1 | -6 | -5 | 8 | 9 | 11 |

    In step 2, we move the elements in [first, n_first)
    to the new position where they are supposed to end up.

    In step 3, we move the elements from the tmp view to
    the range starting at first.
   */

  namespace KE                     = ::Kokkos::Experimental;
  const auto num_elements_on_left  = KE::distance(first, n_first);
  const auto num_elements_on_right = KE::distance(n_first, last);

  // create helper tmp view
  using value_type    = typename IteratorType::value_type;
  using tmp_view_type = Kokkos::View<value_type*, ExecutionSpace>;
  tmp_view_type tmp_view("rotate_impl_for_pivot_in_left_half_impl",
                         num_elements_on_right);
  using tmp_readwrite_iterator_type = decltype(begin(tmp_view));

  // index_type is the same and needed in all steps
  using index_type = typename IteratorType::difference_type;

  // stage 1
  using step1_func_type =
      StdMoveFunctor<index_type, IteratorType, tmp_readwrite_iterator_type>;
  ::Kokkos::parallel_for(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements_on_right),
      step1_func_type(n_first, begin(tmp_view)));

  // stage 2
  using step2_func_type =
      StdMoveFunctor<index_type, IteratorType, IteratorType>;
  ::Kokkos::parallel_for(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements_on_left),
      step2_func_type(first, first + num_elements_on_right));

  // step 3
  using step3_func_type =
      StdMoveFunctor<index_type, tmp_readwrite_iterator_type, IteratorType>;
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, tmp_view.extent(0)),
                         step3_func_type(begin(tmp_view), first));

  ex.fence("Kokkos::rotate: fence after operation");
  return first + (last - n_first);
}

template <class ExecutionSpace, class IteratorType>
IteratorType rotate_with_pivot_in_right_half(const std::string& label,
                                             const ExecutionSpace& ex,
                                             IteratorType first,
                                             IteratorType n_first,
                                             IteratorType last) {
  /*
    This impl is specific for when the n_first iterator points to
    an element that is after the middle of the range.

    If we have:

    | 0 | 1 | 2 | 1 | 4 | 5 | 2 | 2 | 10 | -3 | 1 | -6 | -5 | 8 | 9 | 11 | *
      ^                          mid            ^                          ^
    first                                    n_first			  last

    In step 1, we create a temporary view with extent = distance(first, n_first)
    and *move* the elements from [first, n_first) to tmp view,
    such that tmp view becomes:

    | 0 | 1 | 2 | 1 | 4 | 5 | 2 | 2 | 10 | -3 | 1 |

    In step 2, we move the elements in [n_first, last)
    to the beginning where they are supposed to end up.

    In step 3, we move the elements from the tmp view to
    the range starting at first.
   */

  namespace KE                     = ::Kokkos::Experimental;
  const auto num_elements_on_left  = KE::distance(first, n_first);
  const auto num_elements_on_right = KE::distance(n_first, last);

  // create helper tmp view
  using value_type    = typename IteratorType::value_type;
  using tmp_view_type = Kokkos::View<value_type*, ExecutionSpace>;
  tmp_view_type tmp_view("rotate_impl_for_pivot_in_left_half_impl",
                         num_elements_on_left);
  using tmp_readwrite_iterator_type = decltype(begin(tmp_view));

  // index_type is the same and needed in all steps
  using index_type = typename IteratorType::difference_type;

  // stage 1
  using step1_func_type =
      StdMoveFunctor<index_type, IteratorType, tmp_readwrite_iterator_type>;
  ::Kokkos::parallel_for(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements_on_left),
      step1_func_type(first, begin(tmp_view)));

  // stage 2
  using step2_func_type =
      StdMoveFunctor<index_type, IteratorType, IteratorType>;
  ::Kokkos::parallel_for(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements_on_right),
      step2_func_type(n_first, first));

  // step 3:
  using step3_func_type =
      StdMoveFunctor<index_type, tmp_readwrite_iterator_type, IteratorType>;
  ::Kokkos::parallel_for(
      label, RangePolicy<ExecutionSpace>(ex, 0, tmp_view.extent(0)),
      step3_func_type(begin(tmp_view), first + num_elements_on_right));

  ex.fence("Kokkos::rotate: fence after operation");
  return first + (last - n_first);
}

template <class ExecutionSpace, class IteratorType>
IteratorType rotate_impl(const std::string& label, const ExecutionSpace& ex,
                         IteratorType first, IteratorType n_first,
                         IteratorType last) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);
  Impl::expect_valid_range(first, n_first);
  Impl::expect_valid_range(n_first, last);

  namespace KE                     = ::Kokkos::Experimental;
  const auto num_elements          = KE::distance(first, last);
  const auto n_distance_from_first = KE::distance(first, n_first);
  if (n_distance_from_first <= num_elements / 2) {
    return rotate_with_pivot_in_left_half(label, ex, first, n_first, last);
  } else {
    return rotate_with_pivot_in_right_half(label, ex, first, n_first, last);
  }
}

// ------------------------------------------
// remove_if_impl
// ------------------------------------------
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

// ------------------------------------------
// remove_impl
// ------------------------------------------
template <class ExecutionSpace, class IteratorType, class ValueType>
auto remove_impl(const std::string& label, const ExecutionSpace& ex,
                 IteratorType first, IteratorType last,
                 const ValueType& value) {
  using predicate_type = StdAlgoEqualsValUnaryPredicate<ValueType>;
  return remove_if_impl(label, ex, first, last, predicate_type(value));
}

// ------------------------------------------
// remove_copy_impl
// ------------------------------------------
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

template <class ExecutionSpace, class IteratorType>
IteratorType shift_left_impl(const std::string& label, const ExecutionSpace& ex,
                             IteratorType first, IteratorType last,
                             typename IteratorType::difference_type n) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);
  KOKKOS_EXPECTS(n >= 0);

  // handle trivial cases
  if (n == 0) {
    return last;
  }

  if (n >= Kokkos::Experimental::distance(first, last)) {
    return first;
  }

  /*
    Suppose that n = 5, and our [first,last) spans:

    | 0  | 1  |  2 | 1  | 2  | 1  | 2  | 2  | 10 | -3 | 1  | -6 | *
      ^                         				  ^
    first							 last

    shift_left modifies the range such that we have this data:
    | 1  | 2  | 2  | 10  | -3 | 1  | -6 | x | x  | x  | x  |  x | *
                                          ^
                                   return it pointing here


    and returns an iterator pointing to one past the new end.
    Note: elements marked x are in undefined state because have been moved.

    We implement this in two steps:
    step 1:
      we create a temporary view with extent = distance(first+n, last)
      and *move* assign the elements from [first+n, last) to tmp view, such that
      tmp view becomes:

      | 1  | 2  | 2  | 10  | -3 | 1  | -6 |

    step 2:
      move elements of tmp view back to range starting at first.
   */

  const auto num_elements_to_move =
      ::Kokkos::Experimental::distance(first + n, last);

  // create tmp view
  using value_type    = typename IteratorType::value_type;
  using tmp_view_type = Kokkos::View<value_type*, ExecutionSpace>;
  tmp_view_type tmp_view("shift_left_impl", num_elements_to_move);
  using tmp_readwrite_iterator_type = decltype(begin(tmp_view));

  using index_type = typename IteratorType::difference_type;

  // step 1
  using step1_func_type =
      StdMoveFunctor<index_type, IteratorType, tmp_readwrite_iterator_type>;
  ::Kokkos::parallel_for(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements_to_move),
      step1_func_type(first + n, begin(tmp_view)));

  // step 2
  using step2_func_type =
      StdMoveFunctor<index_type, tmp_readwrite_iterator_type, IteratorType>;
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, tmp_view.extent(0)),
                         step2_func_type(begin(tmp_view), first));

  ex.fence("Kokkos::shift_left: fence after operation");

  return last - n;
}

template <class ExecutionSpace, class IteratorType>
IteratorType shift_right_impl(const std::string& label,
                              const ExecutionSpace& ex, IteratorType first,
                              IteratorType last,
                              typename IteratorType::difference_type n) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);
  KOKKOS_EXPECTS(n >= 0);

  // handle trivial cases
  if (n == 0) {
    return first;
  }

  if (n >= Kokkos::Experimental::distance(first, last)) {
    return last;
  }

  /*
    Suppose that n = 3, and [first,last) spans:

    | 0  | 1  |  2 | 1  | 2  | 1  | 2  | 2  | 10 | -3 | 1  | -6 | *
      ^                         				  ^
    first							 last

    shift_right modifies the range such that we have this data:
    |  x | x  | x  | 0  | 1  |  2 | 1  | 2  | 1  | 2  | 2  | 10 | *
                     ^
             return it points here

    and returns an iterator pointing to the new beginning.
    Note: elements marked x are in undefined state because have been moved.

    We implement this in two steps:
    step 1:
      we create a temporary view with extent = distance(first, last-n)
      and *move* assign the elements from [first, last-n) to tmp view, such that
      tmp view becomes:

      | 0  | 1  |  2 | 1  | 2  | 1  | 2  | 2  | 10 |

    step 2:
      move elements of tmp view back to range starting at first+n.
   */

  const auto num_elements_to_move =
      ::Kokkos::Experimental::distance(first, last - n);

  // create tmp view
  using value_type    = typename IteratorType::value_type;
  using tmp_view_type = Kokkos::View<value_type*, ExecutionSpace>;
  tmp_view_type tmp_view("shift_right_impl", num_elements_to_move);
  using tmp_readwrite_iterator_type = decltype(begin(tmp_view));

  using index_type = typename IteratorType::difference_type;

  // step 1
  using step1_func_type =
      StdMoveFunctor<index_type, IteratorType, tmp_readwrite_iterator_type>;
  ::Kokkos::parallel_for(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements_to_move),
      step1_func_type(first, begin(tmp_view)));

  // step 2
  using step2_func_type =
      StdMoveFunctor<index_type, tmp_readwrite_iterator_type, IteratorType>;
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, tmp_view.extent(0)),
                         step2_func_type(begin(tmp_view), first + n));

  ex.fence("Kokkos::shift_right: fence after operation");

  return first + n;
}

}  // namespace Impl

// -------------------
// reverse_copy
// -------------------
template <class ExecutionSpace, class InputIterator, class OutputIterator>
OutputIterator reverse_copy(const ExecutionSpace& ex, InputIterator first,
                            InputIterator last, OutputIterator d_first) {
  return Impl::reverse_copy_impl("Kokkos::reverse_copy_iterator_api_default",
                                 ex, first, last, d_first);
}

template <class ExecutionSpace, class InputIterator, class OutputIterator>
OutputIterator reverse_copy(const std::string& label, const ExecutionSpace& ex,
                            InputIterator first, InputIterator last,
                            OutputIterator d_first) {
  return Impl::reverse_copy_impl(label, ex, first, last, d_first);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto reverse_copy(const ExecutionSpace& ex,
                  const ::Kokkos::View<DataType1, Properties1...>& source,
                  ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::reverse_copy_impl("Kokkos::reverse_copy_view_api_default", ex,
                                 cbegin(source), cend(source), begin(dest));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto reverse_copy(const std::string& label, const ExecutionSpace& ex,
                  const ::Kokkos::View<DataType1, Properties1...>& source,
                  ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::reverse_copy_impl(label, ex, cbegin(source), cend(source),
                                 begin(dest));
}

// -------------------
// reverse
// -------------------
template <class ExecutionSpace, class InputIterator>
void reverse(const ExecutionSpace& ex, InputIterator first,
             InputIterator last) {
  return Impl::reverse_impl("Kokkos::reverse_iterator_api_default", ex, first,
                            last);
}

template <class ExecutionSpace, class InputIterator>
void reverse(const std::string& label, const ExecutionSpace& ex,
             InputIterator first, InputIterator last) {
  return Impl::reverse_impl(label, ex, first, last);
}

template <class ExecutionSpace, class DataType, class... Properties>
void reverse(const ExecutionSpace& ex,
             const ::Kokkos::View<DataType, Properties...>& view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  namespace KE = ::Kokkos::Experimental;
  return Impl::reverse_impl("Kokkos::reverse_view_api_default", ex,
                            KE::begin(view), KE::end(view));
}

template <class ExecutionSpace, class DataType, class... Properties>
void reverse(const std::string& label, const ExecutionSpace& ex,
             const ::Kokkos::View<DataType, Properties...>& view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  namespace KE = ::Kokkos::Experimental;
  return Impl::reverse_impl(label, ex, KE::begin(view), KE::end(view));
}

// ----------------------
// move
// ----------------------
template <class ExecutionSpace, class InputIterator, class OutputIterator>
OutputIterator move(const ExecutionSpace& ex, InputIterator first,
                    InputIterator last, OutputIterator d_first) {
  return Impl::move_impl("Kokkos::move_iterator_api_default", ex, first, last,
                         d_first);
}

template <class ExecutionSpace, class InputIterator, class OutputIterator>
OutputIterator move(const std::string& label, const ExecutionSpace& ex,
                    InputIterator first, InputIterator last,
                    OutputIterator d_first) {
  return Impl::move_impl(label, ex, first, last, d_first);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto move(const ExecutionSpace& ex,
          const ::Kokkos::View<DataType1, Properties1...>& source,
          ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::move_impl("Kokkos::move_view_api_default", ex, begin(source),
                         end(source), begin(dest));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto move(const std::string& label, const ExecutionSpace& ex,
          const ::Kokkos::View<DataType1, Properties1...>& source,
          ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::move_impl(label, ex, begin(source), end(source), begin(dest));
}

// -------------------
// move_backward
// -------------------
template <class ExecutionSpace, class IteratorType1, class IteratorType2>
IteratorType2 move_backward(const ExecutionSpace& ex, IteratorType1 first,
                            IteratorType1 last, IteratorType2 d_last) {
  return Impl::move_backward_impl("Kokkos::move_backward_iterator_api_default",
                                  ex, first, last, d_last);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto move_backward(const ExecutionSpace& ex,
                   const ::Kokkos::View<DataType1, Properties1...>& source,
                   ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::move_backward_impl("Kokkos::move_backward_view_api_default", ex,
                                  begin(source), end(source), end(dest));
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
IteratorType2 move_backward(const std::string& label, const ExecutionSpace& ex,
                            IteratorType1 first, IteratorType1 last,
                            IteratorType2 d_last) {
  return Impl::move_backward_impl(label, ex, first, last, d_last);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto move_backward(const std::string& label, const ExecutionSpace& ex,
                   const ::Kokkos::View<DataType1, Properties1...>& source,
                   ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::move_backward_impl(label, ex, begin(source), end(source),
                                  end(dest));
}

// ----------------------
// swap_ranges
// ----------------------
template <class ExecutionSpace, class IteratorType1, class IteratorType2>
IteratorType2 swap_ranges(const ExecutionSpace& ex, IteratorType1 first1,
                          IteratorType1 last1, IteratorType2 first2) {
  return Impl::swap_ranges_impl("Kokkos::swap_ranges_iterator_api_default", ex,
                                first1, last1, first2);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto swap_ranges(const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType1, Properties1...>& source,
                 ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  assert(source.extent(0) == dest.extent(0));
  return Impl::swap_ranges_impl("Kokkos::swap_ranges_view_api_default", ex,
                                begin(source), end(source), begin(dest));
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
IteratorType2 swap_ranges(const std::string& label, const ExecutionSpace& ex,
                          IteratorType1 first1, IteratorType1 last1,
                          IteratorType2 first2) {
  return Impl::swap_ranges_impl(label, ex, first1, last1, first2);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto swap_ranges(const std::string& label, const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType1, Properties1...>& source,
                 ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  assert(source.extent(0) == dest.extent(0));
  return Impl::swap_ranges_impl(label, ex, begin(source), end(source),
                                begin(dest));
}

// -------------------
// unique
// -------------------
// note: the enable_if below is to avoid "call to ... is ambiguous"
// for example in the unit test when using a variadic function

// overload set1
template <class ExecutionSpace, class IteratorType>
std::enable_if_t<!::Kokkos::is_view<IteratorType>::value, IteratorType> unique(
    const ExecutionSpace& ex, IteratorType first, IteratorType last) {
  return Impl::unique_impl("Kokkos::unique_iterator_api_default", ex, first,
                           last);
}

template <class ExecutionSpace, class IteratorType>
std::enable_if_t<!::Kokkos::is_view<IteratorType>::value, IteratorType> unique(
    const std::string& label, const ExecutionSpace& ex, IteratorType first,
    IteratorType last) {
  return Impl::unique_impl(label, ex, first, last);
}

template <class ExecutionSpace, class DataType, class... Properties>
auto unique(const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return ::Kokkos::Experimental::unique("Kokkos::unique_view_api_default", ex,
                                        begin(view), end(view));
}

template <class ExecutionSpace, class DataType, class... Properties>
auto unique(const std::string& label, const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return ::Kokkos::Experimental::unique(label, ex, begin(view), end(view));
}

// overload set2
template <class ExecutionSpace, class IteratorType, class BinaryPredicate>
IteratorType unique(const ExecutionSpace& ex, IteratorType first,
                    IteratorType last, BinaryPredicate pred) {
  return Impl::unique_impl("Kokkos::unique_iterator_api_default", ex, first,
                           last, pred);
}

template <class ExecutionSpace, class IteratorType, class BinaryPredicate>
IteratorType unique(const std::string& label, const ExecutionSpace& ex,
                    IteratorType first, IteratorType last,
                    BinaryPredicate pred) {
  return Impl::unique_impl(label, ex, first, last, pred);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class BinaryPredicate>
auto unique(const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view,
            BinaryPredicate pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::unique_impl("Kokkos::unique_view_api_default", ex, begin(view),
                           end(view), std::move(pred));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class BinaryPredicate>
auto unique(const std::string& label, const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view,
            BinaryPredicate pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::unique_impl(label, ex, begin(view), end(view), std::move(pred));
}

// -------------------
// unique_copy
// -------------------
// note: the enable_if below is to avoid "call to ... is ambiguous"
// for example in the unit test when using a variadic function

// overload set1
template <class ExecutionSpace, class InputIterator, class OutputIterator>
std::enable_if_t<!::Kokkos::is_view<InputIterator>::value, OutputIterator>
unique_copy(const ExecutionSpace& ex, InputIterator first, InputIterator last,
            OutputIterator d_first) {
  return Impl::unique_copy_impl("Kokkos::unique_copy_iterator_api_default", ex,
                                first, last, d_first);
}

template <class ExecutionSpace, class InputIterator, class OutputIterator>
std::enable_if_t<!::Kokkos::is_view<InputIterator>::value, OutputIterator>
unique_copy(const std::string& label, const ExecutionSpace& ex,
            InputIterator first, InputIterator last, OutputIterator d_first) {
  return Impl::unique_copy_impl(label, ex, first, last, d_first);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto unique_copy(const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType1, Properties1...>& source,
                 const ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return ::Kokkos::Experimental::unique_copy(
      "Kokkos::unique_copy_view_api_default", ex, cbegin(source), cend(source),
      begin(dest));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto unique_copy(const std::string& label, const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType1, Properties1...>& source,
                 const ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return ::Kokkos::Experimental::unique_copy(label, ex, cbegin(source),
                                             cend(source), begin(dest));
}

// overload set2
template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class BinaryPredicate>
OutputIterator unique_copy(const ExecutionSpace& ex, InputIterator first,
                           InputIterator last, OutputIterator d_first,
                           BinaryPredicate pred) {
  return Impl::unique_copy_impl("Kokkos::unique_copy_iterator_api_default", ex,
                                first, last, d_first, pred);
}

template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class BinaryPredicate>
OutputIterator unique_copy(const std::string& label, const ExecutionSpace& ex,
                           InputIterator first, InputIterator last,
                           OutputIterator d_first, BinaryPredicate pred) {
  return Impl::unique_copy_impl(label, ex, first, last, d_first, pred);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class BinaryPredicate>
auto unique_copy(const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType1, Properties1...>& source,
                 const ::Kokkos::View<DataType2, Properties2...>& dest,
                 BinaryPredicate pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::unique_copy_impl("Kokkos::unique_copy_view_api_default", ex,
                                cbegin(source), cend(source), begin(dest),
                                std::move(pred));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class BinaryPredicate>
auto unique_copy(const std::string& label, const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType1, Properties1...>& source,
                 const ::Kokkos::View<DataType2, Properties2...>& dest,
                 BinaryPredicate pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::unique_copy_impl(label, ex, cbegin(source), cend(source),
                                begin(dest), std::move(pred));
}

// -------------------
// rotate
// -------------------

template <class ExecutionSpace, class IteratorType>
IteratorType rotate(const ExecutionSpace& ex, IteratorType first,
                    IteratorType n_first, IteratorType last) {
  return Impl::rotate_impl("Kokkos::rotate_iterator_api_default", ex, first,
                           n_first, last);
}

template <class ExecutionSpace, class IteratorType>
IteratorType rotate(const std::string& label, const ExecutionSpace& ex,
                    IteratorType first, IteratorType n_first,
                    IteratorType last) {
  return Impl::rotate_impl(label, ex, first, n_first, last);
}

template <class ExecutionSpace, class DataType, class... Properties>
auto rotate(const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view,
            std::size_t n_location) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::rotate_impl("Kokkos::rotate_view_api_default", ex, begin(view),
                           begin(view) + n_location, end(view));
}

template <class ExecutionSpace, class DataType, class... Properties>
auto rotate(const std::string& label, const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view,
            std::size_t n_location) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::rotate_impl(label, ex, begin(view), begin(view) + n_location,
                           end(view));
}

// -------------------
// rotate_copy
// -------------------

template <class ExecutionSpace, class InputIterator, class OutputIterator>
OutputIterator rotate_copy(const ExecutionSpace& ex, InputIterator first,
                           InputIterator n_first, InputIterator last,
                           OutputIterator d_first) {
  return Impl::rotate_copy_impl("Kokkos::rotate_copy_iterator_api_default", ex,
                                first, n_first, last, d_first);
}

template <class ExecutionSpace, class InputIterator, class OutputIterator>
OutputIterator rotate_copy(const std::string& label, const ExecutionSpace& ex,
                           InputIterator first, InputIterator n_first,
                           InputIterator last, OutputIterator d_first) {
  return Impl::rotate_copy_impl(label, ex, first, n_first, last, d_first);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto rotate_copy(const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType1, Properties1...>& source,
                 std::size_t n_location,
                 const ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::rotate_copy_impl("Kokkos::rotate_copy_view_api_default", ex,
                                cbegin(source), cbegin(source) + n_location,
                                cend(source), begin(dest));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto rotate_copy(const std::string& label, const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType1, Properties1...>& source,
                 std::size_t n_location,
                 const ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::rotate_copy_impl(label, ex, cbegin(source),
                                cbegin(source) + n_location, cend(source),
                                begin(dest));
}

// -------------------
// remove_if
// -------------------
template <class ExecutionSpace, class Iterator, class UnaryPredicate>
Iterator remove_if(const ExecutionSpace& ex, Iterator first, Iterator last,
                   UnaryPredicate pred) {
  return Impl::remove_if_impl("Kokkos::remove_if_iterator_api_default", ex,
                              first, last, pred);
}

template <class ExecutionSpace, class Iterator, class UnaryPredicate>
Iterator remove_if(const std::string& label, const ExecutionSpace& ex,
                   Iterator first, Iterator last, UnaryPredicate pred) {
  return Impl::remove_if_impl(label, ex, first, last, pred);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class UnaryPredicate>
auto remove_if(const ExecutionSpace& ex,
               const ::Kokkos::View<DataType, Properties...>& view,
               UnaryPredicate pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  return Impl::remove_if_impl("Kokkos::remove_if_iterator_api_default", ex,
                              ::Kokkos::Experimental::begin(view),
                              ::Kokkos::Experimental::end(view), pred);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class UnaryPredicate>
auto remove_if(const std::string& label, const ExecutionSpace& ex,
               const ::Kokkos::View<DataType, Properties...>& view,
               UnaryPredicate pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::remove_if_impl(label, ex, ::Kokkos::Experimental::begin(view),
                              ::Kokkos::Experimental::end(view), pred);
}

// -------------------
// remove
// -------------------
template <class ExecutionSpace, class Iterator, class ValueType>
Iterator remove(const ExecutionSpace& ex, Iterator first, Iterator last,
                const ValueType& value) {
  return Impl::remove_impl("Kokkos::remove_iterator_api_default", ex, first,
                           last, value);
}

template <class ExecutionSpace, class Iterator, class ValueType>
Iterator remove(const std::string& label, const ExecutionSpace& ex,
                Iterator first, Iterator last, const ValueType& value) {
  return Impl::remove_impl(label, ex, first, last, value);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class ValueType>
auto remove(const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view,
            const ValueType& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::remove_impl("Kokkos::remove_iterator_api_default", ex,
                           ::Kokkos::Experimental::begin(view),
                           ::Kokkos::Experimental::end(view), value);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class ValueType>
auto remove(const std::string& label, const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view,
            const ValueType& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::remove_impl(label, ex, ::Kokkos::Experimental::begin(view),
                           ::Kokkos::Experimental::end(view), value);
}

// -------------------
// remove_copy
// -------------------
template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class ValueType>
OutputIterator remove_copy(const ExecutionSpace& ex, InputIterator first_from,
                           InputIterator last_from, OutputIterator first_dest,
                           const ValueType& value) {
  return Impl::remove_copy_impl("Kokkos::remove_copy_iterator_api_default", ex,
                                first_from, last_from, first_dest, value);
}

template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class ValueType>
OutputIterator remove_copy(const std::string& label, const ExecutionSpace& ex,
                           InputIterator first_from, InputIterator last_from,
                           OutputIterator first_dest, const ValueType& value) {
  return Impl::remove_copy_impl(label, ex, first_from, last_from, first_dest,
                                value);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ValueType>
auto remove_copy(const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType1, Properties1...>& view_from,
                 const ::Kokkos::View<DataType2, Properties2...>& view_dest,
                 const ValueType& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);

  return Impl::remove_copy_impl("Kokkos::remove_copy_iterator_api_default", ex,
                                ::Kokkos::Experimental::cbegin(view_from),
                                ::Kokkos::Experimental::cend(view_from),
                                ::Kokkos::Experimental::begin(view_dest),
                                value);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ValueType>
auto remove_copy(const std::string& label, const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType1, Properties1...>& view_from,
                 const ::Kokkos::View<DataType2, Properties2...>& view_dest,
                 const ValueType& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);

  return Impl::remove_copy_impl(
      label, ex, ::Kokkos::Experimental::cbegin(view_from),
      ::Kokkos::Experimental::cend(view_from),
      ::Kokkos::Experimental::begin(view_dest), value);
}

// -------------------
// remove_copy_if
// -------------------
template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class UnaryPredicate>
OutputIterator remove_copy_if(const ExecutionSpace& ex,
                              InputIterator first_from, InputIterator last_from,
                              OutputIterator first_dest,
                              const UnaryPredicate& pred) {
  return Impl::remove_copy_if_impl(
      "Kokkos::remove_copy_if_iterator_api_default", ex, first_from, last_from,
      first_dest, pred);
}

template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class UnaryPredicate>
OutputIterator remove_copy_if(const std::string& label,
                              const ExecutionSpace& ex,
                              InputIterator first_from, InputIterator last_from,
                              OutputIterator first_dest,
                              const UnaryPredicate& pred) {
  return Impl::remove_copy_if_impl(label, ex, first_from, last_from, first_dest,
                                   pred);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class UnaryPredicate>
auto remove_copy_if(const ExecutionSpace& ex,
                    const ::Kokkos::View<DataType1, Properties1...>& view_from,
                    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
                    const UnaryPredicate& pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);

  return Impl::remove_copy_if_impl(
      "Kokkos::remove_copy_if_iterator_api_default", ex,
      ::Kokkos::Experimental::cbegin(view_from),
      ::Kokkos::Experimental::cend(view_from),
      ::Kokkos::Experimental::begin(view_dest), pred);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class UnaryPredicate>
auto remove_copy_if(const std::string& label, const ExecutionSpace& ex,
                    const ::Kokkos::View<DataType1, Properties1...>& view_from,
                    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
                    const UnaryPredicate& pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);

  return Impl::remove_copy_if_impl(
      label, ex, ::Kokkos::Experimental::cbegin(view_from),
      ::Kokkos::Experimental::cend(view_from),
      ::Kokkos::Experimental::begin(view_dest), pred);
}

// -------------------
// shift_left
// -------------------
template <class ExecutionSpace, class IteratorType>
IteratorType shift_left(const ExecutionSpace& ex, IteratorType first,
                        IteratorType last,
                        typename IteratorType::difference_type n) {
  return Impl::shift_left_impl("Kokkos::shift_left_iterator_api_default", ex,
                               first, last, n);
}

template <class ExecutionSpace, class IteratorType>
IteratorType shift_left(const std::string& label, const ExecutionSpace& ex,
                        IteratorType first, IteratorType last,
                        typename IteratorType::difference_type n) {
  return Impl::shift_left_impl(label, ex, first, last, n);
}

template <class ExecutionSpace, class DataType, class... Properties>
auto shift_left(const ExecutionSpace& ex,
                const ::Kokkos::View<DataType, Properties...>& view,
                typename decltype(begin(view))::difference_type n) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::shift_left_impl("Kokkos::shift_left_view_api_default", ex,
                               begin(view), end(view), n);
}

template <class ExecutionSpace, class DataType, class... Properties>
auto shift_left(const std::string& label, const ExecutionSpace& ex,
                const ::Kokkos::View<DataType, Properties...>& view,
                typename decltype(begin(view))::difference_type n) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::shift_left_impl(label, ex, begin(view), end(view), n);
}

// -------------------
// shift_right
// -------------------
template <class ExecutionSpace, class IteratorType>
IteratorType shift_right(const ExecutionSpace& ex, IteratorType first,
                         IteratorType last,
                         typename IteratorType::difference_type n) {
  return Impl::shift_right_impl("Kokkos::shift_right_iterator_api_default", ex,
                                first, last, n);
}

template <class ExecutionSpace, class IteratorType>
IteratorType shift_right(const std::string& label, const ExecutionSpace& ex,
                         IteratorType first, IteratorType last,
                         typename IteratorType::difference_type n) {
  return Impl::shift_right_impl(label, ex, first, last, n);
}

template <class ExecutionSpace, class DataType, class... Properties>
auto shift_right(const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& view,
                 typename decltype(begin(view))::difference_type n) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::shift_right_impl("Kokkos::shift_right_view_api_default", ex,
                                begin(view), end(view), n);
}

template <class ExecutionSpace, class DataType, class... Properties>
auto shift_right(const std::string& label, const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& view,
                 typename decltype(begin(view))::difference_type n) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::shift_right_impl(label, ex, begin(view), end(view), n);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
