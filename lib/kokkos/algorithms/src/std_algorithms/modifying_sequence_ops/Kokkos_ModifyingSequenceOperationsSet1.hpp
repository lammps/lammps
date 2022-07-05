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

#ifndef KOKKOS_MODIFYING_SEQUENCE_OPERATIONS_SET1_HPP
#define KOKKOS_MODIFYING_SEQUENCE_OPERATIONS_SET1_HPP

#include <Kokkos_Core.hpp>
#include "../Kokkos_BeginEnd.hpp"
#include "../Kokkos_Constraints.hpp"
#include "../Kokkos_ModifyingOperations.hpp"
#include "../Kokkos_NonModifyingSequenceOperations.hpp"

namespace Kokkos {
namespace Experimental {
namespace Impl {

//---------------------------
//
// functors
//
//---------------------------
template <class IndexType, class InputIterator, class OutputIterator>
struct StdCopyFunctor {
  InputIterator m_first;
  OutputIterator m_dest_first;

  KOKKOS_FUNCTION
  void operator()(IndexType i) const { m_dest_first[i] = m_first[i]; }

  KOKKOS_FUNCTION
  StdCopyFunctor(InputIterator _first, OutputIterator _dest_first)
      : m_first(std::move(_first)), m_dest_first(std::move(_dest_first)) {}
};

template <class IndexType, class IteratorType1, class IteratorType2>
struct StdCopyBackwardFunctor {
  static_assert(std::is_signed<IndexType>::value,
                "Kokkos: StdCopyBackwardFunctor requires signed index type");

  IteratorType1 m_last;
  IteratorType2 m_dest_last;

  KOKKOS_FUNCTION
  void operator()(IndexType i) const { m_dest_last[-i - 1] = m_last[-i - 1]; }

  KOKKOS_FUNCTION
  StdCopyBackwardFunctor(IteratorType1 _last, IteratorType2 _dest_last)
      : m_last(std::move(_last)), m_dest_last(std::move(_dest_last)) {}
};

template <class IndexType, class FirstFrom, class FirstDest, class PredType>
struct StdCopyIfFunctor {
  FirstFrom m_first_from;
  FirstDest m_first_dest;
  PredType m_pred;

  KOKKOS_FUNCTION
  StdCopyIfFunctor(FirstFrom first_from, FirstDest first_dest, PredType pred)
      : m_first_from(std::move(first_from)),
        m_first_dest(std::move(first_dest)),
        m_pred(std::move(pred)) {}

  KOKKOS_FUNCTION
  void operator()(const IndexType i, IndexType& update,
                  const bool final_pass) const {
    const auto& myval = m_first_from[i];
    if (final_pass) {
      if (m_pred(myval)) {
        m_first_dest[update] = myval;
      }
    }

    if (m_pred(myval)) {
      update += 1;
    }
  }
};

template <class InputIterator, class T>
struct StdFillFunctor {
  using index_type = typename InputIterator::difference_type;
  InputIterator m_first;
  T m_value;

  KOKKOS_FUNCTION
  void operator()(index_type i) const { m_first[i] = m_value; }

  KOKKOS_FUNCTION
  StdFillFunctor(InputIterator _first, T _value)
      : m_first(std::move(_first)), m_value(std::move(_value)) {}
};

template <class IndexType, class InputIterator, class OutputIterator,
          class UnaryFunctorType>
struct StdTransformFunctor {
  InputIterator m_first;
  OutputIterator m_d_first;
  UnaryFunctorType m_unary_op;

  KOKKOS_FUNCTION
  void operator()(IndexType i) const { m_d_first[i] = m_unary_op(m_first[i]); }

  KOKKOS_FUNCTION
  StdTransformFunctor(InputIterator _first, OutputIterator _m_d_first,
                      UnaryFunctorType _functor)
      : m_first(std::move(_first)),
        m_d_first(std::move(_m_d_first)),
        m_unary_op(std::move(_functor)) {}
};

template <class IndexType, class InputIterator1, class InputIterator2,
          class OutputIterator, class BinaryFunctorType>
struct StdTransformBinaryFunctor {
  InputIterator1 m_first1;
  InputIterator2 m_first2;
  OutputIterator m_d_first;
  BinaryFunctorType m_binary_op;

  KOKKOS_FUNCTION
  void operator()(IndexType i) const {
    m_d_first[i] = m_binary_op(m_first1[i], m_first2[i]);
  }

  KOKKOS_FUNCTION
  StdTransformBinaryFunctor(InputIterator1 _first1, InputIterator2 _first2,
                            OutputIterator _m_d_first,
                            BinaryFunctorType _functor)
      : m_first1(std::move(_first1)),
        m_first2(std::move(_first2)),
        m_d_first(std::move(_m_d_first)),
        m_binary_op(std::move(_functor)) {}
};

template <class IteratorType, class Generator>
struct StdGenerateFunctor {
  using index_type = typename IteratorType::difference_type;
  IteratorType m_first;
  Generator m_generator;

  KOKKOS_FUNCTION
  void operator()(index_type i) const { m_first[i] = m_generator(); }

  KOKKOS_FUNCTION
  StdGenerateFunctor(IteratorType _first, Generator _g)
      : m_first(std::move(_first)), m_generator(std::move(_g)) {}
};

template <class InputIterator, class PredicateType, class NewValueType>
struct StdReplaceIfFunctor {
  using index_type = typename InputIterator::difference_type;

  InputIterator m_first;
  PredicateType m_predicate;
  NewValueType m_new_value;

  KOKKOS_FUNCTION
  void operator()(index_type i) const {
    if (m_predicate(m_first[i])) {
      m_first[i] = m_new_value;
    }
  }

  KOKKOS_FUNCTION
  StdReplaceIfFunctor(InputIterator first, PredicateType pred,
                      NewValueType new_value)
      : m_first(std::move(first)),
        m_predicate(std::move(pred)),
        m_new_value(std::move(new_value)) {}
};

template <class InputIterator, class ValueType>
struct StdReplaceFunctor {
  using index_type = typename InputIterator::difference_type;
  InputIterator m_first;
  ValueType m_old_value;
  ValueType m_new_value;

  KOKKOS_FUNCTION
  void operator()(index_type i) const {
    if (m_first[i] == m_old_value) {
      m_first[i] = m_new_value;
    }
  }

  KOKKOS_FUNCTION
  StdReplaceFunctor(InputIterator first, ValueType old_value,
                    ValueType new_value)
      : m_first(std::move(first)),
        m_old_value(std::move(old_value)),
        m_new_value(std::move(new_value)) {}
};

template <class InputIterator, class OutputIterator, class ValueType>
struct StdReplaceCopyFunctor {
  using index_type = typename InputIterator::difference_type;

  InputIterator m_first_from;
  OutputIterator m_first_dest;
  ValueType m_old_value;
  ValueType m_new_value;

  KOKKOS_FUNCTION
  void operator()(index_type i) const {
    const auto& myvalue_from = m_first_from[i];

    if (myvalue_from == m_old_value) {
      m_first_dest[i] = m_new_value;
    } else {
      m_first_dest[i] = myvalue_from;
    }
  }

  KOKKOS_FUNCTION
  StdReplaceCopyFunctor(InputIterator first_from, OutputIterator first_dest,
                        ValueType old_value, ValueType new_value)
      : m_first_from(std::move(first_from)),
        m_first_dest(std::move(first_dest)),
        m_old_value(std::move(old_value)),
        m_new_value(std::move(new_value)) {}
};

template <class IndexType, class InputIterator, class OutputIterator,
          class PredicateType, class ValueType>
struct StdReplaceIfCopyFunctor {
  InputIterator m_first_from;
  OutputIterator m_first_dest;
  PredicateType m_pred;
  ValueType m_new_value;

  KOKKOS_FUNCTION
  void operator()(IndexType i) const {
    const auto& myvalue_from = m_first_from[i];

    if (m_pred(myvalue_from)) {
      m_first_dest[i] = m_new_value;
    } else {
      m_first_dest[i] = myvalue_from;
    }
  }

  KOKKOS_FUNCTION
  StdReplaceIfCopyFunctor(InputIterator first_from, OutputIterator first_dest,
                          PredicateType pred, ValueType new_value)
      : m_first_from(std::move(first_from)),
        m_first_dest(std::move(first_dest)),
        m_pred(std::move(pred)),
        m_new_value(std::move(new_value)) {}
};

// ------------------------------------------
// copy_impl
// ------------------------------------------
template <class ExecutionSpace, class InputIterator, class OutputIterator>
OutputIterator copy_impl(const std::string& label, const ExecutionSpace& ex,
                         InputIterator first, InputIterator last,
                         OutputIterator d_first) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first, d_first);
  Impl::static_assert_iterators_have_matching_difference_type(first, d_first);
  Impl::expect_valid_range(first, last);

  // aliases
  using index_type = typename InputIterator::difference_type;
  using func_t     = StdCopyFunctor<index_type, InputIterator, OutputIterator>;

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                         func_t(first, d_first));
  ex.fence("Kokkos::copy: fence after operation");

  // return
  return d_first + num_elements;
}

// ------------------------------------------
// copy_n_impl
// ------------------------------------------
template <class ExecutionSpace, class InputIterator, class Size,
          class OutputIterator>
OutputIterator copy_n_impl(const std::string& label, const ExecutionSpace& ex,
                           InputIterator first_from, Size count,
                           OutputIterator first_dest) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first_from, first_dest);
  Impl::static_assert_iterators_have_matching_difference_type(first_from,
                                                              first_dest);

  if (count > 0) {
    return copy_impl(label, ex, first_from, first_from + count, first_dest);
  } else {
    return first_dest;
  }
}

// ------------------------------------------
// copy_backward_impl
// ------------------------------------------
template <class ExecutionSpace, class IteratorType1, class IteratorType2>
IteratorType2 copy_backward_impl(const std::string& label,
                                 const ExecutionSpace& ex, IteratorType1 first,
                                 IteratorType1 last, IteratorType2 d_last) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first, d_last);
  Impl::static_assert_iterators_have_matching_difference_type(first, d_last);
  Impl::expect_valid_range(first, last);

  // aliases
  using index_type = typename IteratorType1::difference_type;
  using func_t =
      StdCopyBackwardFunctor<index_type, IteratorType1, IteratorType2>;

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                         func_t(last, d_last));
  ex.fence("Kokkos::copy_backward: fence after operation");

  // return
  return d_last - num_elements;
}

// ------------------------------------------
// copy_if_impl
// ------------------------------------------
template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class PredicateType>
OutputIterator copy_if_impl(const std::string& label, const ExecutionSpace& ex,
                            InputIterator first, InputIterator last,
                            OutputIterator d_first, PredicateType pred) {
  /*
    To explain the impl, suppose that our data is:

    | 1 | 1 | 2 | 2 | 3 | -2 | 4 | 4 | 4 | 5 | 7 | -10 |

    and we want to copy only the even entries,
    We can use an exclusive scan where the "update"
    is incremented only for the elements that satisfy the predicate.
    This way, the update allows us to track where in the destination
    we need to copy the elements:

    In this case, counting only the even entries, the exlusive scan
    during the final pass would yield:

    | 0 | 0 | 0 | 1 | 2 | 2 | 3 | 4 | 5 | 6 | 6 | 6 |
              *   *       *   *   *   *           *

    which provides the indexing in the destination where
    each starred (*) element needs to be copied to since
    the starred elements are those that satisfy the predicate.
   */

  // checks
  Impl::static_assert_random_access_and_accessible(ex, first, d_first);
  Impl::static_assert_iterators_have_matching_difference_type(first, d_first);
  Impl::expect_valid_range(first, last);

  if (first == last) {
    return d_first;
  } else {
    // aliases
    using index_type = typename InputIterator::difference_type;
    using func_type  = StdCopyIfFunctor<index_type, InputIterator,
                                       OutputIterator, PredicateType>;

    // run
    const auto num_elements = Kokkos::Experimental::distance(first, last);
    index_type count        = 0;
    ::Kokkos::parallel_scan(label,
                            RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                            func_type(first, d_first, pred), count);

    // fence not needed because of the scan accumulating into count
    return d_first + count;
  }
}

// ------------------------------------------
// fill_impl
// ------------------------------------------
template <class ExecutionSpace, class IteratorType, class T>
void fill_impl(const std::string& label, const ExecutionSpace& ex,
               IteratorType first, IteratorType last, const T& value) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                         StdFillFunctor<IteratorType, T>(first, value));
  ex.fence("Kokkos::fill: fence after operation");
}

template <class ExecutionSpace, class IteratorType, class SizeType, class T>
IteratorType fill_n_impl(const std::string& label, const ExecutionSpace& ex,
                         IteratorType first, SizeType n, const T& value) {
  auto last = first + n;
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  if (n <= 0) {
    return first;
  }

  fill_impl(label, ex, first, last, value);
  return last;
}

// ------------------------------------------
// transform_impl
// ------------------------------------------
template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class UnaryOperation>
OutputIterator transform_impl(const std::string& label,
                              const ExecutionSpace& ex, InputIterator first1,
                              InputIterator last1, OutputIterator d_first,
                              UnaryOperation unary_op) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first1, d_first);
  Impl::static_assert_iterators_have_matching_difference_type(first1, d_first);
  Impl::expect_valid_range(first1, last1);

  // aliases
  using index_type = typename InputIterator::difference_type;
  using func_t = StdTransformFunctor<index_type, InputIterator, OutputIterator,
                                     UnaryOperation>;

  // run
  const auto num_elements = Kokkos::Experimental::distance(first1, last1);
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                         func_t(first1, d_first, unary_op));
  ex.fence("Kokkos::transform: fence after operation");

  // return
  return d_first + num_elements;
}

template <class ExecutionSpace, class InputIterator1, class InputIterator2,
          class OutputIterator, class BinaryOperation>
OutputIterator transform_impl(const std::string& label,
                              const ExecutionSpace& ex, InputIterator1 first1,
                              InputIterator1 last1, InputIterator2 first2,
                              OutputIterator d_first,
                              BinaryOperation binary_op) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first1, first2, d_first);
  Impl::static_assert_iterators_have_matching_difference_type(first1, first2,
                                                              d_first);
  Impl::expect_valid_range(first1, last1);

  // aliases
  using index_type = typename InputIterator1::difference_type;
  using func_t =
      StdTransformBinaryFunctor<index_type, InputIterator1, InputIterator2,
                                OutputIterator, BinaryOperation>;

  // run
  const auto num_elements = Kokkos::Experimental::distance(first1, last1);
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                         func_t(first1, first2, d_first, binary_op));
  ex.fence("Kokkos::transform: fence after operation");
  return d_first + num_elements;
}

// ------------------------------------------
// generate_impl
// ------------------------------------------
template <class ExecutionSpace, class IteratorType, class Generator>
void generate_impl(const std::string& label, const ExecutionSpace& ex,
                   IteratorType first, IteratorType last, Generator g) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  // aliases
  using func_t = StdGenerateFunctor<IteratorType, Generator>;

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                         func_t(first, g));
  ex.fence("Kokkos::generate: fence after operation");
}

template <class ExecutionSpace, class IteratorType, class Size, class Generator>
IteratorType generate_n_impl(const std::string& label, const ExecutionSpace& ex,
                             IteratorType first, Size count, Generator g) {
  if (count <= 0) {
    return first;
  }

  generate_impl(label, ex, first, first + count, g);
  return first + count;
}

// ------------------------------------------
// replace_if_impl
// ------------------------------------------
template <class ExecutionSpace, class IteratorType, class PredicateType,
          class ValueType>
void replace_if_impl(const std::string& label, const ExecutionSpace& ex,
                     IteratorType first, IteratorType last, PredicateType pred,
                     const ValueType& new_value) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  // aliases
  using func_t = StdReplaceIfFunctor<IteratorType, PredicateType, ValueType>;

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                         func_t(first, std::move(pred), new_value));
  ex.fence("Kokkos::replace_if: fence after operation");
}

// ------------------------------------------
// replace_impl
// ------------------------------------------
template <class ExecutionSpace, class IteratorType, class ValueType>
void replace_impl(const std::string& label, const ExecutionSpace& ex,
                  IteratorType first, IteratorType last,
                  const ValueType& old_value, const ValueType& new_value) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  // aliases
  using func_t = StdReplaceFunctor<IteratorType, ValueType>;

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                         func_t(first, old_value, new_value));
  ex.fence("Kokkos::replace: fence after operation");
}

// ------------------------------------------
// replace_copy_impl
// ------------------------------------------
template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class ValueType>
OutputIteratorType replace_copy_impl(const std::string& label,
                                     const ExecutionSpace& ex,
                                     InputIteratorType first_from,
                                     InputIteratorType last_from,
                                     OutputIteratorType first_dest,
                                     const ValueType& old_value,
                                     const ValueType& new_value) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first_from, first_dest);
  Impl::static_assert_iterators_have_matching_difference_type(first_from,
                                                              first_dest);
  Impl::expect_valid_range(first_from, last_from);

  // aliases
  using func_t =
      StdReplaceCopyFunctor<InputIteratorType, OutputIteratorType, ValueType>;

  // run
  const auto num_elements =
      Kokkos::Experimental::distance(first_from, last_from);
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                         func_t(first_from, first_dest, old_value, new_value));
  ex.fence("Kokkos::replace_copy: fence after operation");

  // return
  return first_dest + num_elements;
}

// ------------------------------------------
// replace_copy_if_impl
// ------------------------------------------
template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class PredicateType, class ValueType>
OutputIteratorType replace_copy_if_impl(const std::string& label,
                                        const ExecutionSpace& ex,
                                        InputIteratorType first_from,
                                        InputIteratorType last_from,
                                        OutputIteratorType first_dest,
                                        PredicateType pred,
                                        const ValueType& new_value) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first_from, first_dest);
  Impl::static_assert_iterators_have_matching_difference_type(first_from,
                                                              first_dest);
  Impl::expect_valid_range(first_from, last_from);

  // aliases
  using index_type = typename InputIteratorType::difference_type;
  using func_t =
      StdReplaceIfCopyFunctor<index_type, InputIteratorType, OutputIteratorType,
                              PredicateType, ValueType>;

  // run
  const auto num_elements =
      Kokkos::Experimental::distance(first_from, last_from);
  ::Kokkos::parallel_for(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements),
      func_t(first_from, first_dest, std::move(pred), new_value));
  ex.fence("Kokkos::replace_copy_if: fence after operation");

  // return
  return first_dest + num_elements;
}

}  // namespace Impl

// -------------------
// replace copy
// -------------------
template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class ValueType>
OutputIterator replace_copy(const ExecutionSpace& ex, InputIterator first_from,
                            InputIterator last_from, OutputIterator first_dest,
                            const ValueType& old_value,
                            const ValueType& new_value) {
  return Impl::replace_copy_impl("Kokkos::replace_copy_iterator_api", ex,
                                 first_from, last_from, first_dest, old_value,
                                 new_value);
}

template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class ValueType>
OutputIterator replace_copy(const std::string& label, const ExecutionSpace& ex,
                            InputIterator first_from, InputIterator last_from,
                            OutputIterator first_dest,
                            const ValueType& old_value,
                            const ValueType& new_value) {
  return Impl::replace_copy_impl(label, ex, first_from, last_from, first_dest,
                                 old_value, new_value);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ValueType>
auto replace_copy(const ExecutionSpace& ex,
                  const ::Kokkos::View<DataType1, Properties1...>& view_from,
                  const ::Kokkos::View<DataType2, Properties2...>& view_dest,
                  const ValueType& old_value, const ValueType& new_value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  namespace KE = ::Kokkos::Experimental;
  return Impl::replace_copy_impl("Kokkos::replace_copy_view_api", ex,
                                 KE::cbegin(view_from), KE::cend(view_from),
                                 KE::begin(view_dest), old_value, new_value);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ValueType>
auto replace_copy(const std::string& label, const ExecutionSpace& ex,
                  const ::Kokkos::View<DataType1, Properties1...>& view_from,
                  const ::Kokkos::View<DataType2, Properties2...>& view_dest,
                  const ValueType& old_value, const ValueType& new_value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  namespace KE = ::Kokkos::Experimental;
  return Impl::replace_copy_impl(label, ex, KE::cbegin(view_from),
                                 KE::cend(view_from), KE::begin(view_dest),
                                 old_value, new_value);
}

// -------------------
// replace_copy_if
// -------------------
template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class PredicateType, class ValueType>
OutputIterator replace_copy_if(const ExecutionSpace& ex,
                               InputIterator first_from,
                               InputIterator last_from,
                               OutputIterator first_dest, PredicateType pred,
                               const ValueType& new_value) {
  return Impl::replace_copy_if_impl("Kokkos::replace_copy_if_iterator_api", ex,
                                    first_from, last_from, first_dest, pred,
                                    new_value);
}

template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class PredicateType, class ValueType>
OutputIterator replace_copy_if(const std::string& label,
                               const ExecutionSpace& ex,
                               InputIterator first_from,
                               InputIterator last_from,
                               OutputIterator first_dest, PredicateType pred,
                               const ValueType& new_value) {
  return Impl::replace_copy_if_impl(label, ex, first_from, last_from,
                                    first_dest, pred, new_value);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class PredicateType,
          class ValueType>
auto replace_copy_if(const ExecutionSpace& ex,
                     const ::Kokkos::View<DataType1, Properties1...>& view_from,
                     const ::Kokkos::View<DataType2, Properties2...>& view_dest,
                     PredicateType pred, const ValueType& new_value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  namespace KE = ::Kokkos::Experimental;
  return Impl::replace_copy_if_impl("Kokkos::replace_copy_if_view_api", ex,
                                    KE::cbegin(view_from), KE::cend(view_from),
                                    KE::begin(view_dest), pred, new_value);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class PredicateType,
          class ValueType>
auto replace_copy_if(const std::string& label, const ExecutionSpace& ex,
                     const ::Kokkos::View<DataType1, Properties1...>& view_from,
                     const ::Kokkos::View<DataType2, Properties2...>& view_dest,
                     PredicateType pred, const ValueType& new_value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  namespace KE = ::Kokkos::Experimental;
  return Impl::replace_copy_if_impl(label, ex, KE::cbegin(view_from),
                                    KE::cend(view_from), KE::begin(view_dest),
                                    pred, new_value);
}

// -------------------
// replace
// -------------------
template <class ExecutionSpace, class Iterator, class ValueType>
void replace(const ExecutionSpace& ex, Iterator first, Iterator last,
             const ValueType& old_value, const ValueType& new_value) {
  return Impl::replace_impl("Kokkos::replace_iterator_api", ex, first, last,
                            old_value, new_value);
}

template <class ExecutionSpace, class Iterator, class ValueType>
void replace(const std::string& label, const ExecutionSpace& ex, Iterator first,
             Iterator last, const ValueType& old_value,
             const ValueType& new_value) {
  return Impl::replace_impl(label, ex, first, last, old_value, new_value);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class ValueType>
void replace(const ExecutionSpace& ex,
             const ::Kokkos::View<DataType1, Properties1...>& view,
             const ValueType& old_value, const ValueType& new_value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  namespace KE = ::Kokkos::Experimental;
  return Impl::replace_impl("Kokkos::replace_view_api", ex, KE::begin(view),
                            KE::end(view), old_value, new_value);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class ValueType>
void replace(const std::string& label, const ExecutionSpace& ex,
             const ::Kokkos::View<DataType1, Properties1...>& view,
             const ValueType& old_value, const ValueType& new_value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  namespace KE = ::Kokkos::Experimental;
  return Impl::replace_impl(label, ex, KE::begin(view), KE::end(view),
                            old_value, new_value);
}

// -------------------
// replace_if
// -------------------
template <class ExecutionSpace, class InputIterator, class Predicate,
          class ValueType>
void replace_if(const ExecutionSpace& ex, InputIterator first,
                InputIterator last, Predicate pred,
                const ValueType& new_value) {
  return Impl::replace_if_impl("Kokkos::replace_if_iterator_api", ex, first,
                               last, pred, new_value);
}

template <class ExecutionSpace, class InputIterator, class Predicate,
          class ValueType>
void replace_if(const std::string& label, const ExecutionSpace& ex,
                InputIterator first, InputIterator last, Predicate pred,
                const ValueType& new_value) {
  return Impl::replace_if_impl(label, ex, first, last, pred, new_value);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class Predicate, class ValueType>
void replace_if(const ExecutionSpace& ex,
                const ::Kokkos::View<DataType1, Properties1...>& view,
                Predicate pred, const ValueType& new_value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  namespace KE = ::Kokkos::Experimental;
  return Impl::replace_if_impl("Kokkos::replace_if_view_api", ex,
                               KE::begin(view), KE::end(view), pred, new_value);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class Predicate, class ValueType>
void replace_if(const std::string& label, const ExecutionSpace& ex,
                const ::Kokkos::View<DataType1, Properties1...>& view,
                Predicate pred, const ValueType& new_value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  namespace KE = ::Kokkos::Experimental;
  return Impl::replace_if_impl(label, ex, KE::begin(view), KE::end(view), pred,
                               new_value);
}

// -------------------
// copy
// -------------------
template <class ExecutionSpace, class InputIterator, class OutputIterator>
OutputIterator copy(const ExecutionSpace& ex, InputIterator first,
                    InputIterator last, OutputIterator d_first) {
  return Impl::copy_impl("Kokkos::copy_iterator_api_default", ex, first, last,
                         d_first);
}

template <class ExecutionSpace, class InputIterator, class OutputIterator>
OutputIterator copy(const std::string& label, const ExecutionSpace& ex,
                    InputIterator first, InputIterator last,
                    OutputIterator d_first) {
  return Impl::copy_impl(label, ex, first, last, d_first);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto copy(const ExecutionSpace& ex,
          const ::Kokkos::View<DataType1, Properties1...>& source,
          ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  namespace KE = ::Kokkos::Experimental;
  return Impl::copy_impl("Kokkos::copy_view_api_default", ex,
                         KE::cbegin(source), KE::cend(source), KE::begin(dest));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto copy(const std::string& label, const ExecutionSpace& ex,
          const ::Kokkos::View<DataType1, Properties1...>& source,
          ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  namespace KE = ::Kokkos::Experimental;
  return Impl::copy_impl(label, ex, KE::cbegin(source), KE::cend(source),
                         KE::begin(dest));
}

// -------------------
// copy_n
// -------------------
template <class ExecutionSpace, class InputIterator, class Size,
          class OutputIterator>
OutputIterator copy_n(const ExecutionSpace& ex, InputIterator first, Size count,
                      OutputIterator result) {
  return Impl::copy_n_impl("Kokkos::copy_n_iterator_api_default", ex, first,
                           count, result);
}

template <class ExecutionSpace, class InputIterator, class Size,
          class OutputIterator>
OutputIterator copy_n(const std::string& label, const ExecutionSpace& ex,
                      InputIterator first, Size count, OutputIterator result) {
  return Impl::copy_n_impl(label, ex, first, count, result);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class Size, class DataType2, class... Properties2>
auto copy_n(const ExecutionSpace& ex,
            const ::Kokkos::View<DataType1, Properties1...>& source, Size count,
            ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  namespace KE = ::Kokkos::Experimental;
  return Impl::copy_n_impl("Kokkos::copy_n_view_api_default", ex,
                           KE::cbegin(source), count, KE::begin(dest));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class Size, class DataType2, class... Properties2>
auto copy_n(const std::string& label, const ExecutionSpace& ex,
            const ::Kokkos::View<DataType1, Properties1...>& source, Size count,
            ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  namespace KE = ::Kokkos::Experimental;
  return Impl::copy_n_impl(label, ex, KE::cbegin(source), count,
                           KE::begin(dest));
}

// -------------------
// copy_backward
// -------------------
template <class ExecutionSpace, class IteratorType1, class IteratorType2>
IteratorType2 copy_backward(const ExecutionSpace& ex, IteratorType1 first,
                            IteratorType1 last, IteratorType2 d_last) {
  return Impl::copy_backward_impl("Kokkos::copy_backward_iterator_api_default",
                                  ex, first, last, d_last);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
IteratorType2 copy_backward(const std::string& label, const ExecutionSpace& ex,
                            IteratorType1 first, IteratorType1 last,
                            IteratorType2 d_last) {
  return Impl::copy_backward_impl(label, ex, first, last, d_last);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto copy_backward(const ExecutionSpace& ex,
                   const ::Kokkos::View<DataType1, Properties1...>& source,
                   ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::copy_backward_impl("Kokkos::copy_backward_view_api_default", ex,
                                  cbegin(source), cend(source), end(dest));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto copy_backward(const std::string& label, const ExecutionSpace& ex,
                   const ::Kokkos::View<DataType1, Properties1...>& source,
                   ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::copy_backward_impl(label, ex, cbegin(source), cend(source),
                                  end(dest));
}

// -------------------
// copy_if
// -------------------
template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class Predicate>
OutputIterator copy_if(const ExecutionSpace& ex, InputIterator first,
                       InputIterator last, OutputIterator d_first,
                       Predicate pred) {
  return Impl::copy_if_impl("Kokkos::copy_if_iterator_api_default", ex, first,
                            last, d_first, std::move(pred));
}

template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class Predicate>
OutputIterator copy_if(const std::string& label, const ExecutionSpace& ex,
                       InputIterator first, InputIterator last,
                       OutputIterator d_first, Predicate pred) {
  return Impl::copy_if_impl(label, ex, first, last, d_first, std::move(pred));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class Predicate>
auto copy_if(const ExecutionSpace& ex,
             const ::Kokkos::View<DataType1, Properties1...>& source,
             ::Kokkos::View<DataType2, Properties2...>& dest, Predicate pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::copy_if_impl("Kokkos::copy_if_view_api_default", ex,
                            cbegin(source), cend(source), begin(dest),
                            std::move(pred));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class Predicate>
auto copy_if(const std::string& label, const ExecutionSpace& ex,
             const ::Kokkos::View<DataType1, Properties1...>& source,
             ::Kokkos::View<DataType2, Properties2...>& dest, Predicate pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::copy_if_impl(label, ex, cbegin(source), cend(source),
                            begin(dest), std::move(pred));
}

// -------------------
// fill
// -------------------
template <class ExecutionSpace, class IteratorType, class T>
void fill(const ExecutionSpace& ex, IteratorType first, IteratorType last,
          const T& value) {
  Impl::fill_impl("Kokkos::fill_iterator_api_default", ex, first, last, value);
}

template <class ExecutionSpace, class IteratorType, class T>
void fill(const std::string& label, const ExecutionSpace& ex,
          IteratorType first, IteratorType last, const T& value) {
  Impl::fill_impl(label, ex, first, last, value);
}

template <class ExecutionSpace, class DataType, class... Properties, class T>
void fill(const ExecutionSpace& ex,
          const ::Kokkos::View<DataType, Properties...>& view, const T& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  Impl::fill_impl("Kokkos::fill_view_api_default", ex, begin(view), end(view),
                  value);
}

template <class ExecutionSpace, class DataType, class... Properties, class T>
void fill(const std::string& label, const ExecutionSpace& ex,
          const ::Kokkos::View<DataType, Properties...>& view, const T& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  Impl::fill_impl(label, ex, begin(view), end(view), value);
}

// -------------------
// fill_n
// -------------------
template <class ExecutionSpace, class IteratorType, class SizeType, class T>
IteratorType fill_n(const ExecutionSpace& ex, IteratorType first, SizeType n,
                    const T& value) {
  return Impl::fill_n_impl("Kokkos::fill_n_iterator_api_default", ex, first, n,
                           value);
}

template <class ExecutionSpace, class IteratorType, class SizeType, class T>
IteratorType fill_n(const std::string& label, const ExecutionSpace& ex,
                    IteratorType first, SizeType n, const T& value) {
  return Impl::fill_n_impl(label, ex, first, n, value);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class SizeType, class T>
auto fill_n(const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view, SizeType n,
            const T& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  return Impl::fill_n_impl("Kokkos::fill_n_view_api_default", ex, begin(view),
                           n, value);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class SizeType, class T>
auto fill_n(const std::string& label, const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view, SizeType n,
            const T& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  return Impl::fill_n_impl(label, ex, begin(view), n, value);
}

// -------------------
// transform
// -------------------
template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class UnaryOperation>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      InputIterator, OutputIterator>::value,
                  OutputIterator>
transform(const ExecutionSpace& ex, InputIterator first1, InputIterator last1,
          OutputIterator d_first, UnaryOperation unary_op) {
  return Impl::transform_impl("Kokkos::transform_iterator_api_default", ex,
                              first1, last1, d_first, std::move(unary_op));
}

template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class UnaryOperation>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      InputIterator, OutputIterator>::value,
                  OutputIterator>
transform(const std::string& label, const ExecutionSpace& ex,
          InputIterator first1, InputIterator last1, OutputIterator d_first,
          UnaryOperation unary_op) {
  return Impl::transform_impl(label, ex, first1, last1, d_first,
                              std::move(unary_op));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class UnaryOperation>
auto transform(const ExecutionSpace& ex,
               const ::Kokkos::View<DataType1, Properties1...>& source,
               ::Kokkos::View<DataType2, Properties2...>& dest,
               UnaryOperation unary_op) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::transform_impl("Kokkos::transform_view_api_default", ex,
                              begin(source), end(source), begin(dest),
                              std::move(unary_op));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class UnaryOperation>
auto transform(const std::string& label, const ExecutionSpace& ex,
               const ::Kokkos::View<DataType1, Properties1...>& source,
               ::Kokkos::View<DataType2, Properties2...>& dest,
               UnaryOperation unary_op) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::transform_impl(label, ex, begin(source), end(source),
                              begin(dest), std::move(unary_op));
}

template <class ExecutionSpace, class InputIterator1, class InputIterator2,
          class OutputIterator, class BinaryOperation>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      InputIterator1, InputIterator2, OutputIterator>::value,
                  OutputIterator>
transform(const ExecutionSpace& ex, InputIterator1 first1, InputIterator1 last1,
          InputIterator2 first2, OutputIterator d_first,
          BinaryOperation binary_op) {
  return Impl::transform_impl("Kokkos::transform_iterator_api_default", ex,
                              first1, last1, first2, d_first,
                              std::move(binary_op));
}

template <class ExecutionSpace, class InputIterator1, class InputIterator2,
          class OutputIterator, class BinaryOperation>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      InputIterator1, InputIterator2, OutputIterator>::value,
                  OutputIterator>
transform(const std::string& label, const ExecutionSpace& ex,
          InputIterator1 first1, InputIterator1 last1, InputIterator2 first2,
          OutputIterator d_first, BinaryOperation binary_op) {
  return Impl::transform_impl(label, ex, first1, last1, first2, d_first,
                              std::move(binary_op));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class DataType3,
          class... Properties3, class BinaryOperation>
auto transform(const ExecutionSpace& ex,
               const ::Kokkos::View<DataType1, Properties1...>& source1,
               const ::Kokkos::View<DataType2, Properties2...>& source2,
               ::Kokkos::View<DataType3, Properties3...>& dest,
               BinaryOperation binary_op) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source2);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::transform_impl("Kokkos::transform_view_api_default", ex,
                              begin(source1), end(source1), begin(source2),
                              begin(dest), std::move(binary_op));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class DataType3,
          class... Properties3, class BinaryOperation>
auto transform(const std::string& label, const ExecutionSpace& ex,
               const ::Kokkos::View<DataType1, Properties1...>& source1,
               const ::Kokkos::View<DataType2, Properties2...>& source2,
               ::Kokkos::View<DataType3, Properties3...>& dest,
               BinaryOperation binary_op) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source2);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::transform_impl(label, ex, begin(source1), end(source1),
                              begin(source2), begin(dest),
                              std::move(binary_op));
}

// -------------------
// generate
// -------------------
template <class ExecutionSpace, class IteratorType, class Generator>
void generate(const ExecutionSpace& ex, IteratorType first, IteratorType last,
              Generator g) {
  Impl::generate_impl("Kokkos::generate_iterator_api_default", ex, first, last,
                      std::move(g));
}

template <class ExecutionSpace, class IteratorType, class Generator>
void generate(const std::string& label, const ExecutionSpace& ex,
              IteratorType first, IteratorType last, Generator g) {
  Impl::generate_impl(label, ex, first, last, std::move(g));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class Generator>
void generate(const ExecutionSpace& ex,
              const ::Kokkos::View<DataType, Properties...>& view,
              Generator g) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  Impl::generate_impl("Kokkos::generate_view_api_default", ex, begin(view),
                      end(view), std::move(g));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class Generator>
void generate(const std::string& label, const ExecutionSpace& ex,
              const ::Kokkos::View<DataType, Properties...>& view,
              Generator g) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  Impl::generate_impl(label, ex, begin(view), end(view), std::move(g));
}

// -------------------
// generate_n
// -------------------
template <class ExecutionSpace, class IteratorType, class Size, class Generator>
IteratorType generate_n(const ExecutionSpace& ex, IteratorType first,
                        Size count, Generator g) {
  Impl::generate_n_impl("Kokkos::generate_n_iterator_api_default", ex, first,
                        count, std::move(g));
  return first + count;
}

template <class ExecutionSpace, class IteratorType, class Size, class Generator>
IteratorType generate_n(const std::string& label, const ExecutionSpace& ex,
                        IteratorType first, Size count, Generator g) {
  Impl::generate_n_impl(label, ex, first, count, std::move(g));
  return first + count;
}

template <class ExecutionSpace, class DataType, class... Properties, class Size,
          class Generator>
auto generate_n(const ExecutionSpace& ex,
                const ::Kokkos::View<DataType, Properties...>& view, Size count,
                Generator g) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  return Impl::generate_n_impl("Kokkos::generate_n_view_api_default", ex,
                               begin(view), count, std::move(g));
}

template <class ExecutionSpace, class DataType, class... Properties, class Size,
          class Generator>
auto generate_n(const std::string& label, const ExecutionSpace& ex,
                const ::Kokkos::View<DataType, Properties...>& view, Size count,
                Generator g) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  return Impl::generate_n_impl(label, ex, begin(view), count, std::move(g));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
