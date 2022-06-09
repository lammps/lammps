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

#ifndef KOKKOS_STD_NUMERICS_TRANSFORM_REDUCE_HPP
#define KOKKOS_STD_NUMERICS_TRANSFORM_REDUCE_HPP

#include <Kokkos_Core.hpp>
#include "../Kokkos_Constraints.hpp"
#include "../Kokkos_Distance.hpp"
#include "../Kokkos_ModifyingOperations.hpp"
#include "../Kokkos_BeginEnd.hpp"
#include "../Kokkos_ReducerWithArbitraryJoinerNoNeutralElement.hpp"

namespace Kokkos {
namespace Experimental {
namespace Impl {

//
// helper functors
//
template <class ValueType>
struct StdTranformReduceDefaultBinaryTransformFunctor {
  KOKKOS_FUNCTION
  constexpr ValueType operator()(const ValueType& a, const ValueType& b) const {
    return (a * b);
  }
};

template <class ValueType>
struct StdTranformReduceDefaultJoinFunctor {
  KOKKOS_FUNCTION
  constexpr ValueType operator()(const ValueType& a, const ValueType& b) const {
    return a + b;
  }

  KOKKOS_FUNCTION
  constexpr ValueType operator()(const volatile ValueType& a,
                                 const volatile ValueType& b) const {
    return a + b;
  }
};

template <class IteratorType, class ReducerType, class TransformType>
struct StdTransformReduceSingleIntervalFunctor {
  using red_value_type = typename ReducerType::value_type;
  using index_type     = typename IteratorType::difference_type;

  const IteratorType m_first;
  const ReducerType m_reducer;
  const TransformType m_transform;

  KOKKOS_FUNCTION
  void operator()(const index_type i, red_value_type& red_value) const {
    auto tmp_wrapped_value = red_value_type{m_transform(m_first[i]), false};
    if (red_value.is_initial) {
      red_value = tmp_wrapped_value;
    } else {
      m_reducer.join(red_value, tmp_wrapped_value);
    }
  }

  KOKKOS_FUNCTION
  StdTransformReduceSingleIntervalFunctor(IteratorType first,
                                          ReducerType reducer,
                                          TransformType transform)
      : m_first(std::move(first)),
        m_reducer(std::move(reducer)),
        m_transform(std::move(transform)) {}
};

template <class IndexType, class IteratorType1, class IteratorType2,
          class ReducerType, class TransformType>
struct StdTransformReduceTwoIntervalsFunctor {
  using red_value_type = typename ReducerType::value_type;

  const IteratorType1 m_first1;
  const IteratorType2 m_first2;
  const ReducerType m_reducer;
  const TransformType m_transform;

  KOKKOS_FUNCTION
  void operator()(const IndexType i, red_value_type& red_value) const {
    auto tmp_wrapped_value =
        red_value_type{m_transform(m_first1[i], m_first2[i]), false};

    if (red_value.is_initial) {
      red_value = tmp_wrapped_value;
    } else {
      m_reducer.join(red_value, tmp_wrapped_value);
    }
  }

  KOKKOS_FUNCTION
  StdTransformReduceTwoIntervalsFunctor(IteratorType1 first1,
                                        IteratorType2 first2,
                                        ReducerType reducer,
                                        TransformType transform)
      : m_first1(std::move(first1)),
        m_first2(std::move(first2)),
        m_reducer(std::move(reducer)),
        m_transform(std::move(transform)) {}
};

//------------------------------
//
// impl functions
//
//------------------------------

template <class ExecutionSpace, class IteratorType, class ValueType,
          class JoinerType, class UnaryTransformerType>
ValueType transform_reduce_custom_functors_impl(
    const std::string& label, const ExecutionSpace& ex, IteratorType first,
    IteratorType last, ValueType init_reduction_value, JoinerType joiner,
    UnaryTransformerType transformer) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::static_assert_is_not_openmptarget(ex);
  Impl::expect_valid_range(first, last);

  if (first == last) {
    // init is returned, unmodified
    return init_reduction_value;
  }

  // aliases
  using reducer_type =
      ReducerWithArbitraryJoinerNoNeutralElement<ValueType, JoinerType>;
  using functor_type =
      StdTransformReduceSingleIntervalFunctor<IteratorType, reducer_type,
                                              UnaryTransformerType>;
  using reduction_value_type = typename reducer_type::value_type;

  // run
  reduction_value_type result;
  reducer_type reducer(result, joiner);
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_reduce(label,
                            RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                            functor_type(first, reducer, transformer), reducer);

  // fence not needed since reducing into scalar

  // as per standard, transform is not applied to the init value
  // https://en.cppreference.com/w/cpp/algorithm/transform_reduce
  return joiner(result.val, init_reduction_value);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class ValueType, class JoinerType, class BinaryTransformerType>
ValueType transform_reduce_custom_functors_impl(
    const std::string& label, const ExecutionSpace& ex, IteratorType1 first1,
    IteratorType1 last1, IteratorType2 first2, ValueType init_reduction_value,
    JoinerType joiner, BinaryTransformerType transformer) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first1, first2);
  Impl::static_assert_is_not_openmptarget(ex);
  Impl::static_assert_iterators_have_matching_difference_type(first1, first2);
  Impl::expect_valid_range(first1, last1);

  if (first1 == last1) {
    // init is returned, unmodified
    return init_reduction_value;
  }

  // aliases
  using index_type = typename IteratorType1::difference_type;
  using reducer_type =
      ReducerWithArbitraryJoinerNoNeutralElement<ValueType, JoinerType>;
  using functor_type =
      StdTransformReduceTwoIntervalsFunctor<index_type, IteratorType1,
                                            IteratorType2, reducer_type,
                                            BinaryTransformerType>;
  using reduction_value_type = typename reducer_type::value_type;

  // run
  reduction_value_type result;
  reducer_type reducer(result, joiner);

  const auto num_elements = Kokkos::Experimental::distance(first1, last1);
  ::Kokkos::parallel_reduce(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements),
      functor_type(first1, first2, reducer, transformer), reducer);

  // fence not needed since reducing into scalar
  return joiner(result.val, init_reduction_value);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class ValueType>
ValueType transform_reduce_default_functors_impl(
    const std::string& label, const ExecutionSpace& ex, IteratorType1 first1,
    IteratorType1 last1, IteratorType2 first2, ValueType init_reduction_value) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first1, first2);
  Impl::static_assert_is_not_openmptarget(ex);
  Impl::static_assert_iterators_have_matching_difference_type(first1, first2);
  Impl::expect_valid_range(first1, last1);

  // aliases
  using transformer_type =
      Impl::StdTranformReduceDefaultBinaryTransformFunctor<ValueType>;
  using joiner_type = Impl::StdTranformReduceDefaultJoinFunctor<ValueType>;

  return transform_reduce_custom_functors_impl(
      label, ex, first1, last1, first2, std::move(init_reduction_value),
      joiner_type(), transformer_type());
}

}  // end namespace Impl

///////////////////////////////
//
// transform_reduce public API
//
///////////////////////////////

// ----------------------------
// overload set1:
// no custom functors passed, so equivalent to
// transform_reduce(first1, last1, first2, init, plus<>(), multiplies<>());
// ----------------------------
template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class ValueType>
ValueType transform_reduce(const ExecutionSpace& ex, IteratorType1 first1,
                           IteratorType1 last1, IteratorType2 first2,
                           ValueType init_reduction_value) {
  return Impl::transform_reduce_default_functors_impl(
      "Kokkos::transform_reduce_default_functors_iterator_api", ex, first1,
      last1, first2, std::move(init_reduction_value));
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class ValueType>
ValueType transform_reduce(const std::string& label, const ExecutionSpace& ex,
                           IteratorType1 first1, IteratorType1 last1,
                           IteratorType2 first2,
                           ValueType init_reduction_value) {
  return Impl::transform_reduce_default_functors_impl(
      label, ex, first1, last1, first2, std::move(init_reduction_value));
}

// overload1 accepting views
template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ValueType>
ValueType transform_reduce(
    const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& first_view,
    const ::Kokkos::View<DataType2, Properties2...>& second_view,
    ValueType init_reduction_value) {
  namespace KE = ::Kokkos::Experimental;
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(first_view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(second_view);

  return Impl::transform_reduce_default_functors_impl(
      "Kokkos::transform_reduce_default_functors_iterator_api", ex,
      KE::cbegin(first_view), KE::cend(first_view), KE::cbegin(second_view),
      std::move(init_reduction_value));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ValueType>
ValueType transform_reduce(
    const std::string& label, const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& first_view,
    const ::Kokkos::View<DataType2, Properties2...>& second_view,
    ValueType init_reduction_value) {
  namespace KE = ::Kokkos::Experimental;
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(first_view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(second_view);

  return Impl::transform_reduce_default_functors_impl(
      label, ex, KE::cbegin(first_view), KE::cend(first_view),
      KE::cbegin(second_view), std::move(init_reduction_value));
}

//
// overload set2:
// accepts a custom transform and joiner functor
//

// Note the std refers to the arg BinaryReductionOp
// but in the Kokkos naming convention, it corresponds
// to a "joiner" that knows how to join two values
// NOTE: "joiner/transformer" need to be commutative.

// https://en.cppreference.com/w/cpp/algorithm/transform_reduce

// api accepting iterators
template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class ValueType, class BinaryJoinerType, class BinaryTransform>
ValueType transform_reduce(const ExecutionSpace& ex, IteratorType1 first1,
                           IteratorType1 last1, IteratorType2 first2,
                           ValueType init_reduction_value,
                           BinaryJoinerType joiner,
                           BinaryTransform transformer) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  return Impl::transform_reduce_custom_functors_impl(
      "Kokkos::transform_reduce_custom_functors_iterator_api", ex, first1,
      last1, first2, std::move(init_reduction_value), std::move(joiner),
      std::move(transformer));
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class ValueType, class BinaryJoinerType, class BinaryTransform>
ValueType transform_reduce(const std::string& label, const ExecutionSpace& ex,
                           IteratorType1 first1, IteratorType1 last1,
                           IteratorType2 first2, ValueType init_reduction_value,
                           BinaryJoinerType joiner,
                           BinaryTransform transformer) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  return Impl::transform_reduce_custom_functors_impl(
      label, ex, first1, last1, first2, std::move(init_reduction_value),
      std::move(joiner), std::move(transformer));
}

// accepting views
template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ValueType,
          class BinaryJoinerType, class BinaryTransform>
ValueType transform_reduce(
    const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& first_view,
    const ::Kokkos::View<DataType2, Properties2...>& second_view,
    ValueType init_reduction_value, BinaryJoinerType joiner,
    BinaryTransform transformer) {
  namespace KE = ::Kokkos::Experimental;
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(first_view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(second_view);

  return Impl::transform_reduce_custom_functors_impl(
      "Kokkos::transform_reduce_custom_functors_view_api", ex,
      KE::cbegin(first_view), KE::cend(first_view), KE::cbegin(second_view),
      std::move(init_reduction_value), std::move(joiner),
      std::move(transformer));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ValueType,
          class BinaryJoinerType, class BinaryTransform>
ValueType transform_reduce(
    const std::string& label, const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& first_view,
    const ::Kokkos::View<DataType2, Properties2...>& second_view,
    ValueType init_reduction_value, BinaryJoinerType joiner,
    BinaryTransform transformer) {
  namespace KE = ::Kokkos::Experimental;
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(first_view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(second_view);

  return Impl::transform_reduce_custom_functors_impl(
      label, ex, KE::cbegin(first_view), KE::cend(first_view),
      KE::cbegin(second_view), std::move(init_reduction_value),
      std::move(joiner), std::move(transformer));
}

//
// overload set3:
//
// accepting iterators
template <class ExecutionSpace, class IteratorType, class ValueType,
          class BinaryJoinerType, class UnaryTransform>
// need this to avoid ambiguous call
std::enable_if_t<
    ::Kokkos::Experimental::Impl::are_iterators<IteratorType>::value, ValueType>
transform_reduce(const ExecutionSpace& ex, IteratorType first1,
                 IteratorType last1, ValueType init_reduction_value,
                 BinaryJoinerType joiner, UnaryTransform transformer) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  return Impl::transform_reduce_custom_functors_impl(
      "Kokkos::transform_reduce_custom_functors_iterator_api", ex, first1,
      last1, std::move(init_reduction_value), std::move(joiner),
      std::move(transformer));
}

template <class ExecutionSpace, class IteratorType, class ValueType,
          class BinaryJoinerType, class UnaryTransform>
// need this to avoid ambiguous call
std::enable_if_t<
    ::Kokkos::Experimental::Impl::are_iterators<IteratorType>::value, ValueType>
transform_reduce(const std::string& label, const ExecutionSpace& ex,
                 IteratorType first1, IteratorType last1,
                 ValueType init_reduction_value, BinaryJoinerType joiner,
                 UnaryTransform transformer) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  return Impl::transform_reduce_custom_functors_impl(
      label, ex, first1, last1, std::move(init_reduction_value),
      std::move(joiner), std::move(transformer));
}

// accepting views
template <class ExecutionSpace, class DataType, class... Properties,
          class ValueType, class BinaryJoinerType, class UnaryTransform>
ValueType transform_reduce(const ExecutionSpace& ex,
                           const ::Kokkos::View<DataType, Properties...>& view,
                           ValueType init_reduction_value,
                           BinaryJoinerType joiner,
                           UnaryTransform transformer) {
  namespace KE = ::Kokkos::Experimental;
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  return Impl::transform_reduce_custom_functors_impl(
      "Kokkos::transform_reduce_custom_functors_view_api", ex, KE::cbegin(view),
      KE::cend(view), std::move(init_reduction_value), std::move(joiner),
      std::move(transformer));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class ValueType, class BinaryJoinerType, class UnaryTransform>
ValueType transform_reduce(const std::string& label, const ExecutionSpace& ex,
                           const ::Kokkos::View<DataType, Properties...>& view,
                           ValueType init_reduction_value,
                           BinaryJoinerType joiner,
                           UnaryTransform transformer) {
  namespace KE = ::Kokkos::Experimental;
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  return Impl::transform_reduce_custom_functors_impl(
      label, ex, KE::cbegin(view), KE::cend(view),
      std::move(init_reduction_value), std::move(joiner),
      std::move(transformer));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
