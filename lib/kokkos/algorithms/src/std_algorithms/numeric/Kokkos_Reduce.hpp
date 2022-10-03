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

#ifndef KOKKOS_STD_NUMERICS_REDUCE_HPP
#define KOKKOS_STD_NUMERICS_REDUCE_HPP

#include <Kokkos_Core.hpp>
#include "../Kokkos_BeginEnd.hpp"
#include "../Kokkos_Constraints.hpp"
#include "../Kokkos_Distance.hpp"
#include "../Kokkos_ModifyingOperations.hpp"
#include "../Kokkos_ReducerWithArbitraryJoinerNoNeutralElement.hpp"

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class ValueType>
struct StdReduceDefaultJoinFunctor {
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

template <class IteratorType, class ReducerType>
struct StdReduceFunctor {
  using red_value_type = typename ReducerType::value_type;
  using index_type     = typename IteratorType::difference_type;

  const IteratorType m_first;
  const ReducerType m_reducer;

  KOKKOS_FUNCTION
  void operator()(const index_type i, red_value_type& red_value) const {
    auto tmp_wrapped_value = red_value_type{m_first[i], false};

    if (red_value.is_initial) {
      red_value = tmp_wrapped_value;
    } else {
      m_reducer.join(red_value, tmp_wrapped_value);
    }
  }

  KOKKOS_FUNCTION
  StdReduceFunctor(IteratorType first, ReducerType reducer)
      : m_first(std::move(first)), m_reducer(std::move(reducer)) {}
};

//------------------------------
// reduce_custom_functors_impl
//------------------------------
template <class ExecutionSpace, class IteratorType, class ValueType,
          class JoinerType>
ValueType reduce_custom_functors_impl(const std::string& label,
                                      const ExecutionSpace& ex,
                                      IteratorType first, IteratorType last,
                                      ValueType init_reduction_value,
                                      JoinerType joiner) {
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
  using functor_type         = StdReduceFunctor<IteratorType, reducer_type>;
  using reduction_value_type = typename reducer_type::value_type;

  // run
  reduction_value_type result;
  reducer_type reducer(result, joiner);
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_reduce(label,
                            RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                            functor_type(first, reducer), reducer);

  // fence not needed since reducing into scalar
  return joiner(result.val, init_reduction_value);
}

template <class ExecutionSpace, class IteratorType, class ValueType>
ValueType reduce_default_functors_impl(const std::string& label,
                                       const ExecutionSpace& ex,
                                       IteratorType first, IteratorType last,
                                       ValueType init_reduction_value) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::static_assert_is_not_openmptarget(ex);
  Impl::expect_valid_range(first, last);

  using value_type  = Kokkos::Impl::remove_cvref_t<ValueType>;
  using joiner_type = Impl::StdReduceDefaultJoinFunctor<value_type>;
  return reduce_custom_functors_impl(
      label, ex, first, last, std::move(init_reduction_value), joiner_type());
}

}  // end namespace Impl

///////////////////////////////
//
// reduce public API
//
///////////////////////////////

//
// overload set 1
//
template <class ExecutionSpace, class IteratorType>
typename IteratorType::value_type reduce(const ExecutionSpace& ex,
                                         IteratorType first,
                                         IteratorType last) {
  return Impl::reduce_default_functors_impl(
      "Kokkos::reduce_default_functors_iterator_api", ex, first, last,
      typename IteratorType::value_type());
}

template <class ExecutionSpace, class IteratorType>
typename IteratorType::value_type reduce(const std::string& label,
                                         const ExecutionSpace& ex,
                                         IteratorType first,
                                         IteratorType last) {
  return Impl::reduce_default_functors_impl(
      label, ex, first, last, typename IteratorType::value_type());
}

template <class ExecutionSpace, class DataType, class... Properties>
auto reduce(const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view) {
  namespace KE = ::Kokkos::Experimental;
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  using view_type  = ::Kokkos::View<DataType, Properties...>;
  using value_type = typename view_type::value_type;

  return Impl::reduce_default_functors_impl(
      "Kokkos::reduce_default_functors_view_api", ex, KE::cbegin(view),
      KE::cend(view), value_type());
}

template <class ExecutionSpace, class DataType, class... Properties>
auto reduce(const std::string& label, const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view) {
  namespace KE = ::Kokkos::Experimental;
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  using view_type  = ::Kokkos::View<DataType, Properties...>;
  using value_type = typename view_type::value_type;

  return Impl::reduce_default_functors_impl(label, ex, KE::cbegin(view),
                                            KE::cend(view), value_type());
}

//
// overload set2:
//
template <class ExecutionSpace, class IteratorType, class ValueType>
ValueType reduce(const ExecutionSpace& ex, IteratorType first,
                 IteratorType last, ValueType init_reduction_value) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  return Impl::reduce_default_functors_impl(
      "Kokkos::reduce_default_functors_iterator_api", ex, first, last,
      init_reduction_value);
}

template <class ExecutionSpace, class IteratorType, class ValueType>
ValueType reduce(const std::string& label, const ExecutionSpace& ex,
                 IteratorType first, IteratorType last,
                 ValueType init_reduction_value) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  return Impl::reduce_default_functors_impl(label, ex, first, last,
                                            init_reduction_value);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class ValueType>
ValueType reduce(const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& view,
                 ValueType init_reduction_value) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  namespace KE = ::Kokkos::Experimental;
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  return Impl::reduce_default_functors_impl(
      "Kokkos::reduce_default_functors_view_api", ex, KE::cbegin(view),
      KE::cend(view), init_reduction_value);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class ValueType>
ValueType reduce(const std::string& label, const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& view,
                 ValueType init_reduction_value) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  namespace KE = ::Kokkos::Experimental;
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  return Impl::reduce_default_functors_impl(
      label, ex, KE::cbegin(view), KE::cend(view), init_reduction_value);
}

//
// overload set 3
//
template <class ExecutionSpace, class IteratorType, class ValueType,
          class BinaryOp>
ValueType reduce(const ExecutionSpace& ex, IteratorType first,
                 IteratorType last, ValueType init_reduction_value,
                 BinaryOp joiner) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  return Impl::reduce_custom_functors_impl(
      "Kokkos::reduce_default_functors_iterator_api", ex, first, last,
      init_reduction_value, joiner);
}

template <class ExecutionSpace, class IteratorType, class ValueType,
          class BinaryOp>
ValueType reduce(const std::string& label, const ExecutionSpace& ex,
                 IteratorType first, IteratorType last,
                 ValueType init_reduction_value, BinaryOp joiner) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  return Impl::reduce_custom_functors_impl(label, ex, first, last,
                                           init_reduction_value, joiner);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class ValueType, class BinaryOp>
ValueType reduce(const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& view,
                 ValueType init_reduction_value, BinaryOp joiner) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  namespace KE = ::Kokkos::Experimental;
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  return Impl::reduce_custom_functors_impl(
      "Kokkos::reduce_custom_functors_view_api", ex, KE::cbegin(view),
      KE::cend(view), init_reduction_value, joiner);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class ValueType, class BinaryOp>
ValueType reduce(const std::string& label, const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& view,
                 ValueType init_reduction_value, BinaryOp joiner) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  namespace KE = ::Kokkos::Experimental;
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  return Impl::reduce_custom_functors_impl(label, ex, KE::cbegin(view),
                                           KE::cend(view), init_reduction_value,
                                           joiner);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
