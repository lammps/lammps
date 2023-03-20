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

#ifndef KOKKOS_STD_ALGORITHMS_REDUCE_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_REDUCE_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include "Kokkos_ReducerWithArbitraryJoinerNoNeutralElement.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IteratorType, class ValueType>
struct StdReduceDefaultFunctor {
  using index_type = typename IteratorType::difference_type;

  const IteratorType m_first;

  KOKKOS_FUNCTION
  void operator()(const index_type i, ValueType& update) const {
    update += m_first[i];
  }
};

template <class ValueType>
struct StdReduceDefaultJoinFunctor {
  KOKKOS_FUNCTION
  constexpr ValueType operator()(const ValueType& a, const ValueType& b) const {
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

template <typename ValueType>
using has_reduction_identity_sum_t =
    decltype(Kokkos::reduction_identity<ValueType>::sum());

template <class ExecutionSpace, class IteratorType, class ValueType>
ValueType reduce_default_functors_impl(const std::string& label,
                                       const ExecutionSpace& ex,
                                       IteratorType first, IteratorType last,
                                       ValueType init_reduction_value) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::static_assert_is_not_openmptarget(ex);
  Impl::expect_valid_range(first, last);

  using value_type = Kokkos::Impl::remove_cvref_t<ValueType>;

  if (::Kokkos::is_detected<has_reduction_identity_sum_t, value_type>::value) {
    if (first == last) {
      // init is returned, unmodified
      return init_reduction_value;
    }

    using functor_type =
        Impl::StdReduceDefaultFunctor<IteratorType, value_type>;

    // run
    value_type tmp;
    const auto num_elements = Kokkos::Experimental::distance(first, last);
    ::Kokkos::parallel_reduce(label,
                              RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                              functor_type{first}, tmp);
    // fence not needed since reducing into scalar
    tmp += init_reduction_value;
    return tmp;
  } else {
    using joiner_type = Impl::StdReduceDefaultJoinFunctor<value_type>;
    return reduce_custom_functors_impl(
        label, ex, first, last, std::move(init_reduction_value), joiner_type());
  }
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
