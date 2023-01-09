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

#ifndef KOKKOS_STD_ALGORITHMS_INCLUSIVE_SCAN_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_INCLUSIVE_SCAN_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_TransformInclusiveScan.hpp>
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <typename ValueType>
using in_scan_has_reduction_identity_sum_t =
    decltype(Kokkos::reduction_identity<ValueType>::sum());

template <class ExeSpace, class IndexType, class ValueType, class FirstFrom,
          class FirstDest>
struct InclusiveScanDefaultFunctorForKnownIdentityElement {
  using execution_space = ExeSpace;

  FirstFrom m_first_from;
  FirstDest m_first_dest;

  KOKKOS_FUNCTION
  InclusiveScanDefaultFunctorForKnownIdentityElement(FirstFrom first_from,
                                                     FirstDest first_dest)
      : m_first_from(std::move(first_from)),
        m_first_dest(std::move(first_dest)) {}

  KOKKOS_FUNCTION
  void operator()(const IndexType i, ValueType& update,
                  const bool final_pass) const {
    update += m_first_from[i];

    if (final_pass) {
      m_first_dest[i] = update;
    }
  }
};

template <class ExeSpace, class IndexType, class ValueType, class FirstFrom,
          class FirstDest>
struct InclusiveScanDefaultFunctor {
  using execution_space = ExeSpace;
  using value_type      = ValueWrapperForNoNeutralElement<ValueType>;

  FirstFrom m_first_from;
  FirstDest m_first_dest;

  KOKKOS_FUNCTION
  InclusiveScanDefaultFunctor(FirstFrom first_from, FirstDest first_dest)
      : m_first_from(std::move(first_from)),
        m_first_dest(std::move(first_dest)) {}

  KOKKOS_FUNCTION
  void operator()(const IndexType i, value_type& update,
                  const bool final_pass) const {
    const auto tmp = value_type{m_first_from[i], false};
    this->join(update, tmp);

    if (final_pass) {
      m_first_dest[i] = update.val;
    }
  }

  KOKKOS_FUNCTION
  void init(value_type& update) const {
    update.val        = {};
    update.is_initial = true;
  }

  KOKKOS_FUNCTION
  void join(value_type& update, const value_type& input) const {
    if (update.is_initial) {
      update.val = input.val;
    } else {
      update.val = update.val + input.val;
    }
    update.is_initial = false;
  }
};

template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType>
OutputIteratorType inclusive_scan_default_op_impl(
    const std::string& label, const ExecutionSpace& ex,
    InputIteratorType first_from, InputIteratorType last_from,
    OutputIteratorType first_dest) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first_from, first_dest);
  Impl::static_assert_iterators_have_matching_difference_type(first_from,
                                                              first_dest);
  Impl::expect_valid_range(first_from, last_from);

  // aliases
  using index_type = typename InputIteratorType::difference_type;
  using value_type =
      std::remove_const_t<typename InputIteratorType::value_type>;
  using func_type = std::conditional_t<
      ::Kokkos::is_detected<in_scan_has_reduction_identity_sum_t,
                            value_type>::value,
      InclusiveScanDefaultFunctorForKnownIdentityElement<
          ExecutionSpace, index_type, value_type, InputIteratorType,
          OutputIteratorType>,
      InclusiveScanDefaultFunctor<ExecutionSpace, index_type, value_type,
                                  InputIteratorType, OutputIteratorType>>;

  // run
  const auto num_elements =
      Kokkos::Experimental::distance(first_from, last_from);
  ::Kokkos::parallel_scan(label,
                          RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                          func_type(first_from, first_dest));
  ex.fence("Kokkos::inclusive_scan_default_op: fence after operation");

  // return
  return first_dest + num_elements;
}

// -------------------------------------------------------------
// inclusive_scan_custom_binary_op_impl
// -------------------------------------------------------------
template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class BinaryOpType>
OutputIteratorType inclusive_scan_custom_binary_op_impl(
    const std::string& label, const ExecutionSpace& ex,
    InputIteratorType first_from, InputIteratorType last_from,
    OutputIteratorType first_dest, BinaryOpType binary_op) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first_from, first_dest);
  Impl::static_assert_iterators_have_matching_difference_type(first_from,
                                                              first_dest);
  Impl::expect_valid_range(first_from, last_from);

  // aliases
  using index_type = typename InputIteratorType::difference_type;
  using value_type =
      std::remove_const_t<typename InputIteratorType::value_type>;
  using unary_op_type = StdNumericScanIdentityReferenceUnaryFunctor<value_type>;
  using func_type     = TransformInclusiveScanNoInitValueFunctor<
      ExecutionSpace, index_type, value_type, InputIteratorType,
      OutputIteratorType, BinaryOpType, unary_op_type>;

  // run
  const auto num_elements =
      Kokkos::Experimental::distance(first_from, last_from);
  ::Kokkos::parallel_scan(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements),
      func_type(first_from, first_dest, binary_op, unary_op_type()));
  ex.fence("Kokkos::inclusive_scan_custom_binary_op: fence after operation");

  // return
  return first_dest + num_elements;
}

// -------------------------------------------------------------
// inclusive_scan_custom_binary_op_impl with init_value
// -------------------------------------------------------------
template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class BinaryOpType, class ValueType>
OutputIteratorType inclusive_scan_custom_binary_op_impl(
    const std::string& label, const ExecutionSpace& ex,
    InputIteratorType first_from, InputIteratorType last_from,
    OutputIteratorType first_dest, BinaryOpType binary_op,
    ValueType init_value) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first_from, first_dest);
  Impl::static_assert_iterators_have_matching_difference_type(first_from,
                                                              first_dest);
  Impl::expect_valid_range(first_from, last_from);

  // aliases
  using index_type    = typename InputIteratorType::difference_type;
  using unary_op_type = StdNumericScanIdentityReferenceUnaryFunctor<ValueType>;
  using func_type     = TransformInclusiveScanWithInitValueFunctor<
      ExecutionSpace, index_type, ValueType, InputIteratorType,
      OutputIteratorType, BinaryOpType, unary_op_type>;

  // run
  const auto num_elements =
      Kokkos::Experimental::distance(first_from, last_from);
  ::Kokkos::parallel_scan(label,
                          RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                          func_type(first_from, first_dest, binary_op,
                                    unary_op_type(), init_value));
  ex.fence("Kokkos::inclusive_scan_custom_binary_op: fence after operation");

  // return
  return first_dest + num_elements;
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
