//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

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
    if (input.is_initial) return;

    if (update.is_initial) {
      update.val = input.val;
    } else {
      update.val = update.val + input.val;
    }
    update.is_initial = false;
  }
};

//
// exespace impl
//
template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType>
OutputIteratorType inclusive_scan_default_op_exespace_impl(
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
OutputIteratorType inclusive_scan_custom_binary_op_exespace_impl(
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
  using func_type     = ExeSpaceTransformInclusiveScanNoInitValueFunctor<
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
OutputIteratorType inclusive_scan_custom_binary_op_exespace_impl(
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
  using func_type     = ExeSpaceTransformInclusiveScanWithInitValueFunctor<
      ExecutionSpace, index_type, ValueType, InputIteratorType,
      OutputIteratorType, BinaryOpType, unary_op_type>;

  // run
  const auto num_elements =
      Kokkos::Experimental::distance(first_from, last_from);
  ::Kokkos::parallel_scan(label,
                          RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                          func_type(first_from, first_dest, binary_op,
                                    unary_op_type(), std::move(init_value)));
  ex.fence("Kokkos::inclusive_scan_custom_binary_op: fence after operation");

  // return
  return first_dest + num_elements;
}

//
// team impl
//
template <class TeamHandleType, class InputIteratorType,
          class OutputIteratorType>
KOKKOS_FUNCTION OutputIteratorType inclusive_scan_default_op_team_impl(
    const TeamHandleType& teamHandle, InputIteratorType first_from,
    InputIteratorType last_from, OutputIteratorType first_dest) {
  // checks
  Impl::static_assert_random_access_and_accessible(teamHandle, first_from,
                                                   first_dest);
  Impl::static_assert_iterators_have_matching_difference_type(first_from,
                                                              first_dest);
  Impl::expect_valid_range(first_from, last_from);

  using value_type =
      std::remove_const_t<typename InputIteratorType::value_type>;

  // #if defined(KOKKOS_ENABLE_CUDA)

  using exe_space  = typename TeamHandleType::execution_space;
  using index_type = typename InputIteratorType::difference_type;
  using func_type  = std::conditional_t<
      ::Kokkos::is_detected<in_scan_has_reduction_identity_sum_t,
                            value_type>::value,
      InclusiveScanDefaultFunctorForKnownIdentityElement<
          exe_space, index_type, value_type, InputIteratorType,
          OutputIteratorType>,
      InclusiveScanDefaultFunctor<exe_space, index_type, value_type,
                                  InputIteratorType, OutputIteratorType>>;

  // run
  const auto num_elements =
      Kokkos::Experimental::distance(first_from, last_from);
  ::Kokkos::parallel_scan(TeamThreadRange(teamHandle, 0, num_elements),
                          func_type(first_from, first_dest));
  teamHandle.team_barrier();

  // return
  return first_dest + num_elements;
}

// -------------------------------------------------------------
// inclusive_scan_custom_binary_op_impl
// -------------------------------------------------------------
template <class TeamHandleType, class InputIteratorType,
          class OutputIteratorType, class BinaryOpType>
KOKKOS_FUNCTION OutputIteratorType inclusive_scan_custom_binary_op_team_impl(
    const TeamHandleType& teamHandle, InputIteratorType first_from,
    InputIteratorType last_from, OutputIteratorType first_dest,
    BinaryOpType binary_op) {
  // checks
  Impl::static_assert_random_access_and_accessible(teamHandle, first_from,
                                                   first_dest);
  Impl::static_assert_iterators_have_matching_difference_type(first_from,
                                                              first_dest);
  Impl::expect_valid_range(first_from, last_from);

  using value_type =
      std::remove_const_t<typename InputIteratorType::value_type>;

  static_assert(
      ::Kokkos::is_detected_v<ex_scan_has_reduction_identity_sum_t, value_type>,
      "At the moment inclusive_scan doesn't support types without reduction "
      "identity");

  // #if defined(KOKKOS_ENABLE_CUDA)

  // aliases
  using exe_space     = typename TeamHandleType::execution_space;
  using unary_op_type = StdNumericScanIdentityReferenceUnaryFunctor<value_type>;
  using func_type     = TeamTransformInclusiveScanNoInitValueFunctor<
      exe_space, value_type, InputIteratorType, OutputIteratorType,
      BinaryOpType, unary_op_type>;

  // run
  const auto num_elements =
      Kokkos::Experimental::distance(first_from, last_from);

  ::Kokkos::parallel_scan(
      TeamThreadRange(teamHandle, 0, num_elements),
      func_type(first_from, first_dest, binary_op, unary_op_type()));
  teamHandle.team_barrier();

  return first_dest + num_elements;
}

// -------------------------------------------------------------
// inclusive_scan_custom_binary_op_impl with init_value
// -------------------------------------------------------------
template <class TeamHandleType, class InputIteratorType,
          class OutputIteratorType, class BinaryOpType, class ValueType>
KOKKOS_FUNCTION OutputIteratorType inclusive_scan_custom_binary_op_team_impl(
    const TeamHandleType& teamHandle, InputIteratorType first_from,
    InputIteratorType last_from, OutputIteratorType first_dest,
    BinaryOpType binary_op, ValueType init_value) {
  // checks
  Impl::static_assert_random_access_and_accessible(teamHandle, first_from,
                                                   first_dest);
  Impl::static_assert_iterators_have_matching_difference_type(first_from,
                                                              first_dest);
  Impl::expect_valid_range(first_from, last_from);

  static_assert(
      ::Kokkos::is_detected_v<ex_scan_has_reduction_identity_sum_t, ValueType>,
      "At the moment inclusive_scan doesn't support types without reduction "
      "identity");

  // #if defined(KOKKOS_ENABLE_CUDA)

  // aliases
  using exe_space     = typename TeamHandleType::execution_space;
  using unary_op_type = StdNumericScanIdentityReferenceUnaryFunctor<ValueType>;
  using func_type     = TeamTransformInclusiveScanWithInitValueFunctor<
      exe_space, ValueType, InputIteratorType, OutputIteratorType, BinaryOpType,
      unary_op_type>;

  // run
  const auto num_elements =
      Kokkos::Experimental::distance(first_from, last_from);
  ::Kokkos::parallel_scan(TeamThreadRange(teamHandle, 0, num_elements),
                          func_type(first_from, first_dest, binary_op,
                                    unary_op_type(), std::move(init_value)));
  teamHandle.team_barrier();

  // return
  return first_dest + num_elements;
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
