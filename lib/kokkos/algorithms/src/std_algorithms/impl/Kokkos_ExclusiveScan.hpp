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

#ifndef KOKKOS_STD_ALGORITHMS_EXCLUSIVE_SCAN_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_EXCLUSIVE_SCAN_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include "Kokkos_ValueWrapperForNoNeutralElement.hpp"
#include "Kokkos_IdentityReferenceUnaryFunctor.hpp"
#include "Kokkos_FunctorsForExclusiveScan.hpp"
#include <std_algorithms/Kokkos_TransformExclusiveScan.hpp>
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

//
// exespace impl
//
template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class ValueType>
OutputIteratorType exclusive_scan_default_op_exespace_impl(
    const std::string& label, const ExecutionSpace& ex,
    InputIteratorType first_from, InputIteratorType last_from,
    OutputIteratorType first_dest, ValueType init_value) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first_from, first_dest);
  Impl::static_assert_iterators_have_matching_difference_type(first_from,
                                                              first_dest);
  Impl::expect_valid_range(first_from, last_from);

  // does it make sense to do this static_assert too?
  // using input_iterator_value_type = typename InputIteratorType::value_type;
  // static_assert
  //   (std::is_convertible<std::remove_cv_t<input_iterator_value_type>,
  //   ValueType>::value,
  //    "exclusive_scan: InputIteratorType::value_type not convertible to
  //    ValueType");

  // we are unnecessarily duplicating code, but this is on purpose
  // so that we can use the default_op for OpenMPTarget.
  // Originally, I had this implemented as:
  // '''
  // using bop_type   = StdExclusiveScanDefaultJoinFunctor<ValueType>;
  // call exclusive_scan_custom_op_impl(..., bop_type());
  // '''
  // which avoids duplicating the functors, but for OpenMPTarget
  // I cannot use a custom binary op.
  // This is the same problem that occurs for reductions.

  // aliases
  using index_type = typename InputIteratorType::difference_type;
  using func_type  = std::conditional_t<
      ::Kokkos::is_detected<ex_scan_has_reduction_identity_sum_t,
                            ValueType>::value,
      ExclusiveScanDefaultFunctorForKnownNeutralElement<
          ExecutionSpace, index_type, ValueType, InputIteratorType,
          OutputIteratorType>,
      ExclusiveScanDefaultFunctorWithValueWrapper<ExecutionSpace, index_type,
                                                  ValueType, InputIteratorType,
                                                  OutputIteratorType>>;

  // run
  const auto num_elements =
      Kokkos::Experimental::distance(first_from, last_from);
  ::Kokkos::parallel_scan(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements),
      func_type(std::move(init_value), first_from, first_dest));

  ex.fence("Kokkos::exclusive_scan_default_op: fence after operation");

  return first_dest + num_elements;
}

template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class ValueType, class BinaryOpType>
OutputIteratorType exclusive_scan_custom_op_exespace_impl(
    const std::string& label, const ExecutionSpace& ex,
    InputIteratorType first_from, InputIteratorType last_from,
    OutputIteratorType first_dest, ValueType init_value, BinaryOpType bop) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first_from, first_dest);
  Impl::static_assert_iterators_have_matching_difference_type(first_from,
                                                              first_dest);
  Impl::expect_valid_range(first_from, last_from);

  // aliases
  using index_type    = typename InputIteratorType::difference_type;
  using unary_op_type = StdNumericScanIdentityReferenceUnaryFunctor<ValueType>;
  using func_type     = TransformExclusiveScanFunctorWithValueWrapper<
      ExecutionSpace, index_type, ValueType, InputIteratorType,
      OutputIteratorType, BinaryOpType, unary_op_type>;

  // run
  const auto num_elements =
      Kokkos::Experimental::distance(first_from, last_from);
  ::Kokkos::parallel_scan(label,
                          RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                          func_type(std::move(init_value), first_from,
                                    first_dest, bop, unary_op_type()));
  ex.fence("Kokkos::exclusive_scan_custom_op: fence after operation");

  // return
  return first_dest + num_elements;
}

//
// team impl
//
template <class TeamHandleType, class InputIteratorType,
          class OutputIteratorType, class ValueType>
KOKKOS_FUNCTION OutputIteratorType exclusive_scan_default_op_team_impl(
    const TeamHandleType& teamHandle, InputIteratorType first_from,
    InputIteratorType last_from, OutputIteratorType first_dest,
    ValueType init_value) {
  // checks
  Impl::static_assert_random_access_and_accessible(teamHandle, first_from,
                                                   first_dest);
  Impl::static_assert_iterators_have_matching_difference_type(first_from,
                                                              first_dest);
  Impl::expect_valid_range(first_from, last_from);

  static_assert(
      ::Kokkos::is_detected_v<ex_scan_has_reduction_identity_sum_t, ValueType>,
      "The team-level impl of Kokkos::Experimental::exclusive_scan currently "
      "does not support types without reduction identity");

  // aliases
  using exe_space  = typename TeamHandleType::execution_space;
  using index_type = typename InputIteratorType::difference_type;
  using func_type  = ExclusiveScanDefaultFunctorForKnownNeutralElement<
      exe_space, index_type, ValueType, InputIteratorType, OutputIteratorType>;

  const auto num_elements =
      Kokkos::Experimental::distance(first_from, last_from);
  ::Kokkos::parallel_scan(
      TeamThreadRange(teamHandle, 0, num_elements),
      func_type(std::move(init_value), first_from, first_dest));
  teamHandle.team_barrier();
  return first_dest + num_elements;
}

template <class TeamHandleType, class InputIteratorType,
          class OutputIteratorType, class ValueType, class BinaryOpType>
KOKKOS_FUNCTION OutputIteratorType exclusive_scan_custom_op_team_impl(
    const TeamHandleType& teamHandle, InputIteratorType first_from,
    InputIteratorType last_from, OutputIteratorType first_dest,
    ValueType init_value, BinaryOpType bop) {
  // checks
  Impl::static_assert_random_access_and_accessible(teamHandle, first_from,
                                                   first_dest);
  Impl::static_assert_iterators_have_matching_difference_type(first_from,
                                                              first_dest);
  Impl::expect_valid_range(first_from, last_from);

  static_assert(
      ::Kokkos::is_detected_v<ex_scan_has_reduction_identity_sum_t, ValueType>,
      "The team-level impl of Kokkos::Experimental::exclusive_scan currently "
      "does not support types without reduction identity");

  // aliases
  using exe_space     = typename TeamHandleType::execution_space;
  using unary_op_type = StdNumericScanIdentityReferenceUnaryFunctor<ValueType>;
  using index_type    = typename InputIteratorType::difference_type;
  using func_type     = TransformExclusiveScanFunctorWithoutValueWrapper<
      exe_space, index_type, ValueType, InputIteratorType, OutputIteratorType,
      BinaryOpType, unary_op_type>;

  const auto num_elements =
      Kokkos::Experimental::distance(first_from, last_from);
  ::Kokkos::parallel_scan(TeamThreadRange(teamHandle, 0, num_elements),
                          func_type(std::move(init_value), first_from,
                                    first_dest, bop, unary_op_type()));
  teamHandle.team_barrier();

  return first_dest + num_elements;
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
