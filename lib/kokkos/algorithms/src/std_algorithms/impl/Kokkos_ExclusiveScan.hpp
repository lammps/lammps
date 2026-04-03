// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_STD_ALGORITHMS_EXCLUSIVE_SCAN_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_EXCLUSIVE_SCAN_IMPL_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
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

template <class ValueType>
struct StdExclusiveScanDefaultJoinFunctor {
  KOKKOS_FUNCTION
  constexpr ValueType operator()(const ValueType& a, const ValueType& b) const {
    return a + b;
  }
};

template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class ValueType, class BinaryOpType>
OutputIteratorType exclusive_scan_exespace_impl(
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
  using unary_op_type = StdNumericScanIdentityReferenceUnaryFunctor;
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
  ex.fence("Kokkos::exclusive_scan: fence after operation");

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
  using unary_op_type = StdNumericScanIdentityReferenceUnaryFunctor;
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
