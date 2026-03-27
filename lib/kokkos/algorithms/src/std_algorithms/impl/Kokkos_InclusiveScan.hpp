// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_STD_ALGORITHMS_INCLUSIVE_SCAN_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_INCLUSIVE_SCAN_IMPL_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include <Kokkos_Profiling_ScopedRegion.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_TransformInclusiveScan.hpp>
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

#if defined(KOKKOS_ENABLE_CUDA)

// Workaround for `Instruction 'shfl' without '.sync' is not supported on
// .target sm_70 and higher from PTX ISA version 6.4`.
// Also see https://github.com/NVIDIA/cub/pull/170.
#if !defined(CUB_USE_COOPERATIVE_GROUPS)
#define CUB_USE_COOPERATIVE_GROUPS
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wsuggest-override"

#if defined(KOKKOS_COMPILER_CLANG)
// Some versions of Clang fail to compile Thrust, failing with errors like
// this:
//    <snip>/thrust/system/cuda/detail/core/agent_launcher.h:557:11:
//    error: use of undeclared identifier 'va_printf'
// The exact combination of versions for Clang and Thrust (or CUDA) for this
// failure was not investigated, however even very recent version combination
// (Clang 10.0.0 and Cuda 10.0) demonstrated failure.
//
// Defining _CubLog here locally allows us to avoid that code path, however
// disabling some debugging diagnostics
#pragma push_macro("_CubLog")
#ifdef _CubLog
#undef _CubLog
#endif
// NOLINTNEXTLINE(bugprone-reserved-identifier)
#define _CubLog
#include <thrust/distance.h>
#include <thrust/scan.h>
#pragma pop_macro("_CubLog")
#else
#include <thrust/distance.h>
#include <thrust/scan.h>
#endif

#pragma GCC diagnostic pop

#endif

#if defined(KOKKOS_ENABLE_ROCTHRUST)
#include <thrust/distance.h>
#include <thrust/scan.h>
#endif

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

// -------------------------------------------------------------
// inclusive_scan_default_op_exespace_impl
// -------------------------------------------------------------

#if defined(KOKKOS_ENABLE_CUDA)
template <class InputIteratorType, class OutputIteratorType>
OutputIteratorType inclusive_scan_default_op_exespace_impl(
    const std::string& label, const Cuda& ex, InputIteratorType first_from,
    InputIteratorType last_from, OutputIteratorType first_dest) {
  const auto thrust_ex = thrust::cuda::par.on(ex.cuda_stream());

  Kokkos::Profiling::pushRegion(label + " via thrust::inclusive_scan");

  thrust::inclusive_scan(thrust_ex, first_from, last_from, first_dest);

  Kokkos::Profiling::popRegion();

  const auto num_elements = thrust::distance(first_from, last_from);

  return first_dest + num_elements;
}
#endif

#if defined(KOKKOS_ENABLE_ROCTHRUST)
template <class InputIteratorType, class OutputIteratorType>
OutputIteratorType inclusive_scan_default_op_exespace_impl(
    const std::string& label, const HIP& ex, InputIteratorType first_from,
    InputIteratorType last_from, OutputIteratorType first_dest) {
  const auto thrust_ex = thrust::hip::par.on(ex.hip_stream());

  Kokkos::Profiling::pushRegion(label + " via thrust::inclusive_scan");

  thrust::inclusive_scan(thrust_ex, first_from, last_from, first_dest);

  Kokkos::Profiling::popRegion();

  const auto num_elements = thrust::distance(first_from, last_from);

  return first_dest + num_elements;
}
#endif

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

  Kokkos::Profiling::pushRegion(label + " via Kokkos::parallel_scan");

  ::Kokkos::parallel_scan(label,
                          RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                          func_type(first_from, first_dest));
  ex.fence("Kokkos::inclusive_scan_default_op: fence after operation");

  Kokkos::Profiling::popRegion();

  // return
  return first_dest + num_elements;
}

// -------------------------------------------------------------
// inclusive_scan_custom_binary_op_impl
// -------------------------------------------------------------

#if defined(KOKKOS_ENABLE_CUDA)
template <class InputIteratorType, class OutputIteratorType, class BinaryOpType>
OutputIteratorType inclusive_scan_custom_binary_op_exespace_impl(
    const std::string& label, const Cuda& ex, InputIteratorType first_from,
    InputIteratorType last_from, OutputIteratorType first_dest,
    BinaryOpType binary_op) {
  const auto thrust_ex = thrust::cuda::par.on(ex.cuda_stream());

  Kokkos::Profiling::pushRegion(label + " via thrust::inclusive_scan");

  thrust::inclusive_scan(thrust_ex, first_from, last_from, first_dest,
                         binary_op);

  Kokkos::Profiling::popRegion();

  const auto num_elements = thrust::distance(first_from, last_from);

  return first_dest + num_elements;
}
#endif

#if defined(KOKKOS_ENABLE_ROCTHRUST)
template <class InputIteratorType, class OutputIteratorType, class BinaryOpType>
OutputIteratorType inclusive_scan_custom_binary_op_exespace_impl(
    const std::string& label, const HIP& ex, InputIteratorType first_from,
    InputIteratorType last_from, OutputIteratorType first_dest,
    BinaryOpType binary_op) {
  const auto thrust_ex = thrust::hip::par.on(ex.hip_stream());

  Kokkos::Profiling::pushRegion(label + " via thrust::inclusive_scan");

  thrust::inclusive_scan(thrust_ex, first_from, last_from, first_dest,
                         binary_op);

  Kokkos::Profiling::popRegion();

  const auto num_elements = thrust::distance(first_from, last_from);

  return first_dest + num_elements;
}
#endif

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
  using unary_op_type = StdNumericScanIdentityReferenceUnaryFunctor;
  using func_type     = ExeSpaceTransformInclusiveScanNoInitValueFunctor<
      ExecutionSpace, index_type, value_type, InputIteratorType,
      OutputIteratorType, BinaryOpType, unary_op_type>;

  // run
  const auto num_elements =
      Kokkos::Experimental::distance(first_from, last_from);

  Kokkos::Profiling::pushRegion(label + " via Kokkos::parallel_scan");

  ::Kokkos::parallel_scan(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements),
      func_type(first_from, first_dest, binary_op, unary_op_type()));
  ex.fence("Kokkos::inclusive_scan_custom_binary_op: fence after operation");

  Kokkos::Profiling::popRegion();

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
  using unary_op_type = StdNumericScanIdentityReferenceUnaryFunctor;
  using func_type     = ExeSpaceTransformInclusiveScanWithInitValueFunctor<
      ExecutionSpace, index_type, ValueType, InputIteratorType,
      OutputIteratorType, BinaryOpType, unary_op_type>;

  // run
  const auto num_elements =
      Kokkos::Experimental::distance(first_from, last_from);

  Kokkos::Profiling::pushRegion(label + " via Kokkos::parallel_scan");

  ::Kokkos::parallel_scan(label,
                          RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                          func_type(first_from, first_dest, binary_op,
                                    unary_op_type(), std::move(init_value)));
  ex.fence("Kokkos::inclusive_scan_custom_binary_op: fence after operation");

  Kokkos::Profiling::popRegion();

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
  using unary_op_type = StdNumericScanIdentityReferenceUnaryFunctor;
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
  using unary_op_type = StdNumericScanIdentityReferenceUnaryFunctor;
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
