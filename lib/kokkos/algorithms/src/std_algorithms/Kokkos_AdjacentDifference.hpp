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

#ifndef KOKKOS_STD_ALGORITHMS_ADJACENT_DIFFERENCE_HPP
#define KOKKOS_STD_ALGORITHMS_ADJACENT_DIFFERENCE_HPP

#include "impl/Kokkos_AdjacentDifference.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

//
// overload set accepting execution space
//
template <
    typename ExecutionSpace, typename InputIteratorType,
    typename OutputIteratorType,
    std::enable_if_t<!::Kokkos::is_view<InputIteratorType>::value &&
                         ::Kokkos::is_execution_space<ExecutionSpace>::value,
                     int> = 0>
OutputIteratorType adjacent_difference(const ExecutionSpace& ex,
                                       InputIteratorType first_from,
                                       InputIteratorType last_from,
                                       OutputIteratorType first_dest) {
  using value_type1 = typename InputIteratorType::value_type;
  using value_type2 = typename OutputIteratorType::value_type;
  using binary_op =
      Impl::StdAdjacentDifferenceDefaultBinaryOpFunctor<value_type1,
                                                        value_type2>;
  return Impl::adjacent_difference_exespace_impl(
      "Kokkos::adjacent_difference_iterator_api", ex, first_from, last_from,
      first_dest, binary_op());
}

template <
    typename ExecutionSpace, typename InputIteratorType,
    typename OutputIteratorType, typename BinaryOp,
    std::enable_if_t<!::Kokkos::is_view<InputIteratorType>::value &&
                         ::Kokkos::is_execution_space<ExecutionSpace>::value,
                     int> = 0>
OutputIteratorType adjacent_difference(const ExecutionSpace& ex,
                                       InputIteratorType first_from,
                                       InputIteratorType last_from,
                                       OutputIteratorType first_dest,
                                       BinaryOp bin_op) {
  return Impl::adjacent_difference_exespace_impl(
      "Kokkos::adjacent_difference_iterator_api", ex, first_from, last_from,
      first_dest, bin_op);
}

template <
    typename ExecutionSpace, typename InputIteratorType,
    typename OutputIteratorType,
    std::enable_if_t<!::Kokkos::is_view<InputIteratorType>::value &&
                         ::Kokkos::is_execution_space<ExecutionSpace>::value,
                     int> = 0>
OutputIteratorType adjacent_difference(const std::string& label,
                                       const ExecutionSpace& ex,
                                       InputIteratorType first_from,
                                       InputIteratorType last_from,
                                       OutputIteratorType first_dest) {
  using value_type1 = typename InputIteratorType::value_type;
  using value_type2 = typename OutputIteratorType::value_type;
  using binary_op =
      Impl::StdAdjacentDifferenceDefaultBinaryOpFunctor<value_type1,
                                                        value_type2>;
  return Impl::adjacent_difference_exespace_impl(
      label, ex, first_from, last_from, first_dest, binary_op());
}

template <
    typename ExecutionSpace, typename InputIteratorType,
    typename OutputIteratorType, typename BinaryOp,
    std::enable_if_t<!::Kokkos::is_view<InputIteratorType>::value &&
                         ::Kokkos::is_execution_space<ExecutionSpace>::value,
                     int> = 0>
OutputIteratorType adjacent_difference(const std::string& label,
                                       const ExecutionSpace& ex,
                                       InputIteratorType first_from,
                                       InputIteratorType last_from,
                                       OutputIteratorType first_dest,
                                       BinaryOp bin_op) {
  return Impl::adjacent_difference_exespace_impl(label, ex, first_from,
                                                 last_from, first_dest, bin_op);
}

template <typename ExecutionSpace, typename DataType1, typename... Properties1,
          typename DataType2, typename... Properties2,
          std::enable_if_t<::Kokkos::is_execution_space<ExecutionSpace>::value,
                           int> = 0>
auto adjacent_difference(
    const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest) {
  namespace KE = ::Kokkos::Experimental;
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);

  using view_type1  = ::Kokkos::View<DataType1, Properties1...>;
  using view_type2  = ::Kokkos::View<DataType2, Properties2...>;
  using value_type1 = typename view_type1::value_type;
  using value_type2 = typename view_type2::value_type;
  using binary_op =
      Impl::StdAdjacentDifferenceDefaultBinaryOpFunctor<value_type1,
                                                        value_type2>;
  return Impl::adjacent_difference_exespace_impl(
      "Kokkos::adjacent_difference_view_api", ex, KE::cbegin(view_from),
      KE::cend(view_from), KE::begin(view_dest), binary_op());
}

template <typename ExecutionSpace, typename DataType1, typename... Properties1,
          typename DataType2, typename... Properties2, typename BinaryOp,
          std::enable_if_t<::Kokkos::is_execution_space<ExecutionSpace>::value,
                           int> = 0>
auto adjacent_difference(
    const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
    BinaryOp bin_op) {
  namespace KE = ::Kokkos::Experimental;
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  return Impl::adjacent_difference_exespace_impl(
      "Kokkos::adjacent_difference_view_api", ex, KE::cbegin(view_from),
      KE::cend(view_from), KE::begin(view_dest), bin_op);
}

template <typename ExecutionSpace, typename DataType1, typename... Properties1,
          typename DataType2, typename... Properties2,
          std::enable_if_t<::Kokkos::is_execution_space<ExecutionSpace>::value,
                           int> = 0>
auto adjacent_difference(
    const std::string& label, const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest) {
  namespace KE = ::Kokkos::Experimental;
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);

  using view_type1  = ::Kokkos::View<DataType1, Properties1...>;
  using view_type2  = ::Kokkos::View<DataType2, Properties2...>;
  using value_type1 = typename view_type1::value_type;
  using value_type2 = typename view_type2::value_type;
  using binary_op =
      Impl::StdAdjacentDifferenceDefaultBinaryOpFunctor<value_type1,
                                                        value_type2>;

  return Impl::adjacent_difference_exespace_impl(
      label, ex, KE::cbegin(view_from), KE::cend(view_from),
      KE::begin(view_dest), binary_op());
}

template <typename ExecutionSpace, typename DataType1, typename... Properties1,
          typename DataType2, typename... Properties2, typename BinaryOp,
          std::enable_if_t<::Kokkos::is_execution_space<ExecutionSpace>::value,
                           int> = 0>
auto adjacent_difference(
    const std::string& label, const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
    BinaryOp bin_op) {
  namespace KE = ::Kokkos::Experimental;
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  return Impl::adjacent_difference_exespace_impl(
      label, ex, KE::cbegin(view_from), KE::cend(view_from),
      KE::begin(view_dest), bin_op);
}

//
// overload set accepting a team handle
// Note: for now omit the overloads accepting a label
// since they cause issues on device because of the string allocation.
//
template <typename TeamHandleType, typename InputIteratorType,
          typename OutputIteratorType,
          std::enable_if_t<!::Kokkos::is_view<InputIteratorType>::value &&
                               ::Kokkos::is_team_handle<TeamHandleType>::value,
                           int> = 0>
KOKKOS_FUNCTION OutputIteratorType adjacent_difference(
    const TeamHandleType& teamHandle, InputIteratorType first_from,
    InputIteratorType last_from, OutputIteratorType first_dest) {
  using value_type1 = typename InputIteratorType::value_type;
  using value_type2 = typename OutputIteratorType::value_type;
  using binary_op =
      Impl::StdAdjacentDifferenceDefaultBinaryOpFunctor<value_type1,
                                                        value_type2>;
  return Impl::adjacent_difference_team_impl(teamHandle, first_from, last_from,
                                             first_dest, binary_op());
}

template <typename TeamHandleType, typename InputIteratorType,
          typename OutputIteratorType, typename BinaryOp,
          std::enable_if_t<!::Kokkos::is_view<InputIteratorType>::value &&
                               ::Kokkos::is_team_handle<TeamHandleType>::value,
                           int> = 0>
KOKKOS_FUNCTION OutputIteratorType
adjacent_difference(const TeamHandleType& teamHandle,
                    InputIteratorType first_from, InputIteratorType last_from,
                    OutputIteratorType first_dest, BinaryOp bin_op) {
  return Impl::adjacent_difference_team_impl(teamHandle, first_from, last_from,
                                             first_dest, bin_op);
}

template <
    typename TeamHandleType, typename DataType1, typename... Properties1,
    typename DataType2, typename... Properties2,
    std::enable_if_t<::Kokkos::is_team_handle<TeamHandleType>::value, int> = 0>
KOKKOS_FUNCTION auto adjacent_difference(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest) {
  namespace KE = ::Kokkos::Experimental;
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);

  using view_type1  = ::Kokkos::View<DataType1, Properties1...>;
  using view_type2  = ::Kokkos::View<DataType2, Properties2...>;
  using value_type1 = typename view_type1::value_type;
  using value_type2 = typename view_type2::value_type;
  using binary_op =
      Impl::StdAdjacentDifferenceDefaultBinaryOpFunctor<value_type1,
                                                        value_type2>;
  return Impl::adjacent_difference_team_impl(teamHandle, KE::cbegin(view_from),
                                             KE::cend(view_from),
                                             KE::begin(view_dest), binary_op());
}

template <
    typename TeamHandleType, typename DataType1, typename... Properties1,
    typename DataType2, typename... Properties2, typename BinaryOp,
    std::enable_if_t<::Kokkos::is_team_handle<TeamHandleType>::value, int> = 0>
KOKKOS_FUNCTION auto adjacent_difference(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
    BinaryOp bin_op) {
  namespace KE = ::Kokkos::Experimental;
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  return Impl::adjacent_difference_team_impl(teamHandle, KE::cbegin(view_from),
                                             KE::cend(view_from),
                                             KE::begin(view_dest), bin_op);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
