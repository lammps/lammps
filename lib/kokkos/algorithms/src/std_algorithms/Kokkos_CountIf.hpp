// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_STD_ALGORITHMS_COUNT_IF_HPP
#define KOKKOS_STD_ALGORITHMS_COUNT_IF_HPP

#include "impl/Kokkos_CountCountIf.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

//
// overload set accepting execution space
//
template <
    typename ExecutionSpace, typename IteratorType, typename Predicate,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
typename IteratorType::difference_type count_if(const ExecutionSpace& ex,
                                                IteratorType first,
                                                IteratorType last,
                                                Predicate predicate) {
  return Impl::count_if_exespace_impl("Kokkos::count_if_iterator_api_default",
                                      ex, first, last, std::move(predicate));
}

template <
    typename ExecutionSpace, typename IteratorType, typename Predicate,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
typename IteratorType::difference_type count_if(const std::string& label,
                                                const ExecutionSpace& ex,
                                                IteratorType first,
                                                IteratorType last,
                                                Predicate predicate) {
  return Impl::count_if_exespace_impl(label, ex, first, last,
                                      std::move(predicate));
}

template <
    typename ExecutionSpace, typename DataType, typename... Properties,
    typename Predicate,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
auto count_if(const ExecutionSpace& ex,
              const ::Kokkos::View<DataType, Properties...>& v,
              Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::count_if_exespace_impl("Kokkos::count_if_view_api_default", ex,
                                      KE::cbegin(v), KE::cend(v),
                                      std::move(predicate));
}

template <
    typename ExecutionSpace, typename DataType, typename... Properties,
    typename Predicate,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
auto count_if(const std::string& label, const ExecutionSpace& ex,
              const ::Kokkos::View<DataType, Properties...>& v,
              Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::count_if_exespace_impl(label, ex, KE::cbegin(v), KE::cend(v),
                                      std::move(predicate));
}

//
// overload set accepting team handle
//
template <typename TeamHandleType, typename IteratorType, typename Predicate,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION typename IteratorType::difference_type count_if(
    const TeamHandleType& teamHandle, IteratorType first, IteratorType last,
    Predicate predicate) {
  return Impl::count_if_team_impl(teamHandle, first, last,
                                  std::move(predicate));
}

template <typename TeamHandleType, typename DataType, typename... Properties,
          typename Predicate,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION auto count_if(const TeamHandleType& teamHandle,
                              const ::Kokkos::View<DataType, Properties...>& v,
                              Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::count_if_team_impl(teamHandle, KE::cbegin(v), KE::cend(v),
                                  std::move(predicate));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
