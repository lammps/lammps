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

#ifndef KOKKOS_SORT_PUBLIC_API_HPP_
#define KOKKOS_SORT_PUBLIC_API_HPP_

#include "./impl/Kokkos_SortImpl.hpp"
#include <std_algorithms/Kokkos_BeginEnd.hpp>
#include <Kokkos_Core.hpp>
#include <algorithm>

namespace Kokkos {

// ---------------------------------------------------------------
// basic overloads
// ---------------------------------------------------------------

template <class ExecutionSpace, class DataType, class... Properties>
void sort([[maybe_unused]] const ExecutionSpace& exec,
          const Kokkos::View<DataType, Properties...>& view) {
  // constraints
  using ViewType = Kokkos::View<DataType, Properties...>;
  using MemSpace = typename ViewType::memory_space;
  static_assert(
      ViewType::rank == 1 &&
          (std::is_same_v<typename ViewType::array_layout, LayoutRight> ||
           std::is_same_v<typename ViewType::array_layout, LayoutLeft> ||
           std::is_same_v<typename ViewType::array_layout, LayoutStride>),
      "Kokkos::sort without comparator: supports 1D Views with LayoutRight, "
      "LayoutLeft or LayoutStride.");

  static_assert(SpaceAccessibility<ExecutionSpace, MemSpace>::accessible,
                "Kokkos::sort: execution space instance is not able to access "
                "the memory space of the "
                "View argument!");

  if (view.extent(0) <= 1) {
    return;
  }

  if constexpr (Impl::better_off_calling_std_sort_v<ExecutionSpace>) {
    auto first = ::Kokkos::Experimental::begin(view);
    auto last  = ::Kokkos::Experimental::end(view);
    std::sort(first, last);
  } else {
    Impl::sort_device_view_without_comparator(exec, view);
  }
}

template <class DataType, class... Properties>
void sort(const Kokkos::View<DataType, Properties...>& view) {
  using ViewType = Kokkos::View<DataType, Properties...>;
  static_assert(ViewType::rank == 1,
                "Kokkos::sort: currently only supports rank-1 Views.");

  Kokkos::fence("Kokkos::sort: before");

  if (view.extent(0) <= 1) {
    return;
  }

  typename ViewType::execution_space exec;
  sort(exec, view);
  exec.fence("Kokkos::sort: fence after sorting");
}

// ---------------------------------------------------------------
// overloads supporting a custom comparator
// ---------------------------------------------------------------
template <class ExecutionSpace, class ComparatorType, class DataType,
          class... Properties>
void sort([[maybe_unused]] const ExecutionSpace& exec,
          const Kokkos::View<DataType, Properties...>& view,
          const ComparatorType& comparator) {
  // constraints
  using ViewType = Kokkos::View<DataType, Properties...>;
  using MemSpace = typename ViewType::memory_space;
  static_assert(
      ViewType::rank == 1 &&
          (std::is_same_v<typename ViewType::array_layout, LayoutRight> ||
           std::is_same_v<typename ViewType::array_layout, LayoutLeft> ||
           std::is_same_v<typename ViewType::array_layout, LayoutStride>),
      "Kokkos::sort with comparator: supports 1D Views with LayoutRight, "
      "LayoutLeft or LayoutStride.");

  static_assert(SpaceAccessibility<ExecutionSpace, MemSpace>::accessible,
                "Kokkos::sort: execution space instance is not able to access "
                "the memory space of the View argument!");

  if (view.extent(0) <= 1) {
    return;
  }

  if constexpr (Impl::better_off_calling_std_sort_v<ExecutionSpace>) {
    auto first = ::Kokkos::Experimental::begin(view);
    auto last  = ::Kokkos::Experimental::end(view);
    std::sort(first, last, comparator);
  } else {
    Impl::sort_device_view_with_comparator(exec, view, comparator);
  }
}

template <class ComparatorType, class DataType, class... Properties>
void sort(const Kokkos::View<DataType, Properties...>& view,
          const ComparatorType& comparator) {
  using ViewType = Kokkos::View<DataType, Properties...>;
  static_assert(
      ViewType::rank == 1 &&
          (std::is_same_v<typename ViewType::array_layout, LayoutRight> ||
           std::is_same_v<typename ViewType::array_layout, LayoutLeft> ||
           std::is_same_v<typename ViewType::array_layout, LayoutStride>),
      "Kokkos::sort with comparator: supports 1D Views with LayoutRight, "
      "LayoutLeft or LayoutStride.");

  Kokkos::fence("Kokkos::sort with comparator: before");

  if (view.extent(0) <= 1) {
    return;
  }

  typename ViewType::execution_space exec;
  sort(exec, view, comparator);
  exec.fence("Kokkos::sort with comparator: fence after sorting");
}

// ---------------------------------------------------------------
// overloads for sorting a view with a subrange
// specified via integers begin, end
// ---------------------------------------------------------------

template <class ExecutionSpace, class ViewType>
std::enable_if_t<Kokkos::is_execution_space<ExecutionSpace>::value> sort(
    const ExecutionSpace& exec, ViewType view, size_t const begin,
    size_t const end) {
  // view must be rank-1 because the Impl::min_max_functor
  // used below only works for rank-1 views for now
  static_assert(ViewType::rank == 1,
                "Kokkos::sort: currently only supports rank-1 Views.");

  if (view.extent(0) <= 1) {
    return;
  }

  using range_policy = Kokkos::RangePolicy<typename ViewType::execution_space>;
  using CompType     = BinOp1D<ViewType>;

  Kokkos::MinMaxScalar<typename ViewType::non_const_value_type> result;
  Kokkos::MinMax<typename ViewType::non_const_value_type> reducer(result);

  parallel_reduce("Kokkos::Sort::FindExtent", range_policy(exec, begin, end),
                  Impl::min_max_functor<ViewType>(view), reducer);

  if (result.min_val == result.max_val) return;

  BinSort<ViewType, CompType> bin_sort(
      exec, view, begin, end,
      CompType((end - begin) / 2, result.min_val, result.max_val), true);

  bin_sort.create_permute_vector(exec);
  bin_sort.sort(exec, view, begin, end);
}

template <class ViewType>
void sort(ViewType view, size_t const begin, size_t const end) {
  // same constraints as the overload above which this gets dispatched to
  static_assert(ViewType::rank == 1,
                "Kokkos::sort: currently only supports rank-1 Views.");

  Kokkos::fence("Kokkos::sort: before");

  if (view.extent(0) <= 1) {
    return;
  }

  typename ViewType::execution_space exec;
  sort(exec, view, begin, end);
  exec.fence("Kokkos::Sort: fence after sorting");
}

}  // namespace Kokkos
#endif
