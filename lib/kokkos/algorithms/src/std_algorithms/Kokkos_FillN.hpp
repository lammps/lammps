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

#ifndef KOKKOS_STD_ALGORITHMS_FILL_N_HPP
#define KOKKOS_STD_ALGORITHMS_FILL_N_HPP

#include "impl/Kokkos_FillFillN.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

template <class ExecutionSpace, class IteratorType, class SizeType, class T>
IteratorType fill_n(const ExecutionSpace& ex, IteratorType first, SizeType n,
                    const T& value) {
  return Impl::fill_n_impl("Kokkos::fill_n_iterator_api_default", ex, first, n,
                           value);
}

template <class ExecutionSpace, class IteratorType, class SizeType, class T>
IteratorType fill_n(const std::string& label, const ExecutionSpace& ex,
                    IteratorType first, SizeType n, const T& value) {
  return Impl::fill_n_impl(label, ex, first, n, value);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class SizeType, class T>
auto fill_n(const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view, SizeType n,
            const T& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  return Impl::fill_n_impl("Kokkos::fill_n_view_api_default", ex, begin(view),
                           n, value);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class SizeType, class T>
auto fill_n(const std::string& label, const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view, SizeType n,
            const T& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  return Impl::fill_n_impl(label, ex, begin(view), n, value);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
