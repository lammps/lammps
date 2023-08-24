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

#ifndef KOKKOS_STD_ALGORITHMS_COPY_IF_HPP
#define KOKKOS_STD_ALGORITHMS_COPY_IF_HPP

#include "impl/Kokkos_CopyIf.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class Predicate>
OutputIterator copy_if(const ExecutionSpace& ex, InputIterator first,
                       InputIterator last, OutputIterator d_first,
                       Predicate pred) {
  return Impl::copy_if_impl("Kokkos::copy_if_iterator_api_default", ex, first,
                            last, d_first, std::move(pred));
}

template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class Predicate>
OutputIterator copy_if(const std::string& label, const ExecutionSpace& ex,
                       InputIterator first, InputIterator last,
                       OutputIterator d_first, Predicate pred) {
  return Impl::copy_if_impl(label, ex, first, last, d_first, std::move(pred));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class Predicate>
auto copy_if(const ExecutionSpace& ex,
             const ::Kokkos::View<DataType1, Properties1...>& source,
             ::Kokkos::View<DataType2, Properties2...>& dest, Predicate pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::copy_if_impl("Kokkos::copy_if_view_api_default", ex,
                            cbegin(source), cend(source), begin(dest),
                            std::move(pred));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class Predicate>
auto copy_if(const std::string& label, const ExecutionSpace& ex,
             const ::Kokkos::View<DataType1, Properties1...>& source,
             ::Kokkos::View<DataType2, Properties2...>& dest, Predicate pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::copy_if_impl(label, ex, cbegin(source), cend(source),
                            begin(dest), std::move(pred));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
