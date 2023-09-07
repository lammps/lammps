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

#ifndef KOKKOS_STD_ALGORITHMS_COPY_N_HPP
#define KOKKOS_STD_ALGORITHMS_COPY_N_HPP

#include "impl/Kokkos_CopyCopyN.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

template <class ExecutionSpace, class InputIterator, class Size,
          class OutputIterator>
OutputIterator copy_n(const ExecutionSpace& ex, InputIterator first, Size count,
                      OutputIterator result) {
  return Impl::copy_n_impl("Kokkos::copy_n_iterator_api_default", ex, first,
                           count, result);
}

template <class ExecutionSpace, class InputIterator, class Size,
          class OutputIterator>
OutputIterator copy_n(const std::string& label, const ExecutionSpace& ex,
                      InputIterator first, Size count, OutputIterator result) {
  return Impl::copy_n_impl(label, ex, first, count, result);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class Size, class DataType2, class... Properties2>
auto copy_n(const ExecutionSpace& ex,
            const ::Kokkos::View<DataType1, Properties1...>& source, Size count,
            ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  namespace KE = ::Kokkos::Experimental;
  return Impl::copy_n_impl("Kokkos::copy_n_view_api_default", ex,
                           KE::cbegin(source), count, KE::begin(dest));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class Size, class DataType2, class... Properties2>
auto copy_n(const std::string& label, const ExecutionSpace& ex,
            const ::Kokkos::View<DataType1, Properties1...>& source, Size count,
            ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  namespace KE = ::Kokkos::Experimental;
  return Impl::copy_n_impl(label, ex, KE::cbegin(source), count,
                           KE::begin(dest));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
