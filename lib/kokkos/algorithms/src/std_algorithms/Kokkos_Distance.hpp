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

#ifndef KOKKOS_STD_ALGORITHMS_DISTANCE_HPP
#define KOKKOS_STD_ALGORITHMS_DISTANCE_HPP

#include "impl/Kokkos_Constraints.hpp"
#include "impl/Kokkos_RandomAccessIterator.hpp"

namespace Kokkos {
namespace Experimental {

template <class IteratorType>
KOKKOS_INLINE_FUNCTION constexpr typename IteratorType::difference_type
distance(IteratorType first, IteratorType last) {
  static_assert(
      ::Kokkos::Experimental::Impl::are_random_access_iterators<
          IteratorType>::value,
      "Kokkos::Experimental::distance: only implemented for random access "
      "iterators.");

  return last - first;
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
