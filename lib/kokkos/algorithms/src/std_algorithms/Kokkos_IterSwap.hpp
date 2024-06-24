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

#ifndef KOKKOS_STD_ALGORITHMS_ITER_SWAP_HPP
#define KOKKOS_STD_ALGORITHMS_ITER_SWAP_HPP

#include <Kokkos_Core.hpp>
#include "impl/Kokkos_Constraints.hpp"

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IteratorType1, class IteratorType2>
struct StdIterSwapFunctor {
  IteratorType1 m_a;
  IteratorType2 m_b;

  KOKKOS_FUNCTION
  void operator()(int i) const {
    (void)i;
    ::Kokkos::kokkos_swap(*m_a, *m_b);
  }

  KOKKOS_FUNCTION
  StdIterSwapFunctor(IteratorType1 _a, IteratorType2 _b)
      : m_a(std::move(_a)), m_b(std::move(_b)) {}
};

template <class IteratorType1, class IteratorType2>
void iter_swap_impl(IteratorType1 a, IteratorType2 b) {
  // is there a better way to do this maybe?
  ::Kokkos::parallel_for(
      1, StdIterSwapFunctor<IteratorType1, IteratorType2>(a, b));
  Kokkos::DefaultExecutionSpace().fence(
      "Kokkos::iter_swap: fence after operation");
}
}  // namespace Impl
//----------------------------------------------------------------------------

// iter_swap
template <class IteratorType1, class IteratorType2>
void iter_swap(IteratorType1 a, IteratorType2 b) {
  Impl::iter_swap_impl(a, b);
}

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
template <class T>
KOKKOS_DEPRECATED_WITH_COMMENT("Use Kokkos::kokkos_swap instead!")
KOKKOS_FUNCTION
    void swap(T& a, T& b) noexcept(::Kokkos::kokkos_swap(std::declval<T&>(),
                                                         std::declval<T&>())) {
  ::Kokkos::kokkos_swap(a, b);
}
#endif

}  // namespace Experimental
}  // namespace Kokkos

#endif
