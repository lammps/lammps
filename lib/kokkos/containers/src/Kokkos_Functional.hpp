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

#ifndef KOKKOS_FUNCTIONAL_HPP
#define KOKKOS_FUNCTIONAL_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_FUNCTIONAL
#endif

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_Functional_impl.hpp>

namespace Kokkos {

// These should work for most types

template <typename T>
struct pod_hash {
  KOKKOS_FORCEINLINE_FUNCTION
  uint32_t operator()(T const& t) const {
    return Impl::MurmurHash3_x86_32(&t, sizeof(T), 0);
  }

  KOKKOS_FORCEINLINE_FUNCTION
  uint32_t operator()(T const& t, uint32_t seed) const {
    return Impl::MurmurHash3_x86_32(&t, sizeof(T), seed);
  }
};

template <typename T>
struct pod_equal_to {
  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const& a, T const& b) const {
    return Impl::bitwise_equal(&a, &b);
  }
};

template <typename T>
struct pod_not_equal_to {
  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const& a, T const& b) const {
    return !Impl::bitwise_equal(&a, &b);
  }
};

template <typename T>
struct equal_to {
  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const& a, T const& b) const { return a == b; }
};

template <typename T>
struct not_equal_to {
  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const& a, T const& b) const { return a != b; }
};

template <typename T>
struct greater {
  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const& a, T const& b) const { return a > b; }
};

template <typename T>
struct less {
  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const& a, T const& b) const { return a < b; }
};

template <typename T>
struct greater_equal {
  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const& a, T const& b) const { return a >= b; }
};

template <typename T>
struct less_equal {
  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const& a, T const& b) const { return a <= b; }
};

}  // namespace Kokkos

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_FUNCTIONAL
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_FUNCTIONAL
#endif
#endif  // KOKKOS_FUNCTIONAL_HPP
