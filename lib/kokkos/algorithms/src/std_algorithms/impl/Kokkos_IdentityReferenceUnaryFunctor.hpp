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

#ifndef KOKKOS_STD_ALGORITHMS_NUMERIC_IDENTITY_REFERENCE_UNARY_FUNCTOR_HPP
#define KOKKOS_STD_ALGORITHMS_NUMERIC_IDENTITY_REFERENCE_UNARY_FUNCTOR_HPP

#include <Kokkos_Macros.hpp>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class ValueType>
struct StdNumericScanIdentityReferenceUnaryFunctor {
  KOKKOS_FUNCTION
  constexpr const ValueType& operator()(const ValueType& a) const { return a; }
};

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
