// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_STD_ALGORITHMS_NUMERIC_IDENTITY_REFERENCE_UNARY_FUNCTOR_HPP
#define KOKKOS_STD_ALGORITHMS_NUMERIC_IDENTITY_REFERENCE_UNARY_FUNCTOR_HPP

#include <Kokkos_Macros.hpp>

namespace Kokkos {
namespace Experimental {
namespace Impl {

struct StdNumericScanIdentityReferenceUnaryFunctor {
  template <class T>
  KOKKOS_FUNCTION constexpr T&& operator()(T&& t) const {
    return static_cast<T&&>(t);
  }
};

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
