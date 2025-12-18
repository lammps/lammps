// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_STD_ALGORITHMS_VALUE_WRAPPER_FOR_NO_NEUTRAL_ELEMENT_HPP
#define KOKKOS_STD_ALGORITHMS_VALUE_WRAPPER_FOR_NO_NEUTRAL_ELEMENT_HPP

namespace Kokkos {
namespace Experimental {
namespace Impl {

//
// scalar wrapper used for reductions and scans
// when we don't have neutral element
//
template <class Scalar>
struct ValueWrapperForNoNeutralElement {
  Scalar val;
  bool is_initial = true;
};

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
