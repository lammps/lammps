// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

module;

#include <Kokkos_Bitset.hpp>

export module kokkos.bitset;

export {
  namespace Kokkos {
  using ::Kokkos::Bitset;
  using ::Kokkos::ConstBitset;

  using ::Kokkos::deep_copy;
  }  // namespace Kokkos
}
