// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

module;

#include <Kokkos_Random.hpp>

export module kokkos.random;

export {
  namespace Kokkos {
  using ::Kokkos::fill_random;
  using ::Kokkos::rand;
  using ::Kokkos::Random_XorShift1024_Pool;
  using ::Kokkos::Random_XorShift64_Pool;
  }  // namespace Kokkos
}
