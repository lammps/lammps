// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

module;

#include <Kokkos_Functional.hpp>

export module kokkos.functional;

export {
  namespace Kokkos {
  using ::Kokkos::equal_to;
  using ::Kokkos::greater;
  using ::Kokkos::greater_equal;
  using ::Kokkos::less;
  using ::Kokkos::less_equal;
  using ::Kokkos::not_equal_to;
  using ::Kokkos::pod_equal_to;
  using ::Kokkos::pod_hash;
  using ::Kokkos::pod_not_equal_to;
  }  // namespace Kokkos
}
