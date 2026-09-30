// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

module;

#include <Kokkos_UnorderedMap.hpp>

export module kokkos.unordered_map;

export {
  namespace Kokkos {
  using ::Kokkos::UnorderedMap;

  using ::Kokkos::UnorderedMapInsertOpTypes;
  using ::Kokkos::UnorderedMapInsertResult;

  using ::Kokkos::create_mirror;

  using ::Kokkos::deep_copy;
  }  // namespace Kokkos
}
