// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

module;

#include <Kokkos_DualView.hpp>

export module kokkos.dual_view;

export {
  namespace Kokkos {
  using ::Kokkos::DualView;

  using ::Kokkos::is_dual_view;
  using ::Kokkos::is_dual_view_v;

  using ::Kokkos::deep_copy;

  using ::Kokkos::realloc;
  using ::Kokkos::resize;

  using ::Kokkos::subview;
  }  // namespace Kokkos
}
