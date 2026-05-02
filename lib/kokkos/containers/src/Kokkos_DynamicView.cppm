// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

module;

#include <Kokkos_DynamicView.hpp>

export module kokkos.dynamic_view;

export {
  namespace Kokkos {
  namespace Experimental {
  using ::Kokkos::Experimental::DynamicView;
  }

  using ::Kokkos::is_dynamic_view;
  using ::Kokkos::is_dynamic_view_v;

  using ::Kokkos::create_mirror;
  using ::Kokkos::create_mirror_view;
  using ::Kokkos::create_mirror_view_and_copy;

  using ::Kokkos::deep_copy;
  }  // namespace Kokkos
}
