// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

module;

#include <Kokkos_ScatterView.hpp>

export module kokkos.scatter_view;

export {
  namespace Kokkos {

  namespace Experimental {
  using ::Kokkos::Experimental::ScatterView;

  using ::Kokkos::Experimental::contribute;
  using ::Kokkos::Experimental::create_scatter_view;

  using ::Kokkos::Experimental::is_scatter_view;
  using ::Kokkos::Experimental::is_scatter_view_v;

  using ::Kokkos::Experimental::ScatterDuplicated;
  using ::Kokkos::Experimental::ScatterNonDuplicated;

  using ::Kokkos::Experimental::ScatterAccess;

  using ::Kokkos::Experimental::ScatterAtomic;
  using ::Kokkos::Experimental::ScatterNonAtomic;

  using ::Kokkos::Experimental::ScatterMax;
  using ::Kokkos::Experimental::ScatterMin;
  using ::Kokkos::Experimental::ScatterProd;
  using ::Kokkos::Experimental::ScatterSum;
  }  // namespace Experimental

  using ::Kokkos::realloc;
  using ::Kokkos::resize;
  }  // namespace Kokkos
}
