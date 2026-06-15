// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

module;

#include <Kokkos_ScatterView.hpp>

export module kokkos.scatter_view_impl;

export {
  namespace Kokkos::Impl::Experimental {
  using ::Kokkos::Impl::Experimental::DefaultContribution;
  using ::Kokkos::Impl::Experimental::DefaultDuplication;
  }  // namespace Kokkos::Impl::Experimental
}
