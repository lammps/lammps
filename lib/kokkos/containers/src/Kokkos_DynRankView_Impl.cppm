// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

module;

#include <Kokkos_DynRankView.hpp>

export module kokkos.dyn_rank_view_impl;

export {
  namespace Kokkos::Impl {
  using ::Kokkos::Impl::ApplyToViewOfStaticRank;
  using ::Kokkos::Impl::as_view_of_rank_n;
  }  // namespace Kokkos::Impl
}
