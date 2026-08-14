// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

module;

#include <Kokkos_DynRankView.hpp>

export module kokkos.dyn_rank_view;

export {
  namespace Kokkos {
  using ::Kokkos::DynRankView;

  using ::Kokkos::is_dyn_rank_view;
  using ::Kokkos::is_dyn_rank_view_v;

  using ::Kokkos::Subdynrankview;
  using ::Kokkos::subdynrankview;
  using ::Kokkos::subview;

  using ::Kokkos::rank;

  using ::Kokkos::deep_copy;
  using ::Kokkos::realloc;
  using ::Kokkos::resize;

  using ::Kokkos::create_mirror;
  using ::Kokkos::create_mirror_view;
  using ::Kokkos::create_mirror_view_and_copy;

  using ::Kokkos::operator!=;
  using ::Kokkos::operator==;
  }  // namespace Kokkos
}
