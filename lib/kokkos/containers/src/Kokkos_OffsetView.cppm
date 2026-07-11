// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

module;

#include <Kokkos_OffsetView.hpp>

export module kokkos.offset_view;

export {
  namespace Kokkos {
  namespace Experimental {
  using ::Kokkos::Experimental::OffsetView;

  using ::Kokkos::Experimental::is_offset_view;
  using ::Kokkos::Experimental::is_offset_view_v;

  using ::Kokkos::Experimental::index_list_type;
  using ::Kokkos::Experimental::IndexRange;

  using ::Kokkos::Experimental::operator==;
  using ::Kokkos::Experimental::operator!=;
  }  // namespace Experimental

  using ::Kokkos::create_mirror;
  using ::Kokkos::create_mirror_view;
  using ::Kokkos::create_mirror_view_and_copy;

  using ::Kokkos::deep_copy;

  using ::Kokkos::subview;
  }  // namespace Kokkos
}
