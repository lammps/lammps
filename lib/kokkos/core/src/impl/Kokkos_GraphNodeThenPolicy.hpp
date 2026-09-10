// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_KOKKOS_GRAPHNODETHENPOLICY_HPP
#define KOKKOS_IMPL_KOKKOS_GRAPHNODETHENPOLICY_HPP

#include <type_traits>

namespace Kokkos::Experimental {

template <typename WorkTag>
struct ThenPolicy {
  static_assert(std::is_empty_v<WorkTag> || std::is_void_v<WorkTag>);

  using work_tag = WorkTag;
};

}  // namespace Kokkos::Experimental

#endif  // KOKKOS_IMPL_KOKKOS_GRAPHNODETHENPOLICY_HPP
