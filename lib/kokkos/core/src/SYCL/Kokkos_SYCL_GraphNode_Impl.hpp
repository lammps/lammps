// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_SYCL_GRAPHNODE_IMPL_HPP
#define KOKKOS_SYCL_GRAPHNODE_IMPL_HPP

#include <Kokkos_Graph_fwd.hpp>

#include <impl/Kokkos_GraphImpl.hpp>

#include <SYCL/Kokkos_SYCL.hpp>

#include <optional>

namespace Kokkos {
namespace Impl {
template <>
struct GraphNodeBackendSpecificDetails<Kokkos::SYCL> {
  std::optional<sycl::ext::oneapi::experimental::node> node;

  explicit GraphNodeBackendSpecificDetails() = default;

  explicit GraphNodeBackendSpecificDetails(
      _graph_node_is_root_ctor_tag) noexcept {}
};

template <typename Kernel, typename PredecessorRef>
struct GraphNodeBackendDetailsBeforeTypeErasure<Kokkos::SYCL, Kernel,
                                                PredecessorRef> {
 protected:
  GraphNodeBackendDetailsBeforeTypeErasure(
      Kokkos::SYCL const &, Kernel &, PredecessorRef const &,
      GraphNodeBackendSpecificDetails<Kokkos::SYCL> &) noexcept {}

  GraphNodeBackendDetailsBeforeTypeErasure(
      Kokkos::SYCL const &, _graph_node_is_root_ctor_tag,
      GraphNodeBackendSpecificDetails<Kokkos::SYCL> &) noexcept {}
};

}  // namespace Impl
}  // namespace Kokkos

#endif
