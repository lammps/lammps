//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_HIP_GRAPHNODE_IMPL_HPP
#define KOKKOS_HIP_GRAPHNODE_IMPL_HPP

#include <Kokkos_Graph_fwd.hpp>

#include <impl/Kokkos_GraphImpl.hpp>

#include <HIP/Kokkos_HIP.hpp>

namespace Kokkos {
namespace Impl {
template <>
struct GraphNodeBackendSpecificDetails<Kokkos::HIP> {
  hipGraphNode_t node = nullptr;

  explicit GraphNodeBackendSpecificDetails() = default;

  explicit GraphNodeBackendSpecificDetails(
      _graph_node_is_root_ctor_tag) noexcept {}
};

template <typename Kernel, typename PredecessorRef>
struct GraphNodeBackendDetailsBeforeTypeErasure<Kokkos::HIP, Kernel,
                                                PredecessorRef> {
 protected:
  GraphNodeBackendDetailsBeforeTypeErasure(
      Kokkos::HIP const &, Kernel &, PredecessorRef const &,
      GraphNodeBackendSpecificDetails<Kokkos::HIP> &) noexcept {}

  GraphNodeBackendDetailsBeforeTypeErasure(
      Kokkos::HIP const &, _graph_node_is_root_ctor_tag,
      GraphNodeBackendSpecificDetails<Kokkos::HIP> &) noexcept {}
};

}  // namespace Impl
}  // namespace Kokkos

#endif
