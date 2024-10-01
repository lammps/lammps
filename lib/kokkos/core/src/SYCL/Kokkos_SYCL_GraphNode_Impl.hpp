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

#ifndef KOKKOS_SYCL_GRAPHNODE_IMPL_HPP
#define KOKKOS_SYCL_GRAPHNODE_IMPL_HPP

#include <Kokkos_Graph_fwd.hpp>

#include <impl/Kokkos_GraphImpl.hpp>

#include <SYCL/Kokkos_SYCL.hpp>

#include <optional>

namespace Kokkos {
namespace Impl {
template <>
struct GraphNodeBackendSpecificDetails<Kokkos::Experimental::SYCL> {
  std::optional<sycl::ext::oneapi::experimental::node> node;

  explicit GraphNodeBackendSpecificDetails() = default;

  explicit GraphNodeBackendSpecificDetails(
      _graph_node_is_root_ctor_tag) noexcept {}
};

template <typename Kernel, typename PredecessorRef>
struct GraphNodeBackendDetailsBeforeTypeErasure<Kokkos::Experimental::SYCL,
                                                Kernel, PredecessorRef> {
 protected:
  GraphNodeBackendDetailsBeforeTypeErasure(
      Kokkos::Experimental::SYCL const &, Kernel &, PredecessorRef const &,
      GraphNodeBackendSpecificDetails<Kokkos::Experimental::SYCL> &) noexcept {}

  GraphNodeBackendDetailsBeforeTypeErasure(
      Kokkos::Experimental::SYCL const &, _graph_node_is_root_ctor_tag,
      GraphNodeBackendSpecificDetails<Kokkos::Experimental::SYCL> &) noexcept {}
};

}  // namespace Impl
}  // namespace Kokkos

#endif
