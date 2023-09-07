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

#ifndef KOKKOS_IMPL_KOKKOS_GRAPHNODECUSTOMIZATION_HPP
#define KOKKOS_IMPL_KOKKOS_GRAPHNODECUSTOMIZATION_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_Graph_fwd.hpp>
#include <impl/Kokkos_GraphImpl_fwd.hpp>

namespace Kokkos {
namespace Impl {

// Customizable for backends
template <class ExecutionSpace, class Kernel, class PredecessorRef>
struct GraphNodeBackendDetailsBeforeTypeErasure {
 protected:
  //----------------------------------------------------------------------------
  // <editor-fold desc="ctors, destructor, and assignment"> {{{2

  // Required constructors in customizations:
  GraphNodeBackendDetailsBeforeTypeErasure(
      ExecutionSpace const&, Kernel&, PredecessorRef const&,
      GraphNodeBackendSpecificDetails<ExecutionSpace>&
      /* this_as_details */) noexcept {}
  GraphNodeBackendDetailsBeforeTypeErasure(
      ExecutionSpace const&, _graph_node_is_root_ctor_tag,
      GraphNodeBackendSpecificDetails<ExecutionSpace>&
      /* this_as_details */) noexcept {}

  // Not copyable or movable at the concept level, so the default
  // implementation shouldn't be either.
  GraphNodeBackendDetailsBeforeTypeErasure() = delete;

  GraphNodeBackendDetailsBeforeTypeErasure(
      GraphNodeBackendDetailsBeforeTypeErasure const&) = delete;

  GraphNodeBackendDetailsBeforeTypeErasure(
      GraphNodeBackendDetailsBeforeTypeErasure&&) = delete;

  GraphNodeBackendDetailsBeforeTypeErasure& operator   =(
      GraphNodeBackendDetailsBeforeTypeErasure const&) = delete;

  GraphNodeBackendDetailsBeforeTypeErasure& operator=(
      GraphNodeBackendDetailsBeforeTypeErasure&&) = delete;

  ~GraphNodeBackendDetailsBeforeTypeErasure() = default;

  // </editor-fold> end ctors, destructor, and assignment }}}2
  //----------------------------------------------------------------------------
};

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_GRAPHNODECUSTOMIZATION_HPP
