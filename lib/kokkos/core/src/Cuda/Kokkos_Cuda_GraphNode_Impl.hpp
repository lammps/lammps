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

#ifndef KOKKOS_KOKKOS_CUDA_GRAPHNODE_IMPL_HPP
#define KOKKOS_KOKKOS_CUDA_GRAPHNODE_IMPL_HPP

#include <Kokkos_Macros.hpp>

#if defined(KOKKOS_ENABLE_CUDA)

#include <Kokkos_Graph_fwd.hpp>

#include <impl/Kokkos_GraphImpl.hpp>  // GraphAccess needs to be complete

#include <Cuda/Kokkos_Cuda.hpp>

namespace Kokkos {
namespace Impl {

template <>
struct GraphNodeBackendSpecificDetails<Kokkos::Cuda> {
  cudaGraphNode_t node = nullptr;

  //----------------------------------------------------------------------------
  // <editor-fold desc="Ctors, destructor, and assignment"> {{{2

  explicit GraphNodeBackendSpecificDetails() = default;

  explicit GraphNodeBackendSpecificDetails(
      _graph_node_is_root_ctor_tag) noexcept {}

  // </editor-fold> end Ctors, destructor, and assignment }}}2
  //----------------------------------------------------------------------------
};

template <class Kernel, class PredecessorRef>
struct GraphNodeBackendDetailsBeforeTypeErasure<Kokkos::Cuda, Kernel,
                                                PredecessorRef> {
 protected:
  //----------------------------------------------------------------------------
  // <editor-fold desc="ctors, destructor, and assignment"> {{{2

  GraphNodeBackendDetailsBeforeTypeErasure(
      Kokkos::Cuda const&, Kernel&, PredecessorRef const&,
      GraphNodeBackendSpecificDetails<Kokkos::Cuda>&) noexcept {}

  GraphNodeBackendDetailsBeforeTypeErasure(
      Kokkos::Cuda const&, _graph_node_is_root_ctor_tag,
      GraphNodeBackendSpecificDetails<Kokkos::Cuda>&) noexcept {}

  // </editor-fold> end ctors, destructor, and assignment }}}2
  //----------------------------------------------------------------------------
};

}  // end namespace Impl
}  // end namespace Kokkos

#include <Cuda/Kokkos_Cuda_GraphNodeKernel.hpp>

#endif  // defined(KOKKOS_ENABLE_CUDA)
#endif  // KOKKOS_KOKKOS_CUDA_GRAPHNODE_IMPL_HPP
