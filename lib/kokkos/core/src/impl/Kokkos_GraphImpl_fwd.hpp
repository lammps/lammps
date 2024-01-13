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

#ifndef KOKKOS_IMPL_KOKKOS_GRAPHIMPL_FWD_HPP
#define KOKKOS_IMPL_KOKKOS_GRAPHIMPL_FWD_HPP

#include <Kokkos_Macros.hpp>

namespace Kokkos {
namespace Impl {

template <class ExecutionSpace, class Kernel, class Predecessor>
struct GraphNodeImpl;

template <class ExecutionSpace>
struct GraphImpl;

template <class ExecutionSpace, class Policy, class Functor,
          class KernelTypeTag, class... Args>
class GraphNodeKernelImpl;

struct _graph_node_kernel_ctor_tag {};
struct _graph_node_predecessor_ctor_tag {};
struct _graph_node_is_root_ctor_tag {};

struct GraphAccess;

// Customizable for backends
template <class ExecutionSpace>
struct GraphNodeBackendSpecificDetails;

// Customizable for backends
template <class ExecutionSpace, class Kernel, class PredecessorRef>
struct GraphNodeBackendDetailsBeforeTypeErasure;

// TODO move this to a more appropriate place
struct DoNotExplicitlySpecifyThisTemplateParameter;

struct KernelInGraphProperty {};

struct IsGraphKernelTag {};

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_IMPL_KOKKOS_GRAPHIMPL_FWD_HPP
