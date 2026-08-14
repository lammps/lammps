// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_KOKKOS_HOST_GRAPH_FWD_HPP
#define KOKKOS_KOKKOS_HOST_GRAPH_FWD_HPP

#include <Kokkos_Macros.hpp>

namespace Kokkos {
namespace Impl {

template <class ExecutionSpace>
struct GraphNodeKernelDefaultImpl;

template <class ExecutionSpace>
struct GraphNodeAggregateDefaultImpl;

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_HOST_GRAPH_FWD_HPP
