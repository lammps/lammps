// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_KOKKOS_GRAPH_FWD_HPP
#define KOKKOS_KOKKOS_GRAPH_FWD_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_GRAPH_FWD
#endif

#include <Kokkos_Concepts.hpp>

#include <cstdint>

namespace Kokkos {
namespace Experimental {

struct TypeErasedTag {};
struct GraphNodeRootTag {};

template <class ExecutionSpace>
struct Graph;

enum class GraphNodeKind : std::uint8_t {
  TypeErased,
  Root,
  Kernel,
  Aggregate,
  Host,
  Capture,
};

template <Kokkos::ExecutionSpace ExecutionSpace, class Kernel = TypeErasedTag,
          class Predecessor = TypeErasedTag>
class GraphNodeRef;

template <class Worktag = void>
struct ThenPolicy;

}  // end namespace Experimental
}  // end namespace Kokkos

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_GRAPH_FWD
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_GRAPH_FWD
#endif
#endif  // KOKKOS_KOKKOS_GRAPH_FWD_HPP
