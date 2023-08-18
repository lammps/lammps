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

#ifndef KOKKOS_KOKKOS_GRAPH_FWD_HPP
#define KOKKOS_KOKKOS_GRAPH_FWD_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_GRAPH_FWD
#endif

#include <Kokkos_Macros.hpp>

namespace Kokkos {
namespace Experimental {

struct TypeErasedTag {};

template <class ExecutionSpace>
struct Graph;

template <class ExecutionSpace, class Kernel = TypeErasedTag,
          class Predecessor = TypeErasedTag>
class GraphNodeRef;

}  // end namespace Experimental
}  // end namespace Kokkos

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_GRAPH_FWD
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_GRAPH_FWD
#endif
#endif  // KOKKOS_KOKKOS_GRAPH_FWD_HPP
