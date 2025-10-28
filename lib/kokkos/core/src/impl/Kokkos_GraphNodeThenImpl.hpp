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

#ifndef KOKKOS_IMPL_KOKKOS_GRAPHNODETHENIMPL_HPP
#define KOKKOS_IMPL_KOKKOS_GRAPHNODETHENIMPL_HPP

#include <Kokkos_ExecPolicy.hpp>
#include <impl/Kokkos_GraphImpl_fwd.hpp>

namespace Kokkos::Impl {

// Helper for the 'then', such that the user can indeed pass a callable that
// takes no argument.
template <typename Functor>
struct ThenWrapper {
  Functor m_functor;
  template <typename T>
  KOKKOS_FUNCTION void operator()(const T) const {
    m_functor();
  }
};

template <typename ExecutionSpace, typename Functor>
struct GraphNodeThenImpl
    : public GraphNodeKernelImpl<
          ExecutionSpace,
          Kokkos::RangePolicy<ExecutionSpace, IsGraphKernelTag,
                              Kokkos::LaunchBounds<1>>,
          ThenWrapper<Functor>, ParallelForTag> {
  using policy_t  = Kokkos::RangePolicy<ExecutionSpace, IsGraphKernelTag,
                                       Kokkos::LaunchBounds<1>>;
  using base_t    = GraphNodeKernelImpl<ExecutionSpace, policy_t,
                                     ThenWrapper<Functor>, ParallelForTag>;
  using wrapper_t = ThenWrapper<Functor>;

  template <typename Label, typename T>
  GraphNodeThenImpl(Label&& label_, const ExecutionSpace& exec, T&& functor)
      : base_t(std::forward<Label>(label_), exec,
               wrapper_t{std::forward<T>(functor)}, policy_t(exec, 0, 1)) {}
};

}  // namespace Kokkos::Impl

#endif  // KOKKOS_IMPL_KOKKOS_GRAPHNODETHENIMPL_HPP
