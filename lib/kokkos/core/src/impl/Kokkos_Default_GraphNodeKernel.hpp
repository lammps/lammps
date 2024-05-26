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

#ifndef KOKKOS_KOKKOS_HOST_GRAPHNODEKERNEL_HPP
#define KOKKOS_KOKKOS_HOST_GRAPHNODEKERNEL_HPP

#include <Kokkos_Macros.hpp>

#include <impl/Kokkos_Default_Graph_fwd.hpp>

#include <Kokkos_Graph.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_Parallel_Reduce.hpp>

namespace Kokkos {
namespace Impl {

//==============================================================================
// <editor-fold desc="GraphNodeKernelImpl"> {{{1

template <class ExecutionSpace>
struct GraphNodeKernelDefaultImpl {
  // TODO @graphs decide if this should use vtable or intrusive erasure via
  //      function pointers like in the rest of the graph interface
  virtual void execute_kernel() = 0;
};

// TODO Indicate that this kernel specialization is only for the Host somehow?
template <class ExecutionSpace, class PolicyType, class Functor,
          class PatternTag, class... Args>
class GraphNodeKernelImpl
    : public PatternImplSpecializationFromTag<PatternTag, Functor, PolicyType,
                                              Args..., ExecutionSpace>::type,
      public GraphNodeKernelDefaultImpl<ExecutionSpace> {
 public:
  using base_t =
      typename PatternImplSpecializationFromTag<PatternTag, Functor, PolicyType,
                                                Args..., ExecutionSpace>::type;
  using execute_kernel_vtable_base_t =
      GraphNodeKernelDefaultImpl<ExecutionSpace>;
  // We have to use this name here because that's how it was done way back when
  // then implementations of Impl::Parallel*<> were written
  using Policy       = PolicyType;
  using graph_kernel = GraphNodeKernelImpl;

  // TODO @graph kernel name info propagation
  template <class PolicyDeduced, class... ArgsDeduced>
  GraphNodeKernelImpl(std::string const&, ExecutionSpace const&,
                      Functor arg_functor, PolicyDeduced&& arg_policy,
                      ArgsDeduced&&... args)
      : base_t(std::move(arg_functor), (PolicyDeduced &&) arg_policy,
               (ArgsDeduced &&) args...),
        execute_kernel_vtable_base_t() {}

  // FIXME @graph Forward through the instance once that works in the backends
  template <class PolicyDeduced, class... ArgsDeduced>
  GraphNodeKernelImpl(ExecutionSpace const& ex, Functor arg_functor,
                      PolicyDeduced&& arg_policy, ArgsDeduced&&... args)
      : GraphNodeKernelImpl("", ex, std::move(arg_functor),
                            (PolicyDeduced &&) arg_policy,
                            (ArgsDeduced &&) args...) {}

  void execute_kernel() final { this->base_t::execute(); }
};

// </editor-fold> end GraphNodeKernelImpl }}}1
//==============================================================================

template <class ExecutionSpace>
struct GraphNodeAggregateKernelDefaultImpl
    : GraphNodeKernelDefaultImpl<ExecutionSpace> {
  // Aggregates don't need a policy, but for the purposes of checking the static
  // assertions about graph kernels,
  struct Policy {
    using is_graph_kernel = std::true_type;
  };
  using graph_kernel = GraphNodeAggregateKernelDefaultImpl;
  void execute_kernel() final {}
};

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_HOST_GRAPHNODEKERNEL_HPP
