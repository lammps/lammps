// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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

  GraphNodeKernelDefaultImpl() = default;

  explicit GraphNodeKernelDefaultImpl(ExecutionSpace exec)
      : m_execution_space(std::move(exec)) {}

  ExecutionSpace m_execution_space;
};

template <typename ExecutionSpace, typename Functor>
struct GraphNodeThenHostImpl
    : public GraphNodeKernelDefaultImpl<ExecutionSpace> {
  using execute_kernel_vtable_base_t =
      GraphNodeKernelDefaultImpl<ExecutionSpace>;

  explicit GraphNodeThenHostImpl(Functor functor)
      : execute_kernel_vtable_base_t{}, m_functor(std::move(functor)) {}

  void execute_kernel() override final {
    this->m_execution_space.fence(
        "Kokkos::DefaultGraphNode::then_host: fence needed before host "
        "callback");
    m_functor();
  }

  Functor m_functor;
};

// TODO Indicate that this kernel specialization is only for the Host somehow?
template <class ExecutionSpace, class PolicyType, class Functor,
          class PatternTag, class... Args>
class GraphNodeKernelImpl
    : public GraphNodeKernelDefaultImpl<ExecutionSpace>,
      public PatternImplSpecializationFromTag<PatternTag, Functor, PolicyType,
                                              Args..., ExecutionSpace>::type {
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
  GraphNodeKernelImpl(std::string const &, ExecutionSpace const &,
                      Functor arg_functor, PolicyDeduced &&arg_policy,
                      ArgsDeduced &&...args)
      : execute_kernel_vtable_base_t(arg_policy.space()),
        base_t(std::move(arg_functor), (PolicyDeduced &&)arg_policy,
               (ArgsDeduced &&)args...) {}

  // FIXME @graph Forward through the instance once that works in the backends
  template <class PolicyDeduced, class... ArgsDeduced>
  GraphNodeKernelImpl(ExecutionSpace const &ex, Functor arg_functor,
                      PolicyDeduced &&arg_policy, ArgsDeduced &&...args)
      : GraphNodeKernelImpl("", ex, std::move(arg_functor),
                            (PolicyDeduced &&)arg_policy,
                            (ArgsDeduced &&)args...) {
    // FIXME This constructor seem unused.
  }

  void execute_kernel() override final { this->base_t::execute(); }
};

// </editor-fold> end GraphNodeKernelImpl }}}1
//==============================================================================

template <class ExecutionSpace>
struct GraphNodeAggregateDefaultImpl
    : GraphNodeKernelDefaultImpl<ExecutionSpace> {
  void execute_kernel() override final {}
};

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_HOST_GRAPHNODEKERNEL_HPP
