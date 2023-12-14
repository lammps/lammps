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

#ifndef KOKKOS_HOST_GRAPH_IMPL_HPP
#define KOKKOS_HOST_GRAPH_IMPL_HPP

#include <Kokkos_ExecPolicy.hpp>
#include <Kokkos_Graph.hpp>

#include <impl/Kokkos_GraphImpl_fwd.hpp>
#include <impl/Kokkos_Default_Graph_fwd.hpp>

#include <Serial/Kokkos_Serial.hpp>
#include <OpenMP/Kokkos_OpenMP.hpp>
// FIXME @graph other backends?

#include <impl/Kokkos_OptionalRef.hpp>
#include <impl/Kokkos_EBO.hpp>

#include <set>

namespace Kokkos {
namespace Impl {

//==============================================================================
// <editor-fold desc="GraphImpl default implementation"> {{{1

template <class ExecutionSpace>
struct GraphImpl : private ExecutionSpaceInstanceStorage<ExecutionSpace> {
 public:
  using root_node_impl_t =
      GraphNodeImpl<ExecutionSpace, Kokkos::Experimental::TypeErasedTag,
                    Kokkos::Experimental::TypeErasedTag>;

 private:
  using execution_space_instance_storage_base_t =
      ExecutionSpaceInstanceStorage<ExecutionSpace>;

  using node_details_t = GraphNodeBackendSpecificDetails<ExecutionSpace>;
  std::set<std::shared_ptr<node_details_t>> m_sinks;

 public:
  //----------------------------------------------------------------------------
  // <editor-fold desc="Constructors, destructor, and assignment"> {{{2

  // Not moveable or copyable; it spends its whole live as a shared_ptr in the
  // Graph object
  GraphImpl()                 = default;
  GraphImpl(GraphImpl const&) = delete;
  GraphImpl(GraphImpl&&)      = delete;
  GraphImpl& operator=(GraphImpl const&) = delete;
  GraphImpl& operator=(GraphImpl&&) = delete;
  ~GraphImpl()                      = default;

  explicit GraphImpl(ExecutionSpace arg_space)
      : execution_space_instance_storage_base_t(std::move(arg_space)) {}

  // </editor-fold> end Constructors, destructor, and assignment }}}2
  //----------------------------------------------------------------------------

  ExecutionSpace const& get_execution_space() const {
    return this
        ->execution_space_instance_storage_base_t::execution_space_instance();
  }

  //----------------------------------------------------------------------------
  // <editor-fold desc="required customizations"> {{{2

  template <class NodeImpl>
  //  requires NodeImplPtr is a shared_ptr to specialization of GraphNodeImpl
  void add_node(std::shared_ptr<NodeImpl> const& arg_node_ptr) {
    static_assert(
        NodeImpl::kernel_type::Policy::is_graph_kernel::value,
        "Something has gone horribly wrong, but it's too complicated to "
        "explain here.  Buy Daisy a coffee and she'll explain it to you.");
    // Since this is always called before any calls to add_predecessor involving
    // it, we can treat this node as a sink until we discover otherwise.
    arg_node_ptr->node_details_t::set_kernel(arg_node_ptr->get_kernel());
    auto spot = m_sinks.find(arg_node_ptr);
    KOKKOS_ASSERT(spot == m_sinks.end())
    m_sinks.insert(std::move(spot), std::move(arg_node_ptr));
  }

  template <class NodeImplPtr, class PredecessorRef>
  // requires PredecessorRef is a specialization of GraphNodeRef that has
  // already been added to this graph and NodeImpl is a specialization of
  // GraphNodeImpl that has already been added to this graph.
  void add_predecessor(NodeImplPtr arg_node_ptr, PredecessorRef arg_pred_ref) {
    auto node_ptr_spot = m_sinks.find(arg_node_ptr);
    auto pred_ptr      = GraphAccess::get_node_ptr(arg_pred_ref);
    auto pred_ref_spot = m_sinks.find(pred_ptr);
    KOKKOS_ASSERT(node_ptr_spot != m_sinks.end())
    if (pred_ref_spot != m_sinks.end()) {
      // delegate responsibility for executing the predecessor to arg_node
      // and then remove the predecessor from the set of sinks
      (*node_ptr_spot)->set_predecessor(std::move(*pred_ref_spot));
      m_sinks.erase(pred_ref_spot);
    } else {
      // We still want to check that it's executed, even though someone else
      // should have executed it before us
      (*node_ptr_spot)->set_predecessor(std::move(pred_ptr));
    }
  }

  template <class... PredecessorRefs>
  // See requirements/expectations in GraphBuilder
  auto create_aggregate_ptr(PredecessorRefs&&...) {
    // The attachment to predecessors, which is all we really need, happens
    // in the generic layer, which calls through to add_predecessor for
    // each predecessor ref, so all we need to do here is create the (trivial)
    // aggregate node.
    using aggregate_kernel_impl_t =
        GraphNodeAggregateKernelDefaultImpl<ExecutionSpace>;
    using aggregate_node_impl_t =
        GraphNodeImpl<ExecutionSpace, aggregate_kernel_impl_t,
                      Kokkos::Experimental::TypeErasedTag>;
    return GraphAccess::make_node_shared_ptr<aggregate_node_impl_t>(
        this->execution_space_instance(), _graph_node_kernel_ctor_tag{},
        aggregate_kernel_impl_t{});
  }

  auto create_root_node_ptr() {
    auto rv = Kokkos::Impl::GraphAccess::make_node_shared_ptr<root_node_impl_t>(
        get_execution_space(), _graph_node_is_root_ctor_tag{});
    m_sinks.insert(rv);
    return rv;
  }

  void submit() {
    // This reset is gross, but for the purposes of our simple host
    // implementation...
    for (auto& sink : m_sinks) {
      sink->reset_has_executed();
    }
    for (auto& sink : m_sinks) {
      sink->execute_node();
    }
  }

  // </editor-fold> end required customizations }}}2
  //----------------------------------------------------------------------------
};

// </editor-fold> end GraphImpl default implementation }}}1
//==============================================================================

}  // end namespace Impl

}  // end namespace Kokkos

#include <impl/Kokkos_Default_GraphNodeKernel.hpp>
#include <impl/Kokkos_Default_GraphNode_Impl.hpp>

#endif  // KOKKOS_HOST_GRAPH_IMPL_HPP
