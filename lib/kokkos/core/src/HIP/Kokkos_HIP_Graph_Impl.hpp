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

#ifndef KOKKOS_HIP_GRAPH_IMPL_HPP
#define KOKKOS_HIP_GRAPH_IMPL_HPP

#include <Kokkos_Macros.hpp>

#include <Kokkos_Graph_fwd.hpp>

#include <impl/Kokkos_GraphImpl.hpp>
#include <impl/Kokkos_GraphNodeImpl.hpp>

#include <HIP/Kokkos_HIP_GraphNodeKernel.hpp>

namespace Kokkos {
namespace Impl {
template <>
class GraphImpl<Kokkos::HIP> {
 public:
  using node_details_t = GraphNodeBackendSpecificDetails<Kokkos::HIP>;
  using root_node_impl_t =
      GraphNodeImpl<Kokkos::HIP, Kokkos::Experimental::TypeErasedTag,
                    Kokkos::Experimental::TypeErasedTag>;
  using aggregate_kernel_impl_t = HIPGraphNodeAggregateKernel;
  using aggregate_node_impl_t =
      GraphNodeImpl<Kokkos::HIP, aggregate_kernel_impl_t,
                    Kokkos::Experimental::TypeErasedTag>;

  // Not movable or copyable; it spends its whole life as a shared_ptr in the
  // Graph object.
  GraphImpl()                 = delete;
  GraphImpl(GraphImpl const&) = delete;
  GraphImpl(GraphImpl&&)      = delete;
  GraphImpl& operator=(GraphImpl const&) = delete;
  GraphImpl& operator=(GraphImpl&&) = delete;

  ~GraphImpl();

  explicit GraphImpl(Kokkos::HIP instance);

  void add_node(std::shared_ptr<aggregate_node_impl_t> const& arg_node_ptr);

  template <class NodeImpl>
  void add_node(std::shared_ptr<NodeImpl> const& arg_node_ptr);

  template <class NodeImplPtr, class PredecessorRef>
  void add_predecessor(NodeImplPtr arg_node_ptr, PredecessorRef arg_pred_ref);

  void submit();

  Kokkos::HIP const& get_execution_space() const noexcept;

  auto create_root_node_ptr();

  template <class... PredecessorRefs>
  auto create_aggregate_ptr(PredecessorRefs&&...);

 private:
  void instantiate_graph() {
    constexpr size_t error_log_size = 256;
    hipGraphNode_t error_node       = nullptr;
    char error_log[error_log_size];
    KOKKOS_IMPL_HIP_SAFE_CALL(hipGraphInstantiate(
        &m_graph_exec, m_graph, &error_node, error_log, error_log_size));
  }

  Kokkos::HIP m_execution_space;
  hipGraph_t m_graph          = nullptr;
  hipGraphExec_t m_graph_exec = nullptr;
};

inline GraphImpl<Kokkos::HIP>::~GraphImpl() {
  m_execution_space.fence("Kokkos::GraphImpl::~GraphImpl: Graph Destruction");
  KOKKOS_EXPECTS(m_graph);
  if (m_graph_exec) {
    KOKKOS_IMPL_HIP_SAFE_CALL(hipGraphExecDestroy(m_graph_exec));
  }
  KOKKOS_IMPL_HIP_SAFE_CALL(hipGraphDestroy(m_graph));
}

inline GraphImpl<Kokkos::HIP>::GraphImpl(Kokkos::HIP instance)
    : m_execution_space(std::move(instance)) {
  KOKKOS_IMPL_HIP_SAFE_CALL(hipGraphCreate(&m_graph, 0));
}

inline void GraphImpl<Kokkos::HIP>::add_node(
    std::shared_ptr<aggregate_node_impl_t> const& arg_node_ptr) {
  // All of the predecessors are just added as normal, so all we need to
  // do here is add an empty node
  KOKKOS_IMPL_HIP_SAFE_CALL(
      hipGraphAddEmptyNode(&(arg_node_ptr->node_details_t::node), m_graph,
                           /* dependencies = */ nullptr,
                           /* numDependencies = */ 0));
}

// Requires NodeImplPtr is a shared_ptr to specialization of GraphNodeImpl
// Also requires that the kernel has the graph node tag in its policy
template <class NodeImpl>
inline void GraphImpl<Kokkos::HIP>::add_node(
    std::shared_ptr<NodeImpl> const& arg_node_ptr) {
  static_assert(NodeImpl::kernel_type::Policy::is_graph_kernel::value);
  KOKKOS_EXPECTS(arg_node_ptr);
  // The Kernel launch from the execute() method has been shimmed to insert
  // the node into the graph
  auto& kernel = arg_node_ptr->get_kernel();
  auto& node   = static_cast<node_details_t*>(arg_node_ptr.get())->node;
  KOKKOS_EXPECTS(!node);
  kernel.set_hip_graph_ptr(&m_graph);
  kernel.set_hip_graph_node_ptr(&node);
  kernel.execute();
  KOKKOS_ENSURES(node);
}

// Requires PredecessorRef is a specialization of GraphNodeRef that has
// already been added to this graph and NodeImpl is a specialization of
// GraphNodeImpl that has already been added to this graph.
template <class NodeImplPtr, class PredecessorRef>
inline void GraphImpl<Kokkos::HIP>::add_predecessor(
    NodeImplPtr arg_node_ptr, PredecessorRef arg_pred_ref) {
  KOKKOS_EXPECTS(arg_node_ptr);
  auto pred_ptr = GraphAccess::get_node_ptr(arg_pred_ref);
  KOKKOS_EXPECTS(pred_ptr);

  auto const& pred_node = pred_ptr->node_details_t::node;
  KOKKOS_EXPECTS(pred_node);

  auto const& node = arg_node_ptr->node_details_t::node;
  KOKKOS_EXPECTS(node);

  KOKKOS_IMPL_HIP_SAFE_CALL(
      hipGraphAddDependencies(m_graph, &pred_node, &node, 1));
}

inline void GraphImpl<Kokkos::HIP>::submit() {
  if (!m_graph_exec) {
    instantiate_graph();
  }
  KOKKOS_IMPL_HIP_SAFE_CALL(
      hipGraphLaunch(m_graph_exec, m_execution_space.hip_stream()));
}

inline Kokkos::HIP const& GraphImpl<Kokkos::HIP>::get_execution_space() const
    noexcept {
  return m_execution_space;
}

inline auto GraphImpl<Kokkos::HIP>::create_root_node_ptr() {
  KOKKOS_EXPECTS(m_graph);
  KOKKOS_EXPECTS(!m_graph_exec);
  auto rv = std::make_shared<root_node_impl_t>(get_execution_space(),
                                               _graph_node_is_root_ctor_tag{});
  KOKKOS_IMPL_HIP_SAFE_CALL(hipGraphAddEmptyNode(&(rv->node_details_t::node),
                                                 m_graph,
                                                 /* dependencies = */ nullptr,
                                                 /* numDependencies = */ 0));
  KOKKOS_ENSURES(rv->node_details_t::node);
  return rv;
}

template <class... PredecessorRefs>
inline auto GraphImpl<Kokkos::HIP>::create_aggregate_ptr(PredecessorRefs&&...) {
  // The attachment to predecessors, which is all we really need, happens
  // in the generic layer, which calls through to add_predecessor for
  // each predecessor ref, so all we need to do here is create the (trivial)
  // aggregate node.
  return std::make_shared<aggregate_node_impl_t>(m_execution_space,
                                                 _graph_node_kernel_ctor_tag{},
                                                 aggregate_kernel_impl_t{});
}
}  // namespace Impl
}  // namespace Kokkos

#endif
