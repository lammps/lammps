// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_HIP_GRAPH_IMPL_HPP
#define KOKKOS_HIP_GRAPH_IMPL_HPP

#include <Kokkos_Macros.hpp>

#include <Kokkos_Graph_fwd.hpp>

#include <impl/Kokkos_GraphImpl.hpp>
#include <impl/Kokkos_GraphNodeImpl.hpp>

#include <HIP/Kokkos_HIP_GraphNodeKernel.hpp>
#include <HIP/Kokkos_HIP_Instance.hpp>

namespace Kokkos {
namespace Impl {
template <>
class GraphImpl<Kokkos::HIP> {
 public:
  using node_details_t = GraphNodeBackendSpecificDetails<Kokkos::HIP>;
  using root_node_impl_t =
      GraphNodeImpl<Kokkos::HIP, Kokkos::Experimental::TypeErasedTag,
                    Kokkos::Experimental::TypeErasedTag>;
  using aggregate_impl_t = HIPGraphNodeAggregate;
  using aggregate_node_impl_t =
      GraphNodeImpl<Kokkos::HIP, aggregate_impl_t,
                    Kokkos::Experimental::TypeErasedTag>;

  // Not movable or copyable; it spends its whole life as a shared_ptr in the
  // Graph object.
  GraphImpl()                            = delete;
  GraphImpl(GraphImpl const&)            = delete;
  GraphImpl(GraphImpl&&)                 = delete;
  GraphImpl& operator=(GraphImpl const&) = delete;
  GraphImpl& operator=(GraphImpl&&)      = delete;

  ~GraphImpl();

  explicit GraphImpl(Kokkos::HIP instance);

  GraphImpl(Kokkos::HIP instance, hipGraph_t graph);

  void add_node(std::shared_ptr<aggregate_node_impl_t> const& arg_node_ptr);

  template <class NodeImpl>
  std::enable_if_t<
      Kokkos::Impl::is_graph_kernel_v<typename NodeImpl::kernel_type>>
  add_node(std::shared_ptr<NodeImpl> arg_node_ptr);

  template <class NodeImpl>
  std::enable_if_t<
      Kokkos::Impl::is_graph_capture_v<typename NodeImpl::kernel_type>>
  add_node(const Kokkos::HIP& exec, std::shared_ptr<NodeImpl> arg_node_ptr);

  template <class NodeImpl>
  std::enable_if_t<
      Kokkos::Impl::is_graph_then_host_v<typename NodeImpl::kernel_type>>
  add_node(std::shared_ptr<NodeImpl> arg_node_ptr);

  template <class NodeImplPtr, class PredecessorRef>
  void add_predecessor(NodeImplPtr arg_node_ptr, PredecessorRef arg_pred_ref);

  void submit(const Kokkos::HIP& exec);

  Kokkos::HIP const& get_execution_space() const noexcept;

  auto create_root_node_ptr();

  template <class... PredecessorRefs>
  auto create_aggregate_ptr(PredecessorRefs&&...);

  void instantiate() {
    KOKKOS_EXPECTS(!m_graph_exec);
    KOKKOS_IMPL_HIP_SAFE_CALL(
        m_execution_space.impl_internal_space_instance()
            ->hip_graph_instantiate_wrapper(&m_graph_exec, m_graph, nullptr,
                                            nullptr, 0));
    KOKKOS_ENSURES(m_graph_exec);
  }

  hipGraph_t hip_graph() { return m_graph; }
  hipGraphExec_t hip_graph_exec() { return m_graph_exec; }

 private:
  Kokkos::HIP m_execution_space;
  hipGraph_t m_graph          = nullptr;
  hipGraphExec_t m_graph_exec = nullptr;

  bool m_graph_owning = false;

  std::vector<std::shared_ptr<node_details_t>> m_nodes;
};

inline GraphImpl<Kokkos::HIP>::~GraphImpl() {
  m_execution_space.fence("Kokkos::GraphImpl::~GraphImpl: Graph Destruction");
  KOKKOS_EXPECTS(m_graph);
  if (m_graph_exec) {
    KOKKOS_IMPL_HIP_SAFE_CALL(
        m_execution_space.impl_internal_space_instance()
            ->hip_graph_exec_destroy_wrapper(m_graph_exec));
  }
  if (m_graph_owning) {
    KOKKOS_IMPL_HIP_SAFE_CALL(m_execution_space.impl_internal_space_instance()
                                  ->hip_graph_destroy_wrapper(m_graph));
  }
}

inline GraphImpl<Kokkos::HIP>::GraphImpl(Kokkos::HIP instance)
    : m_execution_space(std::move(instance)), m_graph_owning(true) {
  KOKKOS_IMPL_HIP_SAFE_CALL(m_execution_space.impl_internal_space_instance()
                                ->hip_graph_create_wrapper(&m_graph, 0));
}

inline GraphImpl<Kokkos::HIP>::GraphImpl(Kokkos::HIP instance, hipGraph_t graph)
    : m_execution_space(std::move(instance)),
      m_graph(graph),
      m_graph_owning(false) {
  KOKKOS_EXPECTS(graph != nullptr);
}

inline void GraphImpl<Kokkos::HIP>::add_node(
    std::shared_ptr<aggregate_node_impl_t> const& arg_node_ptr) {
  // All of the predecessors are just added as normal, so all we need to
  // do here is add an empty node
  KOKKOS_IMPL_HIP_SAFE_CALL(m_execution_space.impl_internal_space_instance()
                                ->hip_graph_add_empty_node_wrapper(
                                    &(arg_node_ptr->node_details_t::node),
                                    m_graph,
                                    /* dependencies = */ nullptr,
                                    /* numDependencies = */ 0));
}

template <class NodeImpl>
inline std::enable_if_t<
    Kokkos::Impl::is_graph_kernel_v<typename NodeImpl::kernel_type>>
GraphImpl<Kokkos::HIP>::add_node(std::shared_ptr<NodeImpl> arg_node_ptr) {
  static_assert(Kokkos::Impl::is_specialization_of_v<NodeImpl, GraphNodeImpl>);
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
  m_nodes.push_back(std::move(arg_node_ptr));
}

template <class NodeImpl>
inline std::enable_if_t<
    Kokkos::Impl::is_graph_capture_v<typename NodeImpl::kernel_type>>
GraphImpl<Kokkos::HIP>::add_node(const Kokkos::HIP& exec,
                                 std::shared_ptr<NodeImpl> arg_node_ptr) {
  static_assert(Kokkos::Impl::is_specialization_of_v<NodeImpl, GraphNodeImpl>);
  KOKKOS_EXPECTS(bool(arg_node_ptr));

  auto& kernel = arg_node_ptr->get_kernel();
  kernel.capture(exec, m_graph);
  static_cast<node_details_t*>(arg_node_ptr.get())->node = kernel.m_node;

  m_nodes.push_back(std::move(arg_node_ptr));
}

template <class NodeImpl>
inline std::enable_if_t<
    Kokkos::Impl::is_graph_then_host_v<typename NodeImpl::kernel_type>>
GraphImpl<Kokkos::HIP>::add_node(std::shared_ptr<NodeImpl> arg_node_ptr) {
  static_assert(Kokkos::Impl::is_specialization_of_v<NodeImpl, GraphNodeImpl>);
  KOKKOS_EXPECTS(bool(arg_node_ptr));

  auto& kernel = arg_node_ptr->get_kernel();
  kernel.add_to_graph(m_graph);
  static_cast<node_details_t*>(arg_node_ptr.get())->node = kernel.m_node;

  m_nodes.push_back(std::move(arg_node_ptr));
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
      m_execution_space.impl_internal_space_instance()
          ->hip_graph_add_dependencies_wrapper(m_graph, &pred_node, &node, 1));
}

inline void GraphImpl<Kokkos::HIP>::submit(const Kokkos::HIP& exec) {
  if (!m_graph_exec) {
    instantiate();
  }
  KOKKOS_IMPL_HIP_SAFE_CALL(
      exec.impl_internal_space_instance()->hip_graph_launch_wrapper(
          m_graph_exec));
}

inline Kokkos::HIP const& GraphImpl<Kokkos::HIP>::get_execution_space()
    const noexcept {
  return m_execution_space;
}

inline auto GraphImpl<Kokkos::HIP>::create_root_node_ptr() {
  KOKKOS_EXPECTS(m_graph);
  KOKKOS_EXPECTS(!m_graph_exec);
  auto rv = std::make_shared<root_node_impl_t>(get_execution_space(),
                                               _graph_node_is_root_ctor_tag{});
  KOKKOS_IMPL_HIP_SAFE_CALL(m_execution_space.impl_internal_space_instance()
                                ->hip_graph_add_empty_node_wrapper(
                                    &(rv->node_details_t::node), m_graph,
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
  return std::make_shared<aggregate_node_impl_t>(
      m_execution_space, _graph_node_kernel_ctor_tag{}, aggregate_impl_t{});
}
}  // namespace Impl
}  // namespace Kokkos

#endif
