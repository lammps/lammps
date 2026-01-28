// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_KOKKOS_CUDA_GRAPH_IMPL_HPP
#define KOKKOS_KOKKOS_CUDA_GRAPH_IMPL_HPP

#include <Kokkos_Macros.hpp>

#if defined(KOKKOS_ENABLE_CUDA)

#include <Kokkos_Graph_fwd.hpp>

#include <impl/Kokkos_GraphImpl.hpp>  // GraphAccess needs to be complete

// GraphNodeImpl needs to be complete because GraphImpl here is a full
// specialization and not just a partial one
#include <impl/Kokkos_GraphNodeImpl.hpp>
#include <Cuda/Kokkos_Cuda_GraphNode_Impl.hpp>

#include <Cuda/Kokkos_Cuda.hpp>
#include <Cuda/Kokkos_Cuda_Error.hpp>
#include <Cuda/Kokkos_Cuda_Instance.hpp>

namespace Kokkos {
namespace Impl {

template <>
struct GraphImpl<Kokkos::Cuda> {
 public:
  using execution_space = Kokkos::Cuda;

 private:
  execution_space m_execution_space;
  cudaGraph_t m_graph          = nullptr;
  cudaGraphExec_t m_graph_exec = nullptr;

  bool m_graph_owning = false;

  using cuda_graph_flags_t = unsigned int;

  using node_details_t = GraphNodeBackendSpecificDetails<Kokkos::Cuda>;

  std::vector<std::shared_ptr<node_details_t>> m_nodes;

 public:
  void instantiate() {
    KOKKOS_EXPECTS(!m_graph_exec);
    KOKKOS_IMPL_CUDA_SAFE_CALL(
        (m_execution_space.impl_internal_space_instance()
             ->cuda_graph_instantiate_wrapper(&m_graph_exec, m_graph)));
    KOKKOS_ENSURES(m_graph_exec);
    // TODO @graphs print out errors
  }

  using root_node_impl_t =
      GraphNodeImpl<Kokkos::Cuda, Kokkos::Experimental::TypeErasedTag,
                    Kokkos::Experimental::TypeErasedTag>;
  using aggregate_impl_t = CudaGraphNodeAggregate;
  using aggregate_node_impl_t =
      GraphNodeImpl<Kokkos::Cuda, aggregate_impl_t,
                    Kokkos::Experimental::TypeErasedTag>;

  // Not movable or copyable; it spends its whole life as a shared_ptr in the
  // Graph object
  GraphImpl()                            = delete;
  GraphImpl(GraphImpl const&)            = delete;
  GraphImpl(GraphImpl&&)                 = delete;
  GraphImpl& operator=(GraphImpl const&) = delete;
  GraphImpl& operator=(GraphImpl&&)      = delete;
  ~GraphImpl() {
    // TODO @graphs we need to somehow indicate the need for a fence in the
    //              destructor of the GraphImpl object (so that we don't have to
    //              just always do it)
    m_execution_space.fence("Kokkos::GraphImpl::~GraphImpl: Graph Destruction");
    KOKKOS_EXPECTS(bool(m_graph))
    if (bool(m_graph_exec)) {
      KOKKOS_IMPL_CUDA_SAFE_CALL(
          (m_execution_space.impl_internal_space_instance()
               ->cuda_graph_exec_destroy_wrapper(m_graph_exec)));
    }
    if (m_graph_owning) {
      KOKKOS_IMPL_CUDA_SAFE_CALL(
          (m_execution_space.impl_internal_space_instance()
               ->cuda_graph_destroy_wrapper(m_graph)));
    }
  }

  explicit GraphImpl(Kokkos::Cuda arg_instance)
      : m_execution_space(std::move(arg_instance)), m_graph_owning(true) {
    KOKKOS_IMPL_CUDA_SAFE_CALL(
        (m_execution_space.impl_internal_space_instance()
             ->cuda_graph_create_wrapper(&m_graph, cuda_graph_flags_t{0})));
  }

  explicit GraphImpl(Kokkos::Cuda arg_instance, cudaGraph_t graph)
      : m_execution_space(std::move(arg_instance)),
        m_graph(graph),
        m_graph_owning(false) {
    KOKKOS_EXPECTS(graph != nullptr);
  }

  void add_node(std::shared_ptr<aggregate_node_impl_t> const& arg_node_ptr) {
    // All of the predecessors are just added as normal, so all we need to
    // do here is add an empty node
    KOKKOS_IMPL_CUDA_SAFE_CALL(
        (m_execution_space.impl_internal_space_instance()
             ->cuda_graph_add_empty_node_wrapper(
                 &(arg_node_ptr->node_details_t::node), m_graph,
                 /* dependencies = */ nullptr,
                 /* numDependencies = */ 0)));
  }

  template <class NodeImpl>
  std::enable_if_t<
      Kokkos::Impl::is_graph_kernel_v<typename NodeImpl::kernel_type>>
  add_node(std::shared_ptr<NodeImpl> arg_node_ptr) {
    static_assert(
        Kokkos::Impl::is_specialization_of_v<NodeImpl, GraphNodeImpl>);
    KOKKOS_EXPECTS(bool(arg_node_ptr));
    // The Kernel launch from the execute() method has been shimmed to insert
    // the node into the graph
    auto& kernel = arg_node_ptr->get_kernel();
    // note: using arg_node_ptr->node_details_t::node caused an ICE in NVCC 10.1
    auto& cuda_node = static_cast<node_details_t*>(arg_node_ptr.get())->node;
    KOKKOS_EXPECTS(!bool(cuda_node));
    kernel.set_cuda_graph_ptr(&m_graph);
    kernel.set_cuda_graph_node_ptr(&cuda_node);
    kernel.execute();
    KOKKOS_ENSURES(bool(cuda_node));
    m_nodes.push_back(std::move(arg_node_ptr));
  }

  template <class NodeImpl>
  std::enable_if_t<
      Kokkos::Impl::is_graph_capture_v<typename NodeImpl::kernel_type>>
  add_node(const Kokkos::Cuda& exec, std::shared_ptr<NodeImpl> arg_node_ptr) {
    static_assert(
        Kokkos::Impl::is_specialization_of_v<NodeImpl, GraphNodeImpl>);
    KOKKOS_EXPECTS(bool(arg_node_ptr));

    auto& kernel = arg_node_ptr->get_kernel();
    kernel.capture(exec, m_graph);
    static_cast<node_details_t*>(arg_node_ptr.get())->node = kernel.m_node;

    m_nodes.push_back(std::move(arg_node_ptr));
  }

  template <class NodeImpl>
  std::enable_if_t<
      Kokkos::Impl::is_graph_then_host_v<typename NodeImpl::kernel_type>>
  add_node(std::shared_ptr<NodeImpl> arg_node_ptr) {
    static_assert(
        Kokkos::Impl::is_specialization_of_v<NodeImpl, GraphNodeImpl>);
    KOKKOS_EXPECTS(bool(arg_node_ptr));

    auto& kernel = arg_node_ptr->get_kernel();
    kernel.add_to_graph(m_graph);
    static_cast<node_details_t*>(arg_node_ptr.get())->node = kernel.m_node;

    m_nodes.push_back(std::move(arg_node_ptr));
  }

  template <class NodeImplPtr, class PredecessorRef>
  // requires PredecessorRef is a specialization of GraphNodeRef that has
  // already been added to this graph and NodeImpl is a specialization of
  // GraphNodeImpl that has already been added to this graph.
  void add_predecessor(NodeImplPtr arg_node_ptr, PredecessorRef arg_pred_ref) {
    KOKKOS_EXPECTS(bool(arg_node_ptr))
    auto pred_ptr = GraphAccess::get_node_ptr(arg_pred_ref);
    KOKKOS_EXPECTS(bool(pred_ptr))

    // clang-format off
    // NOTE const-qualifiers below are commented out because of an API break
    // from CUDA 10.0 to CUDA 10.1
    // cudaGraphAddDependencies(cudaGraph_t, cudaGraphNode_t*, cudaGraphNode_t*, size_t)
    // cudaGraphAddDependencies(cudaGraph_t, const cudaGraphNode_t*, const cudaGraphNode_t*, size_t)
    // clang-format on
    auto /*const*/& pred_cuda_node = pred_ptr->node_details_t::node;
    KOKKOS_EXPECTS(bool(pred_cuda_node))

    auto /*const*/& cuda_node = arg_node_ptr->node_details_t::node;
    KOKKOS_EXPECTS(bool(cuda_node))

    KOKKOS_IMPL_CUDA_SAFE_CALL(
        (m_execution_space.impl_internal_space_instance()
             ->cuda_graph_add_dependencies_wrapper(m_graph, &pred_cuda_node,
                                                   &cuda_node, 1)));
  }

  void submit(const execution_space& exec) {
    if (!bool(m_graph_exec)) {
      instantiate();
    }
    KOKKOS_IMPL_CUDA_SAFE_CALL(
        (exec.impl_internal_space_instance()->cuda_graph_launch_wrapper(
            m_graph_exec)));
  }

  execution_space const& get_execution_space() const noexcept {
    return m_execution_space;
  }

  auto create_root_node_ptr() {
    KOKKOS_EXPECTS(bool(m_graph))
    KOKKOS_EXPECTS(!bool(m_graph_exec))
    auto rv = std::make_shared<root_node_impl_t>(
        get_execution_space(), _graph_node_is_root_ctor_tag{});
    KOKKOS_IMPL_CUDA_SAFE_CALL(
        (m_execution_space.impl_internal_space_instance()
             ->cuda_graph_add_empty_node_wrapper(&(rv->node_details_t::node),
                                                 m_graph,
                                                 /* dependencies = */ nullptr,
                                                 /* numDependencies = */ 0)));
    KOKKOS_ENSURES(bool(rv->node_details_t::node))
    return rv;
  }

  template <class... PredecessorRefs>
  // See requirements/expectations in GraphBuilder
  auto create_aggregate_ptr(PredecessorRefs&&...) {
    // The attachment to predecessors, which is all we really need, happens
    // in the generic layer, which calls through to add_predecessor for
    // each predecessor ref, so all we need to do here is create the (trivial)
    // aggregate node.
    return std::make_shared<aggregate_node_impl_t>(
        m_execution_space, _graph_node_kernel_ctor_tag{}, aggregate_impl_t{});
  }

  cudaGraph_t cuda_graph() { return m_graph; }
  cudaGraphExec_t cuda_graph_exec() { return m_graph_exec; }
};

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // defined(KOKKOS_ENABLE_CUDA)
#endif  // KOKKOS_KOKKOS_CUDA_GRAPH_IMPL_HPP
