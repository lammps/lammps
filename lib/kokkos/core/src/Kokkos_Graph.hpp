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

#ifndef KOKKOS_GRAPH_HPP
#define KOKKOS_GRAPH_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_GRAPH
#endif

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_Error.hpp>  // KOKKOS_EXPECTS

#include <Kokkos_Graph_fwd.hpp>
#include <impl/Kokkos_GraphImpl_fwd.hpp>

// GraphAccess needs to be defined, not just declared
#include <impl/Kokkos_GraphImpl.hpp>

#include <functional>
#include <memory>

namespace Kokkos {
namespace Experimental {

//==============================================================================
// <editor-fold desc="Graph"> {{{1

template <class ExecutionSpace = DefaultExecutionSpace>
struct [[nodiscard]] Graph {
  static_assert(Kokkos::is_execution_space_v<ExecutionSpace>);

 public:
  //----------------------------------------------------------------------------
  // <editor-fold desc="public member types"> {{{2

  using execution_space = ExecutionSpace;
  using graph           = Graph;
  using root_t          = GraphNodeRef<ExecutionSpace>;

  // </editor-fold> end public member types }}}2
  //----------------------------------------------------------------------------

 private:
  //----------------------------------------------------------------------------
  // <editor-fold desc="friends"> {{{2

  friend struct Kokkos::Impl::GraphAccess;

  // </editor-fold> end friends }}}2
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // <editor-fold desc="private data members"> {{{2

  using impl_t      = Kokkos::Impl::GraphImpl<ExecutionSpace>;
  using root_impl_t = typename impl_t::root_node_impl_t;

  std::shared_ptr<impl_t> m_impl_ptr  = nullptr;
  std::shared_ptr<root_impl_t> m_root = nullptr;

  // </editor-fold> end private data members }}}2
  //----------------------------------------------------------------------------

 public:
  // Construct an empty graph with a root node.
  Graph(ExecutionSpace exec = ExecutionSpace{})
      : m_impl_ptr{std::make_shared<impl_t>(std::move(exec))},
        m_root{m_impl_ptr->create_root_node_ptr()} {}

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || \
    defined(KOKKOS_ENABLE_SYCL)
  // Construct a graph from a native graph, add a root node.
  template <typename T>
#if defined(KOKKOS_ENABLE_CXX20)
    requires std::same_as<ExecutionSpace, Kokkos::DefaultExecutionSpace>
#endif
  Graph(ExecutionSpace exec, T&& native_graph)
      : m_impl_ptr{std::make_shared<impl_t>(std::move(exec),
                                            std::forward<T>(native_graph))},
        m_root{m_impl_ptr->create_root_node_ptr()} {
  }
#endif

  ExecutionSpace const& get_execution_space() const {
    return m_impl_ptr->get_execution_space();
  }

  // Once the graph is instantiated, it is undefined behavior to add nodes.
  // TODO Add a locking mechanism to avoid users shooting themselves
  //      in the foot.
  void instantiate() {
    KOKKOS_EXPECTS(bool(m_impl_ptr))
    (*m_impl_ptr).instantiate();
  }

  auto root_node() const { return root_t{m_impl_ptr, m_root}; }

  void submit(const execution_space& exec) const {
    KOKKOS_EXPECTS(bool(m_impl_ptr))
    (*m_impl_ptr).submit(exec);
  }

  void submit() const { submit(get_execution_space()); }

  decltype(auto) native_graph();

  decltype(auto) native_graph_exec();
};

// </editor-fold> end Graph }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="when_all"> {{{1

template <class... PredecessorRefs>
// constraints (not intended for subsumption, though...)
//   ((remove_cvref_t<PredecessorRefs> is a specialization of
//        GraphNodeRef with get_root().get_graph_impl() as its GraphImpl)
//      && ...)
auto when_all(PredecessorRefs&&... arg_pred_refs) {
  // TODO @graph @desul-integration check the constraints and preconditions
  //                                once we have folded conjunctions from
  //                                desul
  static_assert(sizeof...(PredecessorRefs) > 0,
                "when_all() needs at least one predecessor.");
  auto graph_ptr_impl =
      Kokkos::Impl::GraphAccess::get_graph_weak_ptr(
          std::get<0>(std::forward_as_tuple(arg_pred_refs...)))
          .lock();
  auto node_ptr_impl = graph_ptr_impl->create_aggregate_ptr(arg_pred_refs...);
  graph_ptr_impl->add_node(node_ptr_impl);
  (graph_ptr_impl->add_predecessor(node_ptr_impl, arg_pred_refs), ...);
  return Kokkos::Impl::GraphAccess::make_graph_node_ref(
      std::move(graph_ptr_impl), std::move(node_ptr_impl));
}

// </editor-fold> end when_all }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="create_graph"> {{{1

template <class ExecutionSpace, class Closure>
Graph<ExecutionSpace> create_graph(ExecutionSpace ex, Closure&& arg_closure) {
  // Create a shared pointer to the graph:
  // We need an attorney class here so we have an implementation friend to
  // create a Graph class without graph having public constructors. We can't
  // just make `create_graph` itself a friend because of the way that friend
  // function template injection works.
  Graph<ExecutionSpace> rv{std::move(ex)};
  // Invoke the user's graph construction closure
  ((Closure&&)arg_closure)(rv.root_node());
  // and given them back the graph
  // KOKKOS_ENSURES(rv.m_impl_ptr.use_count() == 1)
  return rv;
}

template <
    class ExecutionSpace = DefaultExecutionSpace,
    class Closure = Kokkos::Impl::DoNotExplicitlySpecifyThisTemplateParameter>
std::enable_if_t<
    !Kokkos::is_execution_space_v<Kokkos::Impl::remove_cvref_t<Closure>>,
    Graph<ExecutionSpace>>
create_graph(Closure&& arg_closure) {
  return create_graph(ExecutionSpace{}, (Closure&&)arg_closure);
}

// </editor-fold> end create_graph }}}1
//==============================================================================

template <class ExecutionSpace>
decltype(auto) Graph<ExecutionSpace>::native_graph() {
  KOKKOS_EXPECTS(bool(m_impl_ptr));
#if defined(KOKKOS_ENABLE_CUDA)
  if constexpr (std::is_same_v<ExecutionSpace, Kokkos::Cuda>) {
    return m_impl_ptr->cuda_graph();
  }
#elif defined(KOKKOS_ENABLE_HIP) && defined(KOKKOS_IMPL_HIP_NATIVE_GRAPH)
  if constexpr (std::is_same_v<ExecutionSpace, Kokkos::HIP>) {
    return m_impl_ptr->hip_graph();
  }
#elif defined(KOKKOS_ENABLE_SYCL) && defined(KOKKOS_IMPL_SYCL_GRAPH_SUPPORT)
  if constexpr (std::is_same_v<ExecutionSpace, Kokkos::SYCL>) {
    return m_impl_ptr->sycl_graph();
  }
#endif
}

template <class ExecutionSpace>
decltype(auto) Graph<ExecutionSpace>::native_graph_exec() {
  KOKKOS_EXPECTS(bool(m_impl_ptr));
#if defined(KOKKOS_ENABLE_CUDA)
  if constexpr (std::is_same_v<ExecutionSpace, Kokkos::Cuda>) {
    return m_impl_ptr->cuda_graph_exec();
  }
#elif defined(KOKKOS_ENABLE_HIP) && defined(KOKKOS_IMPL_HIP_NATIVE_GRAPH)
  if constexpr (std::is_same_v<ExecutionSpace, Kokkos::HIP>) {
    return m_impl_ptr->hip_graph_exec();
  }
#elif defined(KOKKOS_ENABLE_SYCL) && defined(KOKKOS_IMPL_SYCL_GRAPH_SUPPORT)
  if constexpr (std::is_same_v<ExecutionSpace, Kokkos::SYCL>) {
    return m_impl_ptr->sycl_graph_exec();
  }
#endif
}

}  // end namespace Experimental
}  // namespace Kokkos

// Even though these things are separable, include them here for now so that
// the user only needs to include Kokkos_Graph.hpp to get the whole facility.
#include <Kokkos_GraphNode.hpp>

#include <impl/Kokkos_GraphNodeImpl.hpp>
#include <impl/Kokkos_GraphNodeThenImpl.hpp>
#include <impl/Kokkos_Default_Graph_Impl.hpp>
#include <Cuda/Kokkos_Cuda_Graph_Impl.hpp>
#if defined(KOKKOS_ENABLE_HIP)
// The implementation of hipGraph in ROCm 5.2 is bugged, so we cannot use it.
#if defined(KOKKOS_IMPL_HIP_NATIVE_GRAPH)
#include <HIP/Kokkos_HIP_Graph_Impl.hpp>
#endif
#endif
#ifdef KOKKOS_IMPL_SYCL_GRAPH_SUPPORT
#include <SYCL/Kokkos_SYCL_Graph_Impl.hpp>
#endif
#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_GRAPH
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_GRAPH
#endif
#endif  // KOKKOS_GRAPH_HPP
