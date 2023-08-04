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

#ifndef KOKKOS_IMPL_KOKKOS_GRAPHIMPL_HPP
#define KOKKOS_IMPL_KOKKOS_GRAPHIMPL_HPP

#include <Kokkos_Macros.hpp>

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_Graph_fwd.hpp>

#include <Kokkos_Concepts.hpp>  // is_execution_policy
#include <Kokkos_PointerOwnership.hpp>
#include <impl/Kokkos_GraphImpl_fwd.hpp>

#include <memory>  // std::make_shared

namespace Kokkos {
namespace Impl {

struct GraphAccess {
  template <class ExecutionSpace>
  static Kokkos::Experimental::Graph<ExecutionSpace> construct_graph(
      ExecutionSpace ex) {
    //----------------------------------------//
    return Kokkos::Experimental::Graph<ExecutionSpace>{
        std::make_shared<GraphImpl<ExecutionSpace>>(std::move(ex))};
    //----------------------------------------//
  }
  template <class ExecutionSpace>
  static auto create_root_ref(
      Kokkos::Experimental::Graph<ExecutionSpace>& arg_graph) {
    auto const& graph_impl_ptr = arg_graph.m_impl_ptr;

    auto root_ptr = graph_impl_ptr->create_root_node_ptr();

    return Kokkos::Experimental::GraphNodeRef<ExecutionSpace>{
        graph_impl_ptr, std::move(root_ptr)};
  }

  template <class NodeType, class... Args>
  static auto make_node_shared_ptr(Args&&... args) {
    static_assert(
        Kokkos::Impl::is_specialization_of<NodeType, GraphNodeImpl>::value,
        "Kokkos Internal Error in graph interface");
    return std::make_shared<NodeType>((Args &&) args...);
  }

  template <class GraphImplWeakPtr, class ExecutionSpace, class Kernel,
            class Predecessor>
  static auto make_graph_node_ref(
      GraphImplWeakPtr graph_impl,
      std::shared_ptr<
          Kokkos::Impl::GraphNodeImpl<ExecutionSpace, Kernel, Predecessor>>
          pred_impl) {
    //----------------------------------------
    return Kokkos::Experimental::GraphNodeRef<ExecutionSpace, Kernel,
                                              Predecessor>{
        std::move(graph_impl), std::move(pred_impl)};
    //----------------------------------------
  }

  //----------------------------------------------------------------------------
  // <editor-fold desc="accessors for private members of public interface"> {{{2

  template <class NodeRef>
  static auto get_node_ptr(NodeRef&& node_ref) {
    static_assert(
        is_specialization_of<remove_cvref_t<NodeRef>,
                             Kokkos::Experimental::GraphNodeRef>::value,
        "Kokkos Internal Implementation error (bad argument to "
        "`GraphAccess::get_node_ptr()`)");
    return ((NodeRef &&) node_ref).get_node_ptr();
  }

  template <class NodeRef>
  static auto get_graph_weak_ptr(NodeRef&& node_ref) {
    static_assert(
        is_specialization_of<remove_cvref_t<NodeRef>,
                             Kokkos::Experimental::GraphNodeRef>::value,
        "Kokkos Internal Implementation error (bad argument to "
        "`GraphAccess::get_graph_weak_ptr()`)");
    return ((NodeRef &&) node_ref).get_graph_weak_ptr();
  }

  // </editor-fold> end accessors for private members of public interface }}}2
  //----------------------------------------------------------------------------
};

template <class Policy>
struct _add_graph_kernel_tag;

template <template <class...> class PolicyTemplate, class... PolicyTraits>
struct _add_graph_kernel_tag<PolicyTemplate<PolicyTraits...>> {
  using type = PolicyTemplate<PolicyTraits..., IsGraphKernelTag>;
};

}  // end namespace Impl

namespace Experimental {  // but not for users, so...

template <class Policy>
// requires ExecutionPolicy<Policy>
constexpr auto require(Policy const& policy,
                       Kokkos::Impl::KernelInGraphProperty) {
  static_assert(Kokkos::is_execution_policy<Policy>::value,
                "Internal implementation error!");
  return typename Kokkos::Impl::_add_graph_kernel_tag<Policy>::type{policy};
}

}  // end namespace Experimental

}  // end namespace Kokkos

#endif  // KOKKOS_IMPL_KOKKOS_GRAPHIMPL_HPP
