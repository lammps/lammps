/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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
