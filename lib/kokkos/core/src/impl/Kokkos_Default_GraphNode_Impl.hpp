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

#ifndef KOKKOS_KOKKOS_HOST_GRAPHNODE_IMPL_HPP
#define KOKKOS_KOKKOS_HOST_GRAPHNODE_IMPL_HPP

#include <Kokkos_Macros.hpp>

#include <impl/Kokkos_Default_Graph_fwd.hpp>

#include <Kokkos_Graph.hpp>

#include <vector>
#include <memory>

namespace Kokkos {
namespace Impl {

//==============================================================================
// <editor-fold desc="GraphNodeBackendSpecificDetails"> {{{1

template <class ExecutionSpace>
struct GraphNodeBackendSpecificDetails {
 private:
  using execution_space_instance_storage_t =
      ExecutionSpaceInstanceStorage<ExecutionSpace>;
  using default_kernel_impl_t = GraphNodeKernelDefaultImpl<ExecutionSpace>;
  using default_aggregate_kernel_impl_t =
      GraphNodeAggregateKernelDefaultImpl<ExecutionSpace>;

  std::vector<std::shared_ptr<GraphNodeBackendSpecificDetails<ExecutionSpace>>>
      m_predecessors = {};

  Kokkos::ObservingRawPtr<default_kernel_impl_t> m_kernel_ptr = nullptr;

  bool m_has_executed = false;
  bool m_is_aggregate = false;
  bool m_is_root      = false;

  template <class>
  friend struct HostGraphImpl;

 protected:
  //----------------------------------------------------------------------------
  // <editor-fold desc="Ctors, destructor, and assignment"> {{{2

  explicit GraphNodeBackendSpecificDetails() = default;

  explicit GraphNodeBackendSpecificDetails(
      _graph_node_is_root_ctor_tag) noexcept
      : m_has_executed(true), m_is_root(true) {}

  GraphNodeBackendSpecificDetails(GraphNodeBackendSpecificDetails const&) =
      delete;

  GraphNodeBackendSpecificDetails(GraphNodeBackendSpecificDetails&&) noexcept =
      delete;

  GraphNodeBackendSpecificDetails& operator   =(
      GraphNodeBackendSpecificDetails const&) = delete;

  GraphNodeBackendSpecificDetails& operator       =(
      GraphNodeBackendSpecificDetails&&) noexcept = delete;

  ~GraphNodeBackendSpecificDetails() = default;

  // </editor-fold> end Ctors, destructor, and assignment }}}2
  //----------------------------------------------------------------------------

 public:
  void set_kernel(default_kernel_impl_t& arg_kernel) {
    KOKKOS_EXPECTS(m_kernel_ptr == nullptr)
    m_kernel_ptr = &arg_kernel;
  }

  void set_kernel(default_aggregate_kernel_impl_t& arg_kernel) {
    KOKKOS_EXPECTS(m_kernel_ptr == nullptr)
    m_kernel_ptr   = &arg_kernel;
    m_is_aggregate = true;
  }

  void set_predecessor(
      std::shared_ptr<GraphNodeBackendSpecificDetails<ExecutionSpace>>
          arg_pred_impl) {
    // This method delegates responsibility for executing the predecessor to
    // this node.  Each node can have at most one predecessor (which may be an
    // aggregate).
    KOKKOS_EXPECTS(m_predecessors.empty() || m_is_aggregate)
    KOKKOS_EXPECTS(bool(arg_pred_impl))
    KOKKOS_EXPECTS(!m_has_executed)
    m_predecessors.push_back(std::move(arg_pred_impl));
  }

  void execute_node() {
    // This node could have already been executed as the predecessor of some
    // other
    KOKKOS_EXPECTS(bool(m_kernel_ptr) || m_has_executed)
    // Just execute the predecessor here, since calling set_predecessor()
    // delegates the responsibility for running it to us.
    if (!m_has_executed) {
      // I'm pretty sure this doesn't need to be atomic under our current
      // supported semantics, but instinct I have feels like it should be...
      m_has_executed = true;
      for (auto const& predecessor : m_predecessors) {
        predecessor->execute_node();
      }
      m_kernel_ptr->execute_kernel();
    }
    KOKKOS_ENSURES(m_has_executed)
  }

  // This is gross, but for the purposes of our simple default implementation...
  void reset_has_executed() {
    for (auto const& predecessor : m_predecessors) {
      predecessor->reset_has_executed();
    }
    // more readable, probably:
    //   if(!m_is_root) m_has_executed = false;
    m_has_executed = m_is_root;
  }
};

// </editor-fold> end GraphNodeBackendSpecificDetails }}}1
//==============================================================================

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_HOST_GRAPHNODE_IMPL_HPP
