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

#ifndef KOKKOS_SYCL_GRAPHNODEKERNEL_HPP
#define KOKKOS_SYCL_GRAPHNODEKERNEL_HPP

#include <Kokkos_Graph_fwd.hpp>

#include <impl/Kokkos_GraphImpl.hpp>

#include <Kokkos_Parallel.hpp>
#include <Kokkos_Parallel_Reduce.hpp>
#include <Kokkos_PointerOwnership.hpp>

#include <SYCL/Kokkos_SYCL_GraphNode_Impl.hpp>

namespace Kokkos {
namespace Impl {

template <typename PolicyType, typename Functor, typename PatternTag,
          typename... Args>
class GraphNodeKernelImpl<Kokkos::Experimental::SYCL, PolicyType, Functor,
                          PatternTag, Args...>
    : public PatternImplSpecializationFromTag<
          PatternTag, Functor, PolicyType, Args...,
          Kokkos::Experimental::SYCL>::type {
 public:
  using Policy       = PolicyType;
  using graph_kernel = GraphNodeKernelImpl;
  using base_t       = typename PatternImplSpecializationFromTag<
      PatternTag, Functor, Policy, Args..., Kokkos::Experimental::SYCL>::type;

  // TODO use the name and executionspace
  template <typename PolicyDeduced, typename... ArgsDeduced>
  GraphNodeKernelImpl(std::string, Kokkos::Experimental::SYCL const&,
                      Functor arg_functor, PolicyDeduced&& arg_policy,
                      ArgsDeduced&&... args)
      : base_t(std::move(arg_functor), (PolicyDeduced &&) arg_policy,
               (ArgsDeduced &&) args...) {}

  template <typename PolicyDeduced>
  GraphNodeKernelImpl(Kokkos::Experimental::SYCL const& exec_space,
                      Functor arg_functor, PolicyDeduced&& arg_policy)
      : GraphNodeKernelImpl("", exec_space, std::move(arg_functor),
                            (PolicyDeduced &&) arg_policy) {}

  void set_sycl_graph_ptr(
      sycl::ext::oneapi::experimental::command_graph<
          sycl::ext::oneapi::experimental::graph_state::modifiable>*
          arg_graph) {
    m_graph_ptr = arg_graph;
  }

  void set_sycl_graph_node_ptr(
      std::optional<sycl::ext::oneapi::experimental::node>* arg_node) {
    m_graph_node_ptr = arg_node;
  }

  std::optional<sycl::ext::oneapi::experimental::node>& get_sycl_graph_node()
      const {
    return *m_graph_node_ptr;
  }

  sycl::ext::oneapi::experimental::command_graph<
      sycl::ext::oneapi::experimental::graph_state::modifiable>&
  get_sycl_graph() const {
    return *m_graph_ptr;
  }

 private:
  Kokkos::ObservingRawPtr<sycl::ext::oneapi::experimental::command_graph<
      sycl::ext::oneapi::experimental::graph_state::modifiable>>
      m_graph_ptr = nullptr;
  Kokkos::ObservingRawPtr<std::optional<sycl::ext::oneapi::experimental::node>>
      m_graph_node_ptr = nullptr;
};

struct SYCLGraphNodeAggregateKernel {
  using graph_kernel = SYCLGraphNodeAggregateKernel;

  // Aggregates don't need a policy, but for the purposes of checking the static
  // assertions about graph kernels,
  struct Policy {
    using is_graph_kernel = std::true_type;
  };
};

template <typename KernelType,
          typename Tag =
              typename PatternTagFromImplSpecialization<KernelType>::type>
struct get_graph_node_kernel_type
    : type_identity<GraphNodeKernelImpl<
          Kokkos::Experimental::SYCL, typename KernelType::Policy,
          typename KernelType::functor_type, Tag>> {};

template <typename KernelType>
struct get_graph_node_kernel_type<KernelType, Kokkos::ParallelReduceTag>
    : type_identity<GraphNodeKernelImpl<
          Kokkos::Experimental::SYCL, typename KernelType::Policy,
          CombinedFunctorReducer<typename KernelType::FunctorType,
                                 typename KernelType::ReducerType>,
          Kokkos::ParallelReduceTag>> {};

template <typename KernelType>
auto& get_sycl_graph_from_kernel(KernelType const& kernel) {
  using graph_node_kernel_t =
      typename get_graph_node_kernel_type<KernelType>::type;
  auto const& kernel_as_graph_kernel =
      static_cast<graph_node_kernel_t const&>(kernel);
  auto& graph = kernel_as_graph_kernel.get_sycl_graph();

  return graph;
}

template <typename KernelType>
auto& get_sycl_graph_node_from_kernel(KernelType const& kernel) {
  using graph_node_kernel_t =
      typename get_graph_node_kernel_type<KernelType>::type;
  auto const& kernel_as_graph_kernel =
      static_cast<graph_node_kernel_t const&>(kernel);
  auto& graph_node = kernel_as_graph_kernel.get_sycl_graph_node();

  return graph_node;
}

template <typename Kernel, typename Lambda>
void sycl_attach_kernel_to_node(Kernel& kernel, const Lambda& lambda) {
  sycl::ext::oneapi::experimental::command_graph<
      sycl::ext::oneapi::experimental::graph_state::modifiable>& graph =
      Impl::get_sycl_graph_from_kernel(kernel);
  std::optional<sycl::ext::oneapi::experimental::node>& graph_node =
      Impl::get_sycl_graph_node_from_kernel(kernel);
  KOKKOS_ENSURES(!graph_node);
  graph_node = graph.add(lambda);
  KOKKOS_ENSURES(graph_node);
  // FIXME_SYCL_GRAPH not yet implemented in the compiler
  //   KOKKOS_ENSURES(graph_node.get_type() ==
  //   sycl::ext::oneapi::experimental::node_type::kernel)
}

}  // namespace Impl
}  // namespace Kokkos

#endif
