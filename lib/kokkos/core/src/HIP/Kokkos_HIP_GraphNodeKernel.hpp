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

#ifndef KOKKOS_HIP_GRAPHNODEKERNEL_HPP
#define KOKKOS_HIP_GRAPHNODEKERNEL_HPP

#include <Kokkos_Graph_fwd.hpp>

#include <impl/Kokkos_GraphImpl.hpp>
#include <impl/Kokkos_SharedAlloc.hpp>

#include <Kokkos_Parallel.hpp>
#include <Kokkos_Parallel_Reduce.hpp>
#include <Kokkos_PointerOwnership.hpp>

#include <HIP/Kokkos_HIP_SharedAllocationRecord.hpp>
#include <HIP/Kokkos_HIP_GraphNode_Impl.hpp>

namespace Kokkos {
namespace Impl {

template <typename PolicyType, typename Functor, typename PatternTag,
          typename... Args>
class GraphNodeKernelImpl<Kokkos::HIP, PolicyType, Functor, PatternTag, Args...>
    : public PatternImplSpecializationFromTag<PatternTag, Functor, PolicyType,
                                              Args..., Kokkos::HIP>::type {
 public:
  using Policy       = PolicyType;
  using graph_kernel = GraphNodeKernelImpl;
  using base_t =
      typename PatternImplSpecializationFromTag<PatternTag, Functor, Policy,
                                                Args..., Kokkos::HIP>::type;
  using Record = Kokkos::Impl::SharedAllocationRecord<Kokkos::HIPSpace, void>;

  // TODO use the name and executionspace
  template <typename PolicyDeduced, typename... ArgsDeduced>
  GraphNodeKernelImpl(std::string, Kokkos::HIP const&, Functor arg_functor,
                      PolicyDeduced&& arg_policy, ArgsDeduced&&... args)
      : base_t(std::move(arg_functor), (PolicyDeduced &&) arg_policy,
               (ArgsDeduced &&) args...) {}

  template <typename PolicyDeduced>
  GraphNodeKernelImpl(Kokkos::HIP const& exec_space, Functor arg_functor,
                      PolicyDeduced&& arg_policy)
      : GraphNodeKernelImpl("", exec_space, std::move(arg_functor),
                            (PolicyDeduced &&) arg_policy) {}

  ~GraphNodeKernelImpl() {
    if (m_driver_storage) {
      Record::decrement(Record::get_record(m_driver_storage));
    }
  }

  void set_hip_graph_ptr(hipGraph_t* arg_graph_ptr) {
    m_graph_ptr = arg_graph_ptr;
  }

  void set_hip_graph_node_ptr(hipGraphNode_t* arg_node_ptr) {
    m_graph_node_ptr = arg_node_ptr;
  }

  hipGraphNode_t* get_hip_graph_node_ptr() const { return m_graph_node_ptr; }

  hipGraph_t const* get_hip_graph_ptr() const { return m_graph_ptr; }

  Kokkos::ObservingRawPtr<base_t> allocate_driver_memory_buffer() const {
    KOKKOS_EXPECTS(m_driver_storage == nullptr);

    auto* record = Record::allocate(
        Kokkos::HIPSpace{}, "GraphNodeKernel global memory functor storage",
        sizeof(base_t));

    Record::increment(record);
    m_driver_storage = reinterpret_cast<base_t*>(record->data());
    KOKKOS_ENSURES(m_driver_storage != nullptr);

    return m_driver_storage;
  }

 private:
  Kokkos::ObservingRawPtr<const hipGraph_t> m_graph_ptr    = nullptr;
  Kokkos::ObservingRawPtr<hipGraphNode_t> m_graph_node_ptr = nullptr;
  Kokkos::OwningRawPtr<base_t> m_driver_storage            = nullptr;
};

struct HIPGraphNodeAggregateKernel {
  using graph_kernel = HIPGraphNodeAggregateKernel;

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
    : type_identity<
          GraphNodeKernelImpl<Kokkos::HIP, typename KernelType::Policy,
                              typename KernelType::functor_type, Tag>> {};

template <typename KernelType>
struct get_graph_node_kernel_type<KernelType, Kokkos::ParallelReduceTag>
    : type_identity<GraphNodeKernelImpl<
          Kokkos::HIP, typename KernelType::Policy,
          CombinedFunctorReducer<typename KernelType::functor_type,
                                 typename KernelType::reducer_type>,
          Kokkos::ParallelReduceTag>> {};

template <typename KernelType>
auto* allocate_driver_storage_for_kernel(KernelType const& kernel) {
  using graph_node_kernel_t =
      typename get_graph_node_kernel_type<KernelType>::type;
  auto const& kernel_as_graph_kernel =
      static_cast<graph_node_kernel_t const&>(kernel);

  return kernel_as_graph_kernel.allocate_driver_memory_buffer();
}

template <typename KernelType>
auto const& get_hip_graph_from_kernel(KernelType const& kernel) {
  using graph_node_kernel_t =
      typename get_graph_node_kernel_type<KernelType>::type;
  auto const& kernel_as_graph_kernel =
      static_cast<graph_node_kernel_t const&>(kernel);
  hipGraph_t const* graph_ptr = kernel_as_graph_kernel.get_hip_graph_ptr();
  KOKKOS_EXPECTS(graph_ptr != nullptr);

  return *graph_ptr;
}

template <typename KernelType>
auto& get_hip_graph_node_from_kernel(KernelType const& kernel) {
  using graph_node_kernel_t =
      typename get_graph_node_kernel_type<KernelType>::type;
  auto const& kernel_as_graph_kernel =
      static_cast<graph_node_kernel_t const&>(kernel);
  auto* graph_node_ptr = kernel_as_graph_kernel.get_hip_graph_node_ptr();
  KOKKOS_EXPECTS(graph_node_ptr != nullptr);

  return *graph_node_ptr;
}
}  // namespace Impl
}  // namespace Kokkos

#endif
