// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_HIP_GRAPHNODEKERNEL_HPP
#define KOKKOS_HIP_GRAPHNODEKERNEL_HPP

#include <Kokkos_Graph_fwd.hpp>

#include <impl/Kokkos_GraphImpl.hpp>

#include <Kokkos_Parallel.hpp>
#include <Kokkos_Parallel_Reduce.hpp>

#include <HIP/Kokkos_HIP_GraphNode_Impl.hpp>

namespace Kokkos {
namespace Impl {

template <typename Functor>
struct GraphNodeThenHostImpl<Kokkos::HIP, Functor> {
  Functor m_functor;
  hipGraphNode_t m_node = nullptr;

  explicit GraphNodeThenHostImpl(Functor functor)
      : m_functor(std::move(functor)) {}

  static void callback(void* data) {
    reinterpret_cast<Functor*>(data)->operator()();
  }

  void add_to_graph(hipGraph_t graph) {
    hipHostNodeParams params = {};
    params.fn                = callback;
    params.userData          = &m_functor;

    KOKKOS_IMPL_HIP_SAFE_CALL(
        hipGraphAddHostNode(&m_node, graph, nullptr, 0, &params));
  }
};

template <typename Functor>
struct GraphNodeCaptureImpl<Kokkos::HIP, Functor> {
  Functor m_functor;
  hipGraphNode_t m_node = nullptr;

  void capture(const Kokkos::HIP& exec, hipGraph_t graph) {
    KOKKOS_IMPL_HIP_SAFE_CALL(
        hipStreamBeginCapture(exec.hip_stream(), hipStreamCaptureModeGlobal));

    m_functor(exec);

    hipGraph_t captured_subgraph(nullptr);

    KOKKOS_IMPL_HIP_SAFE_CALL(
        hipStreamEndCapture(exec.hip_stream(), &captured_subgraph));

    KOKKOS_IMPL_HIP_SAFE_CALL(hipGraphAddChildGraphNode(&m_node, graph, nullptr,
                                                        0, captured_subgraph));
  }
};

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

 private:
  struct DriverStorageDeleter {
    std::string label;
    HIPSpace mem;

    void operator()(base_t* const ptr) const {
      mem.deallocate(label.c_str(), ptr, sizeof(base_t));
    }
  };

 public:
  template <typename PolicyDeduced, typename... ArgsDeduced>
  GraphNodeKernelImpl(std::string label_, HIP const&, Functor arg_functor,
                      PolicyDeduced&& arg_policy, ArgsDeduced&&... args)
      : base_t(std::move(arg_functor), (PolicyDeduced&&)arg_policy,
               (ArgsDeduced&&)args...),
        label(std::move(label_)) {}

  template <typename PolicyDeduced>
  GraphNodeKernelImpl(Kokkos::HIP const& exec_space, Functor arg_functor,
                      PolicyDeduced&& arg_policy)
      : GraphNodeKernelImpl("[unlabeled]", exec_space, std::move(arg_functor),
                            (PolicyDeduced&&)arg_policy) {}

  void set_hip_graph_ptr(hipGraph_t* arg_graph_ptr) {
    m_graph_ptr = arg_graph_ptr;
  }

  void set_hip_graph_node_ptr(hipGraphNode_t* arg_node_ptr) {
    m_graph_node_ptr = arg_node_ptr;
  }

  hipGraphNode_t* get_hip_graph_node_ptr() const { return m_graph_node_ptr; }

  hipGraph_t const* get_hip_graph_ptr() const { return m_graph_ptr; }

  base_t* allocate_driver_memory_buffer() const {
    KOKKOS_EXPECTS(m_driver_storage == nullptr);

    const auto& exec = this->get_policy().space();

    std::string alloc_label =
        label + " - GraphNodeKernel global memory functor storage";
    auto mem =
        Kokkos::HIPSpace::impl_create(exec.hip_device(), exec.hip_stream());
    auto* ptr = static_cast<base_t*>(
        mem.allocate(exec, alloc_label.c_str(), sizeof(base_t)));
    m_driver_storage = std::unique_ptr<base_t, DriverStorageDeleter>(
        ptr, DriverStorageDeleter{.label = std::move(alloc_label),
                                  .mem   = std::move(mem)});
    KOKKOS_ENSURES(m_driver_storage != nullptr);
    return m_driver_storage.get();
  }

 private:
  hipGraph_t const* m_graph_ptr    = nullptr;
  hipGraphNode_t* m_graph_node_ptr = nullptr;
  mutable std::unique_ptr<base_t, DriverStorageDeleter> m_driver_storage =
      nullptr;
  std::string label;
};

struct HIPGraphNodeAggregate {};

template <typename KernelType,
          typename Tag =
              typename PatternTagFromImplSpecialization<KernelType>::type>
struct get_graph_node_kernel_type
    : std::type_identity<
          GraphNodeKernelImpl<Kokkos::HIP, typename KernelType::Policy,
                              typename KernelType::functor_type, Tag>> {};

template <typename KernelType>
struct get_graph_node_kernel_type<KernelType, Kokkos::ParallelReduceTag>
    : std::type_identity<GraphNodeKernelImpl<
          Kokkos::HIP, typename KernelType::Policy,
          CombinedFunctorReducer<typename KernelType::functor_type,
                                 typename KernelType::reducer_type>,
          Kokkos::ParallelReduceTag>> {};

template <typename KernelType>
auto const& get_graph_node_kernel(KernelType const& kernel) {
  using graph_node_kernel_t =
      typename get_graph_node_kernel_type<KernelType>::type;
  return static_cast<graph_node_kernel_t const&>(kernel);
}

template <typename KernelType>
auto const& get_hip_graph_from_kernel(KernelType const& kernel) {
  hipGraph_t const* graph_ptr =
      get_graph_node_kernel(kernel).get_hip_graph_ptr();
  KOKKOS_EXPECTS(graph_ptr != nullptr);

  return *graph_ptr;
}

template <typename KernelType>
auto& get_hip_graph_node_from_kernel(KernelType const& kernel) {
  auto* graph_node_ptr = get_graph_node_kernel(kernel).get_hip_graph_node_ptr();
  KOKKOS_EXPECTS(graph_node_ptr != nullptr);

  return *graph_node_ptr;
}
}  // namespace Impl
}  // namespace Kokkos

#endif
