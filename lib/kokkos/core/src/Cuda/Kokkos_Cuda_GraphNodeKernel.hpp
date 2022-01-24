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

#ifndef KOKKOS_KOKKOS_CUDA_GRAPHNODEKERNEL_IMPL_HPP
#define KOKKOS_KOKKOS_CUDA_GRAPHNODEKERNEL_IMPL_HPP

#include <Kokkos_Macros.hpp>

#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_CUDA_ENABLE_GRAPHS)

#include <Kokkos_Graph_fwd.hpp>

#include <impl/Kokkos_GraphImpl.hpp>    // GraphAccess needs to be complete
#include <impl/Kokkos_SharedAlloc.hpp>  // SharedAllocationRecord

#include <Kokkos_Parallel.hpp>
#include <Kokkos_Parallel_Reduce.hpp>
#include <Kokkos_PointerOwnership.hpp>

#include <Kokkos_Cuda.hpp>
#include <cuda_runtime_api.h>

namespace Kokkos {
namespace Impl {

template <class PolicyType, class Functor, class PatternTag, class... Args>
class GraphNodeKernelImpl<Kokkos::Cuda, PolicyType, Functor, PatternTag,
                          Args...>
    : public PatternImplSpecializationFromTag<PatternTag, Functor, PolicyType,
                                              Args..., Kokkos::Cuda>::type {
 private:
  using base_t =
      typename PatternImplSpecializationFromTag<PatternTag, Functor, PolicyType,
                                                Args..., Kokkos::Cuda>::type;
  using size_type = Kokkos::Cuda::size_type;
  // These are really functioning as optional references, though I'm not sure
  // that the cudaGraph_t one needs to be since it's a pointer under the
  // covers and we're not modifying it
  Kokkos::ObservingRawPtr<const cudaGraph_t> m_graph_ptr    = nullptr;
  Kokkos::ObservingRawPtr<cudaGraphNode_t> m_graph_node_ptr = nullptr;
  // Note: owned pointer to CudaSpace memory (used for global memory launches),
  // which we're responsible for deallocating, but not responsible for calling
  // its destructor.
  using Record = Kokkos::Impl::SharedAllocationRecord<Kokkos::CudaSpace, void>;
  // Basically, we have to make this mutable for the same reasons that the
  // global kernel buffers in the Cuda instance are mutable...
  mutable Kokkos::OwningRawPtr<base_t> m_driver_storage = nullptr;

 public:
  using Policy       = PolicyType;
  using graph_kernel = GraphNodeKernelImpl;

  // TODO Ensure the execution space of the graph is the same as the one
  //      attached to the policy?
  // TODO @graph kernel name info propagation
  template <class PolicyDeduced, class... ArgsDeduced>
  GraphNodeKernelImpl(std::string, Kokkos::Cuda const&, Functor arg_functor,
                      PolicyDeduced&& arg_policy, ArgsDeduced&&... args)
      // This is super ugly, but it works for now and is the most minimal change
      // to the codebase for now...
      : base_t(std::move(arg_functor), (PolicyDeduced &&) arg_policy,
               (ArgsDeduced &&) args...) {}

  // FIXME @graph Forward through the instance once that works in the backends
  template <class PolicyDeduced>
  GraphNodeKernelImpl(Kokkos::Cuda const& ex, Functor arg_functor,
                      PolicyDeduced&& arg_policy)
      : GraphNodeKernelImpl("", ex, std::move(arg_functor),
                            (PolicyDeduced &&) arg_policy) {}

  ~GraphNodeKernelImpl() {
    if (m_driver_storage) {
      // We should be the only owner, but this is still the easiest way to
      // allocate and deallocate aligned memory for these sorts of things
      Record::decrement(Record::get_record(m_driver_storage));
    }
  }

  void set_cuda_graph_ptr(cudaGraph_t* arg_graph_ptr) {
    m_graph_ptr = arg_graph_ptr;
  }
  void set_cuda_graph_node_ptr(cudaGraphNode_t* arg_node_ptr) {
    m_graph_node_ptr = arg_node_ptr;
  }
  cudaGraphNode_t* get_cuda_graph_node_ptr() const { return m_graph_node_ptr; }
  cudaGraph_t const* get_cuda_graph_ptr() const { return m_graph_ptr; }

  Kokkos::ObservingRawPtr<base_t> allocate_driver_memory_buffer() const {
    KOKKOS_EXPECTS(m_driver_storage == nullptr)

    auto* record = Record::allocate(
        Kokkos::CudaSpace{}, "GraphNodeKernel global memory functor storage",
        sizeof(base_t));

    Record::increment(record);
    m_driver_storage = reinterpret_cast<base_t*>(record->data());
    KOKKOS_ENSURES(m_driver_storage != nullptr)
    return m_driver_storage;
  }
};

struct CudaGraphNodeAggregateKernel {
  using graph_kernel = CudaGraphNodeAggregateKernel;

  // Aggregates don't need a policy, but for the purposes of checking the static
  // assertions about graph kerenls,
  struct Policy {
    using is_graph_kernel = std::true_type;
  };
};

template <class KernelType,
          class Tag =
              typename PatternTagFromImplSpecialization<KernelType>::type>
struct get_graph_node_kernel_type
    : identity<GraphNodeKernelImpl<Kokkos::Cuda, typename KernelType::Policy,
                                   typename KernelType::functor_type, Tag>> {};
template <class KernelType>
struct get_graph_node_kernel_type<KernelType, Kokkos::ParallelReduceTag>
    : identity<GraphNodeKernelImpl<Kokkos::Cuda, typename KernelType::Policy,
                                   typename KernelType::functor_type,
                                   Kokkos::ParallelReduceTag,
                                   typename KernelType::reducer_type>> {};

//==============================================================================
// <editor-fold desc="get_cuda_graph_*() helper functions"> {{{1

template <class KernelType>
auto* allocate_driver_storage_for_kernel(KernelType const& kernel) {
  using graph_node_kernel_t =
      typename get_graph_node_kernel_type<KernelType>::type;
  auto const& kernel_as_graph_kernel =
      static_cast<graph_node_kernel_t const&>(kernel);
  // TODO @graphs we need to somehow indicate the need for a fence in the
  //              destructor of the GraphImpl object (so that we don't have to
  //              just always do it)
  return kernel_as_graph_kernel.allocate_driver_memory_buffer();
}

template <class KernelType>
auto const& get_cuda_graph_from_kernel(KernelType const& kernel) {
  using graph_node_kernel_t =
      typename get_graph_node_kernel_type<KernelType>::type;
  auto const& kernel_as_graph_kernel =
      static_cast<graph_node_kernel_t const&>(kernel);
  cudaGraph_t const* graph_ptr = kernel_as_graph_kernel.get_cuda_graph_ptr();
  KOKKOS_EXPECTS(graph_ptr != nullptr);
  return *graph_ptr;
}

template <class KernelType>
auto& get_cuda_graph_node_from_kernel(KernelType const& kernel) {
  using graph_node_kernel_t =
      typename get_graph_node_kernel_type<KernelType>::type;
  auto const& kernel_as_graph_kernel =
      static_cast<graph_node_kernel_t const&>(kernel);
  auto* graph_node_ptr = kernel_as_graph_kernel.get_cuda_graph_node_ptr();
  KOKKOS_EXPECTS(graph_node_ptr != nullptr);
  return *graph_node_ptr;
}

// </editor-fold> end get_cuda_graph_*() helper functions }}}1
//==============================================================================

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // defined(KOKKOS_ENABLE_CUDA)
#endif  // KOKKOS_KOKKOS_CUDA_GRAPHNODEKERNEL_IMPL_HPP
