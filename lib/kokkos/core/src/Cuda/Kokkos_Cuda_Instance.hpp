// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_CUDA_INSTANCE_HPP_
#define KOKKOS_CUDA_INSTANCE_HPP_

#include <vector>
#include <impl/Kokkos_Tools.hpp>
#include <atomic>
#include <Cuda/Kokkos_Cuda_Error.hpp>
#include <cuda_runtime_api.h>
#include "Kokkos_CudaSpace.hpp"

#include <set>
#include <map>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// These functions fulfill the purpose of allowing to work around
// a suspected system software issue, or to check for race conditions.
// They are not currently a fully officially supported capability.
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
extern "C" void kokkos_impl_cuda_set_serial_execution(bool);
extern "C" bool kokkos_impl_cuda_use_serial_execution();
#endif

#if defined(KOKKOS_COMPILER_NVCC) && !defined(KOKKOS_ARCH_MAXWELL) && \
    !defined(KOKKOS_ARCH_PASCAL)
#define KOKKOS_IMPL_CUDA_USE_GRID_CONSTANT
#endif

namespace Kokkos {
namespace Impl {

struct CudaTraits {
  static constexpr CudaSpace::size_type WarpSize = 32 /* 0x0020 */;
  static constexpr CudaSpace::size_type WarpIndexMask =
      0x001f; /* Mask for warpindex */
  static constexpr CudaSpace::size_type WarpIndexShift =
      5; /* WarpSize == 1 << WarpShift */

  static constexpr CudaSpace::size_type ConstantMemoryUsage =
      0x008000; /* 32k bytes */
  static constexpr CudaSpace::size_type ConstantMemoryCache =
      0x002000; /*  8k bytes */
  static constexpr CudaSpace::size_type KernelArgumentLimit =
#ifdef KOKKOS_IMPL_CUDA_USE_GRID_CONSTANT
      0x008000; /* 32k bytes */
#else
      0x001000; /*  4k bytes */
#endif
  static constexpr CudaSpace::size_type MaxHierarchicalParallelism =
      1024; /* team_size * vector_length */
  using ConstantGlobalBufferType =
      unsigned long[ConstantMemoryUsage / sizeof(unsigned long)];

  static constexpr int ConstantMemoryUseThreshold = 0x000200 /* 512 bytes */;
};

//----------------------------------------------------------------------------

CudaSpace::size_type* cuda_internal_scratch_flags(const Cuda&,
                                                  const std::size_t size);
CudaSpace::size_type* cuda_internal_scratch_space(const Cuda&,
                                                  const std::size_t size);
CudaSpace::size_type* cuda_internal_scratch_unified(const Cuda&,
                                                    const std::size_t size);

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
namespace Kokkos {
namespace Impl {

class CudaInternal {
 private:
  CudaInternal(const CudaInternal&);
  CudaInternal& operator=(const CudaInternal&);
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
  static bool kokkos_impl_cuda_use_serial_execution_v;
#endif

 public:
  using size_type = Cuda::size_type;

  int m_cudaDev = -1;

  // Device Properties
  static int m_cudaArch;
  static int concurrency();

  KOKKOS_IMPL_EXPORT static cudaDeviceProp m_deviceProp;

  // Scratch Spaces for Reductions
  mutable std::size_t m_scratchSpaceCount;
  mutable std::size_t m_scratchFlagsCount;
  mutable std::size_t m_scratchUnifiedCount;
  mutable std::size_t m_scratchFunctorSize;

  mutable size_type* m_scratchSpace;
  mutable size_type* m_scratchFlags;
  mutable size_type* m_scratchUnified;
  mutable size_type* m_scratchFunctor;
  cudaStream_t m_stream;
  uint32_t m_instance_id;

  // Team Scratch Level 1 Space
  int m_n_team_scratch = 10;
  mutable int64_t m_team_scratch_current_size[10];
  mutable void* m_team_scratch_ptr[10];
  mutable std::atomic_int m_team_scratch_pool[10];
  int32_t* m_scratch_locks;
  size_t m_num_scratch_locks;

  bool was_initialized = false;
  bool was_finalized   = false;

  static std::set<int> cuda_devices;
  KOKKOS_IMPL_EXPORT static std::map<int, unsigned long*>
      constantMemHostStagingPerDevice;
  KOKKOS_IMPL_EXPORT static std::map<int, cudaEvent_t>
      constantMemReusablePerDevice;
  KOKKOS_IMPL_EXPORT static std::map<int, std::mutex> constantMemMutexPerDevice;

  static CudaInternal& singleton();

  int verify_is_initialized(const char* const label) const;

  int is_initialized() const {
    return nullptr != m_scratchSpace && nullptr != m_scratchFlags;
  }

  void initialize(cudaStream_t stream);
  void finalize();

  void print_configuration(std::ostream&) const;

#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
  static bool cuda_use_serial_execution();
  static void cuda_set_serial_execution(bool);
#endif

  void fence(const std::string&) const;
  void fence() const;

  ~CudaInternal();

  CudaInternal()
      : m_scratchSpaceCount(0),
        m_scratchFlagsCount(0),
        m_scratchUnifiedCount(0),
        m_scratchFunctorSize(0),
        m_scratchSpace(nullptr),
        m_scratchFlags(nullptr),
        m_scratchUnified(nullptr),
        m_scratchFunctor(nullptr),
        m_stream(nullptr),
        m_instance_id(
            Kokkos::Tools::Experimental::Impl::idForInstance<Kokkos::Cuda>(
                reinterpret_cast<uintptr_t>(this))) {
    for (int i = 0; i < m_n_team_scratch; ++i) {
      m_team_scratch_current_size[i] = 0;
      m_team_scratch_ptr[i]          = nullptr;
      m_team_scratch_pool[i]         = 0;
    }
  }

  // Using CUDA API function/objects will be w.r.t. device 0 unless
  // cudaSetDevice(device_id) is called with the correct device_id.
  // The correct device_id is stored in the variable
  // CudaInternal::m_cudaDev set in Cuda::impl_initialize(). In the case
  // where multiple CUDA instances are used, or threads are launched
  // using non-default CUDA execution space after initialization, all CUDA
  // API calls must follow a call to cudaSetDevice(device_id) when an
  // execution space or CudaInternal object is provided to ensure all
  // computation is done on the correct device.

  // FIXME: Not all CUDA API calls require us to set device. Potential
  // performance gain by selectively setting device.

  // Set the device to the one stored by this instance for CUDA API calls.
  void set_cuda_device() const {
    verify_is_initialized("set_cuda_device");
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaSetDevice(m_cudaDev));
  }

  // CUDA API wrappers

  // C API routines
  cudaError_t cuda_device_get_limit_wrapper(size_t* pValue,
                                            cudaLimit limit) const {
    set_cuda_device();
    return cudaDeviceGetLimit(pValue, limit);
  }

  cudaError_t cuda_device_set_limit_wrapper(cudaLimit limit,
                                            size_t value) const {
    set_cuda_device();
    return cudaDeviceSetLimit(limit, value);
  }

  cudaError_t cuda_event_create_with_flags_wrapper(
      cudaEvent_t* event, const unsigned int flags) const {
    set_cuda_device();
    return cudaEventCreateWithFlags(event, flags);
  }

  cudaError_t cuda_event_record_wrapper(cudaEvent_t event) const {
    set_cuda_device();
    return cudaEventRecord(event, m_stream);
  }

  cudaError_t cuda_free_wrapper(void* devPtr) const {
    set_cuda_device();
    return cudaFree(devPtr);
  }

  cudaError_t cuda_graph_add_dependencies_wrapper(
      cudaGraph_t graph, const cudaGraphNode_t* from, const cudaGraphNode_t* to,
      size_t numDependencies) const {
    set_cuda_device();
#if CUDART_VERSION >= 13000
    return cudaGraphAddDependencies(graph, from, to, NULL, numDependencies);
#else
    return cudaGraphAddDependencies(graph, from, to, numDependencies);
#endif
  }

  cudaError_t cuda_graph_add_empty_node_wrapper(
      cudaGraphNode_t* pGraphNode, cudaGraph_t graph,
      const cudaGraphNode_t* pDependencies, size_t numDependencies) const {
    set_cuda_device();
    return cudaGraphAddEmptyNode(pGraphNode, graph, pDependencies,
                                 numDependencies);
  }

  cudaError_t cuda_graph_add_kernel_node_wrapper(
      cudaGraphNode_t* pGraphNode, cudaGraph_t graph,
      const cudaGraphNode_t* pDependencies, size_t numDependencies,
      const cudaKernelNodeParams* pNodeParams) const {
    set_cuda_device();
    return cudaGraphAddKernelNode(pGraphNode, graph, pDependencies,
                                  numDependencies, pNodeParams);
  }

  cudaError_t cuda_graph_create_wrapper(cudaGraph_t* pGraph,
                                        unsigned int flags) const {
    set_cuda_device();
    return cudaGraphCreate(pGraph, flags);
  }

  cudaError_t cuda_graph_destroy_wrapper(cudaGraph_t graph) const {
    set_cuda_device();
    return cudaGraphDestroy(graph);
  }

  cudaError_t cuda_graph_exec_destroy_wrapper(cudaGraphExec_t graphExec) const {
    set_cuda_device();
    return cudaGraphExecDestroy(graphExec);
  }

  cudaError_t cuda_graph_launch_wrapper(cudaGraphExec_t graphExec) const {
    set_cuda_device();
    return cudaGraphLaunch(graphExec, m_stream);
  }

  cudaError_t cuda_malloc_wrapper(void** devPtr, size_t size) const {
    set_cuda_device();
    return cudaMalloc(devPtr, size);
  }

  cudaError_t cuda_malloc_host_wrapper(void** ptr, size_t size) const {
    set_cuda_device();
    return cudaMallocHost(ptr, size);
  }

  cudaError_t cuda_mem_prefetch_async_wrapper(const void* devPtr, size_t count,
                                              int dstDevice) const {
    set_cuda_device();
#if CUDART_VERSION >= 13000
    cudaMemLocation loc = {cudaMemLocationTypeDevice, dstDevice};
    return cudaMemPrefetchAsync(devPtr, count, loc, 0, m_stream);
#else
    return cudaMemPrefetchAsync(devPtr, count, dstDevice, m_stream);
#endif
  }

  cudaError_t cuda_memcpy_wrapper(void* dst, const void* src, size_t count,
                                  cudaMemcpyKind kind) const {
    set_cuda_device();
    return cudaMemcpy(dst, src, count, kind);
  }

  cudaError_t cuda_memcpy_async_wrapper(void* dst, const void* src,
                                        size_t count,
                                        cudaMemcpyKind kind) const {
    set_cuda_device();
    return cudaMemcpyAsync(dst, src, count, kind, m_stream);
  }

  cudaError_t cuda_memcpy_to_symbol_async_wrapper(const void* symbol,
                                                  const void* src, size_t count,
                                                  size_t offset,
                                                  cudaMemcpyKind kind) const {
    set_cuda_device();
    return cudaMemcpyToSymbolAsync(symbol, src, count, offset, kind, m_stream);
  }

  cudaError_t cuda_memset_wrapper(void* devPtr, int value, size_t count) const {
    set_cuda_device();
    return cudaMemset(devPtr, value, count);
  }

  cudaError_t cuda_memset_async_wrapper(void* devPtr, int value,
                                        size_t count) const {
    set_cuda_device();
    return cudaMemsetAsync(devPtr, value, count, m_stream);
  }

  cudaError_t cuda_pointer_get_attributes_wrapper(
      cudaPointerAttributes* attributes, const void* ptr) const {
    set_cuda_device();
    return cudaPointerGetAttributes(attributes, ptr);
  }

  cudaError_t cuda_stream_create_wrapper(cudaStream_t* pStream) const {
    set_cuda_device();
    return cudaStreamCreate(pStream);
  }

  // C++ API routines
  template <typename T>
  cudaError_t cuda_func_get_attributes_wrapper(cudaFuncAttributes* attr,
                                               T* entry) const {
    set_cuda_device();
    return cudaFuncGetAttributes(attr, entry);
  }

  template <typename T>
  cudaError_t cuda_func_set_attribute_wrapper(T* entry, cudaFuncAttribute attr,
                                              int value) const {
    set_cuda_device();
    return cudaFuncSetAttribute(entry, attr, value);
  }

  cudaError_t cuda_graph_instantiate_wrapper(cudaGraphExec_t* pGraphExec,
                                             cudaGraph_t graph) const {
    set_cuda_device();
#if CUDA_VERSION < 12000
    constexpr size_t error_log_size = 256;
    cudaGraphNode_t error_node      = nullptr;
    char error_log[error_log_size];
    return cudaGraphInstantiate(pGraphExec, graph, &error_node, error_log,
                                error_log_size);
#else
    return cudaGraphInstantiate(pGraphExec, graph);
#endif
  }

  // Resizing of reduction related scratch spaces
  size_type* scratch_space(const std::size_t size) const;
  size_type* scratch_flags(const std::size_t size) const;
  size_type* scratch_unified(const std::size_t size) const;
  size_type* scratch_functor(const std::size_t size) const;
  uint32_t impl_get_instance_id() const;
  int acquire_team_scratch_space();
  // Resizing of team level 1 scratch
  void* resize_team_scratch_space(int scratch_pool_id, std::int64_t bytes,
                                  bool force_shrink = false);
  void release_team_scratch_space(int scratch_pool_id);
};
}  // Namespace Impl

namespace Experimental::Impl {
// For each space in partition, create new cudaStream_t on the same device as
// base_instance, ignoring weights
template <class T>
std::vector<Cuda> impl_partition_space(const Cuda& base_instance,
                                       const std::vector<T>& weights) {
  std::vector<Cuda> instances;
  instances.reserve(weights.size());
  std::generate_n(
      std::back_inserter(instances), weights.size(), [&base_instance]() {
        cudaStream_t stream;
        KOKKOS_IMPL_CUDA_SAFE_CALL(base_instance.impl_internal_space_instance()
                                       ->cuda_stream_create_wrapper(&stream));
        return Cuda(stream, Kokkos::Impl::ManageStream::yes);
      });

  return instances;
}
}  // namespace Experimental::Impl

}  // Namespace Kokkos
#endif
