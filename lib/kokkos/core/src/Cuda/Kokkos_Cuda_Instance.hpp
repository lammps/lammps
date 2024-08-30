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
      0x001000; /*  4k bytes */
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
  inline static int m_cudaArch = -1;
  static int concurrency();

  inline static cudaDeviceProp m_deviceProp;

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

  inline static std::set<int> cuda_devices = {};
  inline static std::map<int, unsigned long*> constantMemHostStagingPerDevice =
      {};
  inline static std::map<int, cudaEvent_t> constantMemReusablePerDevice = {};
  inline static std::map<int, std::mutex> constantMemMutexPerDevice     = {};

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

  // Using cudaAPI function/objects will be w.r.t. device 0 unless
  // cudaSetDevice(device_id) is called with the correct device_id.
  // The correct device_id is stored in the variable
  // CudaInternal::m_cudaDev set in Cuda::impl_initialize(). It is not
  // sufficient to call cudaSetDevice(m_cudaDev) during cuda initialization
  // only, however, since if a user creates a new thread, that thread will be
  // given the default cuda env with device_id=0, causing errors when
  // device_id!=0 is requested by the user. To ensure against this, almost all
  // cudaAPI calls, as well as using cudaStream_t variables, must be proceeded
  // by cudaSetDevice(device_id).

  // This function sets device in cudaAPI to device requested at runtime (set in
  // m_cudaDev).
  void set_cuda_device() const {
    verify_is_initialized("set_cuda_device");
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaSetDevice(m_cudaDev));
  }

  // Return the class stream, optionally setting the device id.
  template <bool setCudaDevice = true>
  cudaStream_t get_stream() const {
    if constexpr (setCudaDevice) set_cuda_device();
    return m_stream;
  }

  // The following are wrappers for cudaAPI functions (C and C++ routines) which
  // set the correct device id directly before the cudaAPI call (unless
  // explicitly disabled by providing setCudaDevice=false template).
  // setCudaDevice=true should be used for all API calls which take a stream
  // unless it is guarenteed to be from a cuda instance with the correct device
  // set already (e.g., back-to-back cudaAPI calls in a single function). For
  // cudaAPI functions that take a stream, an optional input stream is
  // available. If no stream is given, the stream for the CudaInternal instance
  // is used. All cudaAPI calls should be wrapped in these interface functions
  // to ensure safety when using threads.

  // Helper function for selecting the correct input stream
  cudaStream_t get_input_stream(cudaStream_t s) const {
    return s == nullptr ? get_stream<false>() : s;
  }

  // C API routines
  template <bool setCudaDevice = true>
  cudaError_t cuda_device_get_limit_wrapper(size_t* pValue,
                                            cudaLimit limit) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaDeviceGetLimit(pValue, limit);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_device_set_limit_wrapper(cudaLimit limit,
                                            size_t value) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaDeviceSetLimit(limit, value);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_event_create_wrapper(cudaEvent_t* event) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaEventCreate(event);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_event_destroy_wrapper(cudaEvent_t event) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaEventDestroy(event);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_event_record_wrapper(cudaEvent_t event,
                                        cudaStream_t stream = nullptr) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaEventRecord(event, get_input_stream(stream));
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_event_synchronize_wrapper(cudaEvent_t event) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaEventSynchronize(event);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_free_wrapper(void* devPtr) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaFree(devPtr);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_free_host_wrapper(void* ptr) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaFreeHost(ptr);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_graph_add_dependencies_wrapper(
      cudaGraph_t graph, const cudaGraphNode_t* from, const cudaGraphNode_t* to,
      size_t numDependencies) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaGraphAddDependencies(graph, from, to, numDependencies);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_graph_add_empty_node_wrapper(
      cudaGraphNode_t* pGraphNode, cudaGraph_t graph,
      const cudaGraphNode_t* pDependencies, size_t numDependencies) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaGraphAddEmptyNode(pGraphNode, graph, pDependencies,
                                 numDependencies);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_graph_add_kernel_node_wrapper(
      cudaGraphNode_t* pGraphNode, cudaGraph_t graph,
      const cudaGraphNode_t* pDependencies, size_t numDependencies,
      const cudaKernelNodeParams* pNodeParams) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaGraphAddKernelNode(pGraphNode, graph, pDependencies,
                                  numDependencies, pNodeParams);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_graph_create_wrapper(cudaGraph_t* pGraph,
                                        unsigned int flags) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaGraphCreate(pGraph, flags);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_graph_destroy_wrapper(cudaGraph_t graph) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaGraphDestroy(graph);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_graph_exec_destroy_wrapper(cudaGraphExec_t graphExec) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaGraphExecDestroy(graphExec);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_graph_launch_wrapper(cudaGraphExec_t graphExec,
                                        cudaStream_t stream = nullptr) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaGraphLaunch(graphExec, get_input_stream(stream));
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_host_alloc_wrapper(void** pHost, size_t size,
                                      unsigned int flags) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaHostAlloc(pHost, size, flags);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_malloc_wrapper(void** devPtr, size_t size) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaMalloc(devPtr, size);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_malloc_host_wrapper(void** ptr, size_t size) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaMallocHost(ptr, size);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_malloc_managed_wrapper(
      void** devPtr, size_t size,
      unsigned int flags = cudaMemAttachGlobal) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaMallocManaged(devPtr, size, flags);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_mem_advise_wrapper(const void* devPtr, size_t count,
                                      cudaMemoryAdvise advice,
                                      int device) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaMemAdvise(devPtr, count, advice, device);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_mem_prefetch_async_wrapper(
      const void* devPtr, size_t count, int dstDevice,
      cudaStream_t stream = nullptr) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaMemPrefetchAsync(devPtr, count, dstDevice,
                                get_input_stream(stream));
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_memcpy_wrapper(void* dst, const void* src, size_t count,
                                  cudaMemcpyKind kind) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaMemcpy(dst, src, count, kind);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_memcpy_async_wrapper(void* dst, const void* src,
                                        size_t count, cudaMemcpyKind kind,
                                        cudaStream_t stream = nullptr) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaMemcpyAsync(dst, src, count, kind, get_input_stream(stream));
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_memcpy_to_symbol_async_wrapper(
      const void* symbol, const void* src, size_t count, size_t offset,
      cudaMemcpyKind kind, cudaStream_t stream = nullptr) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaMemcpyToSymbolAsync(symbol, src, count, offset, kind,
                                   get_input_stream(stream));
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_memset_wrapper(void* devPtr, int value, size_t count) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaMemset(devPtr, value, count);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_memset_async_wrapper(void* devPtr, int value, size_t count,
                                        cudaStream_t stream = nullptr) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaMemsetAsync(devPtr, value, count, get_input_stream(stream));
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_pointer_get_attributes_wrapper(
      cudaPointerAttributes* attributes, const void* ptr) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaPointerGetAttributes(attributes, ptr);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_stream_create_wrapper(cudaStream_t* pStream) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaStreamCreate(pStream);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_stream_destroy_wrapper(cudaStream_t stream) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaStreamDestroy(stream);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_stream_synchronize_wrapper(cudaStream_t stream) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaStreamSynchronize(stream);
  }

  // The following are only available for cuda 11.2 and greater
#if (defined(KOKKOS_ENABLE_IMPL_CUDA_MALLOC_ASYNC) && CUDART_VERSION >= 11020)
  template <bool setCudaDevice = true>
  cudaError_t cuda_malloc_async_wrapper(void** devPtr, size_t size,
                                        cudaStream_t hStream = nullptr) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaMallocAsync(devPtr, size, get_input_stream(hStream));
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_free_async_wrapper(void* devPtr,
                                      cudaStream_t hStream = nullptr) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaFreeAsync(devPtr, get_input_stream(hStream));
  }
#endif

  // C++ API routines
  template <typename T, bool setCudaDevice = true>
  cudaError_t cuda_func_get_attributes_wrapper(cudaFuncAttributes* attr,
                                               T* entry) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaFuncGetAttributes(attr, entry);
  }

  template <typename T, bool setCudaDevice = true>
  cudaError_t cuda_func_set_attribute_wrapper(T* entry, cudaFuncAttribute attr,
                                              int value) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaFuncSetAttribute(entry, attr, value);
  }

  template <bool setCudaDevice = true>
  cudaError_t cuda_graph_instantiate_wrapper(cudaGraphExec_t* pGraphExec,
                                             cudaGraph_t graph,
                                             cudaGraphNode_t* pErrorNode,
                                             char* pLogBuffer,
                                             size_t bufferSize) const {
    if constexpr (setCudaDevice) set_cuda_device();
    return cudaGraphInstantiate(pGraphExec, graph, pErrorNode, pLogBuffer,
                                bufferSize);
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

void create_Cuda_instances(std::vector<Cuda>& instances);
}  // Namespace Impl

namespace Experimental {
// Partitioning an Execution Space: expects space and integer arguments for
// relative weight
//   Customization point for backends
//   Default behavior is to return the passed in instance

template <class... Args>
std::vector<Cuda> partition_space(const Cuda&, Args...) {
  static_assert(
      (... && std::is_arithmetic_v<Args>),
      "Kokkos Error: partitioning arguments must be integers or floats");
  std::vector<Cuda> instances(sizeof...(Args));
  Kokkos::Impl::create_Cuda_instances(instances);
  return instances;
}

template <class T>
std::vector<Cuda> partition_space(const Cuda&, std::vector<T> const& weights) {
  static_assert(
      std::is_arithmetic<T>::value,
      "Kokkos Error: partitioning arguments must be integers or floats");

  // We only care about the number of instances to create and ignore weights
  // otherwise.
  std::vector<Cuda> instances(weights.size());
  Kokkos::Impl::create_Cuda_instances(instances);
  return instances;
}
}  // namespace Experimental

}  // Namespace Kokkos
#endif
