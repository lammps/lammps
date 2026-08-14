#ifndef COLVAR_GPU_SUPPORT_H
#define COLVAR_GPU_SUPPORT_H

#include <vector>
#include <iostream>
#include <unordered_map>

#define COLVARS_STRINGIFY(s) STRINGIFY_HELPER(s)
#define STRINGIFY_HELPER(s) #s

#if defined(COLVARS_CUDA)
#include <cuda_runtime.h>
#ifdef COLVARS_NVTX_PROFILING
#include <nvtx3/nvToolsExt.h>
#endif
#define COLVARS_SYNC_WARP __syncwarp()
#endif // defined(COLVARS_CUDA)

#if defined(COLVARS_HIP)
#include <hip/hip_runtime.h>
#if defined(__HIP_PLATFORM_AMD__)
  /**
   * @note
   * The following macro is from
   * https://github.com/ROCm/clr/blob/29f513cc028f8c1f727dc6e1a731dbee63609f60/hipamd/include/hip/amd_detail/amd_warp_sync_functions.h#L161-L165
   * It is used for ensuring (i) the memory fence and (ii) the synchronization of threads inside a wavefront.
   */
  #define COLVARS_SYNC_WARP do {\
    __builtin_amdgcn_fence(__ATOMIC_RELEASE, "wavefront"); \
    __builtin_amdgcn_wave_barrier(); \
    __builtin_amdgcn_fence(__ATOMIC_ACQUIRE, "wavefront"); \
  } while (0)
#elif defined(__HIP_PLATFORM_NVIDIA__)
  #define COLVARS_SYNC_WARP __syncwarp()
#else
  #error "Unknown HIP platform"
#endif
#endif // defined(COLVARS_HIP)

#if defined(COLVARS_HIP)
#ifndef cudaError_t
#define cudaError_t hipError_t
#endif // cudaError_t

#ifndef cudaFree
#define cudaFree hipFree
#endif // cudaFree

#ifndef cudaFreeHost
#define cudaFreeHost hipFreeHost
#endif // cudaFreeHost

#ifndef cudaFreeAsync
#define cudaFreeAsync hipFreeAsync
#endif // cudaFreeAsync

#ifndef cudaGetErrorString
#define cudaGetErrorString hipGetErrorString
#endif // cudaGetErrorString

#ifndef cudaGraphAddChildGraphNode
#define cudaGraphAddChildGraphNode hipGraphAddChildGraphNode
#endif // cudaGraphAddChildGraphNode

#ifndef cudaGraphAddKernelNode
#define cudaGraphAddKernelNode hipGraphAddKernelNode
#endif // cudaGraphAddKernelNode

#ifndef cudaGraphAddMemcpyNode
#define cudaGraphAddMemcpyNode hipGraphAddMemcpyNode
#endif // cudaGraphAddMemcpyNode

#ifndef cudaGraphAddMemsetNode
#define cudaGraphAddMemsetNode hipGraphAddMemsetNode
#endif // cudaGraphAddMemsetNode

#ifndef cudaGraphCreate
#define cudaGraphCreate hipGraphCreate
#endif // cudaGraphCreate

#ifndef cudaGraphDestroy
#define cudaGraphDestroy hipGraphDestroy
#endif // cudaGraphDestroy

#ifndef cudaGraphExecDestroy
#define cudaGraphExecDestroy hipGraphExecDestroy
#endif // cudaGraphExecDestroy

#ifndef cudaGraphExecMemcpyNodeSetParams
#define cudaGraphExecMemcpyNodeSetParams hipGraphExecMemcpyNodeSetParams
#endif // cudaGraphExecMemcpyNodeSetParams

#ifndef cudaGraphExec_t
#define cudaGraphExec_t hipGraphExec_t
#endif // cudaGraphExec_t

#ifndef cudaGraphInstantiate
#define cudaGraphInstantiate hipGraphInstantiate
#endif // cudaGraphInstantiate

#ifndef cudaGraphInstantiateWithParams
#define cudaGraphInstantiateWithParams hipGraphInstantiateWithParams
#endif // cudaGraphInstantiateWithParams

#ifndef cudaGraphInstantiateParams
#define cudaGraphInstantiateParams hipGraphInstantiateParams
#endif // cudaGraphInstantiateParams

#ifndef cudaGraphInstantiateFlagUpload
#define cudaGraphInstantiateFlagUpload hipGraphInstantiateFlagUpload
#endif // cudaGraphInstantiateFlagUpload

#ifndef cudaGraphInstantiateSuccess
#define cudaGraphInstantiateSuccess hipGraphInstantiateSuccess
#endif // cudaGraphInstantiateSuccess

#ifndef cudaGraphLaunch
#define cudaGraphLaunch hipGraphLaunch
#endif // cudaGraphLaunch

#ifndef cudaGraphNode_t
#define cudaGraphNode_t hipGraphNode_t
#endif // cudaGraphNode_t

#ifndef cudaGraph_t
#define cudaGraph_t hipGraph_t
#endif // cudaGraph_t

#ifndef cudaGraphDebugDotPrint
#define cudaGraphDebugDotPrint hipGraphDebugDotPrint
#endif // cudaGraphDebugDotPrint

#ifndef cudaGraphDebugDotFlags
#define cudaGraphDebugDotFlags hipGraphDebugDotFlags
#endif // cudaGraphDebugDotFlags

#ifndef cudaGraphDebugDotFlagsVerbose
#define cudaGraphDebugDotFlagsVerbose hipGraphDebugDotFlagsVerbose
#endif // cudaGraphDebugDotFlagsVerbose

#ifndef cudaHostAllocMapped
#define cudaHostAllocMapped hipHostAllocMapped
#endif // cudaHostAllocMapped

#ifndef cudaHostAllocDefault
#define cudaHostAllocDefault hipHostAllocDefault
#endif // cudaHostAllocDefault

#ifndef cudaHostAlloc
#define cudaHostAlloc hipHostAlloc
#endif // cudaHostAlloc

#ifndef cudaLaunchKernel
#define cudaLaunchKernel hipLaunchKernel
#endif // cudaLaunchKernel

#ifndef cudaKernelNodeParams
#define cudaKernelNodeParams hipKernelNodeParams
#endif // cudaKernelNodeParams

#ifndef cudaMalloc
#define cudaMalloc hipMalloc
#endif // cudaMalloc

#ifndef cudaMallocAsync
#define cudaMallocAsync hipMallocAsync
#endif // cudaMallocAsync

#ifndef cudaMallocHost
#define cudaMallocHost hipMallocHost
#endif // cudaMallocHost

#ifndef cudaGraphAddMemcpyNode1D
#define cudaGraphAddMemcpyNode1D hipGraphAddMemcpyNode1D
#endif

#ifndef cudaMemcpy
#define cudaMemcpy hipMemcpy
#endif // cudaMemcpy

#ifndef cudaMemcpy3DParms
#define cudaMemcpy3DParms hipMemcpy3DParms
#endif // cudaMemcpy3DParms

#ifndef cudaMemcpyAsync
#define cudaMemcpyAsync hipMemcpyAsync
#endif // cudaMemcpyAsync

#ifndef cudaMemcpyDeviceToDevice
#define cudaMemcpyDeviceToDevice hipMemcpyDeviceToDevice
#endif // cudaMemcpyDeviceToDevice

#ifndef cudaMemcpyDeviceToHost
#define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost
#endif // cudaMemcpyDeviceToHost

#ifndef cudaMemcpyHostToDevice
#define cudaMemcpyHostToDevice hipMemcpyHostToDevice
#endif // cudaMemcpyHostToDevice

#ifndef cudaMemcpyKind
#define cudaMemcpyKind hipMemcpyKind
#endif // cudaMemcpyKind

#ifndef cudaMemset
#define cudaMemset hipMemset
#endif // cudaMemset

#ifndef cudaMemsetParams
#define cudaMemsetParams hipMemsetParams
#endif // cudaMemsetParams

#ifndef cudaMemsetAsync
#define cudaMemsetAsync hipMemsetAsync
#endif // cudaMemsetAsync

#ifndef cudaStreamCreate
#define cudaStreamCreate hipStreamCreate
#endif // cudaStreamCreate

#ifndef cudaStreamDestroy
#define cudaStreamDestroy hipStreamDestroy
#endif // cudaStreamDestroy

#ifndef cudaStreamSynchronize
#define cudaStreamSynchronize hipStreamSynchronize
#endif // cudaStreamSynchronize

#ifndef cudaStream_t
#define cudaStream_t hipStream_t
#endif // cudaStream_t

#ifndef cudaSuccess
#define cudaSuccess hipSuccess
#endif // cudaSuccess

#ifndef make_cudaExtent
#define make_cudaExtent make_hipExtent
#endif // make_cudaExtent

#ifndef make_cudaPitchedPtr
#define make_cudaPitchedPtr make_hipPitchedPtr
#endif // make_cudaPitchedPtr

#ifndef make_cudaPos
#define make_cudaPos make_hipPos
#endif // make_cudaPos

#endif // defined(COLVARS_HIP)

namespace colvars_gpu {

#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)
/// \brief Default block size for CUDA kernels
constexpr unsigned int default_block_size = 128;
/// \brief Default maximum number of blocks for reduction kernels
static unsigned int default_reduce_max_num_blocks = 64;
// static unsigned int default_atom_wise_num_blocks = 64;
#endif

#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)
#define COLVARS_HOST_DEVICE __device__ __host__
#define COLVARS_DEVICE __device__
#else
#define COLVARS_HOST_DEVICE
#define COLVARS_DEVICE
#endif

// HIP does not have cuda::std::array since libhipcxx is not a part of the ROCm distribution,
// so reinvent the wheel...
#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)
template <typename T, unsigned long N>
class array1d {
public:
  T m_data[N];
  using value_type = T;
  using size_type = decltype(N);
  using reference = value_type&;
  using const_reference = const value_type&;
  using pointer = T*;
  using const_pointer = const T*;
  COLVARS_HOST_DEVICE constexpr size_type size() const {return N;}
  COLVARS_HOST_DEVICE reference operator[](size_type pos) {return m_data[pos];}
  COLVARS_HOST_DEVICE const_reference operator[](size_type pos) const {return m_data[pos];}
};
#endif

// TODO: What about SYCL?
#if ( defined(COLVARS_CUDA) || defined(COLVARS_HIP) )
/**
 * @brief Allocator for pinned host memory using cudaHostAlloc
 *
 * This allocator can be used with STL containers to allocate pinned
 * host memory that is page-locked and directly accessible by the GPU.
 *
 * @tparam T The type of elements to allocate
 */
template <typename T>
class CudaHostAllocator {
public:
  using value_type = T;

  CudaHostAllocator() = default;

  template<typename U>
  constexpr CudaHostAllocator(const CudaHostAllocator<U>&) noexcept {}

  friend bool operator==(const CudaHostAllocator&, const CudaHostAllocator&) { return true; }
  friend bool operator!=(const CudaHostAllocator&, const CudaHostAllocator&) { return false; }

  T* allocate(size_t n) {
    T* ptr;
    if (cudaHostAlloc(&ptr, n * sizeof(T), cudaHostAllocMapped) != cudaSuccess) {
      throw std::bad_alloc();
    }
    return ptr;
  }
  void deallocate(T* ptr, size_t n) noexcept {
    (void)cudaFreeHost(ptr);
  }
  template<typename U, typename... Args>
  void construct(U* p, Args&&... args) {
      new(p) U(std::forward<Args>(args)...);
  }

  template<typename U>
  void destroy(U* p) noexcept {
      p->~U();
  }
};
#endif


#if defined(COLVARS_CUDA) || defined (COLVARS_HIP)
/**
 * @brief Check for CUDA errors and report them
 *
 * @param code The CUDA error code to check
 * @param file The source file where the error occurred
 * @param line The line number in the source file
 * @return COLVARS_OK if no error, otherwise the COLVARS_ERROR
 */
int gpuAssert(cudaError_t code, const char *file, int line);
#endif

} // namespace colvars_gpu

#if defined(COLVARS_CUDA) || defined (COLVARS_HIP)
/// \define checkGPUError(ans) { gpuAssert((ans), __FILE__, __LINE__); }
/// \brief Macro to check for CUDA errors
#define checkGPUError(ans) colvars_gpu::gpuAssert((ans), __FILE__, __LINE__);
#endif

namespace colvars_gpu {
#if defined(COLVARS_CUDA) || defined (COLVARS_HIP)

/**
 * @brief Add a CUDA graph node to clear an array to zero (used by add_clear_array_node)
 */
int add_clear_array_node_impl(
  void* dst, const size_t num_elements, const size_t sizeofT,
  cudaGraphNode_t& node_out, cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies);

/**
 * @brief Add a CUDA graph node to copy an array (used by add_copy_node)
 */
int add_copy_node_impl(
  const void* src, void* dst, const size_t num_elements, const size_t sizeofT,
  cudaMemcpyKind kind, cudaGraphNode_t& node_out, cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies);

/**
 * @brief Add a CUDA graph node to clear an array to zero
 *
 * @tparam T The type of elements in the array
 * @param dst Pointer to the device array to clear
 * @param num_elements Number of elements in the array
 * @param node_out Output parameter to receive the created CUDA graph node
 * @param graph The CUDA graph to which the node will be added
 * @param dependencies A vector of CUDA graph nodes that this node depends on
 * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
 */
template <typename T>
int add_clear_array_node(
  T* dst, const size_t num_elements,
  cudaGraphNode_t& node_out, cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  return add_clear_array_node_impl(
    dst, num_elements, sizeof(T), node_out, graph, dependencies);
}

/**
 * @brief Add a CUDA graph node to copy an array
 *
 * @tparam T The type of elements in the array
 * @param src Pointer to the source array
 * @param dst Pointer to the destination array
 * @param num_elements Number of elements to copy
 * @param kind The type of copy (cudaMemcpyKind)
 * @param node_out Output parameter to receive the created CUDA graph node
 * @param graph The CUDA graph to which the node will be added
 * @param dependencies A vector of CUDA graph nodes that this node depends on
 * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
 */
template <typename T>
int add_copy_node(
  const T* src, T* dst, size_t num_elements,
  cudaMemcpyKind kind, cudaGraphNode_t& node_out, cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  return add_copy_node_impl(src, dst, num_elements, sizeof(T),
                            kind, node_out, graph, dependencies);
}

/**
 * @brief Prepare a list of CUDA graph node dependencies
 *
 * This function looks up the specified node names in the provided map
 * and collects the corresponding CUDA graph nodes into the dependencies vector.
 * If any node name is not found, an error is reported.
 *
 * @param[in] node_names A vector of pairs containing node names and a boolean indicating if they are optional
 * @param[out] dependencies A vector to store the collected CUDA graph nodes
 * @param[in] map A map from node names to CUDA graph nodes
 * @param[in] caller_operation_name Optional name of the calling operation for error reporting
 * @return COLVARS_OK if all required nodes are found, otherwise COLVARS_ERROR
 */
int prepare_dependencies(
  const std::vector<std::pair<std::string, bool>>& node_names,
  std::vector<cudaGraphNode_t>& dependencies,
  const std::unordered_map<std::string, cudaGraphNode_t>& map,
  const std::string& caller_operation_name = "");

// NVTX Profiling
#if defined (COLVARS_NVTX_PROFILING)
/**
 * @brief Class for managing NVTX profiling ranges
 *
 * This class encapsulates the functionality to create and manage NVTX
 * profiling ranges. It allows setting a name and color for the range,
 * and provides methods to start and stop the profiling range.
 */
class colvar_nvtx_prof {
public:
  colvar_nvtx_prof();
  void set_name_color(const std::string& name_in, const uint32_t color_in);
  inline void start() {
    nvtxRangePushEx(&nvtx_event_attr);
  }
  inline void stop() {
    nvtxRangePop();
  }
private:
  std::string nvtx_event_name;
  nvtxEventAttributes_t nvtx_event_attr;
};
#endif // defined (COLVARS_NVTX_PROFILING)
#endif // defined(COLVARS_CUDA) || defined (COLVARS_HIP)
}

#endif // COLVAR_GPU_SUPPORT_H
