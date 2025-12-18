// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

/*--------------------------------------------------------------------------*/

#ifndef KOKKOS_HIP_INSTANCE_HPP
#define KOKKOS_HIP_INSTANCE_HPP

#include <HIP/Kokkos_HIP_Space.hpp>
#include <HIP/Kokkos_HIP_Error.hpp>

#include <hip/hip_runtime_api.h>

#include <atomic>
#include <map>
#include <mutex>
#include <set>

namespace Kokkos {
namespace Impl {

struct HIPTraits {
#if defined(KOKKOS_ARCH_AMD_GFX906) || defined(KOKKOS_ARCH_AMD_GFX908) || \
    defined(KOKKOS_ARCH_AMD_GFX90A) || defined(KOKKOS_ARCH_AMD_GFX940) || \
    defined(KOKKOS_ARCH_AMD_GFX942) || defined(KOKKOS_ARCH_AMD_GFX942_APU)
  static constexpr int WarpSize       = 64;
  static constexpr int WarpIndexMask  = 0x003f; /* hexadecimal for 63 */
  static constexpr int WarpIndexShift = 6;      /* WarpSize == 1 << WarpShift*/
#elif defined(KOKKOS_ARCH_AMD_GFX1030) || defined(KOKKOS_ARCH_AMD_GFX1100) || \
    defined(KOKKOS_ARCH_AMD_GFX1103) || defined(KOKKOS_ARCH_AMD_GFX1201)
  static constexpr int WarpSize       = 32;
  static constexpr int WarpIndexMask  = 0x001f; /* hexadecimal for 31 */
  static constexpr int WarpIndexShift = 5;      /* WarpSize == 1 << WarpShift*/
#endif
  static constexpr int ConservativeThreadsPerBlock =
      256;  // conservative fallback blocksize in case of spills
  static constexpr int MaxThreadsPerBlock =
      1024;  // the maximum we can fit in a block
  static constexpr int ConstantMemoryUsage        = 0x008000; /* 32k bytes */
  static constexpr int KernelArgumentLimit        = 0x001000; /*  4k bytes */
  static constexpr int ConstantMemoryUseThreshold = 0x000200; /* 512 bytes */
};

//----------------------------------------------------------------------------

HIP::size_type *hip_internal_scratch_space(const HIP &instance,
                                           const std::size_t size);
HIP::size_type *hip_internal_scratch_flags(const HIP &instance,
                                           const std::size_t size);

//----------------------------------------------------------------------------

// Helper to protect a shared resource from being used by multiple streams
// simultaneously.
// If used properly, only one stream at a time will be able to use the shared
// resource.
// This helper should be used in a thread-safe way.
//
// Typical usage:
// @code
// auto lock_for_acquisition = shared_resource.acquire();
//
// ... do stuff on the shared resource ...
//
// shared_resource.release(std::move(lock_for_acquisition), stream);
//
// ...
//
// auto lock_for_involvement = shared_resource.lock();
// hipStreamSynchronize(stream);
// shared_resource.check_if_involved_and_unlock(std::move(lock_for_involvement),
//                                              stream);
// hipStreamDestroy(stream);
// @endcode
struct SharedResourceLock {
  bool m_need_sync = false;
  std::mutex m_mutex{};
  hipEvent_t m_event   = nullptr;
  hipStream_t m_stream = nullptr;

  // Acquire the right to interact in a thread-safe way.
  [[nodiscard]] auto lock() { return std::unique_lock<std::mutex>{m_mutex}; }

  // The event is created for the current device. The instance is locked first.
  void initialize() {
    auto lock = this->lock();
    if (!m_event) {
      KOKKOS_IMPL_HIP_SAFE_CALL(
          hipEventCreateWithFlags(&m_event, hipEventDisableTiming));
    }
  }

  // Destroying an event can be done even if it is not bound to the current
  // device.
  void finalize() {
    auto lock = this->lock();
    if (m_event) {
      KOKKOS_IMPL_HIP_SAFE_CALL(hipEventDestroy(m_event));
    }
  }

  SharedResourceLock() = default;

  SharedResourceLock(SharedResourceLock const &other)            = delete;
  SharedResourceLock(SharedResourceLock &&other)                 = delete;
  SharedResourceLock &operator=(SharedResourceLock const &other) = delete;
  SharedResourceLock &operator=(SharedResourceLock &&other)      = delete;

  // Acquire the right to use the shared resource. The instance is locked first.
  [[nodiscard]] auto acquire() {
    auto lock = this->lock();
    if (m_need_sync) KOKKOS_IMPL_HIP_SAFE_CALL(hipEventSynchronize(m_event));
    return lock;
  }

  // Record an event in a stream to signal when it's done with the shared
  // resource.
  void release(std::unique_lock<std::mutex> lock, hipStream_t stream) {
    KOKKOS_ENSURES(lock.owns_lock());
    KOKKOS_ENSURES((lock.mutex() == std::addressof(m_mutex)));
    m_stream = stream;
    KOKKOS_IMPL_HIP_SAFE_CALL(hipEventRecord(m_event, m_stream));
    m_need_sync = true;
    lock.unlock();
  }

  // Check if the stream is the one that was used for the event recording.
  // We assume that this function is called once the stream has been
  // synchronized, and we mark that the shared resource lock does not need to
  // synchronize the next time it is acquired.
  // Doing so allows the current stream to be properly destroyed, while ensuring
  // that the next constant memory launch will work fine.
  // See https://github.com/kokkos/kokkos/issues/8006 for more details.
  void check_if_involved_and_unlock(std::unique_lock<std::mutex> lock,
                                    hipStream_t stream) {
    KOKKOS_ENSURES(lock.owns_lock());
    KOKKOS_ENSURES((lock.mutex() == std::addressof(m_mutex)));
    if (m_stream == stream) m_need_sync = false;
    lock.unlock();
  }
};

class HIPInternal {
 private:
  HIPInternal(const HIPInternal &);
  HIPInternal &operator=(const HIPInternal &);

 public:
  using size_type = ::Kokkos::HIP::size_type;

  int m_hipDev = -1;
  static int m_maxThreadsPerSM;

  static hipDeviceProp_t m_deviceProp;

  static int concurrency();

  // Scratch Spaces for Reductions
  std::size_t m_scratchSpaceCount          = 0;
  std::size_t m_scratchFlagsCount          = 0;
  mutable std::size_t m_scratchFunctorSize = 0;

  size_type *m_scratchSpace               = nullptr;
  size_type *m_scratchFlags               = nullptr;
  mutable size_type *m_scratchFunctor     = nullptr;
  mutable size_type *m_scratchFunctorHost = nullptr;
  static std::mutex scratchFunctorMutex;

  hipStream_t m_stream = nullptr;
  uint32_t m_instance_id =
      Kokkos::Tools::Experimental::Impl::idForInstance<HIP>(
          reinterpret_cast<uintptr_t>(this));

  // Team Scratch Level 1 Space
  int m_n_team_scratch                            = 10;
  mutable int64_t m_team_scratch_current_size[10] = {};
  mutable void *m_team_scratch_ptr[10]            = {};
  mutable std::atomic_int m_team_scratch_pool[10] = {};
  int32_t *m_scratch_locks                        = nullptr;
  size_t m_num_scratch_locks                      = 0;

  bool was_finalized = false;

  static std::set<int> hip_devices;
  static std::map<int, unsigned long *> constantMemHostStaging;
  static std::map<int, SharedResourceLock> constantMemReusable;

  static HIPInternal &singleton();

  int verify_is_initialized(const char *const label) const;

  int is_initialized() const {
    return nullptr != m_scratchSpace && nullptr != m_scratchFlags;
  }

  void initialize(hipStream_t stream);
  void finalize();

  void print_configuration(std::ostream &) const;

  void fence() const;
  void fence(const std::string &) const;

  ~HIPInternal();

  HIPInternal() = default;

  // Using HIP API function/objects will be w.r.t. device 0 unless
  // hipSetDevice(device_id) is called with the correct device_id.
  // The correct device_id is stored in the variable
  // HIPInternal::m_hipDev set in HIP::impl_initialize(). In the case
  // where multiple HIP instances are used, or threads are launched
  // using non-default HIP execution space after initialization, all HIP
  // API calls must follow a call to hipSetDevice(device_id) when an
  // execution space or HIPInternal object is provided to ensure all
  // computation is done on the correct device.

  // FIXME: Not all HIP API calls require us to set device. Potential
  // performance gain by selectively setting device.

  // Set the device in to the device stored by this instance for HIP API calls.
  void set_hip_device() const {
    verify_is_initialized("set_hip_device");
    KOKKOS_IMPL_HIP_SAFE_CALL(hipSetDevice(m_hipDev));
  }

  hipError_t hip_free_wrapper(void *ptr) const {
    set_hip_device();
    return hipFree(ptr);
  }

  hipError_t hip_graph_add_dependencies_wrapper(hipGraph_t graph,
                                                const hipGraphNode_t *from,
                                                const hipGraphNode_t *to,
                                                size_t numDependencies) const {
    set_hip_device();
    return hipGraphAddDependencies(graph, from, to, numDependencies);
  }

  hipError_t hip_graph_add_empty_node_wrapper(
      hipGraphNode_t *pGraphNode, hipGraph_t graph,
      const hipGraphNode_t *pDependencies, size_t numDependencies) const {
    set_hip_device();
    return hipGraphAddEmptyNode(pGraphNode, graph, pDependencies,
                                numDependencies);
  }

  hipError_t hip_graph_add_kernel_node_wrapper(
      hipGraphNode_t *pGraphNode, hipGraph_t graph,
      const hipGraphNode_t *pDependencies, size_t numDependencies,
      const hipKernelNodeParams *pNodeParams) const {
    set_hip_device();
    return hipGraphAddKernelNode(pGraphNode, graph, pDependencies,
                                 numDependencies, pNodeParams);
  }

  hipError_t hip_graph_create_wrapper(hipGraph_t *pGraph,
                                      unsigned int flags) const {
    set_hip_device();
    return hipGraphCreate(pGraph, flags);
  }

  hipError_t hip_graph_destroy_wrapper(hipGraph_t graph) const {
    set_hip_device();
    return hipGraphDestroy(graph);
  }

  hipError_t hip_graph_exec_destroy_wrapper(hipGraphExec_t graphExec) const {
    set_hip_device();
    return hipGraphExecDestroy(graphExec);
  }

  hipError_t hip_graph_instantiate_wrapper(hipGraphExec_t *pGraphExec,
                                           hipGraph_t graph,
                                           hipGraphNode_t *pErrorNode,
                                           char *pLogBuffer,
                                           size_t bufferSize) const {
    set_hip_device();
    return hipGraphInstantiate(pGraphExec, graph, pErrorNode, pLogBuffer,
                               bufferSize);
  }

  hipError_t hip_graph_launch_wrapper(hipGraphExec_t graphExec) const {
    set_hip_device();
    return hipGraphLaunch(graphExec, m_stream);
  }

  hipError_t hip_host_malloc_wrapper(
      void **ptr, size_t size,
      unsigned int flags = hipHostMallocDefault) const {
    set_hip_device();
    return hipHostMalloc(ptr, size, flags);
  }

  hipError_t hip_memcpy_async_wrapper(void *dst, const void *src,
                                      size_t sizeBytes,
                                      hipMemcpyKind kind) const {
    set_hip_device();
    return hipMemcpyAsync(dst, src, sizeBytes, kind, m_stream);
  }

  hipError_t hip_memcpy_to_symbol_async_wrapper(const void *symbol,
                                                const void *src,
                                                size_t sizeBytes, size_t offset,
                                                hipMemcpyKind kind) const {
    set_hip_device();
    return hipMemcpyToSymbolAsync(symbol, src, sizeBytes, offset, kind,
                                  m_stream);
  }

  hipError_t hip_memset_wrapper(void *dst, int value, size_t sizeBytes) const {
    set_hip_device();
    return hipMemset(dst, value, sizeBytes);
  }

  hipError_t hip_memset_async_wrapper(void *dst, int value,
                                      size_t sizeBytes) const {
    set_hip_device();
    return hipMemsetAsync(dst, value, sizeBytes, m_stream);
  }

  hipError_t hip_stream_create_wrapper(hipStream_t *pStream) const {
    set_hip_device();
    return hipStreamCreate(pStream);
  }

  // Resizing of reduction related scratch spaces
  size_type *scratch_space(std::size_t const size);
  size_type *scratch_flags(std::size_t const size);
  size_type *stage_functor_for_execution(void const *driver,
                                         std::size_t const size) const;
  uint32_t impl_get_instance_id() const noexcept;
  int acquire_team_scratch_space();
  // Resizing of team level 1 scratch
  void *resize_team_scratch_space(int scratch_pool_id, std::int64_t bytes,
                                  bool force_shrink = false);
  void release_team_scratch_space(int scratch_pool_id);
};
}  // namespace Impl

namespace Experimental::Impl {
// For each space in partition, create new hipStream_t on the same device as
// base_instance, ignoring weights
template <class T>
std::vector<HIP> impl_partition_space(const HIP &base_instance,
                                      const std::vector<T> &weights) {
  std::vector<HIP> instances;
  instances.reserve(weights.size());
  std::generate_n(
      std::back_inserter(instances), weights.size(), [&base_instance]() {
        hipStream_t stream;
        KOKKOS_IMPL_HIP_SAFE_CALL(base_instance.impl_internal_space_instance()
                                      ->hip_stream_create_wrapper(&stream));
        return HIP(stream, Kokkos::Impl::ManageStream::yes);
      });

  return instances;
}
}  // namespace Experimental::Impl
}  // namespace Kokkos

#endif
