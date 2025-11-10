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

/*--------------------------------------------------------------------------*/
/* Kokkos interfaces */

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <Kokkos_Core.hpp>

#include <HIP/Kokkos_HIP_Instance.hpp>
#include <HIP/Kokkos_HIP.hpp>
#include <HIP/Kokkos_HIP_Space.hpp>
#include <HIP/Kokkos_HIP_IsXnack.hpp>
#include <impl/Kokkos_CheckedIntegerOps.hpp>
#include <impl/Kokkos_DeviceManagement.hpp>
#include <impl/Kokkos_Error.hpp>

/*--------------------------------------------------------------------------*/
/* Standard 'C' libraries */
#include <stdlib.h>

/* Standard 'C++' libraries */
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE
__device__ __constant__ unsigned long kokkos_impl_hip_constant_memory_buffer
    [Kokkos::Impl::HIPTraits::ConstantMemoryUsage / sizeof(unsigned long)];
#endif

namespace Kokkos {
namespace Impl {
Kokkos::View<uint32_t *, HIPSpace> hip_global_unique_token_locks(
    bool deallocate) {
  static Kokkos::View<uint32_t *, HIPSpace> locks =
      Kokkos::View<uint32_t *, HIPSpace>();
  if (!deallocate && locks.extent(0) == 0)
    locks = Kokkos::View<uint32_t *, HIPSpace>(
        "Kokkos::UniqueToken<HIP>::m_locks", HIPInternal::concurrency());
  if (deallocate) locks = Kokkos::View<uint32_t *, HIPSpace>();
  return locks;
}
}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos {
namespace Impl {

namespace {

using ScratchGrain = Kokkos::HIP::size_type[Impl::HIPTraits::WarpSize];
constexpr auto sizeScratchGrain = sizeof(ScratchGrain);

std::size_t scratch_count(const std::size_t size) {
  return (size + sizeScratchGrain - 1) / sizeScratchGrain;
}

}  // namespace

//----------------------------------------------------------------------------

int HIPInternal::concurrency() {
  static int const concurrency =
      m_maxThreadsPerSM * m_deviceProp.multiProcessorCount;

  return concurrency;
}

void HIPInternal::print_configuration(std::ostream &s) const {
  s << "macro  KOKKOS_ENABLE_HIP : defined" << '\n';
#if defined(HIP_VERSION)
  s << "macro  HIP_VERSION = " << HIP_VERSION << " = version "
    << HIP_VERSION_MAJOR << '.' << HIP_VERSION_MINOR << '.' << HIP_VERSION_PATCH
    << '\n';
#endif

  s << "macro KOKKOS_ENABLE_ROCTHRUST : "
#if defined(KOKKOS_ENABLE_ROCTHRUST)
    << "defined\n";
#else
    << "undefined\n";
#endif

  s << "macro KOKKOS_ENABLE_IMPL_HIP_MALLOC_ASYNC: ";
#ifdef KOKKOS_ENABLE_IMPL_HIP_MALLOC_ASYNC
  s << "yes\n";
#else
  s << "no\n";
#endif

  for (int i : get_visible_devices()) {
    hipDeviceProp_t hipProp;
    KOKKOS_IMPL_HIP_SAFE_CALL(hipGetDeviceProperties(&hipProp, i));
    std::string gpu_type = hipProp.integrated == 1 ? "APU" : "dGPU";

    s << "Kokkos::HIP[ " << i << " ] "
      << "gcnArch " << hipProp.gcnArchName;
    if (m_hipDev == i) s << " : Selected";
    s << '\n'
      << "  Total Global Memory: "
      << ::Kokkos::Impl::human_memory_size(hipProp.totalGlobalMem) << '\n'
      << "  Shared Memory per Block: "
      << ::Kokkos::Impl::human_memory_size(hipProp.sharedMemPerBlock) << '\n'
      << "  APU or dGPU: " << gpu_type << '\n'
      << "  Is Large Bar: " << hipProp.isLargeBar << '\n'
      << "  Supports Managed Memory: " << hipProp.managedMemory << '\n'
      << "  Architecture capable of accessing system allocated memory: "
      << gpu_arch_can_access_system_allocations() << '\n'
      << "  System allows accessing system allocated memory on GPU: "
      << (xnack_boot_config_has_hmm_mirror() && xnack_environment_enabled() &&
          gpu_arch_can_access_system_allocations())
      << '\n'
      << "  Wavefront Size: " << hipProp.warpSize << '\n';
  }
}

//----------------------------------------------------------------------------

HIPInternal::~HIPInternal() {
  if (m_scratchSpace || m_scratchFlags) {
    std::cerr << "Kokkos::HIP ERROR: Failed to call "
                 "Kokkos::HIP::finalize()"
              << std::endl;
    std::cerr.flush();
  }

  m_scratchSpaceCount = 0;
  m_scratchFlagsCount = 0;
  m_scratchSpace      = nullptr;
  m_scratchFlags      = nullptr;
  m_stream            = nullptr;
}

int HIPInternal::verify_is_initialized(const char *const label) const {
  if (m_hipDev < 0) {
    Kokkos::abort((std::string("Kokkos::HIP::") + label +
                   " : ERROR device not initialized\n")
                      .c_str());
  }
  return 0 <= m_hipDev;
}

uint32_t HIPInternal::impl_get_instance_id() const noexcept {
  return m_instance_id;
}
HIPInternal &HIPInternal::singleton() {
  static HIPInternal *self = nullptr;
  if (!self) {
    self = new HIPInternal();
  }
  return *self;
}

void HIPInternal::fence() const {
  fence("Kokkos::HIPInternal::fence: Unnamed Internal Fence");
}
void HIPInternal::fence(const std::string &name) const {
  Kokkos::Tools::Experimental::Impl::profile_fence_event<Kokkos::HIP>(
      name,
      Kokkos::Tools::Experimental::Impl::DirectFenceIDHandle{
          impl_get_instance_id()},
      [&]() { KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamSynchronize(m_stream)); });
}

void HIPInternal::initialize(hipStream_t stream) {
  KOKKOS_EXPECTS(!is_initialized());

  if (was_finalized)
    Kokkos::abort("Calling HIP::initialize after HIP::finalize is illegal\n");

    // Get the device ID. If this is ROCm 5.6 or later, we can query this from
    // the provided stream and potentially use multiple GPU devices. For
    // ROCm 5.5 or earlier, we must use the singleton device id and there are no
    // checks possible for the device id matching the device the stream was
    // created on.
#if (HIP_VERSION_MAJOR > 5 || \
     (HIP_VERSION_MAJOR == 5 && HIP_VERSION_MINOR >= 6))
  KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamGetDevice(stream, &m_hipDev));
#else
  m_hipDev = singleton().m_hipDev;
#endif
  KOKKOS_IMPL_HIP_SAFE_CALL(hipSetDevice(m_hipDev));
  hip_devices.insert(m_hipDev);

  m_stream = stream;

  // Allocate a staging buffer for constant mem in pinned host memory.
  if (!constantMemHostStaging[m_hipDev]) {
    void *constant_mem_void_ptr = nullptr;
    KOKKOS_IMPL_HIP_SAFE_CALL(hip_host_malloc_wrapper(
        &constant_mem_void_ptr, Impl::HIPTraits::ConstantMemoryUsage));
    constantMemHostStaging[m_hipDev] =
        static_cast<unsigned long *>(constant_mem_void_ptr);
  }

  // Initialize the shared resource locking to avoid overwriting the driver of
  // the previous kernel launch.
  constantMemReusable[m_hipDev].initialize();

  //----------------------------------
  // Multiblock reduction uses scratch flags for counters
  // and scratch space for partial reduction values.
  // Allocate some initial space.  This will grow as needed.
  {
    // Maximum number of warps,
    // at most one warp per thread in a warp for reduction.
    unsigned int maxWarpCount =
        m_deviceProp.maxThreadsPerBlock / Impl::HIPTraits::WarpSize;
    if (Impl::HIPTraits::WarpSize < maxWarpCount) {
      maxWarpCount = Impl::HIPTraits::WarpSize;
    }

    const unsigned reduce_block_count =
        maxWarpCount * Impl::HIPTraits::WarpSize;

    (void)scratch_flags(static_cast<size_t>(reduce_block_count * 2) *
                        sizeof(size_type));
    (void)scratch_space(static_cast<size_t>(reduce_block_count * 16) *
                        sizeof(size_type));
  }

  m_num_scratch_locks = concurrency();
  KOKKOS_IMPL_HIP_SAFE_CALL(
      hipMalloc(&m_scratch_locks, sizeof(int32_t) * m_num_scratch_locks));
  KOKKOS_IMPL_HIP_SAFE_CALL(
      hipMemset(m_scratch_locks, 0, sizeof(int32_t) * m_num_scratch_locks));
}

//----------------------------------------------------------------------------

Kokkos::HIP::size_type *HIPInternal::scratch_space(const std::size_t size) {
  if (verify_is_initialized("scratch_space") &&
      m_scratchSpaceCount < scratch_count(size)) {
    auto mem_space = Kokkos::HIPSpace::impl_create(m_hipDev, m_stream);

    if (m_scratchSpace) {
      mem_space.deallocate(m_scratchSpace,
                           m_scratchSpaceCount * sizeScratchGrain);
    }

    m_scratchSpaceCount = scratch_count(size);

    std::size_t alloc_size =
        multiply_overflow_abort(m_scratchSpaceCount, sizeScratchGrain);
    m_scratchSpace = static_cast<size_type *>(
        mem_space.allocate("Kokkos::InternalScratchSpace", alloc_size));
  }

  return m_scratchSpace;
}

Kokkos::HIP::size_type *HIPInternal::scratch_flags(const std::size_t size) {
  if (verify_is_initialized("scratch_flags") &&
      m_scratchFlagsCount < scratch_count(size)) {
    auto mem_space = Kokkos::HIPSpace::impl_create(m_hipDev, m_stream);

    if (m_scratchFlags) {
      mem_space.deallocate(m_scratchFlags,
                           m_scratchFlagsCount * sizeScratchGrain);
    }

    m_scratchFlagsCount = scratch_count(size);

    std::size_t alloc_size =
        multiply_overflow_abort(m_scratchFlagsCount, sizeScratchGrain);
    m_scratchFlags = static_cast<size_type *>(
        mem_space.allocate("Kokkos::InternalScratchFlags", alloc_size));

    // We only zero-initialize the allocation when we actually allocate.
    // It's the responsibility of the features using scratch_flags,
    // namely parallel_reduce and parallel_scan, to reset the used values to 0.
    KOKKOS_IMPL_HIP_SAFE_CALL(
        hip_memset_wrapper(m_scratchFlags, 0, alloc_size));
  }

  return m_scratchFlags;
}

Kokkos::HIP::size_type *HIPInternal::stage_functor_for_execution(
    void const *driver, std::size_t const size) const {
  if (verify_is_initialized("scratch_functor") && m_scratchFunctorSize < size) {
    auto device_mem_space = Kokkos::HIPSpace::impl_create(m_hipDev, m_stream);
    auto host_mem_space =
        Kokkos::HIPHostPinnedSpace::impl_create(m_hipDev, m_stream);

    if (m_scratchFunctor) {
      device_mem_space.deallocate(m_scratchFunctor, m_scratchFunctorSize);
      host_mem_space.deallocate(m_scratchFunctorHost, m_scratchFunctorSize);
    }

    m_scratchFunctorSize = size;

    m_scratchFunctor     = static_cast<size_type *>(device_mem_space.allocate(
        "Kokkos::InternalScratchFunctor", m_scratchFunctorSize));
    m_scratchFunctorHost = static_cast<size_type *>(host_mem_space.allocate(
        "Kokkos::InternalScratchFunctorHost", m_scratchFunctorSize));
  }

  // When using HSA_XNACK=1, it is necessary to copy the driver to the host to
  // ensure that the driver is not destroyed before the computation is done.
  // Without this fix, all the atomic tests fail. It is not obvious that this
  // problem is limited to HSA_XNACK=1 even if all the tests pass when
  // HSA_XNACK=0. That's why we always copy the driver.
  KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamSynchronize(m_stream));
  std::memcpy(m_scratchFunctorHost, driver, size);
  KOKKOS_IMPL_HIP_SAFE_CALL(hip_memcpy_async_wrapper(
      m_scratchFunctor, m_scratchFunctorHost, size, hipMemcpyDefault));

  return m_scratchFunctor;
}

int HIPInternal::acquire_team_scratch_space() {
  int current_team_scratch = 0;
  int zero                 = 0;
  while (!m_team_scratch_pool[current_team_scratch].compare_exchange_weak(
      zero, 1, std::memory_order_release, std::memory_order_relaxed)) {
    current_team_scratch = (current_team_scratch + 1) % m_n_team_scratch;
  }

  return current_team_scratch;
}

void *HIPInternal::resize_team_scratch_space(int scratch_pool_id,
                                             std::int64_t bytes,
                                             bool force_shrink) {
  // Multiple ParallelFor/Reduce Teams can call this function at the same time
  // and invalidate the m_team_scratch_ptr. We use a pool to avoid any race
  // condition.
  auto mem_space = Kokkos::HIPSpace::impl_create(m_hipDev, m_stream);
  if (m_team_scratch_current_size[scratch_pool_id] == 0) {
    m_team_scratch_current_size[scratch_pool_id] = bytes;
    m_team_scratch_ptr[scratch_pool_id] =
        mem_space.allocate("Kokkos::HIPSpace::TeamScratchMemory",
                           m_team_scratch_current_size[scratch_pool_id]);
  }
  if ((bytes > m_team_scratch_current_size[scratch_pool_id]) ||
      ((bytes < m_team_scratch_current_size[scratch_pool_id]) &&
       (force_shrink))) {
    mem_space.deallocate("Kokkos::HIPSpace::TeamScratchMemory",
                         m_team_scratch_ptr[scratch_pool_id],
                         m_team_scratch_current_size[scratch_pool_id]);
    m_team_scratch_current_size[scratch_pool_id] = bytes;
    m_team_scratch_ptr[scratch_pool_id] =
        mem_space.allocate("Kokkos::HIPSpace::TeamScratchMemory", bytes);
  }
  return m_team_scratch_ptr[scratch_pool_id];
}

void HIPInternal::release_team_scratch_space(int scratch_pool_id) {
  m_team_scratch_pool[scratch_pool_id] = 0;
}

//----------------------------------------------------------------------------

void HIPInternal::finalize() {
  // First, lock the shared resource locking helper.
  // Then, fence the stream and check if it was involved in the last constant
  // memory launch.
  // Locking is required to avoid a race condition, i.e. it prevents another
  // thread from launching another kernel in-between the fence
  // and the 'check_if_involved_and_unlock'.
  auto lock = HIPInternal::constantMemReusable[m_hipDev].lock();
  this->fence("Kokkos::HIPInternal::finalize: fence on finalization");
  HIPInternal::constantMemReusable[m_hipDev].check_if_involved_and_unlock(
      std::move(lock), m_stream);

  was_finalized = true;

  auto device_mem_space = Kokkos::HIPSpace::impl_create(m_hipDev, m_stream);
  if (nullptr != m_scratchSpace || nullptr != m_scratchFlags) {
    device_mem_space.deallocate(m_scratchFlags,
                                m_scratchSpaceCount * sizeScratchGrain);
    device_mem_space.deallocate(m_scratchSpace,
                                m_scratchFlagsCount * sizeScratchGrain);

    if (m_scratchFunctorSize > 0) {
      device_mem_space.deallocate(m_scratchFunctor, m_scratchFunctorSize);
      auto host_mem_space =
          Kokkos::HIPHostPinnedSpace::impl_create(m_hipDev, m_stream);
      host_mem_space.deallocate(m_scratchFunctorHost, m_scratchFunctorSize);
    }
  }

  for (int i = 0; i < m_n_team_scratch; ++i) {
    if (m_team_scratch_current_size[i] > 0)
      device_mem_space.deallocate(m_team_scratch_ptr[i],
                                  m_team_scratch_current_size[i]);
  }

  m_scratchSpaceCount = 0;
  m_scratchFlagsCount = 0;
  m_scratchSpace      = nullptr;
  m_scratchFlags      = nullptr;
  for (int i = 0; i < m_n_team_scratch; ++i) {
    m_team_scratch_current_size[i] = 0;
    m_team_scratch_ptr[i]          = nullptr;
  }

  KOKKOS_IMPL_HIP_SAFE_CALL(hip_free_wrapper(m_scratch_locks));
  m_scratch_locks     = nullptr;
  m_num_scratch_locks = 0;
  m_hipDev            = -1;
}

int HIPInternal::m_maxThreadsPerSM = 0;

hipDeviceProp_t HIPInternal::m_deviceProp;

std::mutex HIPInternal::scratchFunctorMutex;

std::set<int> HIPInternal::hip_devices                             = {};
std::map<int, unsigned long *> HIPInternal::constantMemHostStaging = {};
std::map<int, SharedResourceLock> HIPInternal::constantMemReusable = {};

//----------------------------------------------------------------------------

Kokkos::HIP::size_type hip_internal_multiprocessor_count() {
  return HIPInternal::singleton().m_deviceProp.multiProcessorCount;
}

Kokkos::HIP::size_type *hip_internal_scratch_space(const HIP &instance,
                                                   const std::size_t size) {
  return instance.impl_internal_space_instance()->scratch_space(size);
}

Kokkos::HIP::size_type *hip_internal_scratch_flags(const HIP &instance,
                                                   const std::size_t size) {
  return instance.impl_internal_space_instance()->scratch_flags(size);
}

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {
void hip_internal_error_throw(hipError_t e, const char *name, const char *file,
                              const int line) {
  std::ostringstream out;
  out << name << " error( " << hipGetErrorName(e)
      << "): " << hipGetErrorString(e);
  if (file) {
    out << " " << file << ":" << line;
  }
  throw_runtime_exception(out.str());
}
}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
