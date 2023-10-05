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

//----------------------------------------------------------------------------

int HIPInternal::concurrency() {
  static int const concurrency = m_deviceProp.maxThreadsPerMultiProcessor *
                                 m_deviceProp.multiProcessorCount;
  return concurrency;
}

void HIPInternal::print_configuration(std::ostream &s) const {
  s << "macro  KOKKOS_ENABLE_HIP : defined" << '\n';
#if defined(HIP_VERSION)
  s << "macro  HIP_VERSION = " << HIP_VERSION << " = version "
    << HIP_VERSION_MAJOR << '.' << HIP_VERSION_MINOR << '.' << HIP_VERSION_PATCH
    << '\n';
#endif

  int hipDevCount;
  KOKKOS_IMPL_HIP_SAFE_CALL(hipGetDeviceCount(&hipDevCount));

  for (int i = 0; i < hipDevCount; ++i) {
    hipDeviceProp_t hipProp;
    KOKKOS_IMPL_HIP_SAFE_CALL(hipGetDeviceProperties(&hipProp, i));

    s << "Kokkos::HIP[ " << i << " ] "
      << "gcnArch " << hipProp.gcnArch << ", Total Global Memory: "
      << ::Kokkos::Impl::human_memory_size(hipProp.totalGlobalMem)
      << ", Shared Memory per Block: "
      << ::Kokkos::Impl::human_memory_size(hipProp.sharedMemPerBlock);
    if (m_hipDev == i) s << " : Selected";
    s << '\n';
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

void HIPInternal::initialize(hipStream_t stream, bool manage_stream) {
  if (was_finalized)
    Kokkos::abort("Calling HIP::initialize after HIP::finalize is illegal\n");

  if (is_initialized()) return;

  if (!HostSpace::execution_space::impl_is_initialized()) {
    const std::string msg(
        "HIP::initialize ERROR : HostSpace::execution_space "
        "is not initialized");
    Kokkos::Impl::throw_runtime_exception(msg);
  }

  const bool ok_init = nullptr == m_scratchSpace || nullptr == m_scratchFlags;

  if (ok_init) {
    m_stream        = stream;
    m_manage_stream = manage_stream;

    //----------------------------------
    // Multiblock reduction uses scratch flags for counters
    // and scratch space for partial reduction values.
    // Allocate some initial space.  This will grow as needed.
    {
      const unsigned reduce_block_count =
          m_maxWarpCount * Impl::HIPTraits::WarpSize;

      (void)scratch_flags(reduce_block_count * 2 * sizeof(size_type));
      (void)scratch_space(reduce_block_count * 16 * sizeof(size_type));
    }
  } else {
    std::ostringstream msg;
    msg << "Kokkos::HIP::initialize(" << m_hipDev
        << ") FAILED : Already initialized";
    Kokkos::Impl::throw_runtime_exception(msg.str());
  }

  m_num_scratch_locks = concurrency();
  KOKKOS_IMPL_HIP_SAFE_CALL(
      hipMalloc(&m_scratch_locks, sizeof(int32_t) * m_num_scratch_locks));
  KOKKOS_IMPL_HIP_SAFE_CALL(
      hipMemset(m_scratch_locks, 0, sizeof(int32_t) * m_num_scratch_locks));
}

//----------------------------------------------------------------------------

using ScratchGrain = Kokkos::HIP::size_type[Impl::HIPTraits::WarpSize];
enum { sizeScratchGrain = sizeof(ScratchGrain) };

Kokkos::HIP::size_type *HIPInternal::scratch_space(const std::size_t size) {
  if (verify_is_initialized("scratch_space") &&
      m_scratchSpaceCount * sizeScratchGrain < size) {
    m_scratchSpaceCount = (size + sizeScratchGrain - 1) / sizeScratchGrain;

    using Record = Kokkos::Impl::SharedAllocationRecord<Kokkos::HIPSpace, void>;

    if (m_scratchSpace) Record::decrement(Record::get_record(m_scratchSpace));

    Record *const r =
        Record::allocate(Kokkos::HIPSpace(), "Kokkos::InternalScratchSpace",
                         (sizeScratchGrain * m_scratchSpaceCount));

    Record::increment(r);

    m_scratchSpace = reinterpret_cast<size_type *>(r->data());
  }

  return m_scratchSpace;
}

Kokkos::HIP::size_type *HIPInternal::scratch_flags(const std::size_t size) {
  if (verify_is_initialized("scratch_flags") &&
      m_scratchFlagsCount * sizeScratchGrain < size) {
    m_scratchFlagsCount = (size + sizeScratchGrain - 1) / sizeScratchGrain;

    using Record = Kokkos::Impl::SharedAllocationRecord<Kokkos::HIPSpace, void>;

    if (m_scratchFlags) Record::decrement(Record::get_record(m_scratchFlags));

    Record *const r =
        Record::allocate(Kokkos::HIPSpace(), "Kokkos::InternalScratchFlags",
                         (sizeScratchGrain * m_scratchFlagsCount));

    Record::increment(r);

    m_scratchFlags = reinterpret_cast<size_type *>(r->data());

    KOKKOS_IMPL_HIP_SAFE_CALL(
        hipMemset(m_scratchFlags, 0, m_scratchFlagsCount * sizeScratchGrain));
  }

  return m_scratchFlags;
}

Kokkos::HIP::size_type *HIPInternal::stage_functor_for_execution(
    void const *driver, std::size_t const size) const {
  if (verify_is_initialized("scratch_functor") && m_scratchFunctorSize < size) {
    m_scratchFunctorSize = size;

    using Record = Kokkos::Impl::SharedAllocationRecord<Kokkos::HIPSpace, void>;
    using RecordHost =
        Kokkos::Impl::SharedAllocationRecord<Kokkos::HIPHostPinnedSpace, void>;

    if (m_scratchFunctor) {
      Record::decrement(Record::get_record(m_scratchFunctor));
      RecordHost::decrement(RecordHost::get_record(m_scratchFunctorHost));
    }

    Record *const r =
        Record::allocate(Kokkos::HIPSpace(), "Kokkos::InternalScratchFunctor",
                         m_scratchFunctorSize);
    RecordHost *const r_host = RecordHost::allocate(
        Kokkos::HIPHostPinnedSpace(), "Kokkos::InternalScratchFunctorHost",
        m_scratchFunctorSize);

    Record::increment(r);
    RecordHost::increment(r_host);

    m_scratchFunctor     = reinterpret_cast<size_type *>(r->data());
    m_scratchFunctorHost = reinterpret_cast<size_type *>(r_host->data());
  }

  // When using HSA_XNACK=1, it is necessary to copy the driver to the host to
  // ensure that the driver is not destroyed before the computation is done.
  // Without this fix, all the atomic tests fail. It is not obvious that this
  // problem is limited to HSA_XNACK=1 even if all the tests pass when
  // HSA_XNACK=0. That's why we always copy the driver.
  KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamSynchronize(m_stream));
  std::memcpy(m_scratchFunctorHost, driver, size);
  KOKKOS_IMPL_HIP_SAFE_CALL(hipMemcpyAsync(m_scratchFunctor,
                                           m_scratchFunctorHost, size,
                                           hipMemcpyDefault, m_stream));

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
  if (m_team_scratch_current_size[scratch_pool_id] == 0) {
    m_team_scratch_current_size[scratch_pool_id] = bytes;
    m_team_scratch_ptr[scratch_pool_id] =
        Kokkos::kokkos_malloc<Kokkos::HIPSpace>(
            "Kokkos::HIPSpace::TeamScratchMemory",
            m_team_scratch_current_size[scratch_pool_id]);
  }
  if ((bytes > m_team_scratch_current_size[scratch_pool_id]) ||
      ((bytes < m_team_scratch_current_size[scratch_pool_id]) &&
       (force_shrink))) {
    m_team_scratch_current_size[scratch_pool_id] = bytes;
    m_team_scratch_ptr[scratch_pool_id] =
        Kokkos::kokkos_realloc<Kokkos::HIPSpace>(
            m_team_scratch_ptr[scratch_pool_id],
            m_team_scratch_current_size[scratch_pool_id]);
  }
  return m_team_scratch_ptr[scratch_pool_id];
}

void HIPInternal::release_team_scratch_space(int scratch_pool_id) {
  m_team_scratch_pool[scratch_pool_id] = 0;
}

//----------------------------------------------------------------------------

void HIPInternal::finalize() {
  this->fence("Kokkos::HIPInternal::finalize: fence on finalization");
  was_finalized = true;

  if (this == &singleton()) {
    (void)Kokkos::Impl::hip_global_unique_token_locks(true);
    desul::Impl::finalize_lock_arrays();  // FIXME

    KOKKOS_IMPL_HIP_SAFE_CALL(hipHostFree(constantMemHostStaging));
    KOKKOS_IMPL_HIP_SAFE_CALL(hipEventDestroy(constantMemReusable));
  }

  if (nullptr != m_scratchSpace || nullptr != m_scratchFlags) {
    using RecordHIP = Kokkos::Impl::SharedAllocationRecord<Kokkos::HIPSpace>;

    RecordHIP::decrement(RecordHIP::get_record(m_scratchFlags));
    RecordHIP::decrement(RecordHIP::get_record(m_scratchSpace));

    if (m_scratchFunctorSize > 0) {
      RecordHIP::decrement(RecordHIP::get_record(m_scratchFunctor));
      RecordHIP::decrement(RecordHIP::get_record(m_scratchFunctorHost));
    }
  }

  for (int i = 0; i < m_n_team_scratch; ++i) {
    if (m_team_scratch_current_size[i] > 0)
      Kokkos::kokkos_free<Kokkos::HIPSpace>(m_team_scratch_ptr[i]);
  }

  if (m_manage_stream && m_stream != nullptr)
    KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamDestroy(m_stream));

  m_scratchSpaceCount = 0;
  m_scratchFlagsCount = 0;
  m_scratchSpace      = nullptr;
  m_scratchFlags      = nullptr;
  m_stream            = nullptr;
  for (int i = 0; i < m_n_team_scratch; ++i) {
    m_team_scratch_current_size[i] = 0;
    m_team_scratch_ptr[i]          = nullptr;
  }

  KOKKOS_IMPL_HIP_SAFE_CALL(hipFree(m_scratch_locks));
  m_scratch_locks     = nullptr;
  m_num_scratch_locks = 0;
}

//----------------------------------------------------------------------------

Kokkos::HIP::size_type hip_internal_multiprocessor_count() {
  return HIPInternal::singleton().m_multiProcCount;
}

Kokkos::HIP::size_type hip_internal_maximum_warp_count() {
  return HIPInternal::singleton().m_maxWarpCount;
}

std::array<Kokkos::HIP::size_type, 3> hip_internal_maximum_grid_count() {
  return HIPInternal::singleton().m_maxBlock;
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

namespace Kokkos {
HIP::size_type HIP::detect_device_count() {
  int hipDevCount;
  KOKKOS_IMPL_HIP_SAFE_CALL(hipGetDeviceCount(&hipDevCount));
  return hipDevCount;
}
}  // namespace Kokkos
