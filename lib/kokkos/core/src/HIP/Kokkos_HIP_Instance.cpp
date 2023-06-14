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

/*--------------------------------------------------------------------------*/
/* Kokkos interfaces */

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <Kokkos_Core.hpp>

#include <HIP/Kokkos_HIP_Instance.hpp>
#include <Kokkos_HIP.hpp>
#include <Kokkos_HIP_Space.hpp>
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
    [Kokkos::Experimental::Impl::HIPTraits::ConstantMemoryUsage /
     sizeof(unsigned long)];
#endif

namespace Kokkos {
namespace Impl {
Kokkos::View<uint32_t *, Kokkos::Experimental::HIPSpace>
hip_global_unique_token_locks(bool deallocate) {
  static Kokkos::View<uint32_t *, Kokkos::Experimental::HIPSpace> locks =
      Kokkos::View<uint32_t *, Kokkos::Experimental::HIPSpace>();
  if (!deallocate && locks.extent(0) == 0)
    locks = Kokkos::View<uint32_t *, Kokkos::Experimental::HIPSpace>(
        "Kokkos::UniqueToken<HIP>::m_locks",
        Kokkos::Experimental::HIP().concurrency());
  if (deallocate)
    locks = Kokkos::View<uint32_t *, Kokkos::Experimental::HIPSpace>();
  return locks;
}
}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos {
namespace Experimental {
namespace {
class HIPInternalDevices {
 public:
  enum { MAXIMUM_DEVICE_COUNT = 64 };
  struct hipDeviceProp_t m_hipProp[MAXIMUM_DEVICE_COUNT];
  int m_hipDevCount;

  HIPInternalDevices();

  static HIPInternalDevices const &singleton();
};

HIPInternalDevices::HIPInternalDevices() {
  KOKKOS_IMPL_HIP_SAFE_CALL(hipGetDeviceCount(&m_hipDevCount));

  if (m_hipDevCount > MAXIMUM_DEVICE_COUNT) {
    Kokkos::abort(
        "Sorry, you have more GPUs per node than we thought anybody would ever "
        "have. Please report this to github.com/kokkos/kokkos.");
  }
  for (int i = 0; i < m_hipDevCount; ++i) {
    KOKKOS_IMPL_HIP_SAFE_CALL(hipGetDeviceProperties(m_hipProp + i, i));
  }
}

const HIPInternalDevices &HIPInternalDevices::singleton() {
  static HIPInternalDevices self;
  return self;
}
}  // namespace

unsigned long *Impl::HIPInternal::constantMemHostStaging = nullptr;
hipEvent_t Impl::HIPInternal::constantMemReusable        = nullptr;
std::mutex Impl::HIPInternal::constantMemMutex;

namespace Impl {

//----------------------------------------------------------------------------

void HIPInternal::print_configuration(std::ostream &s) const {
  const HIPInternalDevices &dev_info = HIPInternalDevices::singleton();

  s << "macro  KOKKOS_ENABLE_HIP : defined" << '\n';
#if defined(HIP_VERSION)
  s << "macro  HIP_VERSION = " << HIP_VERSION << " = version "
    << HIP_VERSION_MAJOR << '.' << HIP_VERSION_MINOR << '.' << HIP_VERSION_PATCH
    << '\n';
#endif

  for (int i = 0; i < dev_info.m_hipDevCount; ++i) {
    s << "Kokkos::Experimental::HIP[ " << i << " ] "
      << dev_info.m_hipProp[i].name << " version "
      << (dev_info.m_hipProp[i].major) << "." << dev_info.m_hipProp[i].minor
      << ", Total Global Memory: "
      << ::Kokkos::Impl::human_memory_size(dev_info.m_hipProp[i].totalGlobalMem)
      << ", Shared Memory per Block: "
      << ::Kokkos::Impl::human_memory_size(
             dev_info.m_hipProp[i].sharedMemPerBlock);
    if (m_hipDev == i) s << " : Selected";
    s << '\n';
  }
}

//----------------------------------------------------------------------------

HIPInternal::~HIPInternal() {
  if (m_scratchSpace || m_scratchFlags) {
    std::cerr << "Kokkos::Experimental::HIP ERROR: Failed to call "
                 "Kokkos::Experimental::HIP::finalize()"
              << std::endl;
    std::cerr.flush();
  }

  m_hipDev            = -1;
  m_hipArch           = -1;
  m_multiProcCount    = 0;
  m_maxWarpCount      = 0;
  m_maxSharedWords    = 0;
  m_maxShmemPerBlock  = 0;
  m_scratchSpaceCount = 0;
  m_scratchFlagsCount = 0;
  m_scratchSpace      = nullptr;
  m_scratchFlags      = nullptr;
  m_stream            = nullptr;
}

int HIPInternal::verify_is_initialized(const char *const label) const {
  if (m_hipDev < 0) {
    Kokkos::abort((std::string("Kokkos::Experimental::HIP::") + label +
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
  Kokkos::Tools::Experimental::Impl::profile_fence_event<
      Kokkos::Experimental::HIP>(
      name,
      Kokkos::Tools::Experimental::Impl::DirectFenceIDHandle{
          impl_get_instance_id()},
      [&]() {
        KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamSynchronize(m_stream));
        // can reset our cycle id now as well
        m_cycleId = 0;
      });
}

void HIPInternal::initialize(int hip_device_id, hipStream_t stream,
                             bool manage_stream) {
  if (was_finalized)
    Kokkos::abort("Calling HIP::initialize after HIP::finalize is illegal\n");

  if (is_initialized()) return;

  int constexpr WordSize = sizeof(size_type);

  if (!HostSpace::execution_space::impl_is_initialized()) {
    const std::string msg(
        "HIP::initialize ERROR : HostSpace::execution_space "
        "is not initialized");
    Kokkos::Impl::throw_runtime_exception(msg);
  }

  const HIPInternalDevices &dev_info = HIPInternalDevices::singleton();

  const bool ok_init = nullptr == m_scratchSpace || nullptr == m_scratchFlags;

  // Need at least a GPU device
  const bool ok_id =
      0 <= hip_device_id && hip_device_id < dev_info.m_hipDevCount;

  if (ok_init && ok_id) {
    const struct hipDeviceProp_t &hipProp = dev_info.m_hipProp[hip_device_id];

    m_hipDev     = hip_device_id;
    m_deviceProp = hipProp;

    KOKKOS_IMPL_HIP_SAFE_CALL(hipSetDevice(m_hipDev));

    m_stream                    = stream;
    m_manage_stream             = manage_stream;
    m_team_scratch_current_size = 0;
    m_team_scratch_ptr          = nullptr;

    // number of multiprocessors
    m_multiProcCount = hipProp.multiProcessorCount;

    //----------------------------------
    // Maximum number of warps,
    // at most one warp per thread in a warp for reduction.
    m_maxWarpCount = hipProp.maxThreadsPerBlock / Impl::HIPTraits::WarpSize;
    if (HIPTraits::WarpSize < m_maxWarpCount) {
      m_maxWarpCount = Impl::HIPTraits::WarpSize;
    }
    m_maxSharedWords = hipProp.sharedMemPerBlock / WordSize;

    //----------------------------------
    // Maximum number of blocks
    m_maxBlock[0] = hipProp.maxGridSize[0];
    m_maxBlock[1] = hipProp.maxGridSize[1];
    m_maxBlock[2] = hipProp.maxGridSize[2];

    // theoretically, we can get 40 WF's / CU, but only can sustain 32
    // see
    // https://github.com/ROCm-Developer-Tools/HIP/blob/a0b5dfd625d99af7e288629747b40dd057183173/vdi/hip_platform.cpp#L742
    m_maxWavesPerCU = 32;
    // FIXME_HIP - Nick to implement this upstream
    //             Register count comes from Sec. 2.2. "Data Sharing" of the
    //             Vega 7nm ISA document (see the diagram)
    //             https://developer.amd.com/wp-content/resources/Vega_7nm_Shader_ISA.pdf
    //             VGPRS = 4 (SIMD/CU) * 256 VGPR/SIMD * 64 registers / VGPR =
    //             65536 VGPR/CU
    m_regsPerSM        = 65536;
    m_shmemPerSM       = hipProp.maxSharedMemoryPerMultiProcessor;
    m_maxShmemPerBlock = hipProp.sharedMemPerBlock;
    m_maxThreadsPerSM  = m_maxWavesPerCU * HIPTraits::WarpSize;
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
    //----------------------------------
    // Concurrent bitset for obtaining unique tokens from within
    // an executing kernel.
    {
      const int32_t buffer_bound =
          Kokkos::Impl::concurrent_bitset::buffer_bound(HIP::concurrency());

      // Allocate and initialize uint32_t[ buffer_bound ]

      using Record =
          Kokkos::Impl::SharedAllocationRecord<Kokkos::Experimental::HIPSpace,
                                               void>;

      Record *const r = Record::allocate(Kokkos::Experimental::HIPSpace(),
                                         "Kokkos::InternalScratchBitset",
                                         sizeof(uint32_t) * buffer_bound);

      Record::increment(r);
    }
    //----------------------------------

  } else {
    std::ostringstream msg;
    msg << "Kokkos::Experimental::HIP::initialize(" << hip_device_id
        << ") FAILED";

    if (!ok_init) {
      msg << " : Already initialized";
    }
    if (!ok_id) {
      msg << " : Device identifier out of range "
          << "[0.." << dev_info.m_hipDevCount - 1 << "]";
    }
    Kokkos::Impl::throw_runtime_exception(msg.str());
  }

  // Init the array for used for arbitrarily sized atomics
  if (m_stream == nullptr) ::Kokkos::Impl::initialize_host_hip_lock_arrays();

  // Allocate a staging buffer for constant mem in pinned host memory
  // and an event to avoid overwriting driver for previous kernel launches
  if (m_stream == nullptr) {
    KOKKOS_IMPL_HIP_SAFE_CALL(hipHostMalloc((void **)&constantMemHostStaging,
                                            HIPTraits::ConstantMemoryUsage));

    KOKKOS_IMPL_HIP_SAFE_CALL(hipEventCreate(&constantMemReusable));
  }

  KOKKOS_IMPL_HIP_SAFE_CALL(
      hipMalloc(&m_scratch_locks, sizeof(int32_t) * HIP::concurrency()));
  KOKKOS_IMPL_HIP_SAFE_CALL(
      hipMemset(m_scratch_locks, 0, sizeof(int32_t) * HIP::concurrency()));
}

//----------------------------------------------------------------------------

using ScratchGrain =
    Kokkos::Experimental::HIP::size_type[Impl::HIPTraits::WarpSize];
enum { sizeScratchGrain = sizeof(ScratchGrain) };

Kokkos::Experimental::HIP::size_type *HIPInternal::scratch_space(
    const std::size_t size) {
  if (verify_is_initialized("scratch_space") &&
      m_scratchSpaceCount * sizeScratchGrain < size) {
    m_scratchSpaceCount = (size + sizeScratchGrain - 1) / sizeScratchGrain;

    using Record =
        Kokkos::Impl::SharedAllocationRecord<Kokkos::Experimental::HIPSpace,
                                             void>;

    if (m_scratchSpace) Record::decrement(Record::get_record(m_scratchSpace));

    Record *const r = Record::allocate(
        Kokkos::Experimental::HIPSpace(), "Kokkos::InternalScratchSpace",
        (sizeScratchGrain * m_scratchSpaceCount));

    Record::increment(r);

    m_scratchSpace = reinterpret_cast<size_type *>(r->data());
  }

  return m_scratchSpace;
}

Kokkos::Experimental::HIP::size_type *HIPInternal::scratch_flags(
    const std::size_t size) {
  if (verify_is_initialized("scratch_flags") &&
      m_scratchFlagsCount * sizeScratchGrain < size) {
    m_scratchFlagsCount = (size + sizeScratchGrain - 1) / sizeScratchGrain;

    using Record =
        Kokkos::Impl::SharedAllocationRecord<Kokkos::Experimental::HIPSpace,
                                             void>;

    if (m_scratchFlags) Record::decrement(Record::get_record(m_scratchFlags));

    Record *const r = Record::allocate(
        Kokkos::Experimental::HIPSpace(), "Kokkos::InternalScratchFlags",
        (sizeScratchGrain * m_scratchFlagsCount));

    Record::increment(r);

    m_scratchFlags = reinterpret_cast<size_type *>(r->data());

    KOKKOS_IMPL_HIP_SAFE_CALL(
        hipMemset(m_scratchFlags, 0, m_scratchFlagsCount * sizeScratchGrain));
  }

  return m_scratchFlags;
}

void *HIPInternal::resize_team_scratch_space(std::int64_t bytes,
                                             bool force_shrink) {
  if (m_team_scratch_current_size == 0) {
    m_team_scratch_current_size = bytes;
    m_team_scratch_ptr = Kokkos::kokkos_malloc<Kokkos::Experimental::HIPSpace>(
        "Kokkos::HIPSpace::TeamScratchMemory", m_team_scratch_current_size);
  }
  if ((bytes > m_team_scratch_current_size) ||
      ((bytes < m_team_scratch_current_size) && (force_shrink))) {
    m_team_scratch_current_size = bytes;
    m_team_scratch_ptr = Kokkos::kokkos_realloc<Kokkos::Experimental::HIPSpace>(
        m_team_scratch_ptr, m_team_scratch_current_size);
  }
  return m_team_scratch_ptr;
}

//----------------------------------------------------------------------------

void HIPInternal::finalize() {
  this->fence("Kokkos::HIPInternal::finalize: fence on finalization");
  was_finalized = true;

  if (this == &singleton()) {
    (void)Kokkos::Impl::hip_global_unique_token_locks(true);
    Kokkos::Impl::finalize_host_hip_lock_arrays();

    KOKKOS_IMPL_HIP_SAFE_CALL(hipHostFree(constantMemHostStaging));
    KOKKOS_IMPL_HIP_SAFE_CALL(hipEventDestroy(constantMemReusable));
  }

  if (nullptr != m_scratchSpace || nullptr != m_scratchFlags) {
    using RecordHIP =
        Kokkos::Impl::SharedAllocationRecord<Kokkos::Experimental::HIPSpace>;

    RecordHIP::decrement(RecordHIP::get_record(m_scratchFlags));
    RecordHIP::decrement(RecordHIP::get_record(m_scratchSpace));

    if (m_team_scratch_current_size > 0)
      Kokkos::kokkos_free<Kokkos::Experimental::HIPSpace>(m_team_scratch_ptr);

    if (m_manage_stream && m_stream != nullptr)
      KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamDestroy(m_stream));
  }

  m_hipDev                    = -1;
  m_hipArch                   = -1;
  m_multiProcCount            = 0;
  m_maxWarpCount              = 0;
  m_maxBlock                  = {0, 0, 0};
  m_maxSharedWords            = 0;
  m_maxShmemPerBlock          = 0;
  m_scratchSpaceCount         = 0;
  m_scratchFlagsCount         = 0;
  m_scratchSpace              = nullptr;
  m_scratchFlags              = nullptr;
  m_stream                    = nullptr;
  m_team_scratch_current_size = 0;
  m_team_scratch_ptr          = nullptr;

  KOKKOS_IMPL_HIP_SAFE_CALL(hipFree(m_scratch_locks));
  m_scratch_locks = nullptr;

  if (nullptr != d_driverWorkArray) {
    KOKKOS_IMPL_HIP_SAFE_CALL(hipHostFree(d_driverWorkArray));
    d_driverWorkArray = nullptr;
  }
}

char *HIPInternal::get_next_driver(size_t driverTypeSize) const {
  if (d_driverWorkArray == nullptr) {
    KOKKOS_IMPL_HIP_SAFE_CALL(
        hipHostMalloc(&d_driverWorkArray,
                      m_maxDriverCycles * m_maxDriverTypeSize * sizeof(char),
                      hipHostMallocNonCoherent));
  }
  if (driverTypeSize > m_maxDriverTypeSize) {
    // fence handles the cycle id reset for us
    fence(
        "Kokkos::HIPInternal::get_next_driver: fence before reallocating "
        "resources");
    KOKKOS_IMPL_HIP_SAFE_CALL(hipHostFree(d_driverWorkArray));
    m_maxDriverTypeSize = driverTypeSize;
    if (m_maxDriverTypeSize % 128 != 0)
      m_maxDriverTypeSize =
          m_maxDriverTypeSize + 128 - m_maxDriverTypeSize % 128;
    KOKKOS_IMPL_HIP_SAFE_CALL(
        hipHostMalloc(&d_driverWorkArray,
                      m_maxDriverCycles * m_maxDriverTypeSize * sizeof(char),
                      hipHostMallocNonCoherent));
  } else {
    m_cycleId = (m_cycleId + 1) % m_maxDriverCycles;
    if (m_cycleId == 0) {
      // ensure any outstanding kernels are completed before we wrap around
      fence(
          "Kokkos::HIPInternal::get_next_driver: fence before reusing first "
          "driver");
    }
  }
  return &d_driverWorkArray[m_maxDriverTypeSize * m_cycleId];
}

//----------------------------------------------------------------------------

Kokkos::Experimental::HIP::size_type hip_internal_multiprocessor_count() {
  return HIPInternal::singleton().m_multiProcCount;
}

Kokkos::Experimental::HIP::size_type hip_internal_maximum_warp_count() {
  return HIPInternal::singleton().m_maxWarpCount;
}

std::array<Kokkos::Experimental::HIP::size_type, 3>
hip_internal_maximum_grid_count() {
  return HIPInternal::singleton().m_maxBlock;
}

Kokkos::Experimental::HIP::size_type *hip_internal_scratch_space(
    const HIP &instance, const std::size_t size) {
  return instance.impl_internal_space_instance()->scratch_space(size);
}

Kokkos::Experimental::HIP::size_type *hip_internal_scratch_flags(
    const HIP &instance, const std::size_t size) {
  return instance.impl_internal_space_instance()->scratch_flags(size);
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {
void hip_device_synchronize(const std::string &name) {
  Kokkos::Tools::Experimental::Impl::profile_fence_event<
      Kokkos::Experimental::HIP>(
      name,
      Kokkos::Tools::Experimental::SpecialSynchronizationCases::
          GlobalDeviceSynchronization,
      [&]() { KOKKOS_IMPL_HIP_SAFE_CALL(hipDeviceSynchronize()); });
}

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
namespace Experimental {
HIP::size_type HIP::detect_device_count() {
  return HIPInternalDevices::singleton().m_hipDevCount;
}
}  // namespace Experimental
}  // namespace Kokkos
