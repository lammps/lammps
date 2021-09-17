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
  HIP_SAFE_CALL(hipGetDeviceCount(&m_hipDevCount));

  if (m_hipDevCount > MAXIMUM_DEVICE_COUNT) {
    Kokkos::abort(
        "Sorry, you have more GPUs per node than we thought anybody would ever "
        "have. Please report this to github.com/kokkos/kokkos.");
  }
  for (int i = 0; i < m_hipDevCount; ++i) {
    HIP_SAFE_CALL(hipGetDeviceProperties(m_hipProp + i, i));
  }
}

const HIPInternalDevices &HIPInternalDevices::singleton() {
  static HIPInternalDevices self;
  return self;
}
}  // namespace

namespace Impl {

//----------------------------------------------------------------------------

void HIPInternal::print_configuration(std::ostream &s) const {
  const HIPInternalDevices &dev_info = HIPInternalDevices::singleton();

  s << "macro  KOKKOS_ENABLE_HIP : defined" << '\n';
#if defined(HIP_VERSION)
  s << "macro  HIP_VERSION = " << HIP_VERSION << " = version "
    << HIP_VERSION / 100 << "." << HIP_VERSION % 100 << '\n';
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
  if (m_scratchSpace || m_scratchFlags || m_scratchConcurrentBitset) {
    std::cerr << "Kokkos::Experimental::HIP ERROR: Failed to call "
                 "Kokkos::Experimental::HIP::finalize()"
              << std::endl;
    std::cerr.flush();
  }

  m_hipDev                  = -1;
  m_hipArch                 = -1;
  m_multiProcCount          = 0;
  m_maxWarpCount            = 0;
  m_maxSharedWords          = 0;
  m_maxShmemPerBlock        = 0;
  m_scratchSpaceCount       = 0;
  m_scratchFlagsCount       = 0;
  m_scratchSpace            = nullptr;
  m_scratchFlags            = nullptr;
  m_scratchConcurrentBitset = nullptr;
  m_stream                  = nullptr;
}

int HIPInternal::verify_is_initialized(const char *const label) const {
  if (m_hipDev < 0) {
    std::cerr << "Kokkos::Experimental::HIP::" << label
              << " : ERROR device not initialized" << std::endl;
  }
  return 0 <= m_hipDev;
}

HIPInternal &HIPInternal::singleton() {
  static HIPInternal *self = nullptr;
  if (!self) {
    self = new HIPInternal();
  }
  return *self;
}

void HIPInternal::fence() const {
  HIP_SAFE_CALL(hipStreamSynchronize(m_stream));
  // can reset our cycle id now as well
  m_cycleId = 0;
}

void HIPInternal::initialize(int hip_device_id, hipStream_t stream) {
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

    HIP_SAFE_CALL(hipSetDevice(m_hipDev));

    m_stream                    = stream;
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
    m_maxBlock = hipProp.maxGridSize[0];

    // theoretically, we can get 40 WF's / CU, but only can sustain 32
    // see
    // https://github.com/ROCm-Developer-Tools/HIP/blob/a0b5dfd625d99af7e288629747b40dd057183173/vdi/hip_platform.cpp#L742
    m_maxBlocksPerSM = 32;
    // FIXME_HIP - Nick to implement this upstream
    //             Register count comes from Sec. 2.2. "Data Sharing" of the
    //             Vega 7nm ISA document (see the diagram)
    //             https://developer.amd.com/wp-content/resources/Vega_7nm_Shader_ISA.pdf
    //             VGPRS = 4 (SIMD/CU) * 256 VGPR/SIMD * 64 registers / VGPR =
    //             65536 VGPR/CU
    m_regsPerSM        = 65536;
    m_shmemPerSM       = hipProp.maxSharedMemoryPerMultiProcessor;
    m_maxShmemPerBlock = hipProp.sharedMemPerBlock;
    m_maxThreadsPerSM  = m_maxBlocksPerSM * HIPTraits::WarpSize;
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

      m_scratchConcurrentBitset = reinterpret_cast<uint32_t *>(r->data());

      HIP_SAFE_CALL(hipMemset(m_scratchConcurrentBitset, 0,
                              sizeof(uint32_t) * buffer_bound));
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
}

//----------------------------------------------------------------------------

using ScratchGrain =
    Kokkos::Experimental::HIP::size_type[Impl::HIPTraits::WarpSize];
enum { sizeScratchGrain = sizeof(ScratchGrain) };

Kokkos::Experimental::HIP::size_type *HIPInternal::scratch_space(
    const Kokkos::Experimental::HIP::size_type size) {
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
    const Kokkos::Experimental::HIP::size_type size) {
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

    HIP_SAFE_CALL(
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
  this->fence();
  was_finalized = true;
  if (nullptr != m_scratchSpace || nullptr != m_scratchFlags) {
    using RecordHIP =
        Kokkos::Impl::SharedAllocationRecord<Kokkos::Experimental::HIPSpace>;

    RecordHIP::decrement(RecordHIP::get_record(m_scratchFlags));
    RecordHIP::decrement(RecordHIP::get_record(m_scratchSpace));
    RecordHIP::decrement(RecordHIP::get_record(m_scratchConcurrentBitset));

    if (m_team_scratch_current_size > 0)
      Kokkos::kokkos_free<Kokkos::Experimental::HIPSpace>(m_team_scratch_ptr);

    m_hipDev                    = -1;
    m_hipArch                   = -1;
    m_multiProcCount            = 0;
    m_maxWarpCount              = 0;
    m_maxBlock                  = 0;
    m_maxSharedWords            = 0;
    m_maxShmemPerBlock          = 0;
    m_scratchSpaceCount         = 0;
    m_scratchFlagsCount         = 0;
    m_scratchSpace              = nullptr;
    m_scratchFlags              = nullptr;
    m_scratchConcurrentBitset   = nullptr;
    m_stream                    = nullptr;
    m_team_scratch_current_size = 0;
    m_team_scratch_ptr          = nullptr;
  }
  if (nullptr != d_driverWorkArray) {
    HIP_SAFE_CALL(hipHostFree(d_driverWorkArray));
    d_driverWorkArray = nullptr;
  }
}

char *HIPInternal::get_next_driver(size_t driverTypeSize) const {
  std::lock_guard<std::mutex> const lock(m_mutexWorkArray);
  if (d_driverWorkArray == nullptr) {
    HIP_SAFE_CALL(
        hipHostMalloc(&d_driverWorkArray,
                      m_maxDriverCycles * m_maxDriverTypeSize * sizeof(char),
                      hipHostMallocNonCoherent));
  }
  if (driverTypeSize > m_maxDriverTypeSize) {
    // fence handles the cycle id reset for us
    fence();
    HIP_SAFE_CALL(hipHostFree(d_driverWorkArray));
    m_maxDriverTypeSize = driverTypeSize;
    if (m_maxDriverTypeSize % 128 != 0)
      m_maxDriverTypeSize =
          m_maxDriverTypeSize + 128 - m_maxDriverTypeSize % 128;
    HIP_SAFE_CALL(
        hipHostMalloc(&d_driverWorkArray,
                      m_maxDriverCycles * m_maxDriverTypeSize * sizeof(char),
                      hipHostMallocNonCoherent));
  } else {
    m_cycleId = (m_cycleId + 1) % m_maxDriverCycles;
    if (m_cycleId == 0) {
      // ensure any outstanding kernels are completed before we wrap around
      fence();
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

Kokkos::Experimental::HIP::size_type hip_internal_maximum_grid_count() {
  return HIPInternal::singleton().m_maxBlock;
}

Kokkos::Experimental::HIP::size_type *hip_internal_scratch_space(
    const Kokkos::Experimental::HIP::size_type size) {
  return HIPInternal::singleton().scratch_space(size);
}

Kokkos::Experimental::HIP::size_type *hip_internal_scratch_flags(
    const Kokkos::Experimental::HIP::size_type size) {
  return HIPInternal::singleton().scratch_flags(size);
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {
void hip_device_synchronize() { HIP_SAFE_CALL(hipDeviceSynchronize()); }

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
