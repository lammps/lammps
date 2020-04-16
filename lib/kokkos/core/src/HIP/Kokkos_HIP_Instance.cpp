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

int HIPInternal::was_initialized = 0;
int HIPInternal::was_finalized   = 0;
//----------------------------------------------------------------------------

void HIPInternal::print_configuration(std::ostream & /*s*/) const {
  // FIXME_HIP
  Kokkos::abort("print_configuration not implemented!\n");
  /*const HIPInternalDevices & dev_info = HIPInternalDevices::singleton();

#if defined( KOKKOS_ENABLE_HIP )
    s << "macro  KOKKOS_ENABLE_HIP      : defined" << std::endl ;
#endif
#if defined( __hcc_version__ )
    s << "macro  __hcc_version__          = " << __hcc_version__
      << std::endl ;
#endif

  for ( int i = 0 ; i < dev_info.m_hipDevCount ; ++i ) {
    s << "Kokkos::Experimental::HIP[ " << i << " ] "
      << dev_info.m_hipProp[i].name
      << " version " << (dev_info.m_hipProp[i].major) << "." <<
dev_info.m_hipProp[i].minor
      << ", Total Global Memory: " <<
human_memory_size(dev_info.m_hipProp[i].totalGlobalMem)
      << ", Shared Memory per Wavefront: " <<
human_memory_size(dev_info.m_hipProp[i].sharedMemPerWavefront); if ( m_hipDev ==
i ) s << " : Selected" ; s << std::endl ;
  }*/
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
  m_scratchSpace      = 0;
  m_scratchFlags      = 0;
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

void HIPInternal::initialize(int hip_device_id) {
  if (was_finalized)
    Kokkos::abort("Calling HIP::initialize after HIP::finalize is illegal\n");

  if (is_initialized()) return;

  enum { WordSize = sizeof(size_type) };

  if (!HostSpace::execution_space::impl_is_initialized()) {
    const std::string msg(
        "HIP::initialize ERROR : HostSpace::execution_space "
        "is not initialized");
    Kokkos::Impl::throw_runtime_exception(msg);
  }

  const HIPInternalDevices &dev_info = HIPInternalDevices::singleton();

  const bool ok_init = 0 == m_scratchSpace || 0 == m_scratchFlags;

  // Need at least a GPU device
  const bool ok_id =
      0 <= hip_device_id && hip_device_id < dev_info.m_hipDevCount;

  if (ok_init && ok_id) {
    const struct hipDeviceProp_t &hipProp = dev_info.m_hipProp[hip_device_id];

    m_hipDev = hip_device_id;

    hipSetDevice(m_hipDev);

    // FIXME_HIP for now always uses default stream
    m_stream = 0;

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

    m_shmemPerSM         = hipProp.maxSharedMemoryPerMultiProcessor;
    m_maxShmemPerBlock   = hipProp.sharedMemPerBlock;
    m_maxThreadsPerSM    = hipProp.maxThreadsPerMultiProcessor;
    m_maxThreadsPerBlock = hipProp.maxThreadsPerBlock;

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
  // FIXME_HIP uncomment this when global variable works
  // if (m_stream == 0) ::Kokkos::Impl::initialize_host_hip_lock_arrays();
}

//----------------------------------------------------------------------------

typedef Kokkos::Experimental::HIP::size_type
    ScratchGrain[Impl::HIPTraits::WarpSize];
enum { sizeScratchGrain = sizeof(ScratchGrain) };

Kokkos::Experimental::HIP::size_type *HIPInternal::scratch_space(
    const Kokkos::Experimental::HIP::size_type size) {
  if (verify_is_initialized("scratch_space") &&
      m_scratchSpaceCount * sizeScratchGrain < size) {
    m_scratchSpaceCount = (size + sizeScratchGrain - 1) / sizeScratchGrain;

    typedef Kokkos::Impl::SharedAllocationRecord<Kokkos::Experimental::HIPSpace,
                                                 void>
        Record;

    static Record *const r = Record::allocate(
        Kokkos::Experimental::HIPSpace(), "InternalScratchSpace",
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

    typedef Kokkos::Impl::SharedAllocationRecord<Kokkos::Experimental::HIPSpace,
                                                 void>
        Record;

    Record *const r = Record::allocate(
        Kokkos::Experimental::HIPSpace(), "InternalScratchFlags",
        (sizeScratchGrain * m_scratchFlagsCount));

    Record::increment(r);

    m_scratchFlags = reinterpret_cast<size_type *>(r->data());

    hipMemset(m_scratchFlags, 0, m_scratchFlagsCount * sizeScratchGrain);
  }

  return m_scratchFlags;
}

//----------------------------------------------------------------------------

void HIPInternal::finalize() {
  HIP().fence();
  was_finalized = 1;
  if (0 != m_scratchSpace || 0 != m_scratchFlags) {
    typedef Kokkos::Impl::SharedAllocationRecord<Kokkos::Experimental::HIPSpace>
        RecordHIP;

    RecordHIP::decrement(RecordHIP::get_record(m_scratchFlags));
    RecordHIP::decrement(RecordHIP::get_record(m_scratchSpace));

    m_hipDev            = -1;
    m_hipArch           = -1;
    m_multiProcCount    = 0;
    m_maxWarpCount      = 0;
    m_maxBlock          = 0;
    m_maxSharedWords    = 0;
    m_maxShmemPerBlock  = 0;
    m_scratchSpaceCount = 0;
    m_scratchFlagsCount = 0;
    m_scratchSpace      = 0;
    m_scratchFlags      = 0;
  }
}

//----------------------------------------------------------------------------

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
