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

/* only compile this file if ROCM is enabled for Kokkos */
#ifdef KOKKOS_ENABLE_ROCM

//#include <ROCm/Kokkos_ROCm_Internal.hpp>
#include <impl/Kokkos_Error.hpp>
#include <Kokkos_ROCmSpace.hpp>
#include <ROCm/Kokkos_ROCm_Exec.hpp>

/*--------------------------------------------------------------------------*/
/* Standard 'C' libraries */
#include <stdlib.h>

/* Standard 'C++' libraries */
#include <vector>
#include <iostream>
#include <sstream>
#include <string>

// KOKKOS_INLINE_FUNCTION
// Kokkos::Impl::ROCmLockArraysStruct kokkos_impl_rocm_lock_arrays ;

/*--------------------------------------------------------------------------*/
namespace Kokkos {
namespace Impl {

#if 0
namespace {
__global__
void query_rocm_kernel_arch( int * d_arch )
{
#if defined(__HCC_ACCELERATOR__)
  *d_arch = OCM_ARCH__ ;
#else
  *d_arch = 0 ;
#endif
}

/** Query what compute capability is actually launched to the device: */
int rocm_kernel_arch()
{
  int * d_arch = 0 ;
  rocmMalloc( (void **) & d_arch , sizeof(int) );
  query_rocm_kernel_arch<<<1,1>>>( d_arch );
  int arch = 0 ;
  rocmMemcpy( & arch , d_arch , sizeof(int) , rocmMemcpyDefault );
  rocmFree( d_arch );
  return arch ;
}
bool rocm_launch_blocking()
{
  const char * env = getenv("ROCM_LAUNCH_BLOCKING");

  if (env == 0) return false;

  return atoi(env);
}

}
#endif

// true device memory allocation, not visible from host
void* rocm_device_allocate(int size) {
  void* ptr;
  hc::accelerator acc;
  ptr = hc::am_alloc(size, acc, 0);
  return ptr;
}

// host pinned allocation
// flag = 1, non-coherent, host resident, but with gpu address space pointer
// flag = 2, coherent, host resident, but with host address space pointer
void* rocm_hostpinned_allocate(int size) {
  void* ptr;
  hc::accelerator acc;
  ptr = hc::am_alloc(size, acc, 2);
  return ptr;
}
// same free used by all rocm memory allocations
void rocm_device_free(void* ptr) { hc::am_free(ptr); }

KOKKOS_INLINE_FUNCTION
void rocm_device_synchronize() {
  hc::accelerator_view av   = hc::accelerator().get_default_view();
  hc::completion_future fut = av.create_marker();
  fut.wait();
}

void rocm_internal_error_throw(const char* name, const char* file,
                               const int line) {
#if 0
  std::ostringstream out ;
  out << name << " error( " << rocmGetErrorName(e) << "): " << rocmGetErrorString(e);
  if (file) {
    out << " " << file << ":" << line;
  }
  throw_runtime_exception( out.str() );
#endif
}

//----------------------------------------------------------------------------
// Some significant rocm device properties:
//
// rocmDeviceProp::name                : Text label for device
// rocmDeviceProp::major               : Device major number
// rocmDeviceProp::minor               : Device minor number
// rocmDeviceProp::workgroupSize       : number of threads per workgroup
// rocmDeviceProp::multiProcessorCount : number of multiprocessors
// rocmDeviceProp::sharedMemPerBlock   : capacity of shared memory per wavefront
// rocmDeviceProp::totalConstMem       : capacity of constant memory
// rocmDeviceProp::totalGlobalMem      : capacity of global memory
// rocmDeviceProp::maxGridSize[3]      : maximum grid size

//
//
// the data we have available from a ROCm accelerator
// std::wstring get_device_path()
// std::wstring get_description()
// unsigned int get_version()
// bool get_has_display()
// size_t get_dedicated_memory()
// bool get_supports_double_precision()
// bool get_supports_limited_double_precision()
// bool get_is_debug()
// bool get_supports_cpu_shared_memory()
// size_t get_max_tile_static_size()
// unsigned int get_cu_count()
// bool has_cpu_accessible_am()
struct rocmDeviceProp {
  char name[256];
  char description[256];
  unsigned int version;
  int device_type;
  int device_ordinal;
  int major;
  int minor;
  size_t totalGlobalMem;
  size_t sharedMemPerWavefront;
  int WavefrontSize;
  int WorkgroupSize;
  int MaxTileCount;
  int maxThreadsPerWorkgroup;
  int multiProcessorCount;
  int canMapHostMemory;
  bool APU;
};

void rocmGetDeviceProperties(struct rocmDeviceProp* devProp, int device) {
  std::wstring s;
  int i, n;
  hc::accelerator acc;
  std::vector<hc::accelerator> accv = acc.get_all();

  hc::accelerator a = accv[device];

  s = a.get_device_path();
  i = 0;
  for (wchar_t c : s)
    if ((n = std::wctomb(&devProp->name[i], c)) > 0) i += n;

  /* assume a CPU */
  devProp->version = a.get_version();
  devProp->major   = a.get_version() >> 16;  // for CPU, these are meaningless
  devProp->minor   = a.get_version() & 0xff;
  devProp->device_ordinal = 0;

  /* is this an AMD graphics card */
  if ((devProp->name[0] == 'g') && (devProp->name[1] == 'f') &&
      (devProp->name[2] == 'x')) {
    /* for AMD cards, the name has the format gfxMmmO */

    devProp->device_type = ((devProp->name[3] - 0x30) << 16) +
                           ((devProp->name[4] - 0x30) << 8) +
                           (devProp->name[5] - 0x30);
    devProp->device_ordinal = devProp->name[6] - 0x30;
    devProp->major          = devProp->name[3] - 0x30;
    devProp->minor          = devProp->name[5] - 0x30;
  }

  s = a.get_description();
  i = 0;
  for (wchar_t c : s)
    if ((n = std::wctomb(&devProp->description[i], c)) > 0) i += n;
  devProp->totalGlobalMem        = a.get_dedicated_memory();
  devProp->sharedMemPerWavefront = a.get_max_tile_static_size();
  devProp->WavefrontSize         = 64;
  devProp->WorkgroupSize         = 256;  // preferred
  devProp->MaxTileCount =
      409600;  // as defined in /opt/rocm/hcc-lc/include/hsa_new.h
  devProp->maxThreadsPerWorkgroup = 1024;
  devProp->multiProcessorCount    = a.get_cu_count();
  devProp->canMapHostMemory       = a.get_supports_cpu_shared_memory();
  // Kaveri has 64KB L2 per CU, 16KB L1, 64KB Vector Regs/SIMD, or 128
  // regs/thread GCN has 64KB LDS per CU

  // Kaveri APU is 7:0:0
  // Carrizo APU is 8:0:1
  devProp->APU = (((devProp->major == 7) && (devProp->minor == 0)) |
                  ((devProp->major == 8) && (devProp->minor == 1)))
                     ? true
                     : false;
}

namespace {

class ROCmInternalDevices {
 public:
  enum { MAXIMUM_DEVICE_COUNT = 64 };
  struct rocmDeviceProp m_rocmProp[MAXIMUM_DEVICE_COUNT];
  int m_rocmDevCount;

  ROCmInternalDevices();

  static const ROCmInternalDevices& singleton();
};

ROCmInternalDevices::ROCmInternalDevices() {
  hc::accelerator acc;
  std::vector<hc::accelerator> accv = acc.get_all();
  m_rocmDevCount                    = accv.size();

  if (m_rocmDevCount > MAXIMUM_DEVICE_COUNT) {
    Kokkos::abort(
        "Sorry, you have more GPUs per node than we thought anybody would ever "
        "have. Please report this to github.com/kokkos/kokkos.");
  }
  for (int i = 0; i < m_rocmDevCount; ++i) {
    rocmGetDeviceProperties(m_rocmProp + i, i);
  }
}

const ROCmInternalDevices& ROCmInternalDevices::singleton() {
  static ROCmInternalDevices* self = nullptr;
  if (!self) {
    self = new ROCmInternalDevices();
  }
  return *self;
}

}  // namespace

//----------------------------------------------------------------------------

class ROCmInternal {
 private:
  ROCmInternal(const ROCmInternal&);
  ROCmInternal& operator=(const ROCmInternal&);

 public:
  typedef Kokkos::Experimental::ROCm::size_type size_type;

  int m_rocmDev;
  int m_rocmArch;
  unsigned m_multiProcCount;
  unsigned m_maxWorkgroup;
  unsigned m_maxSharedWords;
  size_type m_scratchSpaceCount;
  size_type m_scratchFlagsCount;
  size_type* m_scratchSpace;
  size_type* m_scratchFlags;

  static int was_finalized;

  static ROCmInternal& singleton();

  int verify_is_initialized(const char* const label) const;

  int is_initialized() const {
    return 0 != m_scratchSpace && 0 != m_scratchFlags;
  }

  void initialize(int rocm_device_id);
  void finalize();

  void print_configuration(std::ostream&) const;

  ~ROCmInternal();

  ROCmInternal()
      : m_rocmDev(-1),
        m_rocmArch(-1),
        m_multiProcCount(0),
        m_maxWorkgroup(0),
        m_maxSharedWords(0),
        m_scratchSpaceCount(0),
        m_scratchFlagsCount(0),
        m_scratchSpace(0),
        m_scratchFlags(0) {}

  size_type* scratch_space(const size_type size);
  size_type* scratch_flags(const size_type size);
};

int ROCmInternal::was_finalized = 0;
//----------------------------------------------------------------------------

void ROCmInternal::print_configuration(std::ostream& s) const {
  const ROCmInternalDevices& dev_info = ROCmInternalDevices::singleton();

#if defined(KOKKOS_ENABLE_ROCM)
  s << "macro  KOKKOS_ENABLE_ROCM      : defined" << std::endl;
#endif
#if defined(__hcc_version__)
  s << "macro  __hcc_version__          = " << __hcc_version__ << std::endl;
#endif

  for (int i = 0; i < dev_info.m_rocmDevCount; ++i) {
    s << "Kokkos::Experimental::ROCm[ " << i << " ] "
      << dev_info.m_rocmProp[i].name << " version "
      << (dev_info.m_rocmProp[i].major) << "." << dev_info.m_rocmProp[i].minor
      << ", Total Global Memory: "
      << human_memory_size(dev_info.m_rocmProp[i].totalGlobalMem)
      << ", Shared Memory per Wavefront: "
      << human_memory_size(dev_info.m_rocmProp[i].sharedMemPerWavefront);
    if (m_rocmDev == i) s << " : Selected";
    s << std::endl;
  }
}

//----------------------------------------------------------------------------

ROCmInternal::~ROCmInternal() {
  if (m_scratchSpace || m_scratchFlags) {
    std::cerr << "Kokkos::Experimental::ROCm ERROR: Failed to call "
                 "Kokkos::Experimental::ROCm::finalize()"
              << std::endl;
    std::cerr.flush();
  }

  m_rocmDev           = -1;
  m_rocmArch          = -1;
  m_multiProcCount    = 0;
  m_maxWorkgroup      = 0;
  m_maxSharedWords    = 0;
  m_scratchSpaceCount = 0;
  m_scratchFlagsCount = 0;
  m_scratchSpace      = 0;
  m_scratchFlags      = 0;
}

int ROCmInternal::verify_is_initialized(const char* const label) const {
  if (m_rocmDev < 0) {
    std::cerr << "Kokkos::Experimental::ROCm::" << label
              << " : ERROR device not initialized" << std::endl;
  }
  return 0 <= m_rocmDev;
}

ROCmInternal& ROCmInternal::singleton() {
  static ROCmInternal* self = nullptr;
  if (!self) {
    self = new ROCmInternal();
  }
  return *self;
}

void ROCmInternal::initialize(int rocm_device_id) {
  if (was_finalized)
    Kokkos::abort("Calling ROCm::initialize after ROCm::finalize is illegal\n");

  if (is_initialized()) return;

  enum { WordSize = sizeof(size_type) };

  if (!HostSpace::execution_space::is_initialized()) {
    const std::string msg(
        "ROCm::initialize ERROR : HostSpace::execution_space is not "
        "initialized");
    throw_runtime_exception(msg);
  }

  const ROCmInternalDevices& dev_info = ROCmInternalDevices::singleton();

  const bool ok_init = 0 == m_scratchSpace || 0 == m_scratchFlags;

  const bool ok_id =
      1 <= rocm_device_id && rocm_device_id < dev_info.m_rocmDevCount;

  // Need at least a GPU device

  const bool ok_dev =
      ok_id && (1 <= dev_info.m_rocmProp[rocm_device_id].major &&
                0 <= dev_info.m_rocmProp[rocm_device_id].minor);
  if (ok_init && ok_dev) {
    const struct rocmDeviceProp& rocmProp = dev_info.m_rocmProp[rocm_device_id];

    m_rocmDev = rocm_device_id;

    //  rocmSetDevice( m_rocmDev ) );
    Kokkos::Impl::rocm_device_synchronize();

    /*
        // Query what compute capability architecture a kernel executes:
        m_rocmArch = rocm_kernel_arch();
        if ( m_rocmArch != rocmProp.major * 100 + rocmProp.minor * 10 ) {
          std::cerr << "Kokkos::Experimental::ROCm::initialize WARNING: running
       kernels compiled for compute capability "
                    << ( m_rocmArch / 100 ) << "." << ( ( m_rocmArch % 100 ) /
       10 )
                    << " on device with compute capability "
                    << rocmProp.major << "." << rocmProp.minor
                    << " , this will likely reduce potential performance."
                    << std::endl ;
        }
    */
    // number of multiprocessors

    m_multiProcCount = rocmProp.multiProcessorCount;

    //----------------------------------
    // Maximum number of wavefronts,
    // at most one workgroup per thread in a workgroup for reduction.

    m_maxSharedWords = rocmProp.sharedMemPerWavefront / WordSize;

    //----------------------------------
    // Maximum number of Workgroups:

    m_maxWorkgroup =
        5 * rocmProp.multiProcessorCount;  // TODO: confirm usage and value

    //----------------------------------
    // Multiblock reduction uses scratch flags for counters
    // and scratch space for partial reduction values.
    // Allocate some initial space.  This will grow as needed.

    {
      const unsigned reduce_block_count =
          m_maxWorkgroup * Impl::ROCmTraits::WorkgroupSize;

      (void)scratch_flags(reduce_block_count * 2 * sizeof(size_type));
      (void)scratch_space(reduce_block_count * 16 * sizeof(size_type));
    }
    //----------------------------------

  } else {
    std::ostringstream msg;
    msg << "Kokkos::Experimental::ROCm::initialize(" << rocm_device_id
        << ") FAILED";

    if (!ok_init) {
      msg << " : Already initialized";
    }
    if (!ok_id) {
      msg << " : Device identifier out of range "
          << "[0.." << (dev_info.m_rocmDevCount - 1) << "]";
    } else if (!ok_dev) {
      msg << " : Device ";
      msg << dev_info.m_rocmProp[rocm_device_id].major;
      msg << ".";
      msg << dev_info.m_rocmProp[rocm_device_id].minor;
      msg << " Need at least a GPU";
      msg << std::endl;
    }
    Kokkos::Impl::throw_runtime_exception(msg.str());
  }

  // Init the array for used for arbitrarily sized atomics
  Kokkos::Impl::init_lock_arrays_rocm_space();

  //  Kokkos::Impl::ROCmLockArraysStruct locks;
  //  locks.atomic = atomic_lock_array_rocm_space_ptr(false);
  //  locks.scratch = scratch_lock_array_rocm_space_ptr(false);
  //  locks.threadid = threadid_lock_array_rocm_space_ptr(false);
  //  rocmMemcpyToSymbol( kokkos_impl_rocm_lock_arrays , & locks ,
  //  sizeof(ROCmLockArraysStruct) );
}

//----------------------------------------------------------------------------

typedef Kokkos::Experimental::ROCm::size_type
    ScratchGrain[Impl::ROCmTraits::WorkgroupSize];
enum { sizeScratchGrain = sizeof(ScratchGrain) };

void rocmMemset(Kokkos::Experimental::ROCm::size_type* ptr,
                Kokkos::Experimental::ROCm::size_type value,
                Kokkos::Experimental::ROCm::size_type size) {
  char* mptr = (char*)ptr;
  /*   parallel_for_each(hc::extent<1>(size),
                      [=, &ptr]
                      (hc::index<1> idx) __HC__
     {
        int i = idx[0];
        ptr[i] = value;
     }).wait();*/
}

Kokkos::Experimental::ROCm::size_type* ROCmInternal::scratch_flags(
    const Kokkos::Experimental::ROCm::size_type size) {
  if (verify_is_initialized("scratch_flags") &&
      m_scratchFlagsCount * sizeScratchGrain < size) {
    m_scratchFlagsCount = (size + sizeScratchGrain - 1) / sizeScratchGrain;

    typedef Kokkos::Impl::SharedAllocationRecord<
        Kokkos::Experimental::ROCmSpace, void>
        Record;

    Record* const r = Record::allocate(
        Kokkos::Experimental::ROCmSpace(), "InternalScratchFlags",
        (sizeScratchGrain * m_scratchFlagsCount));

    Record::increment(r);

    m_scratchFlags = reinterpret_cast<size_type*>(r->data());

    rocmMemset(m_scratchFlags, 0, m_scratchFlagsCount * sizeScratchGrain);
  }

  return m_scratchFlags;
}

Kokkos::Experimental::ROCm::size_type* ROCmInternal::scratch_space(
    const Kokkos::Experimental::ROCm::size_type size) {
  if (verify_is_initialized("scratch_space") &&
      m_scratchSpaceCount * sizeScratchGrain < size) {
    m_scratchSpaceCount = (size + sizeScratchGrain - 1) / sizeScratchGrain;

    typedef Kokkos::Impl::SharedAllocationRecord<
        Kokkos::Experimental::ROCmSpace, void>
        Record;

    static Record* const r = Record::allocate(
        Kokkos::Experimental::ROCmSpace(), "InternalScratchSpace",
        (sizeScratchGrain * m_scratchSpaceCount));

    Record::increment(r);

    m_scratchSpace = reinterpret_cast<size_type*>(r->data());
  }

  return m_scratchSpace;
}

//----------------------------------------------------------------------------

void ROCmInternal::finalize() {
  Kokkos::Impl::rocm_device_synchronize();
  was_finalized = 1;
  if (0 != m_scratchSpace || 0 != m_scratchFlags) {
    //    atomic_lock_array_rocm_space_ptr(false);
    //    scratch_lock_array_rocm_space_ptr(false);
    //    threadid_lock_array_rocm_space_ptr(false);

    typedef Kokkos::Impl::SharedAllocationRecord<
        Kokkos::Experimental::ROCmSpace>
        RecordROCm;
    typedef Kokkos::Impl::SharedAllocationRecord<
        Kokkos::Experimental::ROCmHostPinnedSpace>
        RecordHost;

    RecordROCm::decrement(RecordROCm::get_record(m_scratchFlags));
    RecordROCm::decrement(RecordROCm::get_record(m_scratchSpace));

    m_rocmDev           = -1;
    m_multiProcCount    = 0;
    m_maxWorkgroup      = 0;
    m_maxSharedWords    = 0;
    m_scratchSpaceCount = 0;
    m_scratchFlagsCount = 0;
    m_scratchSpace      = 0;
    m_scratchFlags      = 0;
  }
}

//----------------------------------------------------------------------------

Kokkos::Experimental::ROCm::size_type rocm_internal_cu_count() {
  return ROCmInternal::singleton().m_multiProcCount;
}

Kokkos::Experimental::ROCm::size_type rocm_internal_maximum_extent_size() {
  return ROCmInternal::singleton().m_maxWorkgroup;
}

Kokkos::Experimental::ROCm::size_type rocm_internal_maximum_shared_words() {
  return ROCmInternal::singleton().m_maxSharedWords;
}

Kokkos::Experimental::ROCm::size_type* rocm_internal_scratch_space(
    const Kokkos::Experimental::ROCm::size_type size) {
  return ROCmInternal::singleton().scratch_space(size);
}

Kokkos::Experimental::ROCm::size_type* rocm_internal_scratch_flags(
    const Kokkos::Experimental::ROCm::size_type size) {
  return ROCmInternal::singleton().scratch_flags(size);
}

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {

// ROCm::size_type ROCm::detect_device_count()
//{ return Impl::ROCmInternalDevices::singleton().m_rocmDevCount ; }

int ROCm::concurrency() {
#if defined(KOKKOS_ARCH_KAVERI)
  return 8 * 64 * 40;  // 20480 kaveri
#else
  return 32 * 8 * 40;  // 81920 fiji and hawaii
#endif
}
int ROCm::is_initialized() {
  return Kokkos::Impl::ROCmInternal::singleton().is_initialized();
}

void ROCm::initialize(const ROCm::SelectDevice config) {
  Kokkos::Impl::ROCmInternal::singleton().initialize(config.rocm_device_id);

#if defined(KOKKOS_ENABLE_PROFILING)
  Kokkos::Profiling::initialize();
#endif
}

#if 0
std::vector<unsigned>
ROCm::detect_device_arch()
{
  const Impl::ROCmInternalDevices & s = Impl::ROCmInternalDevices::singleton();

  std::vector<unsigned> output( s.m_rocmDevCount );

  for ( int i = 0 ; i < s.m_rocmDevCount ; ++i ) {
    output[i] = s.m_rocmProp[i].major * 100 + s.m_rocmProp[i].minor ;
  }

  return output ;
}

ROCm::size_type ROCm::device_arch()
{
  return 1 ;
}
#endif

void ROCm::finalize() {
  Kokkos::Impl::ROCmInternal::singleton().finalize();

#if defined(KOKKOS_ENABLE_PROFILING)
  Kokkos::Profiling::finalize();
#endif
}

ROCm::ROCm() : m_device(Kokkos::Impl::ROCmInternal::singleton().m_rocmDev) {
  Kokkos::Impl::ROCmInternal::singleton().verify_is_initialized(
      "ROCm instance constructor");
}

bool ROCm::isAPU(int device) {
  const Kokkos::Impl::ROCmInternalDevices& dev_info =
      Kokkos::Impl::ROCmInternalDevices::singleton();
  return (dev_info.m_rocmProp[device].APU);
}

bool ROCm::isAPU() { return ROCm::isAPU(rocm_device()); }

// ROCm::ROCm( const int instance_id )
//  : m_device( Impl::ROCmInternal::singleton().m_rocmDev )
//{}

void ROCm::print_configuration(std::ostream& s, const bool) {
  Kokkos::Impl::ROCmInternal::singleton().print_configuration(s);
}

bool ROCm::sleep() { return false; }

bool ROCm::wake() { return true; }

void ROCm::fence() { Kokkos::Impl::rocm_device_synchronize(); }

const char* ROCm::name() { return "ROCm"; }

}  // namespace Experimental
}  // namespace Kokkos

#endif  // KOKKOS_ENABLE_ROCM
//----------------------------------------------------------------------------
