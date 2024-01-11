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

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_CUDA

#include <Kokkos_Core.hpp>

#include <Cuda/Kokkos_Cuda_Error.hpp>
#include <Cuda/Kokkos_Cuda_BlockSize_Deduction.hpp>
#include <Cuda/Kokkos_Cuda_Instance.hpp>
#include <Cuda/Kokkos_Cuda_UniqueToken.hpp>
#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_Tools.hpp>
#include <impl/Kokkos_CheckedIntegerOps.hpp>
#include <impl/Kokkos_DeviceManagement.hpp>
#include <impl/Kokkos_ExecSpaceManager.hpp>

/*--------------------------------------------------------------------------*/
/* Standard 'C' libraries */
#include <cstdlib>

/* Standard 'C++' libraries */
#include <vector>
#include <iostream>
#include <sstream>
#include <string>

#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
namespace Kokkos {
namespace Impl {

bool CudaInternal::kokkos_impl_cuda_use_serial_execution_v = false;

void CudaInternal::cuda_set_serial_execution(bool val) {
  CudaInternal::kokkos_impl_cuda_use_serial_execution_v = val;
}
bool CudaInternal::cuda_use_serial_execution() {
  return CudaInternal::kokkos_impl_cuda_use_serial_execution_v;
}

}  // namespace Impl
}  // namespace Kokkos

void kokkos_impl_cuda_set_serial_execution(bool val) {
  Kokkos::Impl::CudaInternal::cuda_set_serial_execution(val);
}
bool kokkos_impl_cuda_use_serial_execution() {
  return Kokkos::Impl::CudaInternal::cuda_use_serial_execution();
}
#endif

#ifdef KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE

__device__ __constant__ unsigned long kokkos_impl_cuda_constant_memory_buffer
    [Kokkos::Impl::CudaTraits::ConstantMemoryUsage / sizeof(unsigned long)];

#endif

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

namespace {

__global__ void query_cuda_kernel_arch(int *d_arch) {
#ifdef _NVHPC_CUDA
  *d_arch = __builtin_current_device_sm() * 10;
#else
#if defined(__CUDA_ARCH__)
  *d_arch = __CUDA_ARCH__;
#else
  *d_arch = 0;
#endif
#endif
}

/** Query what compute capability is actually launched to the device: */
int cuda_kernel_arch() {
  int arch    = 0;
  int *d_arch = nullptr;

  KOKKOS_IMPL_CUDA_SAFE_CALL((CudaInternal::singleton().cuda_malloc_wrapper(
      reinterpret_cast<void **>(&d_arch), sizeof(int))));
  KOKKOS_IMPL_CUDA_SAFE_CALL((CudaInternal::singleton().cuda_memcpy_wrapper(
      d_arch, &arch, sizeof(int), cudaMemcpyDefault)));

  query_cuda_kernel_arch<<<1, 1>>>(d_arch);

  KOKKOS_IMPL_CUDA_SAFE_CALL((CudaInternal::singleton().cuda_memcpy_wrapper(
      &arch, d_arch, sizeof(int), cudaMemcpyDefault)));
  KOKKOS_IMPL_CUDA_SAFE_CALL(
      (CudaInternal::singleton().cuda_free_wrapper(d_arch)));
  return arch;
}

constexpr auto sizeScratchGrain =
    sizeof(Cuda::size_type[Impl::CudaTraits::WarpSize]);

std::size_t scratch_count(const std::size_t size) {
  return (size + sizeScratchGrain - 1) / sizeScratchGrain;
}

}  // namespace

Kokkos::View<uint32_t *, Kokkos::CudaSpace> cuda_global_unique_token_locks(
    bool deallocate) {
  static Kokkos::View<uint32_t *, Kokkos::CudaSpace> locks =
      Kokkos::View<uint32_t *, Kokkos::CudaSpace>();
  if (!deallocate && locks.extent(0) == 0)
    locks = Kokkos::View<uint32_t *, Kokkos::CudaSpace>(
        "Kokkos::UniqueToken<Cuda>::m_locks", Kokkos::Cuda().concurrency());
  if (deallocate) locks = Kokkos::View<uint32_t *, Kokkos::CudaSpace>();
  return locks;
}

// FIXME_CUDA_MULTIPLE_DEVICES
void cuda_device_synchronize(const std::string &name) {
  Kokkos::Tools::Experimental::Impl::profile_fence_event<Kokkos::Cuda>(
      name,
      Kokkos::Tools::Experimental::SpecialSynchronizationCases::
          GlobalDeviceSynchronization,
#if defined(KOKKOS_COMPILER_CLANG)
      // annotate with __host__ silence a clang warning about using
      // cudaDeviceSynchronize in device code
      [] __host__() {
        KOKKOS_IMPL_CUDA_SAFE_CALL(
            (CudaInternal::singleton().cuda_device_synchronize_wrapper()));
      });
#else
      []() {
        KOKKOS_IMPL_CUDA_SAFE_CALL(
            (CudaInternal::singleton().cuda_device_synchronize_wrapper()));
      });
#endif
}

void cuda_stream_synchronize(const cudaStream_t stream, const CudaInternal *ptr,
                             const std::string &name) {
  Kokkos::Tools::Experimental::Impl::profile_fence_event<Kokkos::Cuda>(
      name,
      Kokkos::Tools::Experimental::Impl::DirectFenceIDHandle{
          ptr->impl_get_instance_id()},
      [&]() {
        KOKKOS_IMPL_CUDA_SAFE_CALL(
            (ptr->cuda_stream_synchronize_wrapper(stream)));
      });
}

void cuda_stream_synchronize(
    const cudaStream_t stream,
    Kokkos::Tools::Experimental::SpecialSynchronizationCases reason,
    const std::string &name) {
  Kokkos::Tools::Experimental::Impl::profile_fence_event<Kokkos::Cuda>(
      name, reason, [&]() {
        KOKKOS_IMPL_CUDA_SAFE_CALL(
            (Impl::CudaInternal::singleton().cuda_stream_synchronize_wrapper(
                stream)));
      });
}

void cuda_internal_error_throw(cudaError e, const char *name, const char *file,
                               const int line) {
  std::ostringstream out;
  out << name << " error( "
      << CudaInternal::singleton().cuda_get_error_name_wrapper<false>(e)
      << "): "
      << CudaInternal::singleton().cuda_get_error_string_wrapper<false>(e);
  if (file) {
    out << " " << file << ":" << line;
  }
  throw_runtime_exception(out.str());
}

void cuda_internal_error_abort(cudaError e, const char *name, const char *file,
                               const int line) {
  std::ostringstream out;
  out << name << " error( "
      << CudaInternal::singleton().cuda_get_error_name_wrapper<false>(e)
      << "): "
      << CudaInternal::singleton().cuda_get_error_string_wrapper<false>(e);
  if (file) {
    out << " " << file << ":" << line;
  }
  // FIXME Call Kokkos::Impl::host_abort instead of Kokkos::abort to avoid a
  // warning about Kokkos::abort returning in some cases.
  host_abort(out.str().c_str());
}

//----------------------------------------------------------------------------
// Some significant cuda device properties:
//
// cudaDeviceProp::name                : Text label for device
// cudaDeviceProp::major               : Device major number
// cudaDeviceProp::minor               : Device minor number
// cudaDeviceProp::warpSize            : number of threads per warp
// cudaDeviceProp::multiProcessorCount : number of multiprocessors
// cudaDeviceProp::sharedMemPerBlock   : capacity of shared memory per block
// cudaDeviceProp::totalConstMem       : capacity of constant memory
// cudaDeviceProp::totalGlobalMem      : capacity of global memory
// cudaDeviceProp::maxGridSize[3]      : maximum grid size

//
//  Section 4.4.2.4 of the CUDA Toolkit Reference Manual
//
// struct cudaDeviceProp {
//   char name[256];
//   size_t totalGlobalMem;
//   size_t sharedMemPerBlock;
//   int regsPerBlock;
//   int warpSize;
//   size_t memPitch;
//   int maxThreadsPerBlock;
//   int maxThreadsDim[3];
//   int maxGridSize[3];
//   size_t totalConstMem;
//   int major;
//   int minor;
//   int clockRate;
//   size_t textureAlignment;
//   int deviceOverlap;
//   int multiProcessorCount;
//   int kernelExecTimeoutEnabled;
//   int integrated;
//   int canMapHostMemory;
//   int computeMode;
//   int concurrentKernels;
//   int ECCEnabled;
//   int pciBusID;
//   int pciDeviceID;
//   int tccDriver;
//   int asyncEngineCount;
//   int unifiedAddressing;
//   int memoryClockRate;
//   int memoryBusWidth;
//   int l2CacheSize;
//   int maxThreadsPerMultiProcessor;
// };

namespace {

class CudaInternalDevices {
 public:
  enum { MAXIMUM_DEVICE_COUNT = 64 };
  struct cudaDeviceProp m_cudaProp[MAXIMUM_DEVICE_COUNT];
  int m_cudaDevCount;

  CudaInternalDevices();

  static const CudaInternalDevices &singleton();
};

CudaInternalDevices::CudaInternalDevices() {
  // See 'cudaSetDeviceFlags' for host-device thread interaction
  // Section 4.4.2.6 of the CUDA Toolkit Reference Manual

  KOKKOS_IMPL_CUDA_SAFE_CALL(
      (CudaInternal::singleton().cuda_get_device_count_wrapper<false>(
          &m_cudaDevCount)));

  if (m_cudaDevCount > MAXIMUM_DEVICE_COUNT) {
    Kokkos::abort(
        "Sorry, you have more GPUs per node than we thought anybody would ever "
        "have. Please report this to github.com/kokkos/kokkos.");
  }
  for (int i = 0; i < m_cudaDevCount; ++i) {
    KOKKOS_IMPL_CUDA_SAFE_CALL(
        (CudaInternal::singleton().cuda_get_device_properties_wrapper<false>(
            m_cudaProp + i, i)));
  }
}

const CudaInternalDevices &CudaInternalDevices::singleton() {
  static CudaInternalDevices self;
  return self;
}

}  // namespace

//----------------------------------------------------------------------------

int Impl::CudaInternal::concurrency() {
  static int const concurrency = m_deviceProp.maxThreadsPerMultiProcessor *
                                 m_deviceProp.multiProcessorCount;
  return concurrency;
}

void CudaInternal::print_configuration(std::ostream &s) const {
  const CudaInternalDevices &dev_info = CudaInternalDevices::singleton();

#if defined(KOKKOS_ENABLE_CUDA)
  s << "macro  KOKKOS_ENABLE_CUDA      : defined\n";
#endif
#if defined(CUDA_VERSION)
  s << "macro  CUDA_VERSION          = " << CUDA_VERSION << " = version "
    << CUDA_VERSION / 1000 << "." << (CUDA_VERSION % 1000) / 10 << '\n';
#endif

  for (int i = 0; i < dev_info.m_cudaDevCount; ++i) {
    s << "Kokkos::Cuda[ " << i << " ] " << dev_info.m_cudaProp[i].name
      << " capability " << dev_info.m_cudaProp[i].major << "."
      << dev_info.m_cudaProp[i].minor << ", Total Global Memory: "
      << human_memory_size(dev_info.m_cudaProp[i].totalGlobalMem)
      << ", Shared Memory per Block: "
      << human_memory_size(dev_info.m_cudaProp[i].sharedMemPerBlock);
    if (m_cudaDev == i) s << " : Selected";
    s << std::endl;
  }
}

//----------------------------------------------------------------------------

CudaInternal::~CudaInternal() {
  if (m_stream || m_scratchSpace || m_scratchFlags || m_scratchUnified) {
    std::cerr << "Kokkos::Cuda ERROR: Failed to call Kokkos::Cuda::finalize()"
              << std::endl;
  }

  m_scratchSpaceCount   = 0;
  m_scratchFlagsCount   = 0;
  m_scratchUnifiedCount = 0;
  m_scratchSpace        = nullptr;
  m_scratchFlags        = nullptr;
  m_scratchUnified      = nullptr;
  m_stream              = nullptr;
  for (int i = 0; i < m_n_team_scratch; ++i) {
    m_team_scratch_current_size[i] = 0;
    m_team_scratch_ptr[i]          = nullptr;
  }
}

int CudaInternal::verify_is_initialized(const char *const label) const {
  if (m_cudaDev < 0) {
    Kokkos::abort((std::string("Kokkos::Cuda::") + label +
                   " : ERROR device not initialized\n")
                      .c_str());
  }
  return 0 <= m_cudaDev;
}
uint32_t CudaInternal::impl_get_instance_id() const { return m_instance_id; }
CudaInternal &CudaInternal::singleton() {
  static CudaInternal self;
  return self;
}
void CudaInternal::fence(const std::string &name) const {
  Impl::cuda_stream_synchronize(get_stream(), this, name);
}
void CudaInternal::fence() const {
  fence("Kokkos::CudaInternal::fence(): Unnamed Instance Fence");
}

void CudaInternal::initialize(cudaStream_t stream, bool manage_stream) {
  KOKKOS_EXPECTS(!is_initialized());

  if (was_finalized)
    Kokkos::abort("Calling Cuda::initialize after Cuda::finalize is illegal\n");
  was_initialized = true;

  //----------------------------------
  // Multiblock reduction uses scratch flags for counters
  // and scratch space for partial reduction values.
  // Allocate some initial space.  This will grow as needed.

  {
    const unsigned reduce_block_count =
        m_maxWarpCount * Impl::CudaTraits::WarpSize;

    (void)scratch_unified(16 * sizeof(size_type));
    (void)scratch_flags(reduce_block_count * 2 * sizeof(size_type));
    (void)scratch_space(reduce_block_count * 16 * sizeof(size_type));
  }

  // Init the array for used for arbitrarily sized atomics
  if (this == &singleton()) {
    desul::Impl::init_lock_arrays();  // FIXME
  }

  // Allocate a staging buffer for constant mem in pinned host memory
  // and an event to avoid overwriting driver for previous kernel launches
  if (this == &singleton()) {
    KOKKOS_IMPL_CUDA_SAFE_CALL((cuda_malloc_host_wrapper(
        reinterpret_cast<void **>(&constantMemHostStaging),
        CudaTraits::ConstantMemoryUsage)));

    KOKKOS_IMPL_CUDA_SAFE_CALL(
        (cuda_event_create_wrapper(&constantMemReusable)));
  }

  m_stream        = stream;
  m_manage_stream = manage_stream;
  for (int i = 0; i < m_n_team_scratch; ++i) {
    m_team_scratch_current_size[i] = 0;
    m_team_scratch_ptr[i]          = nullptr;
  }

  m_num_scratch_locks = concurrency();
  KOKKOS_IMPL_CUDA_SAFE_CALL(
      (cuda_malloc_wrapper(reinterpret_cast<void **>(&m_scratch_locks),
                           sizeof(int32_t) * m_num_scratch_locks)));
  KOKKOS_IMPL_CUDA_SAFE_CALL((cuda_memset_wrapper(
      m_scratch_locks, 0, sizeof(int32_t) * m_num_scratch_locks)));
}

//----------------------------------------------------------------------------

Cuda::size_type *CudaInternal::scratch_flags(const std::size_t size) const {
  if (verify_is_initialized("scratch_flags") &&
      m_scratchFlagsCount < scratch_count(size)) {
    m_scratchFlagsCount = scratch_count(size);

    using Record =
        Kokkos::Impl::SharedAllocationRecord<Kokkos::CudaSpace, void>;

    if (m_scratchFlags) Record::decrement(Record::get_record(m_scratchFlags));

    std::size_t alloc_size =
        multiply_overflow_abort(m_scratchFlagsCount, sizeScratchGrain);
    Record *const r = Record::allocate(
        Kokkos::CudaSpace(), "Kokkos::InternalScratchFlags", alloc_size);

    Record::increment(r);

    m_scratchFlags = reinterpret_cast<size_type *>(r->data());

    KOKKOS_IMPL_CUDA_SAFE_CALL(
        (cuda_memset_wrapper(m_scratchFlags, 0, alloc_size)));
  }

  return m_scratchFlags;
}

Cuda::size_type *CudaInternal::scratch_space(const std::size_t size) const {
  if (verify_is_initialized("scratch_space") &&
      m_scratchSpaceCount < scratch_count(size)) {
    m_scratchSpaceCount = scratch_count(size);

    using Record =
        Kokkos::Impl::SharedAllocationRecord<Kokkos::CudaSpace, void>;

    if (m_scratchSpace) Record::decrement(Record::get_record(m_scratchSpace));

    std::size_t alloc_size =
        multiply_overflow_abort(m_scratchSpaceCount, sizeScratchGrain);
    Record *const r = Record::allocate(
        Kokkos::CudaSpace(), "Kokkos::InternalScratchSpace", alloc_size);

    Record::increment(r);

    m_scratchSpace = reinterpret_cast<size_type *>(r->data());
  }

  return m_scratchSpace;
}

Cuda::size_type *CudaInternal::scratch_unified(const std::size_t size) const {
  if (verify_is_initialized("scratch_unified") &&
      m_scratchUnifiedCount < scratch_count(size)) {
    m_scratchUnifiedCount = scratch_count(size);

    using Record =
        Kokkos::Impl::SharedAllocationRecord<Kokkos::CudaHostPinnedSpace, void>;

    if (m_scratchUnified)
      Record::decrement(Record::get_record(m_scratchUnified));

    std::size_t alloc_size =
        multiply_overflow_abort(m_scratchUnifiedCount, sizeScratchGrain);
    Record *const r =
        Record::allocate(Kokkos::CudaHostPinnedSpace(),
                         "Kokkos::InternalScratchUnified", alloc_size);

    Record::increment(r);

    m_scratchUnified = reinterpret_cast<size_type *>(r->data());
  }

  return m_scratchUnified;
}

Cuda::size_type *CudaInternal::scratch_functor(const std::size_t size) const {
  if (verify_is_initialized("scratch_functor") && m_scratchFunctorSize < size) {
    m_scratchFunctorSize = size;

    using Record =
        Kokkos::Impl::SharedAllocationRecord<Kokkos::CudaSpace, void>;

    if (m_scratchFunctor)
      Record::decrement(Record::get_record(m_scratchFunctor));

    Record *const r =
        Record::allocate(Kokkos::CudaSpace(), "Kokkos::InternalScratchFunctor",
                         m_scratchFunctorSize);

    Record::increment(r);

    m_scratchFunctor = reinterpret_cast<size_type *>(r->data());
  }

  return m_scratchFunctor;
}

int CudaInternal::acquire_team_scratch_space() {
  int current_team_scratch = 0;
  int zero                 = 0;
  while (!m_team_scratch_pool[current_team_scratch].compare_exchange_weak(
      zero, 1, std::memory_order_release, std::memory_order_relaxed)) {
    current_team_scratch = (current_team_scratch + 1) % m_n_team_scratch;
  }

  return current_team_scratch;
}

void *CudaInternal::resize_team_scratch_space(int scratch_pool_id,
                                              std::int64_t bytes,
                                              bool force_shrink) {
  // Multiple ParallelFor/Reduce Teams can call this function at the same time
  // and invalidate the m_team_scratch_ptr. We use a pool to avoid any race
  // condition.
  if (m_team_scratch_current_size[scratch_pool_id] == 0) {
    m_team_scratch_current_size[scratch_pool_id] = bytes;
    m_team_scratch_ptr[scratch_pool_id] =
        Kokkos::kokkos_malloc<Kokkos::CudaSpace>(
            "Kokkos::CudaSpace::TeamScratchMemory",
            m_team_scratch_current_size[scratch_pool_id]);
  }
  if ((bytes > m_team_scratch_current_size[scratch_pool_id]) ||
      ((bytes < m_team_scratch_current_size[scratch_pool_id]) &&
       (force_shrink))) {
    m_team_scratch_current_size[scratch_pool_id] = bytes;
    m_team_scratch_ptr[scratch_pool_id] =
        Kokkos::kokkos_realloc<Kokkos::CudaSpace>(
            m_team_scratch_ptr[scratch_pool_id],
            m_team_scratch_current_size[scratch_pool_id]);
  }
  return m_team_scratch_ptr[scratch_pool_id];
}

void CudaInternal::release_team_scratch_space(int scratch_pool_id) {
  m_team_scratch_pool[scratch_pool_id] = 0;
}

//----------------------------------------------------------------------------

void CudaInternal::finalize() {
  // skip if finalize() has already been called
  if (was_finalized) return;

  was_finalized = true;

  // Only finalize this if we're the singleton
  if (this == &singleton()) {
    (void)Impl::cuda_global_unique_token_locks(true);
    desul::Impl::finalize_lock_arrays();  // FIXME

    KOKKOS_IMPL_CUDA_SAFE_CALL(
        (cuda_free_host_wrapper(constantMemHostStaging)));
    KOKKOS_IMPL_CUDA_SAFE_CALL(
        (cuda_event_destroy_wrapper(constantMemReusable)));
    auto &deep_copy_space =
        Kokkos::Impl::cuda_get_deep_copy_space(/*initialize*/ false);
    if (deep_copy_space)
      deep_copy_space->impl_internal_space_instance()->finalize();
    KOKKOS_IMPL_CUDA_SAFE_CALL(
        (cuda_stream_destroy_wrapper(cuda_get_deep_copy_stream())));
  }

  if (nullptr != m_scratchSpace || nullptr != m_scratchFlags) {
    using RecordCuda = Kokkos::Impl::SharedAllocationRecord<CudaSpace>;
    using RecordHost =
        Kokkos::Impl::SharedAllocationRecord<CudaHostPinnedSpace>;

    RecordCuda::decrement(RecordCuda::get_record(m_scratchFlags));
    RecordCuda::decrement(RecordCuda::get_record(m_scratchSpace));
    RecordHost::decrement(RecordHost::get_record(m_scratchUnified));
    if (m_scratchFunctorSize > 0)
      RecordCuda::decrement(RecordCuda::get_record(m_scratchFunctor));
  }

  for (int i = 0; i < m_n_team_scratch; ++i) {
    if (m_team_scratch_current_size[i] > 0)
      Kokkos::kokkos_free<Kokkos::CudaSpace>(m_team_scratch_ptr[i]);
  }

  if (m_manage_stream && get_stream() != nullptr)
    KOKKOS_IMPL_CUDA_SAFE_CALL((cuda_stream_destroy_wrapper(m_stream)));

  m_scratchSpaceCount   = 0;
  m_scratchFlagsCount   = 0;
  m_scratchUnifiedCount = 0;
  m_scratchSpace        = nullptr;
  m_scratchFlags        = nullptr;
  m_scratchUnified      = nullptr;
  m_stream              = nullptr;
  for (int i = 0; i < m_n_team_scratch; ++i) {
    m_team_scratch_current_size[i] = 0;
    m_team_scratch_ptr[i]          = nullptr;
  }

  KOKKOS_IMPL_CUDA_SAFE_CALL((cuda_free_wrapper(m_scratch_locks)));
  m_scratch_locks     = nullptr;
  m_num_scratch_locks = 0;
}

//----------------------------------------------------------------------------

Cuda::size_type cuda_internal_multiprocessor_count() {
  return CudaInternal::singleton().m_multiProcCount;
}

CudaSpace::size_type cuda_internal_maximum_concurrent_block_count() {
#if defined(KOKKOS_ARCH_KEPLER)
  // Compute capability 3.0 through 3.7
  enum : int { max_resident_blocks_per_multiprocessor = 16 };
#else
  // Compute capability 5.0 through 6.2
  enum : int { max_resident_blocks_per_multiprocessor = 32 };
#endif
  return CudaInternal::singleton().m_multiProcCount *
         max_resident_blocks_per_multiprocessor;
};

Cuda::size_type cuda_internal_maximum_warp_count() {
  return CudaInternal::singleton().m_maxWarpCount;
}

std::array<Cuda::size_type, 3> cuda_internal_maximum_grid_count() {
  return CudaInternal::singleton().m_maxBlock;
}

Cuda::size_type *cuda_internal_scratch_space(const Cuda &instance,
                                             const std::size_t size) {
  return instance.impl_internal_space_instance()->scratch_space(size);
}

Cuda::size_type *cuda_internal_scratch_flags(const Cuda &instance,
                                             const std::size_t size) {
  return instance.impl_internal_space_instance()->scratch_flags(size);
}

Cuda::size_type *cuda_internal_scratch_unified(const Cuda &instance,
                                               const std::size_t size) {
  return instance.impl_internal_space_instance()->scratch_unified(size);
}

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

Cuda::size_type Cuda::detect_device_count() {
  return Impl::CudaInternalDevices::singleton().m_cudaDevCount;
}

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
int Cuda::concurrency() {
#else
int Cuda::concurrency() const {
#endif
  return Impl::CudaInternal::concurrency();
}

int Cuda::impl_is_initialized() {
  return Impl::CudaInternal::singleton().is_initialized();
}

void Cuda::impl_initialize(InitializationSettings const &settings) {
  const int cuda_device_id = Impl::get_gpu(settings);
  const auto &dev_info     = Impl::CudaInternalDevices::singleton();

  const struct cudaDeviceProp &cudaProp = dev_info.m_cudaProp[cuda_device_id];

  Impl::CudaInternal::m_cudaDev    = cuda_device_id;
  Impl::CudaInternal::m_deviceProp = cudaProp;

  Kokkos::Impl::cuda_device_synchronize(
      "Kokkos::CudaInternal::initialize: Fence on space initialization");

  // Query what compute capability architecture a kernel executes:
  Impl::CudaInternal::m_cudaArch = Impl::cuda_kernel_arch();

  if (Impl::CudaInternal::m_cudaArch == 0) {
    std::stringstream ss;
    ss << "Kokkos::Cuda::initialize ERROR: likely mismatch of architecture\n";
    std::string msg = ss.str();
    Kokkos::abort(msg.c_str());
  }

  int compiled_major = Impl::CudaInternal::m_cudaArch / 100;
  int compiled_minor = (Impl::CudaInternal::m_cudaArch % 100) / 10;

  if ((compiled_major > cudaProp.major) ||
      ((compiled_major == cudaProp.major) &&
       (compiled_minor > cudaProp.minor))) {
    std::stringstream ss;
    ss << "Kokkos::Cuda::initialize ERROR: running kernels compiled for "
          "compute capability "
       << compiled_major << "." << compiled_minor
       << " on device with compute capability " << cudaProp.major << "."
       << cudaProp.minor << " is not supported by CUDA!\n";
    std::string msg = ss.str();
    Kokkos::abort(msg.c_str());
  }
  if (Kokkos::show_warnings() &&
      (compiled_major != cudaProp.major || compiled_minor != cudaProp.minor)) {
    std::cerr << "Kokkos::Cuda::initialize WARNING: running kernels compiled "
                 "for compute capability "
              << compiled_major << "." << compiled_minor
              << " on device with compute capability " << cudaProp.major << "."
              << cudaProp.minor
              << " , this will likely reduce potential performance."
              << std::endl;
  }

  //----------------------------------

#ifdef KOKKOS_ENABLE_CUDA_UVM
  const char *env_force_device_alloc =
      getenv("CUDA_MANAGED_FORCE_DEVICE_ALLOC");
  bool force_device_alloc;
  if (env_force_device_alloc == nullptr)
    force_device_alloc = false;
  else
    force_device_alloc = std::stoi(env_force_device_alloc) != 0;

  const char *env_visible_devices = getenv("CUDA_VISIBLE_DEVICES");
  bool visible_devices_one        = true;
  if (env_visible_devices == nullptr) visible_devices_one = false;

  if (Kokkos::show_warnings() &&
      (!visible_devices_one && !force_device_alloc)) {
    std::cerr << R"warning(
Kokkos::Cuda::initialize WARNING: Cuda is allocating into UVMSpace by default
                                  without setting CUDA_MANAGED_FORCE_DEVICE_ALLOC=1 or
                                  setting CUDA_VISIBLE_DEVICES.
                                  This could on multi GPU systems lead to severe performance"
                                  penalties.)warning"
              << std::endl;
  }
#endif

  //----------------------------------
  // number of multiprocessors
  Impl::CudaInternal::m_multiProcCount = cudaProp.multiProcessorCount;

  //----------------------------------
  // Maximum number of warps,
  // at most one warp per thread in a warp for reduction.
  Impl::CudaInternal::m_maxWarpCount =
      cudaProp.maxThreadsPerBlock / Impl::CudaTraits::WarpSize;

  if (Impl::CudaTraits::WarpSize < Impl::CudaInternal::m_maxWarpCount) {
    Impl::CudaInternal::m_maxWarpCount = Impl::CudaTraits::WarpSize;
  }

  //----------------------------------
  // Maximum number of blocks:

  Impl::CudaInternal::m_maxBlock[0] = cudaProp.maxGridSize[0];
  Impl::CudaInternal::m_maxBlock[1] = cudaProp.maxGridSize[1];
  Impl::CudaInternal::m_maxBlock[2] = cudaProp.maxGridSize[2];

  Impl::CudaInternal::m_shmemPerSM       = cudaProp.sharedMemPerMultiprocessor;
  Impl::CudaInternal::m_maxShmemPerBlock = cudaProp.sharedMemPerBlock;
  Impl::CudaInternal::m_maxBlocksPerSM =
      Impl::CudaInternal::m_cudaArch < 500
          ? 16
          : (Impl::CudaInternal::m_cudaArch < 750
                 ? 32
                 : (Impl::CudaInternal::m_cudaArch == 750 ? 16 : 32));
  Impl::CudaInternal::m_maxThreadsPerSM = cudaProp.maxThreadsPerMultiProcessor;
  Impl::CudaInternal::m_maxThreadsPerBlock = cudaProp.maxThreadsPerBlock;

  //----------------------------------

  cudaStream_t singleton_stream;
  KOKKOS_IMPL_CUDA_SAFE_CALL(
      (Impl::CudaInternal::singleton().cuda_stream_create_wrapper(
          &singleton_stream)));

  auto &cuda_singleton = Impl::CudaInternal::singleton();
  cuda_singleton.initialize(singleton_stream, /*manage*/ true);
}

std::vector<unsigned> Cuda::detect_device_arch() {
  const Impl::CudaInternalDevices &s = Impl::CudaInternalDevices::singleton();

  std::vector<unsigned> output(s.m_cudaDevCount);

  for (int i = 0; i < s.m_cudaDevCount; ++i) {
    output[i] = s.m_cudaProp[i].major * 100 + s.m_cudaProp[i].minor;
  }

  return output;
}

Cuda::size_type Cuda::device_arch() {
  const int dev_id = Impl::CudaInternal::singleton().m_cudaDev;

  int dev_arch = 0;

  if (0 <= dev_id) {
    const struct cudaDeviceProp &cudaProp =
        Impl::CudaInternalDevices::singleton().m_cudaProp[dev_id];

    dev_arch = cudaProp.major * 100 + cudaProp.minor;
  }

  return dev_arch;
}

void Cuda::impl_finalize() { Impl::CudaInternal::singleton().finalize(); }

Cuda::Cuda()
    : m_space_instance(&Impl::CudaInternal::singleton(),
                       [](Impl::CudaInternal *) {}) {
  Impl::CudaInternal::singleton().verify_is_initialized(
      "Cuda instance constructor");
}

KOKKOS_DEPRECATED Cuda::Cuda(cudaStream_t stream, bool manage_stream)
    : Cuda(stream,
           manage_stream ? Impl::ManageStream::yes : Impl::ManageStream::no) {}

Cuda::Cuda(cudaStream_t stream, Impl::ManageStream manage_stream)
    : m_space_instance(new Impl::CudaInternal, [](Impl::CudaInternal *ptr) {
        ptr->finalize();
        delete ptr;
      }) {
  Impl::CudaInternal::singleton().verify_is_initialized(
      "Cuda instance constructor");
  m_space_instance->initialize(stream, static_cast<bool>(manage_stream));
}

void Cuda::print_configuration(std::ostream &os, bool /*verbose*/) const {
  os << "Device Execution Space:\n";
  os << "  KOKKOS_ENABLE_CUDA: yes\n";

  os << "Cuda Options:\n";
  os << "  KOKKOS_ENABLE_CUDA_LAMBDA: ";
#ifdef KOKKOS_ENABLE_CUDA_LAMBDA
  os << "yes\n";
#else
  os << "no\n";
#endif
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  os << "  KOKKOS_ENABLE_CUDA_LDG_INTRINSIC: ";
  os << "yes\n";
#endif
  os << "  KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE: ";
#ifdef KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE
  os << "yes\n";
#else
  os << "no\n";
#endif
  os << "  KOKKOS_ENABLE_CUDA_UVM: ";
#ifdef KOKKOS_ENABLE_CUDA_UVM
  os << "yes\n";
#else
  os << "no\n";
#endif
  os << "  KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA: ";
#ifdef KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA
  os << "yes\n";
#else
  os << "no\n";
#endif
  os << "  KOKKOS_ENABLE_IMPL_CUDA_MALLOC_ASYNC: ";
#ifdef KOKKOS_ENABLE_IMPL_CUDA_MALLOC_ASYNC
  os << "yes\n";
#else
  os << "no\n";
#endif

  os << "\nCuda Runtime Configuration:\n";

  m_space_instance->print_configuration(os);
}

void Cuda::impl_static_fence(const std::string &name) {
  Kokkos::Impl::cuda_device_synchronize(name);
}

void Cuda::fence(const std::string &name) const {
  m_space_instance->fence(name);
}

const char *Cuda::name() { return "Cuda"; }
uint32_t Cuda::impl_instance_id() const noexcept {
  return m_space_instance->impl_get_instance_id();
}

cudaStream_t Cuda::cuda_stream() const {
  return m_space_instance->get_stream();
}
int Cuda::cuda_device() const { return m_space_instance->m_cudaDev; }
const cudaDeviceProp &Cuda::cuda_device_prop() const {
  return m_space_instance->m_deviceProp;
}

namespace Impl {

int g_cuda_space_factory_initialized =
    initialize_space_factory<Cuda>("150_Cuda");

}  // namespace Impl

}  // namespace Kokkos

void Kokkos::Impl::create_Cuda_instances(std::vector<Cuda> &instances) {
  for (int s = 0; s < int(instances.size()); s++) {
    cudaStream_t stream;
    KOKKOS_IMPL_CUDA_SAFE_CALL((
        instances[s].impl_internal_space_instance()->cuda_stream_create_wrapper(
            &stream)));
    instances[s] = Cuda(stream, ManageStream::yes);
  }
}

#else

void KOKKOS_CORE_SRC_CUDA_IMPL_PREVENT_LINK_ERROR() {}

#endif  // KOKKOS_ENABLE_CUDA
