#ifndef KOKKOS_CUDA_INSTANCE_HPP_
#define KOKKOS_CUDA_INSTANCE_HPP_

#include <vector>
#include <impl/Kokkos_Tools.hpp>
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
  enum : CudaSpace::size_type { WarpSize = 32 /* 0x0020 */ };
  enum : CudaSpace::size_type {
    WarpIndexMask = 0x001f /* Mask for warpindex */
  };
  enum : CudaSpace::size_type {
    WarpIndexShift = 5 /* WarpSize == 1 << WarpShift */
  };

  enum : CudaSpace::size_type {
    ConstantMemoryUsage = 0x008000 /* 32k bytes */
  };
  enum : CudaSpace::size_type {
    ConstantMemoryCache = 0x002000 /*  8k bytes */
  };
  enum : CudaSpace::size_type {
    KernelArgumentLimit = 0x001000 /*  4k bytes */
  };
  enum : CudaSpace::size_type {
    MaxHierarchicalParallelism = 1024 /* team_size * vector_length */
  };
  using ConstantGlobalBufferType =
      unsigned long[ConstantMemoryUsage / sizeof(unsigned long)];

  enum { ConstantMemoryUseThreshold = 0x000200 /* 512 bytes */ };

  KOKKOS_INLINE_FUNCTION static CudaSpace::size_type warp_count(
      CudaSpace::size_type i) {
    return (i + WarpIndexMask) >> WarpIndexShift;
  }

  KOKKOS_INLINE_FUNCTION static CudaSpace::size_type warp_align(
      CudaSpace::size_type i) {
    constexpr CudaSpace::size_type Mask = ~WarpIndexMask;
    return (i + WarpIndexMask) & Mask;
  }
};

//----------------------------------------------------------------------------

CudaSpace::size_type cuda_internal_multiprocessor_count();
CudaSpace::size_type cuda_internal_maximum_warp_count();
CudaSpace::size_type cuda_internal_maximum_grid_count();
CudaSpace::size_type cuda_internal_maximum_shared_words();

CudaSpace::size_type cuda_internal_maximum_concurrent_block_count();

CudaSpace::size_type* cuda_internal_scratch_flags(
    const Cuda&, const CudaSpace::size_type size);
CudaSpace::size_type* cuda_internal_scratch_space(
    const Cuda&, const CudaSpace::size_type size);
CudaSpace::size_type* cuda_internal_scratch_unified(
    const Cuda&, const CudaSpace::size_type size);

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

  int m_cudaDev;

  // Device Properties
  int m_cudaArch;
  unsigned m_multiProcCount;
  unsigned m_maxWarpCount;
  unsigned m_maxBlock;
  unsigned m_maxSharedWords;
  uint32_t m_maxConcurrency;
  int m_shmemPerSM;
  int m_maxShmemPerBlock;
  int m_regsPerSM;
  int m_maxBlocksPerSM;
  int m_maxThreadsPerSM;
  int m_maxThreadsPerBlock;

  cudaDeviceProp m_deviceProp;

  // Scratch Spaces for Reductions
  mutable size_type m_scratchSpaceCount;
  mutable size_type m_scratchFlagsCount;
  mutable size_type m_scratchUnifiedCount;
  mutable size_type m_scratchFunctorSize;

  size_type m_scratchUnifiedSupported;
  size_type m_streamCount;
  mutable size_type* m_scratchSpace;
  mutable size_type* m_scratchFlags;
  mutable size_type* m_scratchUnified;
  mutable size_type* m_scratchFunctor;
  uint32_t* m_scratchConcurrentBitset;
  cudaStream_t m_stream;

  // Team Scratch Level 1 Space
  mutable int64_t m_team_scratch_current_size;
  mutable void* m_team_scratch_ptr;

  bool was_initialized = false;
  bool was_finalized   = false;

  // FIXME_CUDA: these want to be per-device, not per-stream...  use of 'static'
  //  here will break once there are multiple devices though
  static unsigned long* constantMemHostStaging;
  static cudaEvent_t constantMemReusable;

  static CudaInternal& singleton();

  int verify_is_initialized(const char* const label) const;

  int is_initialized() const {
    return nullptr != m_scratchSpace && nullptr != m_scratchFlags;
  }

  void initialize(int cuda_device_id, cudaStream_t stream = nullptr);
  void finalize();

  void print_configuration(std::ostream&) const;

#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
  static bool cuda_use_serial_execution();
  static void cuda_set_serial_execution(bool);
#endif

  void fence() const;

  ~CudaInternal();

  CudaInternal()
      : m_cudaDev(-1),
        m_cudaArch(-1),
        m_multiProcCount(0),
        m_maxWarpCount(0),
        m_maxBlock(0),
        m_maxSharedWords(0),
        m_maxConcurrency(0),
        m_shmemPerSM(0),
        m_maxShmemPerBlock(0),
        m_regsPerSM(0),
        m_maxBlocksPerSM(0),
        m_maxThreadsPerSM(0),
        m_maxThreadsPerBlock(0),
        m_scratchSpaceCount(0),
        m_scratchFlagsCount(0),
        m_scratchUnifiedCount(0),
        m_scratchFunctorSize(0),
        m_scratchUnifiedSupported(0),
        m_streamCount(0),
        m_scratchSpace(nullptr),
        m_scratchFlags(nullptr),
        m_scratchUnified(nullptr),
        m_scratchFunctor(nullptr),
        m_scratchConcurrentBitset(nullptr),
        m_stream(nullptr),
        m_team_scratch_current_size(0),
        m_team_scratch_ptr(nullptr) {}

  // Resizing of reduction related scratch spaces
  size_type* scratch_space(const size_type size) const;
  size_type* scratch_flags(const size_type size) const;
  size_type* scratch_unified(const size_type size) const;
  size_type* scratch_functor(const size_type size) const;

  // Resizing of team level 1 scratch
  void* resize_team_scratch_space(std::int64_t bytes,
                                  bool force_shrink = false);
};

}  // Namespace Impl
}  // Namespace Kokkos
#endif
