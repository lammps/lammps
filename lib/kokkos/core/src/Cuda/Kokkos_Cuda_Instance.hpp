#ifndef KOKKOS_CUDA_INSTANCE_HPP_
#define KOKKOS_CUDA_INSTANCE_HPP_

#include <vector>
#include <impl/Kokkos_Tools.hpp>
#include <atomic>
#include <Cuda/Kokkos_Cuda_Error.hpp>

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
  static constexpr CudaSpace::size_type WarpSize = 32 /* 0x0020 */;
  static constexpr CudaSpace::size_type WarpIndexMask =
      0x001f; /* Mask for warpindex */
  static constexpr CudaSpace::size_type WarpIndexShift =
      5; /* WarpSize == 1 << WarpShift */

  static constexpr CudaSpace::size_type ConstantMemoryUsage =
      0x008000; /* 32k bytes */
  static constexpr CudaSpace::size_type ConstantMemoryCache =
      0x002000; /*  8k bytes */
  static constexpr CudaSpace::size_type KernelArgumentLimit =
      0x001000; /*  4k bytes */
  static constexpr CudaSpace::size_type MaxHierarchicalParallelism =
      1024; /* team_size * vector_length */
  using ConstantGlobalBufferType =
      unsigned long[ConstantMemoryUsage / sizeof(unsigned long)];

  static constexpr int ConstantMemoryUseThreshold = 0x000200 /* 512 bytes */;

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
std::array<CudaSpace::size_type, 3> cuda_internal_maximum_grid_count();
CudaSpace::size_type cuda_internal_maximum_shared_words();

CudaSpace::size_type cuda_internal_maximum_concurrent_block_count();

CudaSpace::size_type* cuda_internal_scratch_flags(const Cuda&,
                                                  const std::size_t size);
CudaSpace::size_type* cuda_internal_scratch_space(const Cuda&,
                                                  const std::size_t size);
CudaSpace::size_type* cuda_internal_scratch_unified(const Cuda&,
                                                    const std::size_t size);

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
  std::array<size_type, 3> m_maxBlock;
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
  mutable std::size_t m_scratchSpaceCount;
  mutable std::size_t m_scratchFlagsCount;
  mutable std::size_t m_scratchUnifiedCount;
  mutable std::size_t m_scratchFunctorSize;

  size_type m_scratchUnifiedSupported;
  size_type m_streamCount;
  mutable size_type* m_scratchSpace;
  mutable size_type* m_scratchFlags;
  mutable size_type* m_scratchUnified;
  mutable size_type* m_scratchFunctor;
  cudaStream_t m_stream;
  uint32_t m_instance_id;
  bool m_manage_stream;

  // Team Scratch Level 1 Space
  int m_n_team_scratch = 10;
  mutable int64_t m_team_scratch_current_size[10];
  mutable void* m_team_scratch_ptr[10];
  mutable std::atomic_int m_team_scratch_pool[10];
  std::int32_t* m_scratch_locks;

  bool was_initialized = false;
  bool was_finalized   = false;

  // FIXME_CUDA: these want to be per-device, not per-stream...  use of 'static'
  //  here will break once there are multiple devices though
  static unsigned long* constantMemHostStaging;
  static cudaEvent_t constantMemReusable;
  static std::mutex constantMemMutex;

  static CudaInternal& singleton();

  int verify_is_initialized(const char* const label) const;

  int is_initialized() const {
    return nullptr != m_scratchSpace && nullptr != m_scratchFlags;
  }

  void initialize(int cuda_device_id, cudaStream_t stream = nullptr,
                  bool manage_stream = false);
  void finalize();

  void print_configuration(std::ostream&) const;

#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
  static bool cuda_use_serial_execution();
  static void cuda_set_serial_execution(bool);
#endif

  void fence(const std::string&) const;
  void fence() const;

  ~CudaInternal();

  CudaInternal()
      : m_cudaDev(-1),
        m_cudaArch(-1),
        m_multiProcCount(0),
        m_maxWarpCount(0),
        m_maxBlock({0, 0, 0}),
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
        m_stream(nullptr),
        m_instance_id(
            Kokkos::Tools::Experimental::Impl::idForInstance<Kokkos::Cuda>(
                reinterpret_cast<uintptr_t>(this))) {
    for (int i = 0; i < m_n_team_scratch; ++i) {
      m_team_scratch_current_size[i] = 0;
      m_team_scratch_ptr[i]          = nullptr;
      m_team_scratch_pool[i]         = 0;
    }
  }

  // Resizing of reduction related scratch spaces
  size_type* scratch_space(const std::size_t size) const;
  size_type* scratch_flags(const std::size_t size) const;
  size_type* scratch_unified(const std::size_t size) const;
  size_type* scratch_functor(const std::size_t size) const;
  uint32_t impl_get_instance_id() const;
  // Resizing of team level 1 scratch
  std::pair<void*, int> resize_team_scratch_space(std::int64_t bytes,
                                                  bool force_shrink = false);
};

}  // Namespace Impl

namespace Experimental {
// Partitioning an Execution Space: expects space and integer arguments for
// relative weight
//   Customization point for backends
//   Default behavior is to return the passed in instance

namespace Impl {
inline void create_Cuda_instances(std::vector<Cuda>& instances) {
  for (int s = 0; s < int(instances.size()); s++) {
    cudaStream_t stream;
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaStreamCreate(&stream));
    instances[s] = Cuda(stream, true);
  }
}
}  // namespace Impl

template <class... Args>
std::vector<Cuda> partition_space(const Cuda&, Args...) {
#ifdef __cpp_fold_expressions
  static_assert(
      (... && std::is_arithmetic_v<Args>),
      "Kokkos Error: partitioning arguments must be integers or floats");
#endif
  std::vector<Cuda> instances(sizeof...(Args));
  Impl::create_Cuda_instances(instances);
  return instances;
}

template <class T>
std::vector<Cuda> partition_space(const Cuda&, std::vector<T>& weights) {
  static_assert(
      std::is_arithmetic<T>::value,
      "Kokkos Error: partitioning arguments must be integers or floats");

  std::vector<Cuda> instances(weights.size());
  Impl::create_Cuda_instances(instances);
  return instances;
}
}  // namespace Experimental

}  // Namespace Kokkos
#endif
