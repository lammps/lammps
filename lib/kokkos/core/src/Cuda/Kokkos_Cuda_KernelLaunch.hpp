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

#ifndef KOKKOS_CUDAEXEC_HPP
#define KOKKOS_CUDAEXEC_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_CUDA

#include <string>
#include <cstdint>
#include <Kokkos_Parallel.hpp>
#include <impl/Kokkos_Error.hpp>
#include <Cuda/Kokkos_Cuda_abort.hpp>
#include <Cuda/Kokkos_Cuda_Error.hpp>
#include <Cuda/Kokkos_Cuda_Locks.hpp>
#include <Cuda/Kokkos_Cuda_Instance.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined(__CUDACC__)

/** \brief  Access to constant memory on the device */
#ifdef KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE

__device__ __constant__ extern unsigned long
    kokkos_impl_cuda_constant_memory_buffer[];

#else

__device__ __constant__ unsigned long kokkos_impl_cuda_constant_memory_buffer
    [Kokkos::Impl::CudaTraits::ConstantMemoryUsage / sizeof(unsigned long)];

#endif

namespace Kokkos {
namespace Impl {
void* cuda_resize_scratch_space(std::int64_t bytes, bool force_shrink = false);
}
}  // namespace Kokkos

template <typename T>
inline __device__ T* kokkos_impl_cuda_shared_memory() {
  extern __shared__ Kokkos::CudaSpace::size_type sh[];
  return (T*)sh;
}

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
// See section B.17 of Cuda C Programming Guide Version 3.2
// for discussion of
//   __launch_bounds__(maxThreadsPerBlock,minBlocksPerMultiprocessor)
// function qualifier which could be used to improve performance.
//----------------------------------------------------------------------------
// Maximize L1 cache and minimize shared memory:
//   cudaFuncSetCacheConfig(MyKernel, cudaFuncCachePreferL1 );
// For 2.0 capability: 48 KB L1 and 16 KB shared
//----------------------------------------------------------------------------

template <class DriverType>
__global__ static void cuda_parallel_launch_constant_memory() {
  const DriverType& driver =
      *((const DriverType*)kokkos_impl_cuda_constant_memory_buffer);

  driver();
}

template <class DriverType, unsigned int maxTperB, unsigned int minBperSM>
__global__ __launch_bounds__(
    maxTperB, minBperSM) static void cuda_parallel_launch_constant_memory() {
  const DriverType& driver =
      *((const DriverType*)kokkos_impl_cuda_constant_memory_buffer);

  driver();
}

template <class DriverType>
__global__ static void cuda_parallel_launch_local_memory(
    const DriverType driver) {
  driver();
}

template <class DriverType, unsigned int maxTperB, unsigned int minBperSM>
__global__ __launch_bounds__(
    maxTperB,
    minBperSM) static void cuda_parallel_launch_local_memory(const DriverType
                                                                 driver) {
  driver();
}

template <class DriverType>
__global__ static void cuda_parallel_launch_global_memory(
    const DriverType* driver) {
  driver->operator()();
}

template <class DriverType, unsigned int maxTperB, unsigned int minBperSM>
__global__ __launch_bounds__(
    maxTperB,
    minBperSM) static void cuda_parallel_launch_global_memory(const DriverType*
                                                                  driver) {
  driver->operator()();
}

template <class DriverType>
__global__ static void cuda_parallel_launch_constant_or_global_memory(
    const DriverType* driver_ptr) {
  const DriverType& driver =
      driver_ptr != NULL
          ? *driver_ptr
          : *((const DriverType*)kokkos_impl_cuda_constant_memory_buffer);

  driver();
}

template <class DriverType, unsigned int maxTperB, unsigned int minBperSM>
__global__
__launch_bounds__(maxTperB, minBperSM) static void cuda_parallel_launch_constant_or_global_memory(
    const DriverType* driver_ptr) {
  const DriverType& driver =
      driver_ptr != NULL
          ? *driver_ptr
          : *((const DriverType*)kokkos_impl_cuda_constant_memory_buffer);

  driver();
}

template <class DriverType>
struct DeduceCudaLaunchMechanism {
  constexpr static const Kokkos::Experimental::WorkItemProperty::
      HintLightWeight_t light_weight =
          Kokkos::Experimental::WorkItemProperty::HintLightWeight;
  constexpr static const Kokkos::Experimental::WorkItemProperty::
      HintHeavyWeight_t heavy_weight =
          Kokkos::Experimental::WorkItemProperty::HintHeavyWeight;
  constexpr static const typename DriverType::Policy::work_item_property
      property = typename DriverType::Policy::work_item_property();

  static constexpr const Experimental::CudaLaunchMechanism
      valid_launch_mechanism =
          // BuildValidMask
      (sizeof(DriverType) < CudaTraits::KernelArgumentLimit
           ? Experimental::CudaLaunchMechanism::LocalMemory
           : Experimental::CudaLaunchMechanism::Default) |
      (sizeof(DriverType) < CudaTraits::ConstantMemoryUsage
           ? Experimental::CudaLaunchMechanism::ConstantMemory
           : Experimental::CudaLaunchMechanism::Default) |
      Experimental::CudaLaunchMechanism::GlobalMemory;

  static constexpr const Experimental::CudaLaunchMechanism
      requested_launch_mechanism =
          (((property & light_weight) == light_weight)
               ? Experimental::CudaLaunchMechanism::LocalMemory
               : Experimental::CudaLaunchMechanism::ConstantMemory) |
          Experimental::CudaLaunchMechanism::GlobalMemory;

  static constexpr const Experimental::CudaLaunchMechanism
      default_launch_mechanism =
          // BuildValidMask
      (sizeof(DriverType) < CudaTraits::ConstantMemoryUseThreshold)
          ? Experimental::CudaLaunchMechanism::LocalMemory
          : ((sizeof(DriverType) < CudaTraits::ConstantMemoryUsage)
                 ? Experimental::CudaLaunchMechanism::ConstantMemory
                 : Experimental::CudaLaunchMechanism::GlobalMemory);

  //              None                LightWeight    HeavyWeight
  // F<UseT       LCG LCG L  L        LCG  LG L  L    LCG  CG L  C
  // UseT<F<KAL   LCG LCG C  C        LCG  LG C  L    LCG  CG C  C
  // Kal<F<CMU     CG LCG C  C         CG  LG C  G     CG  CG C  C
  // CMU<F          G LCG G  G          G  LG G  G      G  CG G  G
  static constexpr const Experimental::CudaLaunchMechanism launch_mechanism =
      ((property & light_weight) == light_weight)
          ? (sizeof(DriverType) < CudaTraits::KernelArgumentLimit
                 ? Experimental::CudaLaunchMechanism::LocalMemory
                 : Experimental::CudaLaunchMechanism::GlobalMemory)
          : (((property & heavy_weight) == heavy_weight)
                 ? (sizeof(DriverType) < CudaTraits::ConstantMemoryUsage
                        ? Experimental::CudaLaunchMechanism::ConstantMemory
                        : Experimental::CudaLaunchMechanism::GlobalMemory)
                 : (default_launch_mechanism));
};
// Use local memory up to ConstantMemoryUseThreshold
// Use global memory above ConstantMemoryUsage
// In between use ConstantMemory
template <class DriverType, class LaunchBounds = Kokkos::LaunchBounds<>,
          Experimental::CudaLaunchMechanism LaunchMechanism =
              DeduceCudaLaunchMechanism<DriverType>::launch_mechanism>
struct CudaParallelLaunch;

template <class DriverType, unsigned int MaxThreadsPerBlock,
          unsigned int MinBlocksPerSM>
struct CudaParallelLaunch<
    DriverType, Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>,
    Experimental::CudaLaunchMechanism::ConstantMemory> {
  static_assert(sizeof(DriverType) < CudaTraits::ConstantMemoryUsage,
                "Kokkos Error: Requested CudaLaunchConstantMemory with a "
                "Functor larger than 32kB.");
  inline CudaParallelLaunch(const DriverType& driver, const dim3& grid,
                            const dim3& block, const int shmem,
                            const CudaInternal* cuda_instance,
                            const bool prefer_shmem) {
    if ((grid.x != 0) && ((block.x * block.y * block.z) != 0)) {
      // Fence before changing settings and copying closure
      Kokkos::Cuda().fence();

      if (cuda_instance->m_maxShmemPerBlock < shmem) {
        Kokkos::Impl::throw_runtime_exception(std::string(
            "CudaParallelLaunch FAILED: shared memory request is too large"));
      }
#ifndef KOKKOS_ARCH_KEPLER
      // On Kepler the L1 has no benefit since it doesn't cache reads
      else {
        CUDA_SAFE_CALL(cudaFuncSetCacheConfig(
            cuda_parallel_launch_constant_memory<DriverType, MaxThreadsPerBlock,
                                                 MinBlocksPerSM>,
            (prefer_shmem ? cudaFuncCachePreferShared
                          : cudaFuncCachePreferL1)));
      }
#endif

      // Copy functor to constant memory on the device
      cudaMemcpyToSymbolAsync(kokkos_impl_cuda_constant_memory_buffer, &driver,
                              sizeof(DriverType), 0, cudaMemcpyHostToDevice,
                              cudaStream_t(cuda_instance->m_stream));

      KOKKOS_ENSURE_CUDA_LOCK_ARRAYS_ON_DEVICE();

      // Invoke the driver function on the device
      cuda_parallel_launch_constant_memory<DriverType, MaxThreadsPerBlock,
                                           MinBlocksPerSM>
          <<<grid, block, shmem, cuda_instance->m_stream>>>();

#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
      CUDA_SAFE_CALL(cudaGetLastError());
      Kokkos::Cuda().fence();
#endif
    }
  }

  static cudaFuncAttributes get_cuda_func_attributes() {
    cudaFuncAttributes attr;
    CUDA_SAFE_CALL(cudaFuncGetAttributes(
        &attr,
        cuda_parallel_launch_constant_memory<DriverType, MaxThreadsPerBlock,
                                             MinBlocksPerSM>));
    return attr;
  }
};

template <class DriverType>
struct CudaParallelLaunch<DriverType, Kokkos::LaunchBounds<0, 0>,
                          Experimental::CudaLaunchMechanism::ConstantMemory> {
  static_assert(sizeof(DriverType) < CudaTraits::ConstantMemoryUsage,
                "Kokkos Error: Requested CudaLaunchConstantMemory with a "
                "Functor larger than 32kB.");
  inline CudaParallelLaunch(const DriverType& driver, const dim3& grid,
                            const dim3& block, const int shmem,
                            const CudaInternal* cuda_instance,
                            const bool prefer_shmem) {
    if ((grid.x != 0) && ((block.x * block.y * block.z) != 0)) {
      // Fence before changing settings and copying closure
      Kokkos::Cuda().fence();

      if (cuda_instance->m_maxShmemPerBlock < shmem) {
        Kokkos::Impl::throw_runtime_exception(std::string(
            "CudaParallelLaunch FAILED: shared memory request is too large"));
      }
#ifndef KOKKOS_ARCH_KEPLER
      // On Kepler the L1 has no benefit since it doesn't cache reads
      else {
        CUDA_SAFE_CALL(cudaFuncSetCacheConfig(
            cuda_parallel_launch_constant_memory<DriverType>,
            (prefer_shmem ? cudaFuncCachePreferShared
                          : cudaFuncCachePreferL1)));
      }
#endif

      // Copy functor to constant memory on the device
      cudaMemcpyToSymbolAsync(kokkos_impl_cuda_constant_memory_buffer, &driver,
                              sizeof(DriverType), 0, cudaMemcpyHostToDevice,
                              cudaStream_t(cuda_instance->m_stream));

      KOKKOS_ENSURE_CUDA_LOCK_ARRAYS_ON_DEVICE();

      // Invoke the driver function on the device
      cuda_parallel_launch_constant_memory<DriverType>
          <<<grid, block, shmem, cuda_instance->m_stream>>>();

#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
      CUDA_SAFE_CALL(cudaGetLastError());
      Kokkos::Cuda().fence();
#endif
    }
  }

  static cudaFuncAttributes get_cuda_func_attributes() {
    cudaFuncAttributes attr;
    CUDA_SAFE_CALL(cudaFuncGetAttributes(
        &attr, cuda_parallel_launch_constant_memory<DriverType>));
    return attr;
  }
};

template <class DriverType, unsigned int MaxThreadsPerBlock,
          unsigned int MinBlocksPerSM>
struct CudaParallelLaunch<
    DriverType, Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>,
    Experimental::CudaLaunchMechanism::LocalMemory> {
  static_assert(sizeof(DriverType) < CudaTraits::KernelArgumentLimit,
                "Kokkos Error: Requested CudaLaunchLocalMemory with a Functor "
                "larger than 4096 bytes.");
  inline CudaParallelLaunch(const DriverType& driver, const dim3& grid,
                            const dim3& block, const int shmem,
                            const CudaInternal* cuda_instance,
                            const bool prefer_shmem) {
    if ((grid.x != 0) && ((block.x * block.y * block.z) != 0)) {
      if (cuda_instance->m_maxShmemPerBlock < shmem) {
        Kokkos::Impl::throw_runtime_exception(std::string(
            "CudaParallelLaunch FAILED: shared memory request is too large"));
      }
#ifndef KOKKOS_ARCH_KEPLER
      // On Kepler the L1 has no benefit since it doesn't cache reads
      else {
        CUDA_SAFE_CALL(cudaFuncSetCacheConfig(
            cuda_parallel_launch_local_memory<DriverType, MaxThreadsPerBlock,
                                              MinBlocksPerSM>,
            (prefer_shmem ? cudaFuncCachePreferShared
                          : cudaFuncCachePreferL1)));
      }
#endif

      KOKKOS_ENSURE_CUDA_LOCK_ARRAYS_ON_DEVICE();

      // Invoke the driver function on the device
      cuda_parallel_launch_local_memory<DriverType, MaxThreadsPerBlock,
                                        MinBlocksPerSM>
          <<<grid, block, shmem, cuda_instance->m_stream>>>(driver);

#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
      CUDA_SAFE_CALL(cudaGetLastError());
      Kokkos::Cuda().fence();
#endif
    }
  }

  static cudaFuncAttributes get_cuda_func_attributes() {
    cudaFuncAttributes attr;
    CUDA_SAFE_CALL(cudaFuncGetAttributes(
        &attr, cuda_parallel_launch_local_memory<DriverType, MaxThreadsPerBlock,
                                                 MinBlocksPerSM>));
    return attr;
  }
};

template <class DriverType>
struct CudaParallelLaunch<DriverType, Kokkos::LaunchBounds<0, 0>,
                          Experimental::CudaLaunchMechanism::LocalMemory> {
  static_assert(sizeof(DriverType) < CudaTraits::KernelArgumentLimit,
                "Kokkos Error: Requested CudaLaunchLocalMemory with a Functor "
                "larger than 4096 bytes.");
  inline CudaParallelLaunch(const DriverType& driver, const dim3& grid,
                            const dim3& block, const int shmem,
                            const CudaInternal* cuda_instance,
                            const bool prefer_shmem) {
    if ((grid.x != 0) && ((block.x * block.y * block.z) != 0)) {
      if (cuda_instance->m_maxShmemPerBlock < shmem) {
        Kokkos::Impl::throw_runtime_exception(std::string(
            "CudaParallelLaunch FAILED: shared memory request is too large"));
      }
#ifndef KOKKOS_ARCH_KEPLER
      // On Kepler the L1 has no benefit since it doesn't cache reads
      else {
        CUDA_SAFE_CALL(cudaFuncSetCacheConfig(
            cuda_parallel_launch_local_memory<DriverType>,
            (prefer_shmem ? cudaFuncCachePreferShared
                          : cudaFuncCachePreferL1)));
      }
#endif

      KOKKOS_ENSURE_CUDA_LOCK_ARRAYS_ON_DEVICE();

      // Invoke the driver function on the device
      cuda_parallel_launch_local_memory<DriverType>
          <<<grid, block, shmem, cuda_instance->m_stream>>>(driver);

#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
      CUDA_SAFE_CALL(cudaGetLastError());
      Kokkos::Cuda().fence();
#endif
    }
  }

  static cudaFuncAttributes get_cuda_func_attributes() {
    cudaFuncAttributes attr;
    CUDA_SAFE_CALL(cudaFuncGetAttributes(
        &attr, cuda_parallel_launch_local_memory<DriverType>));
    return attr;
  }
};

template <class DriverType, unsigned int MaxThreadsPerBlock,
          unsigned int MinBlocksPerSM>
struct CudaParallelLaunch<
    DriverType, Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>,
    Experimental::CudaLaunchMechanism::GlobalMemory> {
  inline CudaParallelLaunch(const DriverType& driver, const dim3& grid,
                            const dim3& block, const int shmem,
                            CudaInternal* cuda_instance,
                            const bool prefer_shmem) {
    if ((grid.x != 0) && ((block.x * block.y * block.z) != 0)) {
      if (cuda_instance->m_maxShmemPerBlock < shmem) {
        Kokkos::Impl::throw_runtime_exception(std::string(
            "CudaParallelLaunch FAILED: shared memory request is too large"));
      }
#ifndef KOKKOS_ARCH_KEPLER
      // On Kepler the L1 has no benefit since it doesn't cache reads
      else {
        CUDA_SAFE_CALL(cudaFuncSetCacheConfig(
            cuda_parallel_launch_global_memory<DriverType, MaxThreadsPerBlock,
                                               MinBlocksPerSM>,
            (prefer_shmem ? cudaFuncCachePreferShared
                          : cudaFuncCachePreferL1)));
      }
#endif

      KOKKOS_ENSURE_CUDA_LOCK_ARRAYS_ON_DEVICE();

      DriverType* driver_ptr = NULL;
      driver_ptr             = reinterpret_cast<DriverType*>(
          cuda_instance->scratch_functor(sizeof(DriverType)));
      cudaMemcpyAsync(driver_ptr, &driver, sizeof(DriverType),
                      cudaMemcpyDefault, cuda_instance->m_stream);

      // Invoke the driver function on the device
      cuda_parallel_launch_global_memory<DriverType, MaxThreadsPerBlock,
                                         MinBlocksPerSM>
          <<<grid, block, shmem, cuda_instance->m_stream>>>(driver_ptr);

#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
      CUDA_SAFE_CALL(cudaGetLastError());
      Kokkos::Cuda().fence();
#endif
    }
  }
  static cudaFuncAttributes get_cuda_func_attributes() {
    cudaFuncAttributes attr;
    CUDA_SAFE_CALL(cudaFuncGetAttributes(
        &attr,
        cuda_parallel_launch_global_memory<DriverType, MaxThreadsPerBlock,
                                           MinBlocksPerSM>));
    return attr;
  }
};

template <class DriverType>
struct CudaParallelLaunch<DriverType, Kokkos::LaunchBounds<0, 0>,
                          Experimental::CudaLaunchMechanism::GlobalMemory> {
  inline CudaParallelLaunch(const DriverType& driver, const dim3& grid,
                            const dim3& block, const int shmem,
                            CudaInternal* cuda_instance,
                            const bool prefer_shmem) {
    if ((grid.x != 0) && ((block.x * block.y * block.z) != 0)) {
      if (cuda_instance->m_maxShmemPerBlock < shmem) {
        Kokkos::Impl::throw_runtime_exception(std::string(
            "CudaParallelLaunch FAILED: shared memory request is too large"));
      }
#ifndef KOKKOS_ARCH_KEPLER
      // On Kepler the L1 has no benefit since it doesn't cache reads
      else {
        CUDA_SAFE_CALL(cudaFuncSetCacheConfig(
            cuda_parallel_launch_global_memory<DriverType>,
            (prefer_shmem ? cudaFuncCachePreferShared
                          : cudaFuncCachePreferL1)));
      }
#endif

      KOKKOS_ENSURE_CUDA_LOCK_ARRAYS_ON_DEVICE();

      DriverType* driver_ptr = NULL;
      driver_ptr             = reinterpret_cast<DriverType*>(
          cuda_instance->scratch_functor(sizeof(DriverType)));
      cudaMemcpyAsync(driver_ptr, &driver, sizeof(DriverType),
                      cudaMemcpyDefault, cuda_instance->m_stream);

      cuda_parallel_launch_global_memory<DriverType>
          <<<grid, block, shmem, cuda_instance->m_stream>>>(driver_ptr);

#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
      CUDA_SAFE_CALL(cudaGetLastError());
      Kokkos::Cuda().fence();
#endif
    }
  }

  static cudaFuncAttributes get_cuda_func_attributes() {
    cudaFuncAttributes attr;
    CUDA_SAFE_CALL(cudaFuncGetAttributes(
        &attr, cuda_parallel_launch_global_memory<DriverType>));
    return attr;
  }
};
//----------------------------------------------------------------------------

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* defined( __CUDACC__ ) */
#endif /* defined( KOKKOS_ENABLE_CUDA ) */
#endif /* #ifndef KOKKOS_CUDAEXEC_HPP */
