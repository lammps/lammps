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

#ifndef KOKKOS_HIP_KERNEL_LAUNCH_HPP
#define KOKKOS_HIP_KERNEL_LAUNCH_HPP

#include <Kokkos_Macros.hpp>

#if defined(__HIPCC__)

#include <HIP/Kokkos_HIP_Error.hpp>
#include <HIP/Kokkos_HIP_Instance.hpp>
#include <Kokkos_HIP_Space.hpp>
#include <HIP/Kokkos_HIP_Locks.hpp>

// Must use global variable on the device with HIP-Clang
#ifdef __HIP__
#ifdef KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE
__device__ __constant__ extern unsigned long
    kokkos_impl_hip_constant_memory_buffer[];
#else
__device__ __constant__ unsigned long kokkos_impl_hip_constant_memory_buffer
    [Kokkos::Experimental::Impl::HIPTraits::ConstantMemoryUsage /
     sizeof(unsigned long)];
#endif
#endif

namespace Kokkos {
namespace Experimental {
template <typename T>
inline __device__ T *kokkos_impl_hip_shared_memory() {
  extern __shared__ Kokkos::Experimental::HIPSpace::size_type sh[];
  return (T *)sh;
}
}  // namespace Experimental
}  // namespace Kokkos

namespace Kokkos {
namespace Experimental {
namespace Impl {

// The hip_parallel_launch_*_memory code is identical to the cuda code
template <typename DriverType>
__global__ static void hip_parallel_launch_constant_memory() {
  const DriverType &driver = *(reinterpret_cast<const DriverType *>(
      kokkos_impl_hip_constant_memory_buffer));

  driver();
}

template <typename DriverType, unsigned int maxTperB, unsigned int minBperSM>
__global__ __launch_bounds__(
    maxTperB, minBperSM) static void hip_parallel_launch_constant_memory() {
  const DriverType &driver = *(reinterpret_cast<const DriverType *>(
      kokkos_impl_hip_constant_memory_buffer));

  driver();
}

template <class DriverType>
__global__ static void hip_parallel_launch_local_memory(
    const DriverType *driver) {
  // FIXME_HIP driver() pass by copy
  driver->operator()();
}

template <class DriverType, unsigned int maxTperB, unsigned int minBperSM>
__global__ __launch_bounds__(
    maxTperB,
    minBperSM) static void hip_parallel_launch_local_memory(const DriverType
                                                                *driver) {
  // FIXME_HIP driver() pass by copy
  driver->operator()();
}

template <typename DriverType>
__global__ static void hip_parallel_launch_global_memory(
    const DriverType *driver) {
  driver->operator()();
}

template <typename DriverType, unsigned int maxTperB, unsigned int minBperSM>
__global__ __launch_bounds__(
    maxTperB,
    minBperSM) static void hip_parallel_launch_global_memory(const DriverType
                                                                 *driver) {
  driver->operator()();
}

enum class HIPLaunchMechanism : unsigned {
  Default        = 0,
  ConstantMemory = 1,
  GlobalMemory   = 2,
  LocalMemory    = 4
};

constexpr inline HIPLaunchMechanism operator|(HIPLaunchMechanism p1,
                                              HIPLaunchMechanism p2) {
  return static_cast<HIPLaunchMechanism>(static_cast<unsigned>(p1) |
                                         static_cast<unsigned>(p2));
}
constexpr inline HIPLaunchMechanism operator&(HIPLaunchMechanism p1,
                                              HIPLaunchMechanism p2) {
  return static_cast<HIPLaunchMechanism>(static_cast<unsigned>(p1) &
                                         static_cast<unsigned>(p2));
}

template <HIPLaunchMechanism l>
struct HIPDispatchProperties {
  HIPLaunchMechanism launch_mechanism = l;
};

// Use local memory up to ConstantMemoryUseThreshold
// Use global memory above ConstantMemoryUsage
// In between use ConstantMemory
// The following code is identical to the cuda code
template <typename DriverType>
struct DeduceHIPLaunchMechanism {
  static constexpr Kokkos::Experimental::WorkItemProperty::HintLightWeight_t
      light_weight = Kokkos::Experimental::WorkItemProperty::HintLightWeight;
  static constexpr Kokkos::Experimental::WorkItemProperty::HintHeavyWeight_t
      heavy_weight = Kokkos::Experimental::WorkItemProperty::HintHeavyWeight;
  static constexpr typename DriverType::Policy::work_item_property property =
      typename DriverType::Policy::work_item_property();

  static constexpr HIPLaunchMechanism valid_launch_mechanism =
      // BuildValidMask
      (sizeof(DriverType) < HIPTraits::KernelArgumentLimit
           ? HIPLaunchMechanism::LocalMemory
           : HIPLaunchMechanism::Default) |
      (sizeof(DriverType) < HIPTraits::ConstantMemoryUsage
           ? HIPLaunchMechanism::ConstantMemory
           : HIPLaunchMechanism::Default) |
      HIPLaunchMechanism::GlobalMemory;

  static constexpr HIPLaunchMechanism requested_launch_mechanism =
      (((property & light_weight) == light_weight)
           ? HIPLaunchMechanism::LocalMemory
           : HIPLaunchMechanism::ConstantMemory) |
      HIPLaunchMechanism::GlobalMemory;

  static constexpr HIPLaunchMechanism default_launch_mechanism =
      // BuildValidMask
      (sizeof(DriverType) < HIPTraits::ConstantMemoryUseThreshold)
          ? HIPLaunchMechanism::LocalMemory
          : ((sizeof(DriverType) < HIPTraits::ConstantMemoryUsage)
                 ? HIPLaunchMechanism::ConstantMemory
                 : HIPLaunchMechanism::GlobalMemory);

  //              None                LightWeight    HeavyWeight
  // F<UseT       LCG  LCG L  L       LCG  LG L  L   LCG  CG L  C
  // UseT<F<KAL   LCG  LCG C  C       LCG  LG C  L   LCG  CG C  C
  // Kal<F<CMU     CG  LCG C  C        CG  LG C  G    CG  CG C  C
  // CMU<F          G  LCG G  G         G  LG G  G     G  CG G  G
  static constexpr HIPLaunchMechanism launch_mechanism =
      ((property & light_weight) == light_weight)
          ? (sizeof(DriverType) < HIPTraits::KernelArgumentLimit
                 ? HIPLaunchMechanism::LocalMemory
                 : HIPLaunchMechanism::GlobalMemory)
          : (((property & heavy_weight) == heavy_weight)
                 ? (sizeof(DriverType) < HIPTraits::ConstantMemoryUsage
                        ? HIPLaunchMechanism::ConstantMemory
                        : HIPLaunchMechanism::GlobalMemory)
                 : (default_launch_mechanism));
};

template <typename DriverType, typename LaunchBounds,
          HIPLaunchMechanism LaunchMechanism>
struct HIPParallelLaunchKernelFuncData {
  static unsigned int get_scratch_size(
      hipFuncAttributes const &hip_func_attributes) {
    return hip_func_attributes.localSizeBytes;
  }

  static hipFuncAttributes get_hip_func_attributes(void const *kernel_func) {
    static hipFuncAttributes attr = [=]() {
      hipFuncAttributes attr;
      KOKKOS_IMPL_HIP_SAFE_CALL(hipFuncGetAttributes(&attr, kernel_func));
      return attr;
    }();
    return attr;
  }
};

//---------------------------------------------------------------//
// HIPParallelLaunchKernelFunc structure and its specializations //
//---------------------------------------------------------------//
template <typename DriverType, typename LaunchBounds,
          HIPLaunchMechanism LaunchMechanism>
struct HIPParallelLaunchKernelFunc;

// HIPLaunchMechanism::LocalMemory specializations
template <typename DriverType, unsigned int MaxThreadsPerBlock,
          unsigned int MinBlocksPerSM>
struct HIPParallelLaunchKernelFunc<
    DriverType, Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>,
    HIPLaunchMechanism::LocalMemory> {
  using funcdata_t = HIPParallelLaunchKernelFuncData<
      DriverType, Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>,
      HIPLaunchMechanism::LocalMemory>;
  static auto get_kernel_func() {
    return hip_parallel_launch_local_memory<DriverType, MaxThreadsPerBlock,
                                            MinBlocksPerSM>;
  }

  static constexpr auto default_launchbounds() { return false; }

  static auto get_scratch_size() {
    return funcdata_t::get_scratch_size(get_hip_func_attributes());
  }

  static hipFuncAttributes get_hip_func_attributes() {
    return funcdata_t::get_hip_func_attributes(
        reinterpret_cast<void const *>(get_kernel_func()));
  }
};

template <typename DriverType>
struct HIPParallelLaunchKernelFunc<DriverType, Kokkos::LaunchBounds<0, 0>,
                                   HIPLaunchMechanism::LocalMemory> {
  using funcdata_t =
      HIPParallelLaunchKernelFuncData<DriverType, Kokkos::LaunchBounds<0, 0>,
                                      HIPLaunchMechanism::LocalMemory>;
  static auto get_kernel_func() {
    return HIPParallelLaunchKernelFunc<
        DriverType, Kokkos::LaunchBounds<HIPTraits::MaxThreadsPerBlock, 1>,
        HIPLaunchMechanism::LocalMemory>::get_kernel_func();
  }

  static constexpr auto default_launchbounds() { return true; }

  static auto get_scratch_size() {
    return funcdata_t::get_scratch_size(get_hip_func_attributes());
  }

  static hipFuncAttributes get_hip_func_attributes() {
    return funcdata_t::get_hip_func_attributes(
        reinterpret_cast<void const *>(get_kernel_func()));
  }
};

// HIPLaunchMechanism::GlobalMemory specializations
template <typename DriverType, unsigned int MaxThreadsPerBlock,
          unsigned int MinBlocksPerSM>
struct HIPParallelLaunchKernelFunc<
    DriverType, Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>,
    HIPLaunchMechanism::GlobalMemory> {
  using funcdata_t = HIPParallelLaunchKernelFuncData<
      DriverType, Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>,
      HIPLaunchMechanism::GlobalMemory>;
  static auto get_kernel_func() {
    return hip_parallel_launch_global_memory<DriverType, MaxThreadsPerBlock,
                                             MinBlocksPerSM>;
  }

  static constexpr auto default_launchbounds() { return false; }

  static auto get_scratch_size() {
    return funcdata_t::get_scratch_size(get_hip_func_attributes());
  }

  static hipFuncAttributes get_hip_func_attributes() {
    return funcdata_t::get_hip_func_attributes(
        reinterpret_cast<void const *>(get_kernel_func()));
  }
};

template <typename DriverType>
struct HIPParallelLaunchKernelFunc<DriverType, Kokkos::LaunchBounds<0, 0>,
                                   HIPLaunchMechanism::GlobalMemory> {
  using funcdata_t =
      HIPParallelLaunchKernelFuncData<DriverType, Kokkos::LaunchBounds<0, 0>,
                                      HIPLaunchMechanism::GlobalMemory>;
  static auto get_kernel_func() {
    return hip_parallel_launch_global_memory<DriverType>;
  }

  static constexpr auto default_launchbounds() { return true; }

  static auto get_scratch_size() {
    return funcdata_t::get_scratch_size(get_hip_func_attributes());
  }

  static hipFuncAttributes get_hip_func_attributes() {
    return funcdata_t::get_hip_func_attributes(
        reinterpret_cast<void const *>(get_kernel_func()));
  }
};

// HIPLaunchMechanism::ConstantMemory specializations
template <typename DriverType, unsigned int MaxThreadsPerBlock,
          unsigned int MinBlocksPerSM>
struct HIPParallelLaunchKernelFunc<
    DriverType, Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>,
    HIPLaunchMechanism::ConstantMemory> {
  using funcdata_t = HIPParallelLaunchKernelFuncData<
      DriverType, Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>,
      HIPLaunchMechanism::ConstantMemory>;
  static auto get_kernel_func() {
    return hip_parallel_launch_constant_memory<DriverType, MaxThreadsPerBlock,
                                               MinBlocksPerSM>;
  }

  static constexpr auto default_launchbounds() { return false; }

  static auto get_scratch_size() {
    return funcdata_t::get_scratch_size(get_hip_func_attributes());
  }

  static hipFuncAttributes get_hip_func_attributes() {
    return funcdata_t::get_hip_func_attributes(
        reinterpret_cast<void const *>(get_kernel_func()));
  }
};

template <typename DriverType>
struct HIPParallelLaunchKernelFunc<DriverType, Kokkos::LaunchBounds<0, 0>,
                                   HIPLaunchMechanism::ConstantMemory> {
  using funcdata_t =
      HIPParallelLaunchKernelFuncData<DriverType, Kokkos::LaunchBounds<0, 0>,
                                      HIPLaunchMechanism::ConstantMemory>;
  static auto get_kernel_func() {
    return hip_parallel_launch_constant_memory<DriverType>;
  }
  static constexpr auto default_launchbounds() { return true; }

  static auto get_scratch_size() {
    return funcdata_t::get_scratch_size(get_hip_func_attributes());
  }

  static hipFuncAttributes get_hip_func_attributes() {
    return funcdata_t::get_hip_func_attributes(
        reinterpret_cast<void const *>(get_kernel_func()));
  }
};

//------------------------------------------------------------------//
// HIPParallelLaunchKernelInvoker structure and its specializations //
//------------------------------------------------------------------//
template <typename DriverType, typename LaunchBounds,
          HIPLaunchMechanism LaunchMechanism>
struct HIPParallelLaunchKernelInvoker;

// HIPLaunchMechanism::LocalMemory specialization
template <typename DriverType, typename LaunchBounds>
struct HIPParallelLaunchKernelInvoker<DriverType, LaunchBounds,
                                      HIPLaunchMechanism::LocalMemory>
    : HIPParallelLaunchKernelFunc<DriverType, LaunchBounds,
                                  HIPLaunchMechanism::LocalMemory> {
  using base_t = HIPParallelLaunchKernelFunc<DriverType, LaunchBounds,
                                             HIPLaunchMechanism::LocalMemory>;

  static void invoke_kernel(DriverType const *driver, dim3 const &grid,
                            dim3 const &block, int shmem,
                            HIPInternal const *hip_instance) {
    (base_t::get_kernel_func())<<<grid, block, shmem, hip_instance->m_stream>>>(
        driver);
  }
};

// HIPLaunchMechanism::GlobalMemory specialization
template <typename DriverType, typename LaunchBounds>
struct HIPParallelLaunchKernelInvoker<DriverType, LaunchBounds,
                                      HIPLaunchMechanism::GlobalMemory>
    : HIPParallelLaunchKernelFunc<DriverType, LaunchBounds,
                                  HIPLaunchMechanism::GlobalMemory> {
  using base_t = HIPParallelLaunchKernelFunc<DriverType, LaunchBounds,
                                             HIPLaunchMechanism::GlobalMemory>;

  // FIXME_HIP the code is different than cuda because driver cannot be passed
  // by copy
  static void invoke_kernel(DriverType const *driver, dim3 const &grid,
                            dim3 const &block, int shmem,
                            HIPInternal const *hip_instance) {
    (base_t::get_kernel_func())<<<grid, block, shmem, hip_instance->m_stream>>>(
        driver);
  }
};

// HIPLaunchMechanism::ConstantMemory specializations
template <typename DriverType, typename LaunchBounds>
struct HIPParallelLaunchKernelInvoker<DriverType, LaunchBounds,
                                      HIPLaunchMechanism::ConstantMemory>
    : HIPParallelLaunchKernelFunc<DriverType, LaunchBounds,
                                  HIPLaunchMechanism::ConstantMemory> {
  using base_t =
      HIPParallelLaunchKernelFunc<DriverType, LaunchBounds,
                                  HIPLaunchMechanism::ConstantMemory>;
  static_assert(sizeof(DriverType) < HIPTraits::ConstantMemoryUsage,
                "Kokkos Error: Requested HIPLaunchConstantMemory with a "
                "Functor larger than 32kB.");

  static void invoke_kernel(DriverType const *driver, dim3 const &grid,
                            dim3 const &block, int shmem,
                            HIPInternal const *hip_instance) {
    // Wait until the previous kernel that uses the constant buffer is done
    std::lock_guard<std::mutex> lock(HIPInternal::constantMemMutex);
    KOKKOS_IMPL_HIP_SAFE_CALL(
        hipEventSynchronize(HIPInternal::constantMemReusable));

    // Copy functor (synchronously) to staging buffer in pinned host memory
    unsigned long *staging = hip_instance->constantMemHostStaging;
    std::memcpy((void *)staging, (void *)driver, sizeof(DriverType));

    // Copy functor asynchronously from there to constant memory on the device
    KOKKOS_IMPL_HIP_SAFE_CALL(hipMemcpyToSymbolAsync(
        HIP_SYMBOL(kokkos_impl_hip_constant_memory_buffer), staging,
        sizeof(DriverType), 0, hipMemcpyHostToDevice, hip_instance->m_stream));

    // Invoke the driver function on the device
    (base_t::
         get_kernel_func())<<<grid, block, shmem, hip_instance->m_stream>>>();

    // Record an event that says when the constant buffer can be reused
    KOKKOS_IMPL_HIP_SAFE_CALL(hipEventRecord(HIPInternal::constantMemReusable,
                                             hip_instance->m_stream));
  }
};

//-----------------------------//
// HIPParallelLaunch structure //
//-----------------------------//
template <typename DriverType, typename LaunchBounds = Kokkos::LaunchBounds<>,
          HIPLaunchMechanism LaunchMechanism =
              DeduceHIPLaunchMechanism<DriverType>::launch_mechanism>
struct HIPParallelLaunch;

template <typename DriverType, unsigned int MaxThreadsPerBlock,
          unsigned int MinBlocksPerSM, HIPLaunchMechanism LaunchMechanism>
struct HIPParallelLaunch<
    DriverType, Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>,
    LaunchMechanism>
    : HIPParallelLaunchKernelInvoker<
          DriverType, Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>,
          LaunchMechanism> {
  using base_t = HIPParallelLaunchKernelInvoker<
      DriverType, Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>,
      LaunchMechanism>;

  HIPParallelLaunch(const DriverType &driver, const dim3 &grid,
                    const dim3 &block, const int shmem,
                    const HIPInternal *hip_instance,
                    const bool /*prefer_shmem*/) {
    if ((grid.x != 0) && ((block.x * block.y * block.z) != 0)) {
      if (hip_instance->m_maxShmemPerBlock < shmem) {
        Kokkos::Impl::throw_runtime_exception(
            "HIPParallelLaunch FAILED: shared memory request is too large");
      }

      KOKKOS_ENSURE_HIP_LOCK_ARRAYS_ON_DEVICE();

      // Invoke the driver function on the device
      DriverType *d_driver = reinterpret_cast<DriverType *>(
          hip_instance->get_next_driver(sizeof(DriverType)));
      std::memcpy((void *)d_driver, (void *)&driver, sizeof(DriverType));
      base_t::invoke_kernel(d_driver, grid, block, shmem, hip_instance);

#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
      KOKKOS_IMPL_HIP_SAFE_CALL(hipGetLastError());
      hip_instance->fence(
          "Kokkos::Experimental::Impl::HIParallelLaunch: Debug Only Check for "
          "Execution Error");
#endif
    }
  }
};

// convenience method to launch the correct kernel given the launch bounds et
// al.
template <typename DriverType, typename LaunchBounds = Kokkos::LaunchBounds<>,
          HIPLaunchMechanism LaunchMechanism =
              DeduceHIPLaunchMechanism<DriverType>::launch_mechanism>
void hip_parallel_launch(const DriverType &driver, const dim3 &grid,
                         const dim3 &block, const int shmem,
                         const HIPInternal *hip_instance,
                         const bool prefer_shmem) {
#ifndef KOKKOS_ENABLE_HIP_MULTIPLE_KERNEL_INSTANTIATIONS
  HIPParallelLaunch<DriverType, LaunchBounds, LaunchMechanism>(
      driver, grid, block, shmem, hip_instance, prefer_shmem);
#else
  // FIXME_HIP - could be if constexpr for c++17
  if (!HIPParallelLaunch<DriverType, LaunchBounds,
                         LaunchMechanism>::default_launchbounds()) {
    // for user defined, we *always* honor the request
    HIPParallelLaunch<DriverType, LaunchBounds, LaunchMechanism>(
        driver, grid, block, shmem, hip_instance, prefer_shmem);
  } else {
    // we can do what we like
    const unsigned flat_block_size = block.x * block.y * block.z;
    if (flat_block_size <= HIPTraits::ConservativeThreadsPerBlock) {
      // we have to use the large blocksize
      HIPParallelLaunch<
          DriverType,
          Kokkos::LaunchBounds<HIPTraits::ConservativeThreadsPerBlock, 1>,
          LaunchMechanism>(driver, grid, block, shmem, hip_instance,
                           prefer_shmem);
    } else {
      HIPParallelLaunch<DriverType,
                        Kokkos::LaunchBounds<HIPTraits::MaxThreadsPerBlock, 1>,
                        LaunchMechanism>(driver, grid, block, shmem,
                                         hip_instance, prefer_shmem);
    }
  }
#endif
}
}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif

#endif
