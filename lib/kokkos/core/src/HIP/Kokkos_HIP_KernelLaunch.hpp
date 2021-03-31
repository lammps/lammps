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

#include <Kokkos_HIP_Space.hpp>
#include <HIP/Kokkos_HIP_Error.hpp>
#include <HIP/Kokkos_HIP_Instance.hpp>

// Must use global variable on the device with HIP-Clang
#ifdef __HIP__
__device__ __constant__ unsigned long kokkos_impl_hip_constant_memory_buffer
    [Kokkos::Experimental::Impl::HIPTraits::ConstantMemoryUsage /
     sizeof(unsigned long)];
#endif

namespace Kokkos {
namespace Experimental {
template <typename T>
inline __device__ T *kokkos_impl_hip_shared_memory() {
  HIP_DYNAMIC_SHARED(HIPSpace::size_type, sh);
  return (T *)sh;
}
}  // namespace Experimental
}  // namespace Kokkos

namespace Kokkos {
namespace Experimental {
namespace Impl {

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

  driver->operator()();
}

template <class DriverType>
__global__ static void hip_parallel_launch_local_memory(
    const DriverType *driver) {
  driver->operator()();
}

template <class DriverType, unsigned int maxTperB, unsigned int minBperSM>
__global__ __launch_bounds__(
    maxTperB,
    minBperSM) static void hip_parallel_launch_local_memory(const DriverType
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

template <class DriverType, class LaunchBounds = Kokkos::LaunchBounds<>,
          HIPLaunchMechanism LaunchMechanism = HIPLaunchMechanism::LocalMemory>
struct HIPParallelLaunch;

template <class DriverType, unsigned int MaxThreadsPerBlock,
          unsigned int MinBlocksPerSM>
struct HIPParallelLaunch<
    DriverType, Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>,
    HIPLaunchMechanism::LocalMemory> {
  inline HIPParallelLaunch(const DriverType &driver, const dim3 &grid,
                           const dim3 &block, const int shmem,
                           const HIPInternal *hip_instance,
                           const bool /*prefer_shmem*/) {
    if ((grid.x != 0) && ((block.x * block.y * block.z) != 0)) {
      if (hip_instance->m_maxShmemPerBlock < shmem) {
        Kokkos::Impl::throw_runtime_exception(
            "HIPParallelLaunch FAILED: shared memory request is too large");
      }

      KOKKOS_ENSURE_HIP_LOCK_ARRAYS_ON_DEVICE();

      // FIXME_HIP -- there is currently an error copying (some) structs
      // by value to the device in HIP-Clang / VDI
      // As a workaround, we can malloc the DriverType and explictly copy over.
      // To remove once solved in HIP
      DriverType *d_driver;
      HIP_SAFE_CALL(hipMalloc(&d_driver, sizeof(DriverType)));
      HIP_SAFE_CALL(hipMemcpyAsync(d_driver, &driver, sizeof(DriverType),
                                   hipMemcpyHostToDevice,
                                   hip_instance->m_stream));
      hip_parallel_launch_local_memory<DriverType, MaxThreadsPerBlock,
                                       MinBlocksPerSM>
          <<<grid, block, shmem, hip_instance->m_stream>>>(d_driver);

#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
      HIP_SAFE_CALL(hipGetLastError());
      hip_instance->fence();
#endif
      HIP_SAFE_CALL(hipFree(d_driver));
    }
  }

  static hipFuncAttributes get_hip_func_attributes() {
    static hipFuncAttributes attr = []() {
      hipFuncAttributes attr;
      HIP_SAFE_CALL(hipFuncGetAttributes(
          &attr,
          reinterpret_cast<void const *>(
              hip_parallel_launch_local_memory<DriverType, MaxThreadsPerBlock,
                                               MinBlocksPerSM>)));
      return attr;
    }();
    return attr;
  }
};

template <class DriverType>
struct HIPParallelLaunch<DriverType, Kokkos::LaunchBounds<0, 0>,
                         HIPLaunchMechanism::LocalMemory> {
  inline HIPParallelLaunch(const DriverType &driver, const dim3 &grid,
                           const dim3 &block, const int shmem,
                           const HIPInternal *hip_instance,
                           const bool /*prefer_shmem*/) {
    if ((grid.x != 0) && ((block.x * block.y * block.z) != 0)) {
      if (hip_instance->m_maxShmemPerBlock < shmem) {
        Kokkos::Impl::throw_runtime_exception(std::string(
            "HIPParallelLaunch FAILED: shared memory request is too large"));
      }

      KOKKOS_ENSURE_HIP_LOCK_ARRAYS_ON_DEVICE();

      // Invoke the driver function on the device

      // FIXME_HIP -- see note about struct copy by value above
      DriverType *d_driver;
      HIP_SAFE_CALL(hipMalloc(&d_driver, sizeof(DriverType)));
      HIP_SAFE_CALL(hipMemcpyAsync(d_driver, &driver, sizeof(DriverType),
                                   hipMemcpyHostToDevice,
                                   hip_instance->m_stream));
      hip_parallel_launch_local_memory<DriverType, 1024, 1>
          <<<grid, block, shmem, hip_instance->m_stream>>>(d_driver);

#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
      HIP_SAFE_CALL(hipGetLastError());
      hip_instance->fence();
#endif
      HIP_SAFE_CALL(hipFree(d_driver));
    }
  }

  static hipFuncAttributes get_hip_func_attributes() {
    static hipFuncAttributes attr = []() {
      hipFuncAttributes attr;
      HIP_SAFE_CALL(hipFuncGetAttributes(
          &attr, reinterpret_cast<void const *>(
                     hip_parallel_launch_local_memory<DriverType, 1024, 1>)));
      return attr;
    }();
    return attr;
  }
};
}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif

#endif
