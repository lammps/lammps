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

// FIXME_HIP cannot use global variable on the device with ROCm 2.9
//__device__ __constant__ unsigned long kokkos_impl_hip_constant_memory_buffer
//    [Kokkos::Experimental::Impl::HIPTraits::ConstantMemoryUsage /
//     sizeof(unsigned long)];

namespace Kokkos {
namespace Experimental {
template <typename T>
inline __device__ T *kokkos_impl_hip_shared_memory() {
  extern __shared__ HIPSpace::size_type sh[];
  return (T *)sh;
}
}  // namespace Experimental
}  // namespace Kokkos

namespace Kokkos {
namespace Experimental {
namespace Impl {

void *hip_resize_scratch_space(std::int64_t bytes, bool force_shrink = false);

template <typename DriverType>
__global__ static void hip_parallel_launch_constant_memory() {
  __device__ __constant__ unsigned long kokkos_impl_hip_constant_memory_buffer
      [Kokkos::Experimental::Impl::HIPTraits::ConstantMemoryUsage /
       sizeof(unsigned long)];

  const DriverType &driver = *(reinterpret_cast<const DriverType *>(
      kokkos_impl_hip_constant_memory_buffer));

  driver();
}

template <class DriverType>
__global__ static void hip_parallel_launch_local_memory(
    const DriverType driver) {
  driver();
}

template <class DriverType, unsigned int maxTperB, unsigned int minBperSM>
__global__ __launch_bounds__(
    maxTperB,
    minBperSM) static void hip_parallel_launch_local_memory(const DriverType
                                                                driver) {
  driver();
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

      // Invoke the driver function on the device
      printf("%i %i %i | %i %i %i | %i\n", grid.x, grid.y, grid.z, block.x,
             block.y, block.z, shmem);
      printf("Pre Launch Error: %s\n", hipGetErrorName(hipGetLastError()));

      hipLaunchKernelGGL(
          (hip_parallel_launch_local_memory<DriverType, MaxThreadsPerBlock,
                                            MinBlocksPerSM>),
          grid, block, shmem, hip_instance->m_stream, driver);

      Kokkos::Experimental::HIP().fence();
      printf("Post Launch Error: %s\n", hipGetErrorName(hipGetLastError()));
#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
      HIP_SAFE_CALL(hipGetLastError());
      Kokkos::Experimental::HIP().fence();
#endif
    }
  }

  static hipFuncAttributes get_hip_func_attributes() {
    hipFuncAttributes attr;
    hipFuncGetAttributes(
        &attr, hip_parallel_launch_local_memory<DriverType, MaxThreadsPerBlock,
                                                MinBlocksPerSM>);
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

      // Invoke the driver function on the device
      hipLaunchKernelGGL(hip_parallel_launch_local_memory<DriverType>, grid,
                         block, shmem, hip_instance->m_stream, driver);

      Kokkos::Experimental::HIP().fence();
#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
      HIP_SAFE_CALL(hipGetLastError());
      Kokkos::Experimental::HIP().fence();
#endif
    }
  }

  static hipFuncAttributes get_hip_func_attributes() {
    hipFuncAttributes attr;
    hipFuncGetAttributes(&attr,
                         reinterpret_cast<void *>(
                             &hip_parallel_launch_local_memory<DriverType>));
    return attr;
  }
};
}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif

#endif
