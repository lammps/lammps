// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_HIP_BLOCKSIZE_DEDUCTION_HPP
#define KOKKOS_HIP_BLOCKSIZE_DEDUCTION_HPP

#include <functional>
#include <Kokkos_Macros.hpp>
#include <Kokkos_BitManipulation.hpp>

#if defined(__HIPCC__)

#include <HIP/Kokkos_HIP_Instance.hpp>
#include <HIP/Kokkos_HIP_KernelLaunch.hpp>

namespace Kokkos {
namespace Impl {

enum class BlockType { Max, Preferred };

template <typename DriverType, typename LaunchBounds = Kokkos::LaunchBounds<>,
          HIPLaunchMechanism LaunchMechanism =
              DeduceHIPLaunchMechanism<DriverType>::launch_mechanism>
unsigned get_preferred_blocksize_impl(const int hip_device) {
  if constexpr (!HIPParallelLaunch<DriverType, LaunchBounds,
                                   LaunchMechanism>::default_launchbounds()) {
    // use the user specified value
    return LaunchBounds::maxTperB;
  } else {
    if (HIPParallelLaunch<DriverType, LaunchBounds,
                          LaunchMechanism>::get_scratch_size(hip_device) > 0) {
      return HIPTraits::ConservativeThreadsPerBlock;
    }
    return HIPTraits::MaxThreadsPerBlock;
  }
}

template <typename DriverType, typename LaunchBounds = Kokkos::LaunchBounds<>,
          HIPLaunchMechanism LaunchMechanism =
              DeduceHIPLaunchMechanism<DriverType>::launch_mechanism>
constexpr unsigned get_max_blocksize_impl() {
  if constexpr (!HIPParallelLaunch<DriverType, LaunchBounds,
                                   LaunchMechanism>::default_launchbounds()) {
    // use the user specified value
    return LaunchBounds::maxTperB;
  } else {
    // we can always fit 1024 threads blocks if we only care about registers
    // ... and don't mind spilling
    return HIPTraits::MaxThreadsPerBlock;
  }
}

// convenience method to select and return the proper function attributes
// for a kernel, given the launch bounds et al.
template <typename DriverType, typename LaunchBounds = Kokkos::LaunchBounds<>,
          BlockType BlockSize = BlockType::Max,
          HIPLaunchMechanism LaunchMechanism =
              DeduceHIPLaunchMechanism<DriverType>::launch_mechanism>
hipFuncAttributes get_hip_func_attributes_impl(const int hip_device) {
#ifndef KOKKOS_ENABLE_HIP_MULTIPLE_KERNEL_INSTANTIATIONS
  return HIPParallelLaunch<DriverType, LaunchBounds, LaunchMechanism>::
      get_hip_func_attributes(hip_device);
#else
  if constexpr (!HIPParallelLaunch<DriverType, LaunchBounds,
                                   LaunchMechanism>::default_launchbounds()) {
    // for user defined, we *always* honor the request
    return HIPParallelLaunch<DriverType, LaunchBounds, LaunchMechanism>::
        get_hip_func_attributes(hip_device);
  } else {
    if constexpr (BlockSize == BlockType::Max) {
      return HIPParallelLaunch<
          DriverType, Kokkos::LaunchBounds<HIPTraits::MaxThreadsPerBlock, 1>,
          LaunchMechanism>::get_hip_func_attributes(hip_device);
    } else {
      const int blocksize =
          get_preferred_blocksize_impl<DriverType, LaunchBounds,
                                       LaunchMechanism>(hip_device);
      if (blocksize == HIPTraits::MaxThreadsPerBlock) {
        return HIPParallelLaunch<
            DriverType, Kokkos::LaunchBounds<HIPTraits::MaxThreadsPerBlock, 1>,
            LaunchMechanism>::get_hip_func_attributes(hip_device);
      } else {
        return HIPParallelLaunch<
            DriverType,
            Kokkos::LaunchBounds<HIPTraits::ConservativeThreadsPerBlock, 1>,
            LaunchMechanism>::get_hip_func_attributes(hip_device);
      }
    }
  }
#endif
}

// Given an initial block-size limitation based on register usage
// determine the block size to select based on LDS limitation
template <BlockType BlockSize, class DriverType, class LaunchBounds,
          typename ShmemFunctor>
unsigned hip_internal_get_block_size(const HIPInternal *hip_instance,
                                     const ShmemFunctor &f,
                                     const unsigned tperb_reg) {
  // translate LB from CUDA to HIP
  const unsigned min_waves_per_eu =
      LaunchBounds::minBperSM ? LaunchBounds::minBperSM : 1;
  const unsigned shmem_per_sm =
      hip_instance->m_deviceProp.maxSharedMemoryPerMultiProcessor;
  unsigned block_size     = tperb_reg;
  unsigned min_block_size = 0;
  do {
    unsigned total_shmem = f(block_size);
    // find how many threads we can fit with this blocksize based on LDS usage
    unsigned tperb_shmem = total_shmem > shmem_per_sm ? 0 : block_size;

    if constexpr (BlockSize == BlockType::Max) {
      // we want the maximum blocksize possible
      // just wait until we get a case where we can fit the LDS per SM
      if (tperb_shmem) return block_size;
    } else {
      // If total_shmem is zero, we set blocks_per_cu_shmem to a number greater
      // than min_waves_per_eu.
      const unsigned blocks_per_cu_shmem =
          total_shmem == 0 ? min_waves_per_eu + 1 : shmem_per_sm / total_shmem;
      const unsigned tperb = tperb_shmem < tperb_reg ? tperb_shmem : tperb_reg;

      // The logic prefers smaller blocks sizes over larger ones to give more
      // flexibility to the scheduler and to decrease the number of threads
      // launched when using Kokkos::AUTO in TeamPolicy. If the block size is
      // smaller than 256, fall back to BlockType::Max condition.
      if (blocks_per_cu_shmem > min_waves_per_eu &&
          tperb >= HIPTraits::ConservativeThreadsPerBlock) {
        min_block_size = block_size;
      } else if ((min_block_size == 0) && (tperb_shmem)) {
        return block_size;
      }
    }
    block_size >>= 1;
  } while (block_size >= HIPTraits::WarpSize);

  return min_block_size;
}

// Standardized blocksize deduction for parallel constructs with no LDS usage
// Returns the preferred blocksize as dictated by register usage
//
// Note: a returned block_size of zero indicates that the algorithm could not
//       find a valid block size.  The caller is responsible for error handling.
template <typename DriverType, typename LaunchBounds>
unsigned hip_get_preferred_blocksize(const int hip_device) {
  return get_preferred_blocksize_impl<DriverType, LaunchBounds>(hip_device);
}

// Heuristic to compute the block size for non-team parallelism
template <typename DriverType, typename LaunchBounds = Kokkos::LaunchBounds<>,
          HIPLaunchMechanism LaunchMechanism =
              DeduceHIPLaunchMechanism<DriverType>::launch_mechanism>
unsigned get_preferred_blocksize_for_range(HIPInternal const *hip_instance,
                                           size_t requested_parallelism) {
  /* General approach, if the user did not make a launch bounds request
  - If the requested parallelism is less than the available concurrency, get the
  largest block size that would result in at least 1 block per PE, while also:
    - at least 256
    - power of 2
    - no more than 1024
  */

  if constexpr (HIPParallelLaunch<DriverType, LaunchBounds,
                                  LaunchMechanism>::default_launchbounds()) {
    if (requested_parallelism &&
        requested_parallelism < size_t(hip_instance->concurrency())) {
      const unsigned eus = hip_instance->m_deviceProp.multiProcessorCount;
      const unsigned requestedPerEU = (requested_parallelism + eus - 1) / eus;
      // round up to power of 2
      unsigned threadsPerEU = Kokkos::bit_ceil(requestedPerEU);
      threadsPerEU          = std::max(threadsPerEU,
                                       unsigned(HIPTraits::ConservativeThreadsPerBlock));
      threadsPerEU =
          std::min(threadsPerEU, unsigned(HIPTraits::MaxThreadsPerBlock));
      return threadsPerEU;
    }
  }
  const int hip_device = hip_instance->m_hipDev;
  return get_preferred_blocksize_impl<DriverType, LaunchBounds>(hip_device);
}

// Standardized blocksize deduction for parallel constructs with no LDS usage
// Returns the max blocksize as dictated by register usage
//
// Note: a returned block_size of zero indicates that the algorithm could not
//       find a valid block size.  The caller is responsible for error handling.
template <typename DriverType, typename LaunchBounds>
unsigned hip_get_max_blocksize() {
  return get_max_blocksize_impl<DriverType, LaunchBounds>();
}

// Standardized blocksize deduction for non-teams parallel constructs with LDS
// usage Returns the 'preferred' blocksize, as determined by the heuristics in
// hip_internal_get_block_size
//
// The ShmemFunctor takes a single argument of the current blocksize under
// consideration, and returns the LDS usage
//
// requested_parallelism is a hint about how much parallelism was requested
// in the parallel construct
//
// Note: a returned block_size of zero indicates that the algorithm could not
//       find a valid block size.  The caller is responsible for error handling.
template <typename DriverType, typename LaunchBounds, typename ShmemFunctor>
unsigned hip_get_preferred_blocksize(HIPInternal const *hip_instance,
                                     ShmemFunctor const &f) {
  // get preferred blocksize limited by register usage
  const unsigned tperb_reg =
      hip_get_preferred_blocksize<DriverType, LaunchBounds>(
          hip_instance->m_hipDev);
  return hip_internal_get_block_size<BlockType::Preferred, DriverType,
                                     LaunchBounds>(hip_instance, f, tperb_reg);
}

// Standardized blocksize deduction for teams-based parallel constructs with LDS
// usage Returns the 'preferred' blocksize, as determined by the heuristics in
// hip_internal_get_block_size
//
// The ShmemTeamsFunctor takes two arguments: the hipFunctionAttributes and
//  the current blocksize under consideration, and returns the LDS usage
//
// Note: a returned block_size of zero indicates that the algorithm could not
//       find a valid block size.  The caller is responsible for error handling.
template <typename DriverType, typename LaunchBounds,
          typename ShmemTeamsFunctor>
unsigned hip_get_preferred_team_blocksize(HIPInternal const *hip_instance,
                                          ShmemTeamsFunctor const &f) {
  hipFuncAttributes attr = get_hip_func_attributes_impl<
      DriverType, LaunchBounds, BlockType::Preferred>(hip_instance->m_hipDev);
  // get preferred blocksize limited by register usage
  const unsigned tperb_reg =
      hip_get_preferred_blocksize<DriverType, LaunchBounds>(
          hip_instance->m_hipDev);
  return hip_internal_get_block_size<BlockType::Preferred, DriverType,
                                     LaunchBounds>(
      hip_instance, std::bind(f, attr, std::placeholders::_1), tperb_reg);
}

// Standardized blocksize deduction for non-teams parallel constructs with LDS
// usage Returns the maximum possible blocksize, as determined by the heuristics
// in hip_internal_get_block_size
//
// The ShmemFunctor takes a single argument of the current blocksize under
// consideration, and returns the LDS usage
//
// Note: a returned block_size of zero indicates that the algorithm could not
//       find a valid block size.  The caller is responsible for error handling.
template <typename DriverType, typename LaunchBounds, typename ShmemFunctor>
unsigned hip_get_max_blocksize(HIPInternal const *hip_instance,
                               ShmemFunctor const &f) {
  // get max blocksize limited by register usage
  const unsigned tperb_reg = hip_get_max_blocksize<DriverType, LaunchBounds>();
  return hip_internal_get_block_size<BlockType::Max, DriverType, LaunchBounds>(
      hip_instance, f, tperb_reg);
}

// Standardized blocksize deduction for teams-based parallel constructs with LDS
// usage Returns the maximum possible blocksize, as determined by the heuristics
// in hip_internal_get_block_size
//
// The ShmemTeamsFunctor takes two arguments: the hipFunctionAttributes and
//  the current blocksize under consideration, and returns the LDS usage
//
// Note: a returned block_size of zero indicates that the algorithm could not
//       find a valid block size.  The caller is responsible for error handling.
template <typename DriverType, typename LaunchBounds,
          typename ShmemTeamsFunctor>
unsigned hip_get_max_team_blocksize(HIPInternal const *hip_instance,
                                    ShmemTeamsFunctor const &f) {
  hipFuncAttributes attr =
      get_hip_func_attributes_impl<DriverType, LaunchBounds, BlockType::Max>(
          hip_instance->m_hipDev);
  // get max blocksize
  const unsigned tperb_reg = hip_get_max_blocksize<DriverType, LaunchBounds>();
  return hip_internal_get_block_size<BlockType::Max, DriverType, LaunchBounds>(
      hip_instance, std::bind(f, attr, std::placeholders::_1), tperb_reg);
}

}  // namespace Impl
}  // namespace Kokkos

#endif

#endif
