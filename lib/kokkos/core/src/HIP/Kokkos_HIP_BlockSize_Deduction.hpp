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

#ifndef KOKKOS_HIP_BLOCKSIZE_DEDUCTION_HPP
#define KOKKOS_HIP_BLOCKSIZE_DEDUCTION_HPP

#include <Kokkos_Macros.hpp>

#if defined(__HIPCC__)

#include <HIP/Kokkos_HIP_Instance.hpp>
#include <HIP/Kokkos_HIP_KernelLaunch.hpp>

namespace Kokkos {
namespace Experimental {
namespace Impl {
template <typename DriverType, typename LaunchBounds, bool Large>
struct HIPGetMaxBlockSize;

template <typename DriverType, typename LaunchBounds>
int hip_get_max_block_size(typename DriverType::functor_type const &f,
                           size_t const vector_length,
                           size_t const shmem_extra_block,
                           size_t const shmem_extra_thread) {
  return HIPGetMaxBlockSize<DriverType, LaunchBounds, true>::get_block_size(
      f, vector_length, shmem_extra_block, shmem_extra_thread);
}

template <class FunctorType, class LaunchBounds>
int hip_get_max_block_size(const HIPInternal * /*hip_instance*/,
                           const hipFuncAttributes &attr,
                           const FunctorType & /*f*/,
                           const size_t /*vector_length*/,
                           const size_t /*shmem_block*/,
                           const size_t /*shmem_thread*/) {
  // FIXME_HIP find a better algorithm. Be aware that
  // maxThreadsPerMultiProcessor, regsPerBlock, and l2CacheSize are bugged and
  // always return zero
  // https://github.com/ROCm-Developer-Tools/HIP/blob/6c5fa32815650cc20a4f783d09b013610348a4d5/include/hip/hcc_detail/hip_runtime_api.h#L438-L440
  // and we don't have access to the same information than we do for CUDA

  int const max_threads_per_block_mi60 = 1024;
  int const max_threads_per_block      = LaunchBounds::maxTperB == 0
                                        ? max_threads_per_block_mi60
                                        : LaunchBounds::maxTperB;
  return std::min(attr.maxThreadsPerBlock, max_threads_per_block);
}

template <typename DriverType>
struct HIPGetMaxBlockSize<DriverType, Kokkos::LaunchBounds<>, true> {
  static int get_block_size(typename DriverType::functor_type const &f,
                            size_t const vector_length,
                            size_t const shmem_extra_block,
                            size_t const shmem_extra_thread) {
    unsigned int numBlocks = 0;
    int blockSize          = 1024;
    int sharedmem =
        shmem_extra_block + shmem_extra_thread * (blockSize / vector_length) +
        ::Kokkos::Impl::FunctorTeamShmemSize<
            typename DriverType::functor_type>::value(f, blockSize /
                                                             vector_length);
    hipOccupancyMaxActiveBlocksPerMultiprocessor(
        &numBlocks, hip_parallel_launch_constant_memory<DriverType>, blockSize,
        sharedmem);

    if (numBlocks > 0) return blockSize;
    while (blockSize > HIPTraits::WarpSize && numBlocks == 0) {
      blockSize /= 2;
      sharedmem =
          shmem_extra_block + shmem_extra_thread * (blockSize / vector_length) +
          ::Kokkos::Impl::FunctorTeamShmemSize<
              typename DriverType::functor_type>::value(f, blockSize /
                                                               vector_length);

      hipOccupancyMaxActiveBlocksPerMultiprocessor(
          &numBlocks, hip_parallel_launch_constant_memory<DriverType>,
          blockSize, sharedmem);
    }
    int blockSizeUpperBound = blockSize * 2;
    while (blockSize < blockSizeUpperBound && numBlocks > 0) {
      blockSize += HIPTraits::WarpSize;
      sharedmem =
          shmem_extra_block + shmem_extra_thread * (blockSize / vector_length) +
          ::Kokkos::Impl::FunctorTeamShmemSize<
              typename DriverType::functor_type>::value(f, blockSize /
                                                               vector_length);

      hipOccupancyMaxActiveBlocksPerMultiprocessor(
          &numBlocks, hip_parallel_launch_constant_memory<DriverType>,
          blockSize, sharedmem);
    }
    return blockSize - HIPTraits::WarpSize;
  }
};

template <typename DriverType, typename LaunchBounds, bool Large>
struct HIPGetOptBlockSize;

template <typename DriverType, typename LaunchBounds>
int hip_get_opt_block_size(typename DriverType::functor_type const &f,
                           size_t const vector_length,
                           size_t const shmem_extra_block,
                           size_t const shmem_extra_thread) {
  return HIPGetOptBlockSize<
      DriverType, LaunchBounds,
      (HIPTraits::ConstantMemoryUseThreshold <
       sizeof(DriverType))>::get_block_size(f, vector_length, shmem_extra_block,
                                            shmem_extra_thread);
}

template <typename FunctorType, typename LaunchBounds>
int hip_get_opt_block_size(HIPInternal const * /*hip_instance*/,
                           hipFuncAttributes const &attr,
                           FunctorType const & /*f*/,
                           size_t const /*vector_length*/,
                           size_t const /*shmem_block*/,
                           size_t const /*shmem_thread*/) {
  // FIXME_HIP find a better algorithm. Be aware that
  // maxThreadsPerMultiProcessor, regsPerBlock, and l2CacheSize are bugged and
  // always return zero
  // https://github.com/ROCm-Developer-Tools/HIP/blob/6c5fa32815650cc20a4f783d09b013610348a4d5/include/hip/hcc_detail/hip_runtime_api.h#L438-L440
  // and we don't have access to the same information than we do for CUDA

  int const max_threads_per_block_mi60 = 1024;
  int const max_threads_per_block      = LaunchBounds::maxTperB == 0
                                        ? max_threads_per_block_mi60
                                        : LaunchBounds::maxTperB;
  return std::min(attr.maxThreadsPerBlock, max_threads_per_block);
}

// FIXME_HIP the code is identical to the false struct except for
// hip_parallel_launch_constant_memory
template <typename DriverType>
struct HIPGetOptBlockSize<DriverType, Kokkos::LaunchBounds<0, 0>, true> {
  static int get_block_size(typename DriverType::functor_type const &f,
                            size_t const vector_length,
                            size_t const shmem_extra_block,
                            size_t const shmem_extra_thread) {
    int blockSize = HIPTraits::WarpSize / 2;
    int numBlocks;
    int sharedmem;
    int maxOccupancy  = 0;
    int bestBlockSize = 0;

    while (blockSize < 1024) {
      blockSize *= 2;

      // calculate the occupancy with that optBlockSize and check whether its
      // larger than the largest one found so far
      sharedmem =
          shmem_extra_block + shmem_extra_thread * (blockSize / vector_length) +
          ::Kokkos::Impl::FunctorTeamShmemSize<
              typename DriverType::functor_type>::value(f, blockSize /
                                                               vector_length);
      hipOccupancyMaxActiveBlocksPerMultiprocessor(
          &numBlocks, hip_parallel_launch_constant_memory<DriverType>,
          blockSize, sharedmem);
      if (maxOccupancy < numBlocks * blockSize) {
        maxOccupancy  = numBlocks * blockSize;
        bestBlockSize = blockSize;
      }
    }
    return bestBlockSize;
  }
};

template <typename DriverType>
struct HIPGetOptBlockSize<DriverType, Kokkos::LaunchBounds<0, 0>, false> {
  static int get_block_size(const typename DriverType::functor_type &f,
                            const size_t vector_length,
                            const size_t shmem_extra_block,
                            const size_t shmem_extra_thread) {
    int blockSize = HIPTraits::WarpSize / 2;
    int numBlocks;
    int sharedmem;
    int maxOccupancy  = 0;
    int bestBlockSize = 0;

    while (blockSize < 1024) {
      blockSize *= 2;
      sharedmem =
          shmem_extra_block + shmem_extra_thread * (blockSize / vector_length) +
          ::Kokkos::Impl::FunctorTeamShmemSize<
              typename DriverType::functor_type>::value(f, blockSize /
                                                               vector_length);

      hipOccupancyMaxActiveBlocksPerMultiprocessor(
          &numBlocks, hip_parallel_launch_local_memory<DriverType>, blockSize,
          sharedmem);

      if (maxOccupancy < numBlocks * blockSize) {
        maxOccupancy  = numBlocks * blockSize;
        bestBlockSize = blockSize;
      }
    }
    return bestBlockSize;
  }
};

// FIXME_HIP the code is identical to the false struct except for
// hip_parallel_launch_constant_memory
template <typename DriverType, unsigned int MaxThreadsPerBlock,
          unsigned int MinBlocksPerSM>
struct HIPGetOptBlockSize<
    DriverType, Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>,
    true> {
  static int get_block_size(const typename DriverType::functor_type &f,
                            const size_t vector_length,
                            const size_t shmem_extra_block,
                            const size_t shmem_extra_thread) {
    int blockSize = HIPTraits::WarpSize / 2;
    int numBlocks;
    int sharedmem;
    int maxOccupancy  = 0;
    int bestBlockSize = 0;
    int max_threads_per_block =
        std::min(MaxThreadsPerBlock,
                 hip_internal_maximum_warp_count() * HIPTraits::WarpSize);

    while (blockSize < max_threads_per_block) {
      blockSize *= 2;

      // calculate the occupancy with that optBlockSize and check whether its
      // larger than the largest one found so far
      sharedmem =
          shmem_extra_block + shmem_extra_thread * (blockSize / vector_length) +
          ::Kokkos::Impl::FunctorTeamShmemSize<
              typename DriverType::functor_type>::value(f, blockSize /
                                                               vector_length);
      hipOccupancyMaxActiveBlocksPerMultiprocessor(
          &numBlocks,
          hip_parallel_launch_constant_memory<DriverType, MaxThreadsPerBlock,
                                              MinBlocksPerSM>,
          blockSize, sharedmem);
      if (numBlocks >= static_cast<int>(MinBlocksPerSM) &&
          blockSize <= static_cast<int>(MaxThreadsPerBlock)) {
        if (maxOccupancy < numBlocks * blockSize) {
          maxOccupancy  = numBlocks * blockSize;
          bestBlockSize = blockSize;
        }
      }
    }
    if (maxOccupancy > 0) return bestBlockSize;
    return -1;
  }
};

template <typename DriverType, unsigned int MaxThreadsPerBlock,
          unsigned int MinBlocksPerSM>
struct HIPGetOptBlockSize<
    DriverType, Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>,
    false> {
  static int get_block_size(const typename DriverType::functor_type &f,
                            const size_t vector_length,
                            const size_t shmem_extra_block,
                            const size_t shmem_extra_thread) {
    int blockSize = HIPTraits::WarpSize / 2;
    int numBlocks;
    int sharedmem;
    int maxOccupancy  = 0;
    int bestBlockSize = 0;
    int max_threads_per_block =
        std::min(MaxThreadsPerBlock,
                 hip_internal_maximum_warp_count() * HIPTraits::WarpSize);

    while (blockSize < max_threads_per_block) {
      blockSize *= 2;
      sharedmem =
          shmem_extra_block + shmem_extra_thread * (blockSize / vector_length) +
          ::Kokkos::Impl::FunctorTeamShmemSize<
              typename DriverType::functor_type>::value(f, blockSize /
                                                               vector_length);

      hipOccupancyMaxActiveBlocksPerMultiprocessor(
          &numBlocks,
          hip_parallel_launch_local_memory<DriverType, MaxThreadsPerBlock,
                                           MinBlocksPerSM>,
          blockSize, sharedmem);
      if (numBlocks >= int(MinBlocksPerSM) &&
          blockSize <= int(MaxThreadsPerBlock)) {
        if (maxOccupancy < numBlocks * blockSize) {
          maxOccupancy  = numBlocks * blockSize;
          bestBlockSize = blockSize;
        }
      }
    }
    if (maxOccupancy > 0) return bestBlockSize;
    return -1;
  }
};

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif

#endif
