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

template <typename DriverType, bool, int MaxThreadsPerBlock, int MinBlocksPerSM>
void hipOccupancy(int *numBlocks, int blockSize, int sharedmem) {
  // FIXME_HIP - currently the "constant" path is unimplemented.
  //             we should look at whether it's functional, and
  //             perform some simple scaling studies to see when /
  //             if the constant launcher outperforms the current
  //             pass by pointer shared launcher
  HIP_SAFE_CALL(hipOccupancyMaxActiveBlocksPerMultiprocessor(
      numBlocks,
      hip_parallel_launch_local_memory<DriverType, MaxThreadsPerBlock,
                                       MinBlocksPerSM>,
      blockSize, sharedmem));
}

template <typename DriverType, bool constant>
void hipOccupancy(int *numBlocks, int blockSize, int sharedmem) {
  hipOccupancy<DriverType, constant, HIPTraits::MaxThreadsPerBlock, 1>(
      numBlocks, blockSize, sharedmem);
}
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

template <class FunctorType, class LaunchBounds, typename F>
int hip_internal_get_block_size(const F &condition_check,
                                const HIPInternal *hip_instance,
                                const hipFuncAttributes &attr,
                                const FunctorType &f,
                                const size_t vector_length,
                                const size_t shmem_block,
                                const size_t shmem_thread) {
  const int min_blocks_per_sm =
      LaunchBounds::minBperSM == 0 ? 1 : LaunchBounds::minBperSM;
  const int max_threads_per_block = LaunchBounds::maxTperB == 0
                                        ? HIPTraits::MaxThreadsPerBlock
                                        : LaunchBounds::maxTperB;

  const int regs_per_wavefront  = std::max(attr.numRegs, 1);
  const int regs_per_sm         = hip_instance->m_regsPerSM;
  const int shmem_per_sm        = hip_instance->m_shmemPerSM;
  const int max_shmem_per_block = hip_instance->m_maxShmemPerBlock;
  const int max_blocks_per_sm   = hip_instance->m_maxBlocksPerSM;
  const int max_threads_per_sm  = hip_instance->m_maxThreadsPerSM;

  int block_size = max_threads_per_block;
  KOKKOS_ASSERT(block_size > 0);
  const int blocks_per_warp =
      (block_size + HIPTraits::WarpSize - 1) / HIPTraits::WarpSize;

  int functor_shmem = ::Kokkos::Impl::FunctorTeamShmemSize<FunctorType>::value(
      f, block_size / vector_length);
  int total_shmem = shmem_block + shmem_thread * (block_size / vector_length) +
                    functor_shmem + attr.sharedSizeBytes;
  int max_blocks_regs = regs_per_sm / (regs_per_wavefront * blocks_per_warp);
  int max_blocks_shmem =
      (total_shmem < max_shmem_per_block)
          ? (total_shmem > 0 ? shmem_per_sm / total_shmem : max_blocks_regs)
          : 0;
  int blocks_per_sm  = std::min(max_blocks_regs, max_blocks_shmem);
  int threads_per_sm = blocks_per_sm * block_size;
  if (threads_per_sm > max_threads_per_sm) {
    blocks_per_sm  = max_threads_per_sm / block_size;
    threads_per_sm = blocks_per_sm * block_size;
  }
  int opt_block_size =
      (blocks_per_sm >= min_blocks_per_sm) ? block_size : min_blocks_per_sm;
  int opt_threads_per_sm = threads_per_sm;
  // printf("BlockSizeMax: %i Shmem: %i %i %i %i Regs: %i %i Blocks: %i %i
  // Achieved: %i %i Opt: %i %i\n",block_size,
  //   shmem_per_sm,max_shmem_per_block,functor_shmem,total_shmem,
  //   regs_per_sm,regs_per_wavefront,max_blocks_shmem,max_blocks_regs,blocks_per_sm,threads_per_sm,opt_block_size,opt_threads_per_sm);
  block_size -= HIPTraits::WarpSize;
  while (condition_check(blocks_per_sm) &&
         (block_size >= HIPTraits::WarpSize)) {
    functor_shmem = ::Kokkos::Impl::FunctorTeamShmemSize<FunctorType>::value(
        f, block_size / vector_length);
    total_shmem = shmem_block + shmem_thread * (block_size / vector_length) +
                  functor_shmem + attr.sharedSizeBytes;
    max_blocks_regs = regs_per_sm / (regs_per_wavefront * blocks_per_warp);
    max_blocks_shmem =
        (total_shmem < max_shmem_per_block)
            ? (total_shmem > 0 ? shmem_per_sm / total_shmem : max_blocks_regs)
            : 0;
    blocks_per_sm  = std::min(max_blocks_regs, max_blocks_shmem);
    threads_per_sm = blocks_per_sm * block_size;
    if (threads_per_sm > max_threads_per_sm) {
      blocks_per_sm  = max_threads_per_sm / block_size;
      threads_per_sm = blocks_per_sm * block_size;
    }
    if ((blocks_per_sm >= min_blocks_per_sm) &&
        (blocks_per_sm <= max_blocks_per_sm)) {
      if (threads_per_sm >= opt_threads_per_sm) {
        opt_block_size     = block_size;
        opt_threads_per_sm = threads_per_sm;
      }
    }
    // printf("BlockSizeMax: %i Shmem: %i %i %i %i Regs: %i %i Blocks: %i %i
    // Achieved: %i %i Opt: %i %i\n",block_size,
    //   shmem_per_sm,max_shmem_per_block,functor_shmem,total_shmem,
    //   regs_per_sm,regs_per_wavefront,max_blocks_shmem,max_blocks_regs,blocks_per_sm,threads_per_sm,opt_block_size,opt_threads_per_sm);
    block_size -= HIPTraits::WarpSize;
  }
  return opt_block_size;
}

template <class FunctorType, class LaunchBounds>
int hip_get_max_block_size(const HIPInternal *hip_instance,
                           const hipFuncAttributes &attr, const FunctorType &f,
                           const size_t vector_length, const size_t shmem_block,
                           const size_t shmem_thread) {
  return hip_internal_get_block_size<FunctorType, LaunchBounds>(
      [](int x) { return x == 0; }, hip_instance, attr, f, vector_length,
      shmem_block, shmem_thread);
}
template <typename DriverType, class LaunchBounds>
struct HIPGetMaxBlockSize<DriverType, LaunchBounds, true> {
  static int get_block_size(typename DriverType::functor_type const &f,
                            size_t const vector_length,
                            size_t const shmem_extra_block,
                            size_t const shmem_extra_thread) {
    int numBlocks = 0;
    int blockSize = LaunchBounds::maxTperB == 0 ? 1024 : LaunchBounds::maxTperB;
    int sharedmem =
        shmem_extra_block + shmem_extra_thread * (blockSize / vector_length) +
        ::Kokkos::Impl::FunctorTeamShmemSize<
            typename DriverType::functor_type>::value(f, blockSize /
                                                             vector_length);

    hipOccupancy<DriverType, true>(&numBlocks, blockSize, sharedmem);

    if (numBlocks > 0) return blockSize;
    while (blockSize > HIPTraits::WarpSize && numBlocks == 0) {
      blockSize /= 2;
      sharedmem =
          shmem_extra_block + shmem_extra_thread * (blockSize / vector_length) +
          ::Kokkos::Impl::FunctorTeamShmemSize<
              typename DriverType::functor_type>::value(f, blockSize /
                                                               vector_length);

      hipOccupancy<DriverType, true>(&numBlocks, blockSize, sharedmem);
    }
    int blockSizeUpperBound = blockSize * 2;
    while (blockSize < blockSizeUpperBound && numBlocks > 0) {
      blockSize += HIPTraits::WarpSize;
      sharedmem =
          shmem_extra_block + shmem_extra_thread * (blockSize / vector_length) +
          ::Kokkos::Impl::FunctorTeamShmemSize<
              typename DriverType::functor_type>::value(f, blockSize /
                                                               vector_length);

      hipOccupancy<DriverType, true>(&numBlocks, blockSize, sharedmem);
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
int hip_get_opt_block_size(HIPInternal const *hip_instance,
                           hipFuncAttributes const &attr, FunctorType const &f,
                           size_t const vector_length, size_t const shmem_block,
                           size_t const shmem_thread) {
  return hip_internal_get_block_size<FunctorType, LaunchBounds>(
      [](int) { return true; }, hip_instance, attr, f, vector_length,
      shmem_block, shmem_thread);
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

    while (blockSize < HIPTraits::MaxThreadsPerBlock) {
      blockSize *= 2;

      // calculate the occupancy with that optBlockSize and check whether its
      // larger than the largest one found so far
      sharedmem =
          shmem_extra_block + shmem_extra_thread * (blockSize / vector_length) +
          ::Kokkos::Impl::FunctorTeamShmemSize<
              typename DriverType::functor_type>::value(f, blockSize /
                                                               vector_length);
      hipOccupancy<DriverType, true>(&numBlocks, blockSize, sharedmem);
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

    while (blockSize < HIPTraits::MaxThreadsPerBlock) {
      blockSize *= 2;
      sharedmem =
          shmem_extra_block + shmem_extra_thread * (blockSize / vector_length) +
          ::Kokkos::Impl::FunctorTeamShmemSize<
              typename DriverType::functor_type>::value(f, blockSize /
                                                               vector_length);

      hipOccupancy<DriverType, false>(&numBlocks, blockSize, sharedmem);

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
      hipOccupancy<DriverType, true, MaxThreadsPerBlock, MinBlocksPerSM>(
          &numBlocks, blockSize, sharedmem);
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

      hipOccupancy<DriverType, false, MaxThreadsPerBlock, MinBlocksPerSM>(
          &numBlocks, blockSize, sharedmem);
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
