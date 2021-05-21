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

#ifndef KOKKOS_CUDA_INTERNAL_HPP
#define KOKKOS_CUDA_INTERNAL_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_CUDA

#include <Cuda/Kokkos_Cuda_Error.hpp>

namespace Kokkos {
namespace Impl {

inline int cuda_max_active_blocks_per_sm(cudaDeviceProp const& properties,
                                         cudaFuncAttributes const& attributes,
                                         int block_size, size_t dynamic_shmem) {
  // Limits due do registers/SM
  int const regs_per_sm     = properties.regsPerMultiprocessor;
  int const regs_per_thread = attributes.numRegs;
  int const max_blocks_regs = regs_per_sm / (regs_per_thread * block_size);

  // Limits due to shared memory/SM
  size_t const shmem_per_sm            = properties.sharedMemPerMultiprocessor;
  size_t const shmem_per_block         = properties.sharedMemPerBlock;
  size_t const static_shmem            = attributes.sharedSizeBytes;
  size_t const dynamic_shmem_per_block = attributes.maxDynamicSharedSizeBytes;
  size_t const total_shmem             = static_shmem + dynamic_shmem;

  int const max_blocks_shmem =
      total_shmem > shmem_per_block || dynamic_shmem > dynamic_shmem_per_block
          ? 0
          : (total_shmem > 0 ? (int)shmem_per_sm / total_shmem
                             : max_blocks_regs);

  // Limits due to blocks/SM
#if CUDA_VERSION >= 11000
  int const max_blocks_per_sm = properties.maxBlocksPerMultiProcessor;
#else
  int const max_blocks_per_sm = [&properties]() {
    switch (properties.major) {
      case 3: return 16;
      case 5:
      case 6: return 32;
      case 7: {
        int isTuring = properties.minor == 5;
        return (isTuring) ? 16 : 32;
      }
      default:
        throw_runtime_exception("Unknown device in cuda block size deduction");
        return 0;
    }
  }();
#endif

  // Overall occupancy in blocks
  return std::min({max_blocks_regs, max_blocks_shmem, max_blocks_per_sm});
}

template <typename UnaryFunction, typename LaunchBounds>
inline int cuda_deduce_block_size(bool early_termination,
                                  cudaDeviceProp const& properties,
                                  cudaFuncAttributes const& attributes,
                                  UnaryFunction block_size_to_dynamic_shmem,
                                  LaunchBounds) {
  // Limits
  int const max_threads_per_sm = properties.maxThreadsPerMultiProcessor;
  // unsure if I need to do that or if this is already accounted for in the
  // functor attributes
  int const max_threads_per_block =
      std::min(LaunchBounds::maxTperB == 0 ? (int)properties.maxThreadsPerBlock
                                           : (int)LaunchBounds::maxTperB,
               attributes.maxThreadsPerBlock);
  int const min_blocks_per_sm =
      LaunchBounds::minBperSM == 0 ? 1 : LaunchBounds::minBperSM;

  // Recorded maximum
  int opt_block_size     = 0;
  int opt_threads_per_sm = 0;

  for (int block_size = max_threads_per_block; block_size > 0;
       block_size -= 32) {
    size_t const dynamic_shmem = block_size_to_dynamic_shmem(block_size);

    int blocks_per_sm = cuda_max_active_blocks_per_sm(
        properties, attributes, block_size, dynamic_shmem);

    int threads_per_sm = blocks_per_sm * block_size;

    if (threads_per_sm > max_threads_per_sm) {
      blocks_per_sm  = max_threads_per_sm / block_size;
      threads_per_sm = blocks_per_sm * block_size;
    }

    if (blocks_per_sm >= min_blocks_per_sm) {
      if (threads_per_sm >= opt_threads_per_sm) {
        opt_block_size     = block_size;
        opt_threads_per_sm = threads_per_sm;
      }
    }

    if (early_termination && blocks_per_sm != 0) break;
  }

  return opt_block_size;
}

template <class FunctorType, class LaunchBounds>
int cuda_get_max_block_size(const CudaInternal* cuda_instance,
                            const cudaFuncAttributes& attr,
                            const FunctorType& f, const size_t vector_length,
                            const size_t shmem_block,
                            const size_t shmem_thread) {
  (void)cuda_instance;

  auto const& prop = Kokkos::Cuda().cuda_device_prop();

  auto const block_size_to_dynamic_shmem = [&f, vector_length, shmem_block,
                                            shmem_thread](int block_size) {
    size_t const functor_shmem =
        Kokkos::Impl::FunctorTeamShmemSize<FunctorType>::value(
            f, block_size / vector_length);

    size_t const dynamic_shmem = shmem_block +
                                 shmem_thread * (block_size / vector_length) +
                                 functor_shmem;
    return dynamic_shmem;
  };

  return cuda_deduce_block_size(true, prop, attr, block_size_to_dynamic_shmem,
                                LaunchBounds{});
}

template <class FunctorType, class LaunchBounds>
int cuda_get_opt_block_size(const CudaInternal* cuda_instance,
                            const cudaFuncAttributes& attr,
                            const FunctorType& f, const size_t vector_length,
                            const size_t shmem_block,
                            const size_t shmem_thread) {
  (void)cuda_instance;

  auto const& prop = Kokkos::Cuda().cuda_device_prop();

  auto const block_size_to_dynamic_shmem = [&f, vector_length, shmem_block,
                                            shmem_thread](int block_size) {
    size_t const functor_shmem =
        Kokkos::Impl::FunctorTeamShmemSize<FunctorType>::value(
            f, block_size / vector_length);

    size_t const dynamic_shmem = shmem_block +
                                 shmem_thread * (block_size / vector_length) +
                                 functor_shmem;
    return dynamic_shmem;
  };

  return cuda_deduce_block_size(false, prop, attr, block_size_to_dynamic_shmem,
                                LaunchBounds{});
}

// Assuming cudaFuncSetCacheConfig(MyKernel, cudaFuncCachePreferL1)
// NOTE these number can be obtained several ways:
// * One option is to download the CUDA Occupancy Calculator spreadsheet, select
// "Compute Capability" first and check what is the smallest "Shared Memory
// Size Config" that is available.  The "Shared Memory Per Multiprocessor" in
// bytes is then to be found below in the summary.
// * Another option would be to look for the information in the "Tuning
// Guide(s)" of the CUDA Toolkit Documentation for each GPU architecture, in
// the "Shared Memory" section (more tedious)
inline size_t get_shmem_per_sm_prefer_l1(cudaDeviceProp const& properties) {
  int const compute_capability = properties.major * 10 + properties.minor;
  return [compute_capability]() {
    switch (compute_capability) {
      case 30:
      case 32:
      case 35: return 16;
      case 37: return 80;
      case 50:
      case 53:
      case 60:
      case 62: return 64;
      case 52:
      case 61: return 96;
      case 70:
      case 80: return 8;
      case 75: return 32;
      default:
        Kokkos::Impl::throw_runtime_exception(
            "Unknown device in cuda block size deduction");
    }
    return 0;
  }() * 1024;
}
}  // namespace Impl
}  // namespace Kokkos

#endif  // KOKKOS_ENABLE_CUDA
#endif  /* #ifndef KOKKOS_CUDA_INTERNAL_HPP */
