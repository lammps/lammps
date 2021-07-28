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

#ifndef KOKKOS_HIP_PARALLEL_MDRANGE_HPP
#define KOKKOS_HIP_PARALLEL_MDRANGE_HPP

#include <HIP/Kokkos_HIP_BlockSize_Deduction.hpp>
#include <HIP/Kokkos_HIP_KernelLaunch.hpp>
#include <HIP/Kokkos_HIP_ReduceScan.hpp>
#include <KokkosExp_MDRangePolicy.hpp>
#include <impl/KokkosExp_IterateTileGPU.hpp>
#include <Kokkos_Parallel.hpp>

namespace Kokkos {
namespace Impl {
// ParallelFor
template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::MDRangePolicy<Traits...>,
                  Kokkos::Experimental::HIP> {
 public:
  using Policy = Kokkos::MDRangePolicy<Traits...>;

 private:
  using array_index_type = typename Policy::array_index_type;
  using index_type       = typename Policy::index_type;
  using LaunchBounds     = typename Policy::launch_bounds;

  const FunctorType m_functor;
  const Policy m_policy;

  ParallelFor()        = delete;
  ParallelFor& operator=(ParallelFor const&) = delete;

 public:
  inline __device__ void operator()() const {
    Kokkos::Impl::DeviceIterateTile<Policy::rank, Policy, FunctorType,
                                    typename Policy::work_tag>(m_policy,
                                                               m_functor)
        .exec_range();
  }

  inline void execute() const {
    if (m_policy.m_num_tiles == 0) return;
    array_index_type const maxblocks = static_cast<array_index_type>(
        m_policy.space().impl_internal_space_instance()->m_maxBlock);
    if (Policy::rank == 2) {
      dim3 const block(m_policy.m_tile[0], m_policy.m_tile[1], 1);
      dim3 const grid(
          std::min((m_policy.m_upper[0] - m_policy.m_lower[0] + block.x - 1) /
                       block.x,
                   maxblocks),
          std::min((m_policy.m_upper[1] - m_policy.m_lower[1] + block.y - 1) /
                       block.y,
                   maxblocks),
          1);
      Kokkos::Experimental::Impl::HIPParallelLaunch<ParallelFor, LaunchBounds>(
          *this, grid, block, 0,
          m_policy.space().impl_internal_space_instance(), false);
    } else if (Policy::rank == 3) {
      dim3 const block(m_policy.m_tile[0], m_policy.m_tile[1],
                       m_policy.m_tile[2]);
      dim3 const grid(
          std::min((m_policy.m_upper[0] - m_policy.m_lower[0] + block.x - 1) /
                       block.x,
                   maxblocks),
          std::min((m_policy.m_upper[1] - m_policy.m_lower[1] + block.y - 1) /
                       block.y,
                   maxblocks),
          std::min((m_policy.m_upper[2] - m_policy.m_lower[2] + block.z - 1) /
                       block.z,
                   maxblocks));
      Kokkos::Experimental::Impl::HIPParallelLaunch<ParallelFor, LaunchBounds>(
          *this, grid, block, 0,
          m_policy.space().impl_internal_space_instance(), false);
    } else if (Policy::rank == 4) {
      // id0,id1 encoded within threadIdx.x; id2 to threadIdx.y; id3 to
      // threadIdx.z
      dim3 const block(m_policy.m_tile[0] * m_policy.m_tile[1],
                       m_policy.m_tile[2], m_policy.m_tile[3]);
      dim3 const grid(
          std::min(static_cast<uint32_t>(m_policy.m_tile_end[0] *
                                         m_policy.m_tile_end[1]),
                   static_cast<uint32_t>(maxblocks)),
          std::min((m_policy.m_upper[2] - m_policy.m_lower[2] + block.y - 1) /
                       block.y,
                   maxblocks),
          std::min((m_policy.m_upper[3] - m_policy.m_lower[3] + block.z - 1) /
                       block.z,
                   maxblocks));
      Kokkos::Experimental::Impl::HIPParallelLaunch<ParallelFor, LaunchBounds>(
          *this, grid, block, 0,
          m_policy.space().impl_internal_space_instance(), false);
    } else if (Policy::rank == 5) {
      // id0,id1 encoded within threadIdx.x; id2,id3 to threadIdx.y; id4
      // to threadIdx.z
      dim3 const block(m_policy.m_tile[0] * m_policy.m_tile[1],
                       m_policy.m_tile[2] * m_policy.m_tile[3],
                       m_policy.m_tile[4]);
      dim3 const grid(
          std::min(static_cast<index_type>(m_policy.m_tile_end[0] *
                                           m_policy.m_tile_end[1]),
                   static_cast<index_type>(maxblocks)),
          std::min(static_cast<index_type>(m_policy.m_tile_end[2] *
                                           m_policy.m_tile_end[3]),
                   static_cast<index_type>(maxblocks)),
          std::min((m_policy.m_upper[4] - m_policy.m_lower[4] + block.z - 1) /
                       block.z,
                   maxblocks));
      Kokkos::Experimental::Impl::HIPParallelLaunch<ParallelFor, LaunchBounds>(
          *this, grid, block, 0,
          m_policy.space().impl_internal_space_instance(), false);
    } else if (Policy::rank == 6) {
      // id0,id1 encoded within threadIdx.x; id2,id3 to threadIdx.y;
      // id4,id5 to threadIdx.z
      dim3 const block(m_policy.m_tile[0] * m_policy.m_tile[1],
                       m_policy.m_tile[2] * m_policy.m_tile[3],
                       m_policy.m_tile[4] * m_policy.m_tile[5]);
      dim3 const grid(std::min(static_cast<index_type>(m_policy.m_tile_end[0] *
                                                       m_policy.m_tile_end[1]),
                               static_cast<index_type>(maxblocks)),
                      std::min(static_cast<index_type>(m_policy.m_tile_end[2] *
                                                       m_policy.m_tile_end[3]),
                               static_cast<index_type>(maxblocks)),
                      std::min(static_cast<index_type>(m_policy.m_tile_end[4] *
                                                       m_policy.m_tile_end[5]),
                               static_cast<index_type>(maxblocks)));
      Kokkos::Experimental::Impl::HIPParallelLaunch<ParallelFor, LaunchBounds>(
          *this, grid, block, 0,
          m_policy.space().impl_internal_space_instance(), false);
    } else {
      Kokkos::abort("Kokkos::MDRange Error: Exceeded rank bounds with HIP\n");
    }

  }  // end execute

  ParallelFor(FunctorType const& arg_functor, Policy const& arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}

  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy& pol, const Functor&) {
    using closure_type =
        ParallelFor<FunctorType, Kokkos::MDRangePolicy<Traits...>,
                    Kokkos::Experimental::HIP>;
    hipFuncAttributes attr = Kokkos::Experimental::Impl::HIPParallelLaunch<
        closure_type, LaunchBounds>::get_hip_func_attributes();
    auto const& prop = pol.space().hip_device_prop();
    // Limits due to registers/SM, MDRange doesn't have
    // shared memory constraints
    int const regs_per_sm        = prop.regsPerMultiprocessor;
    int const regs_per_thread    = attr.numRegs;
    int const max_threads_per_sm = regs_per_sm / regs_per_thread;
    return std::min(
        max_threads_per_sm,
        static_cast<int>(
            Kokkos::Experimental::Impl::HIPTraits::MaxThreadsPerBlock));
  }
};

// ParallelReduce
template <class FunctorType, class ReducerType, class... Traits>
class ParallelReduce<FunctorType, Kokkos::MDRangePolicy<Traits...>, ReducerType,
                     Kokkos::Experimental::HIP> {
 public:
  using Policy = Kokkos::MDRangePolicy<Traits...>;

 private:
  using array_index_type = typename Policy::array_index_type;
  using index_type       = typename Policy::index_type;

  using WorkTag      = typename Policy::work_tag;
  using Member       = typename Policy::member_type;
  using LaunchBounds = typename Policy::launch_bounds;

  using ReducerConditional =
      Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                         FunctorType, ReducerType>;
  using ReducerTypeFwd = typename ReducerConditional::type;
  using WorkTagFwd =
      typename Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                                  WorkTag, void>::type;

  using ValueTraits =
      Kokkos::Impl::FunctorValueTraits<ReducerTypeFwd, WorkTagFwd>;
  using ValueInit = Kokkos::Impl::FunctorValueInit<ReducerTypeFwd, WorkTagFwd>;
  using ValueJoin = Kokkos::Impl::FunctorValueJoin<ReducerTypeFwd, WorkTagFwd>;

 public:
  using pointer_type   = typename ValueTraits::pointer_type;
  using value_type     = typename ValueTraits::value_type;
  using reference_type = typename ValueTraits::reference_type;
  using functor_type   = FunctorType;
  using size_type      = Experimental::HIP::size_type;

  // Algorithmic constraints: blockSize is a power of two AND blockDim.y ==
  // blockDim.z == 1

  const FunctorType m_functor;
  const Policy m_policy;  // used for workrange and nwork
  const ReducerType m_reducer;
  const pointer_type m_result_ptr;
  const bool m_result_ptr_device_accessible;
  size_type* m_scratch_space;
  size_type* m_scratch_flags;

  using DeviceIteratePattern = typename Kokkos::Impl::Reduce::DeviceIterateTile<
      Policy::rank, Policy, FunctorType, WorkTag, reference_type>;

 public:
  inline __device__ void exec_range(reference_type update) const {
    DeviceIteratePattern(m_policy, m_functor, update).exec_range();
  }

  inline __device__ void operator()() const {
    const integral_nonzero_constant<size_type, ValueTraits::StaticValueSize /
                                                   sizeof(size_type)>
        word_count(ValueTraits::value_size(
                       ReducerConditional::select(m_functor, m_reducer)) /
                   sizeof(size_type));

    {
      reference_type value = ValueInit::init(
          ReducerConditional::select(m_functor, m_reducer),
          Experimental::kokkos_impl_hip_shared_memory<size_type>() +
              threadIdx.y * word_count.value);

      // Number of blocks is bounded so that the reduction can be limited to two
      // passes. Each thread block is given an approximately equal amount of
      // work to perform. Accumulate the values for this block. The accumulation
      // ordering does not match the final pass, but is arithmatically
      // equivalent.

      this->exec_range(value);
    }

    // Reduce with final value at blockDim.y - 1 location.
    // Problem: non power-of-two blockDim
    if (::Kokkos::Impl::hip_single_inter_block_reduce_scan<
            false, ReducerTypeFwd, WorkTagFwd>(
            ReducerConditional::select(m_functor, m_reducer), blockIdx.x,
            gridDim.x, Experimental::kokkos_impl_hip_shared_memory<size_type>(),
            m_scratch_space, m_scratch_flags)) {
      // This is the final block with the final result at the final threads'
      // location
      size_type* const shared =
          Experimental::kokkos_impl_hip_shared_memory<size_type>() +
          (blockDim.y - 1) * word_count.value;
      size_type* const global = m_result_ptr_device_accessible
                                    ? reinterpret_cast<size_type*>(m_result_ptr)
                                    : m_scratch_space;

      if (threadIdx.y == 0) {
        Kokkos::Impl::FunctorFinal<ReducerTypeFwd, WorkTagFwd>::final(
            ReducerConditional::select(m_functor, m_reducer), shared);
      }

      if (Experimental::Impl::HIPTraits::WarpSize < word_count.value) {
        __syncthreads();
      }

      for (unsigned i = threadIdx.y; i < word_count.value; i += blockDim.y) {
        global[i] = shared[i];
      }
    }
  }

  // Determine block size constrained by shared memory:
  // This is copy/paste from Kokkos_HIP_Parallel_Range
  inline unsigned local_block_size(const FunctorType& f) {
    unsigned int n =
        ::Kokkos::Experimental::Impl::HIPTraits::MaxThreadsPerBlock;
    int shmem_size = ::Kokkos::Impl::hip_single_inter_block_reduce_scan_shmem<
        false, FunctorType, WorkTag>(f, n);
    using closure_type = Impl::ParallelReduce<FunctorType, Policy, ReducerType>;
    hipFuncAttributes attr = ::Kokkos::Experimental::Impl::HIPParallelLaunch<
        closure_type, LaunchBounds>::get_hip_func_attributes();
    while (
        (n &&
         (m_policy.space().impl_internal_space_instance()->m_maxShmemPerBlock <
          shmem_size)) ||
        (n >
         static_cast<unsigned>(
             ::Kokkos::Experimental::Impl::hip_get_max_block_size<FunctorType,
                                                                  LaunchBounds>(
                 m_policy.space().impl_internal_space_instance(), attr, f, 1,
                 shmem_size, 0)))) {
      n >>= 1;
      shmem_size = ::Kokkos::Impl::hip_single_inter_block_reduce_scan_shmem<
          false, FunctorType, WorkTag>(f, n);
    }
    return n;
  }

  inline void execute() {
    const int nwork = m_policy.m_num_tiles;
    if (nwork) {
      int block_size = m_policy.m_prod_tile_dims;
      // CONSTRAINT: Algorithm requires block_size >= product of tile dimensions
      // Nearest power of two
      int exponent_pow_two    = std::ceil(std::log2(block_size));
      block_size              = std::pow(2, exponent_pow_two);
      int suggested_blocksize = local_block_size(m_functor);

      block_size = (block_size > suggested_blocksize)
                       ? block_size
                       : suggested_blocksize;  // Note: block_size must be less
                                               // than or equal to 512

      m_scratch_space =
          ::Kokkos::Experimental::Impl::hip_internal_scratch_space(
              ValueTraits::value_size(
                  ReducerConditional::select(m_functor, m_reducer)) *
              block_size /* block_size == max block_count */);
      m_scratch_flags =
          ::Kokkos::Experimental::Impl::hip_internal_scratch_flags(
              sizeof(size_type));

      // REQUIRED ( 1 , N , 1 )
      const dim3 block(1, block_size, 1);
      // Required grid.x <= block.y
      const dim3 grid(std::min(static_cast<uint32_t>(block.y),
                               static_cast<uint32_t>(nwork)),
                      1, 1);

      const int shmem =
          ::Kokkos::Impl::hip_single_inter_block_reduce_scan_shmem<
              false, FunctorType, WorkTag>(m_functor, block.y);

      Kokkos::Experimental::Impl::HIPParallelLaunch<ParallelReduce,
                                                    LaunchBounds>(
          *this, grid, block, shmem,
          m_policy.space().impl_internal_space_instance(),
          false);  // copy to device and execute

      if (!m_result_ptr_device_accessible) {
        m_policy.space().fence();

        if (m_result_ptr) {
          const int size = ValueTraits::value_size(
              ReducerConditional::select(m_functor, m_reducer));
          DeepCopy<HostSpace, Experimental::HIPSpace>(m_result_ptr,
                                                      m_scratch_space, size);
        }
      }
    } else {
      if (m_result_ptr) {
        ValueInit::init(ReducerConditional::select(m_functor, m_reducer),
                        m_result_ptr);
      }
    }
  }

  template <class ViewType>
  ParallelReduce(const FunctorType& arg_functor, const Policy& arg_policy,
                 const ViewType& arg_result,
                 typename std::enable_if<Kokkos::is_view<ViewType>::value,
                                         void*>::type = nullptr)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(InvalidType()),
        m_result_ptr(arg_result.data()),
        m_result_ptr_device_accessible(
            MemorySpaceAccess<Kokkos::Experimental::HIPSpace,
                              typename ViewType::memory_space>::accessible),
        m_scratch_space(nullptr),
        m_scratch_flags(nullptr) {}

  ParallelReduce(const FunctorType& arg_functor, const Policy& arg_policy,
                 const ReducerType& reducer)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(reducer),
        m_result_ptr(reducer.view().data()),
        m_result_ptr_device_accessible(
            MemorySpaceAccess<Kokkos::Experimental::HIPSpace,
                              typename ReducerType::result_view_type::
                                  memory_space>::accessible),
        m_scratch_space(nullptr),
        m_scratch_flags(nullptr) {}
  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy& pol, const Functor&) {
    using closure_type =
        ParallelReduce<FunctorType, Kokkos::MDRangePolicy<Traits...>,
                       ReducerType, Kokkos::Experimental::HIP>;
    hipFuncAttributes attr = Kokkos::Experimental::Impl::HIPParallelLaunch<
        closure_type, LaunchBounds>::get_hip_func_attributes();
    auto const& prop = pol.space().hip_device_prop();
    // Limits due do registers/SM
    int const regs_per_sm        = prop.regsPerMultiprocessor;
    int const regs_per_thread    = attr.numRegs;
    int const max_threads_per_sm = regs_per_sm / regs_per_thread;
    return std::min(
        max_threads_per_sm,
        static_cast<int>(
            Kokkos::Experimental::Impl::HIPTraits::MaxThreadsPerBlock));
  }
};
}  // namespace Impl
}  // namespace Kokkos

#endif
