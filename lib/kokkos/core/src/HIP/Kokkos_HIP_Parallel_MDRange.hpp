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
    using ClosureType =
        ParallelFor<FunctorType, Policy, Kokkos::Experimental::HIP>;
    if (m_policy.m_num_tiles == 0) return;
    auto const maxblocks =
        Kokkos::Experimental::Impl::hip_internal_maximum_grid_count();
    if (Policy::rank == 2) {
      dim3 const block(m_policy.m_tile[0], m_policy.m_tile[1], 1);
      dim3 const grid(
          std::min<array_index_type>(
              (m_policy.m_upper[0] - m_policy.m_lower[0] + block.x - 1) /
                  block.x,
              maxblocks[0]),
          std::min<array_index_type>(
              (m_policy.m_upper[1] - m_policy.m_lower[1] + block.y - 1) /
                  block.y,
              maxblocks[1]),
          1);
      Kokkos::Experimental::Impl::hip_parallel_launch<ClosureType,
                                                      LaunchBounds>(
          *this, grid, block, 0,
          m_policy.space().impl_internal_space_instance(), false);
    } else if (Policy::rank == 3) {
      dim3 const block(m_policy.m_tile[0], m_policy.m_tile[1],
                       m_policy.m_tile[2]);
      dim3 const grid(
          std::min<array_index_type>(
              (m_policy.m_upper[0] - m_policy.m_lower[0] + block.x - 1) /
                  block.x,
              maxblocks[0]),
          std::min<array_index_type>(
              (m_policy.m_upper[1] - m_policy.m_lower[1] + block.y - 1) /
                  block.y,
              maxblocks[1]),
          std::min<array_index_type>(
              (m_policy.m_upper[2] - m_policy.m_lower[2] + block.z - 1) /
                  block.z,
              maxblocks[2]));
      Kokkos::Experimental::Impl::hip_parallel_launch<ClosureType,
                                                      LaunchBounds>(
          *this, grid, block, 0,
          m_policy.space().impl_internal_space_instance(), false);
    } else if (Policy::rank == 4) {
      // id0,id1 encoded within threadIdx.x; id2 to threadIdx.y; id3 to
      // threadIdx.z
      dim3 const block(m_policy.m_tile[0] * m_policy.m_tile[1],
                       m_policy.m_tile[2], m_policy.m_tile[3]);
      dim3 const grid(
          std::min<array_index_type>(
              m_policy.m_tile_end[0] * m_policy.m_tile_end[1], maxblocks[0]),
          std::min<array_index_type>(
              (m_policy.m_upper[2] - m_policy.m_lower[2] + block.y - 1) /
                  block.y,
              maxblocks[1]),
          std::min<array_index_type>(
              (m_policy.m_upper[3] - m_policy.m_lower[3] + block.z - 1) /
                  block.z,
              maxblocks[2]));
      Kokkos::Experimental::Impl::hip_parallel_launch<ClosureType,
                                                      LaunchBounds>(
          *this, grid, block, 0,
          m_policy.space().impl_internal_space_instance(), false);
    } else if (Policy::rank == 5) {
      // id0,id1 encoded within threadIdx.x; id2,id3 to threadIdx.y; id4
      // to threadIdx.z
      dim3 const block(m_policy.m_tile[0] * m_policy.m_tile[1],
                       m_policy.m_tile[2] * m_policy.m_tile[3],
                       m_policy.m_tile[4]);
      dim3 const grid(
          std::min<array_index_type>(
              m_policy.m_tile_end[0] * m_policy.m_tile_end[1], maxblocks[0]),
          std::min<array_index_type>(
              m_policy.m_tile_end[2] * m_policy.m_tile_end[3], maxblocks[1]),
          std::min<array_index_type>(
              (m_policy.m_upper[4] - m_policy.m_lower[4] + block.z - 1) /
                  block.z,
              maxblocks[2]));
      Kokkos::Experimental::Impl::hip_parallel_launch<ClosureType,
                                                      LaunchBounds>(
          *this, grid, block, 0,
          m_policy.space().impl_internal_space_instance(), false);
    } else if (Policy::rank == 6) {
      // id0,id1 encoded within threadIdx.x; id2,id3 to threadIdx.y;
      // id4,id5 to threadIdx.z
      dim3 const block(m_policy.m_tile[0] * m_policy.m_tile[1],
                       m_policy.m_tile[2] * m_policy.m_tile[3],
                       m_policy.m_tile[4] * m_policy.m_tile[5]);
      dim3 const grid(
          std::min<array_index_type>(
              m_policy.m_tile_end[0] * m_policy.m_tile_end[1], maxblocks[0]),
          std::min<array_index_type>(
              m_policy.m_tile_end[2] * m_policy.m_tile_end[3], maxblocks[1]),
          std::min<array_index_type>(
              m_policy.m_tile_end[4] * m_policy.m_tile_end[5], maxblocks[2]));
      Kokkos::Experimental::Impl::hip_parallel_launch<ClosureType,
                                                      LaunchBounds>(
          *this, grid, block, 0,
          m_policy.space().impl_internal_space_instance(), false);
    } else {
      Kokkos::abort("Kokkos::MDRange Error: Exceeded rank bounds with HIP\n");
    }

  }  // end execute

  ParallelFor(FunctorType const& arg_functor, Policy const& arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}

  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy&, const Functor&) {
    using closure_type =
        ParallelFor<FunctorType, Kokkos::MDRangePolicy<Traits...>,
                    Kokkos::Experimental::HIP>;
    unsigned block_size =
        Kokkos::Experimental::Impl::hip_get_max_blocksize<closure_type,
                                                          LaunchBounds>();
    if (block_size == 0)
      Kokkos::Impl::throw_runtime_exception(
          std::string("Kokkos::Impl::ParallelFor< HIP > could not find a valid "
                      "tile size."));
    return block_size;
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

  using Analysis =
      Kokkos::Impl::FunctorAnalysis<FunctorPatternInterface::REDUCE, Policy,
                                    ReducerTypeFwd>;

 public:
  using pointer_type   = typename Analysis::pointer_type;
  using value_type     = typename Analysis::value_type;
  using reference_type = typename Analysis::reference_type;
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
  // Only let one Parallel/Scan modify the shared memory. The
  // constructor acquires the mutex which is released in the destructor.
  std::lock_guard<std::mutex> m_shared_memory_lock;

  using DeviceIteratePattern = typename Kokkos::Impl::Reduce::DeviceIterateTile<
      Policy::rank, Policy, FunctorType, WorkTag, reference_type>;

 public:
  inline __device__ void exec_range(reference_type update) const {
    DeviceIteratePattern(m_policy, m_functor, update).exec_range();
  }

  inline __device__ void operator()() const {
    typename Analysis::Reducer final_reducer(
        &ReducerConditional::select(m_functor, m_reducer));

    const integral_nonzero_constant<size_type, Analysis::StaticValueSize /
                                                   sizeof(size_type)>
        word_count(Analysis::value_size(
                       ReducerConditional::select(m_functor, m_reducer)) /
                   sizeof(size_type));

    {
      reference_type value = final_reducer.init(reinterpret_cast<pointer_type>(
          Experimental::kokkos_impl_hip_shared_memory<size_type>() +
          threadIdx.y * word_count.value));

      // Number of blocks is bounded so that the reduction can be limited to two
      // passes. Each thread block is given an approximately equal amount of
      // work to perform. Accumulate the values for this block. The accumulation
      // ordering does not match the final pass, but is arithmatically
      // equivalent.

      this->exec_range(value);
    }

    // Reduce with final value at blockDim.y - 1 location.
    // Problem: non power-of-two blockDim
    if (::Kokkos::Impl::hip_single_inter_block_reduce_scan<false>(
            final_reducer, blockIdx.x, gridDim.x,
            Experimental::kokkos_impl_hip_shared_memory<size_type>(),
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
        final_reducer.final(reinterpret_cast<value_type*>(shared));
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
    const auto& instance = m_policy.space().impl_internal_space_instance();
    auto shmem_functor   = [&f](unsigned n) {
      return hip_single_inter_block_reduce_scan_shmem<false, FunctorType,
                                                      WorkTag>(f, n);
    };
    using closure_type = ParallelReduce<FunctorType, Policy, ReducerType,
                                        Kokkos::Experimental::HIP>;

    unsigned block_size =
        Kokkos::Experimental::Impl::hip_get_preferred_blocksize<closure_type,
                                                                LaunchBounds>(
            instance, shmem_functor);
    if (block_size == 0) {
      Kokkos::Impl::throw_runtime_exception(
          std::string("Kokkos::Impl::ParallelReduce< HIP > could not find a "
                      "valid tile size."));
    }
    return block_size;
  }

  inline void execute() {
    typename Analysis::Reducer final_reducer(
        &ReducerConditional::select(m_functor, m_reducer));

    using ClosureType = ParallelReduce<FunctorType, Policy, ReducerType,
                                       Kokkos::Experimental::HIP>;
    const auto nwork  = m_policy.m_num_tiles;
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
              m_policy.space(),
              Analysis::value_size(
                  ReducerConditional::select(m_functor, m_reducer)) *
                  block_size /* block_size == max block_count */);
      m_scratch_flags =
          ::Kokkos::Experimental::Impl::hip_internal_scratch_flags(
              m_policy.space(), sizeof(size_type));

      // REQUIRED ( 1 , N , 1 )
      const dim3 block(1, block_size, 1);
      // Required grid.x <= block.y
      const dim3 grid(std::min(static_cast<uint32_t>(block.y),
                               static_cast<uint32_t>(nwork)),
                      1, 1);

      const int shmem =
          ::Kokkos::Impl::hip_single_inter_block_reduce_scan_shmem<
              false, FunctorType, WorkTag>(m_functor, block.y);

      Kokkos::Experimental::Impl::hip_parallel_launch<ClosureType,
                                                      LaunchBounds>(
          *this, grid, block, shmem,
          m_policy.space().impl_internal_space_instance(),
          false);  // copy to device and execute

      if (!m_result_ptr_device_accessible && m_result_ptr) {
        const int size = Analysis::value_size(
            ReducerConditional::select(m_functor, m_reducer));
        DeepCopy<HostSpace, Experimental::HIPSpace, Experimental::HIP>(
            m_policy.space(), m_result_ptr, m_scratch_space, size);
      }
    } else {
      if (m_result_ptr) {
        final_reducer.init(m_result_ptr);
      }
    }
  }

  template <class ViewType>
  ParallelReduce(
      const FunctorType& arg_functor, const Policy& arg_policy,
      const ViewType& arg_result,
      std::enable_if_t<Kokkos::is_view<ViewType>::value, void*> = nullptr)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(InvalidType()),
        m_result_ptr(arg_result.data()),
        m_result_ptr_device_accessible(
            MemorySpaceAccess<Kokkos::Experimental::HIPSpace,
                              typename ViewType::memory_space>::accessible),
        m_scratch_space(nullptr),
        m_scratch_flags(nullptr),
        m_shared_memory_lock(m_policy.space()
                                 .impl_internal_space_instance()
                                 ->m_mutexSharedMemory) {}

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
        m_scratch_flags(nullptr),
        m_shared_memory_lock(m_policy.space()
                                 .impl_internal_space_instance()
                                 ->m_mutexSharedMemory) {}

  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy&, const Functor&) {
    using closure_type =
        ParallelReduce<FunctorType, Kokkos::MDRangePolicy<Traits...>,
                       ReducerType, Kokkos::Experimental::HIP>;
    unsigned block_size =
        Kokkos::Experimental::Impl::hip_get_max_blocksize<closure_type,
                                                          LaunchBounds>();
    if (block_size == 0) {
      Kokkos::Impl::throw_runtime_exception(
          std::string("Kokkos::Impl::ParallelReduce< HIP > could not find a "
                      "valid tile size."));
    }
    return block_size;
  }
};
}  // namespace Impl
}  // namespace Kokkos

#endif
