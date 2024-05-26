//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_HIP_PARALLEL_REDUCE_RANGE_HPP
#define KOKKOS_HIP_PARALLEL_REDUCE_RANGE_HPP

#include <Kokkos_Parallel.hpp>

#include <HIP/Kokkos_HIP_BlockSize_Deduction.hpp>
#include <HIP/Kokkos_HIP_KernelLaunch.hpp>
#include <HIP/Kokkos_HIP_ReduceScan.hpp>
#include <HIP/Kokkos_HIP_Shuffle_Reduce.hpp>

namespace Kokkos {
namespace Impl {

template <class CombinedFunctorReducerType, class... Traits>
class ParallelReduce<CombinedFunctorReducerType, Kokkos::RangePolicy<Traits...>,
                     Kokkos::HIP> {
 public:
  using Policy      = Kokkos::RangePolicy<Traits...>;
  using FunctorType = typename CombinedFunctorReducerType::functor_type;
  using ReducerType = typename CombinedFunctorReducerType::reducer_type;

 private:
  using WorkRange    = typename Policy::WorkRange;
  using WorkTag      = typename Policy::work_tag;
  using Member       = typename Policy::member_type;
  using LaunchBounds = typename Policy::launch_bounds;

 public:
  using pointer_type   = typename ReducerType::pointer_type;
  using value_type     = typename ReducerType::value_type;
  using reference_type = typename ReducerType::reference_type;
  using functor_type   = FunctorType;
  using reducer_type   = ReducerType;
  using size_type      = Kokkos::HIP::size_type;
  using index_type     = typename Policy::index_type;
  // Conditionally set word_size_type to int16_t or int8_t if value_type is
  // smaller than int32_t (Kokkos::HIP::size_type)
  // word_size_type is used to determine the word count, shared memory buffer
  // size, and global memory buffer size before the scan is performed.
  // Within the scan, the word count is recomputed based on word_size_type
  // and when calculating indexes into the shared/global memory buffers for
  // performing the scan, word_size_type is used again.
  // For scalars > 4 bytes in size, indexing into shared/global memory relies
  // on the block and grid dimensions to ensure that we index at the correct
  // offset rather than at every 4 byte word; such that, when the join is
  // performed, we have the correct data that was copied over in chunks of 4
  // bytes.
  using word_size_type = std::conditional_t<
      sizeof(value_type) < sizeof(size_type),
      std::conditional_t<sizeof(value_type) == 2, int16_t, int8_t>, size_type>;

  // Algorithmic constraints: blockSize is a power of two AND blockDim.y ==
  // blockDim.z == 1

  const CombinedFunctorReducerType m_functor_reducer;
  const Policy m_policy;
  const pointer_type m_result_ptr;
  const bool m_result_ptr_device_accessible;
  const bool m_result_ptr_host_accessible;
  word_size_type* m_scratch_space = nullptr;
  size_type* m_scratch_flags      = nullptr;

  static constexpr bool UseShflReduction = false;

 private:
  struct ShflReductionTag {};
  struct SHMEMReductionTag {};

  // Make the exec_range calls call to Reduce::DeviceIterateTile
  template <class TagType>
  __device__ inline std::enable_if_t<std::is_void<TagType>::value> exec_range(
      const Member& i, reference_type update) const {
    m_functor_reducer.get_functor()(i, update);
  }

  template <class TagType>
  __device__ inline std::enable_if_t<!std::is_void<TagType>::value> exec_range(
      const Member& i, reference_type update) const {
    m_functor_reducer.get_functor()(TagType(), i, update);
  }

 public:
  __device__ inline void operator()() const {
    using ReductionTag = std::conditional_t<UseShflReduction, ShflReductionTag,
                                            SHMEMReductionTag>;
    run(ReductionTag{});
  }

  __device__ inline void run(SHMEMReductionTag) const {
    const ReducerType& reducer = m_functor_reducer.get_reducer();
    const integral_nonzero_constant<word_size_type,
                                    ReducerType::static_value_size() /
                                        sizeof(word_size_type)>
        word_count(reducer.value_size() / sizeof(word_size_type));

    {
      reference_type value = reducer.init(reinterpret_cast<pointer_type>(
          ::Kokkos::kokkos_impl_hip_shared_memory<word_size_type>() +
          threadIdx.y * word_count.value));

      // Number of blocks is bounded so that the reduction can be limited to two
      // passes. Each thread block is given an approximately equal amount of
      // work to perform. Accumulate the values for this block. The accumulation
      // ordering does not match the final pass, but is arithmetically
      // equivalent.

      const WorkRange range(m_policy, blockIdx.x, gridDim.x);

      for (Member iwork = range.begin() + threadIdx.y, iwork_end = range.end();
           iwork < iwork_end; iwork += blockDim.y) {
        this->template exec_range<WorkTag>(iwork, value);
      }
    }

    // Reduce with final value at blockDim.y - 1 location.
    // Shortcut for length zero reduction
    bool do_final_reduction = m_policy.begin() == m_policy.end();
    if (!do_final_reduction)
      do_final_reduction = hip_single_inter_block_reduce_scan<false>(
          reducer, blockIdx.x, gridDim.x,
          ::Kokkos::kokkos_impl_hip_shared_memory<word_size_type>(),
          m_scratch_space, m_scratch_flags);
    if (do_final_reduction) {
      // This is the final block with the final result at the final threads'
      // location

      word_size_type* const shared =
          ::Kokkos::kokkos_impl_hip_shared_memory<word_size_type>() +
          (blockDim.y - 1) * word_count.value;
      word_size_type* const global =
          m_result_ptr_device_accessible
              ? reinterpret_cast<word_size_type*>(m_result_ptr)
              : m_scratch_space;

      if (threadIdx.y == 0) {
        reducer.final(reinterpret_cast<value_type*>(shared));
      }

      if (::Kokkos::Impl::HIPTraits::WarpSize < word_count.value) {
        __syncthreads();
      }

      for (unsigned i = threadIdx.y; i < word_count.value; i += blockDim.y) {
        global[i] = shared[i];
      }
    }
  }

  __device__ inline void run(ShflReductionTag) const {
    const ReducerType& reducer = m_functor_reducer.get_reducer();

    value_type value;
    reducer.init(&value);
    // Number of blocks is bounded so that the reduction can be limited to two
    // passes. Each thread block is given an approximately equal amount of work
    // to perform. Accumulate the values for this block. The accumulation
    // ordering does not match the final pass, but is arithmetically equivalent.

    WorkRange const range(m_policy, blockIdx.x, gridDim.x);

    for (Member iwork = range.begin() + threadIdx.y, iwork_end = range.end();
         iwork < iwork_end; iwork += blockDim.y) {
      this->template exec_range<WorkTag>(iwork, value);
    }

    pointer_type const result = reinterpret_cast<pointer_type>(m_scratch_space);

    int max_active_thread = static_cast<int>(range.end() - range.begin()) <
                                    static_cast<int>(blockDim.y)
                                ? range.end() - range.begin()
                                : blockDim.y;

    max_active_thread =
        (max_active_thread == 0) ? blockDim.y : max_active_thread;

    value_type init;
    reducer.init(&init);
    if (m_policy.begin() == m_policy.end()) {
      reducer.final(&value);
      pointer_type const final_result =
          m_result_ptr_device_accessible ? m_result_ptr : result;
      *final_result = value;
    } else if (Impl::hip_inter_block_shuffle_reduction<>(
                   value, init, reducer, m_scratch_space, result,
                   m_scratch_flags, max_active_thread)) {
      unsigned int const id = threadIdx.y * blockDim.x + threadIdx.x;
      if (id == 0) {
        reducer.final(&value);
        pointer_type const final_result =
            m_result_ptr_device_accessible ? m_result_ptr : result;
        *final_result = value;
      }
    }
  }

  // Determine block size constrained by shared memory:
  inline unsigned local_block_size(const FunctorType& f) {
    const auto& instance = m_policy.space().impl_internal_space_instance();
    auto shmem_functor   = [&f](unsigned n) {
      return hip_single_inter_block_reduce_scan_shmem<false, WorkTag,
                                                      value_type>(f, n);
    };
    return Kokkos::Impl::hip_get_preferred_blocksize<ParallelReduce,
                                                     LaunchBounds>(
        instance, shmem_functor);
  }

  inline void execute() {
    const ReducerType& reducer = m_functor_reducer.get_reducer();

    const index_type nwork     = m_policy.end() - m_policy.begin();
    const bool need_device_set = ReducerType::has_init_member_function() ||
                                 ReducerType::has_final_member_function() ||
                                 !m_result_ptr_host_accessible ||
                                 !std::is_same<ReducerType, InvalidType>::value;
    if ((nwork > 0) || need_device_set) {
      const int block_size = local_block_size(m_functor_reducer.get_functor());
      if (block_size == 0) {
        Kokkos::Impl::throw_runtime_exception(
            std::string("Kokkos::Impl::ParallelReduce< HIP > could not find a "
                        "valid execution configuration."));
      }

      // REQUIRED ( 1 , N , 1 )
      dim3 block(1, block_size, 1);
      // use a slightly less constrained, but still well bounded limit for
      // scratch
      int nblocks = (nwork + block.y - 1) / block.y;
      // Heuristic deciding the value of nblocks.
      // The general idea here is we want to:
      //    1. Not undersubscribe the device (i.e., we want at least
      //    preferred_block_min blocks)
      //    2. Have each thread reduce > 1 value to minimize overheads
      //    3. Limit the total # of blocks, to avoid unbounded scratch space
      constexpr int block_max           = 4096;
      constexpr int preferred_block_min = 1024;

      if (nblocks < preferred_block_min) {
        // keep blocks as is, already have low parallelism
      } else if (nblocks > block_max) {
        // "large dispatch" -> already have lots of parallelism
        nblocks = block_max;
      } else {
        // in the intermediate range, try to have each thread process multiple
        // items to offset the cost of the reduction (with not enough
        // parallelism to hide it)
        int items_per_thread =
            (nwork + nblocks * block_size - 1) / (nblocks * block_size);
        if (items_per_thread < 4) {
          int ratio = std::min(
              (nblocks + preferred_block_min - 1) / preferred_block_min,
              (4 + items_per_thread - 1) / items_per_thread);
          nblocks /= ratio;
        }
      }

      // TODO: down casting these uses more space than required?
      m_scratch_space =
          (word_size_type*)::Kokkos::Impl::hip_internal_scratch_space(
              m_policy.space(), reducer.value_size() * nblocks);
      // Intentionally do not downcast to word_size_type since we use HIP
      // atomics in Kokkos_HIP_ReduceScan.hpp
      m_scratch_flags = ::Kokkos::Impl::hip_internal_scratch_flags(
          m_policy.space(), sizeof(size_type));
      // Required grid.x <= block.y
      dim3 grid(nblocks, 1, 1);

      if (nwork == 0) {
        block = dim3(1, 1, 1);
        grid  = dim3(1, 1, 1);
      }
      const int shmem =
          UseShflReduction
              ? 0
              : hip_single_inter_block_reduce_scan_shmem<false, WorkTag,
                                                         value_type>(
                    m_functor_reducer.get_functor(), block.y);

      Kokkos::Impl::hip_parallel_launch<ParallelReduce, LaunchBounds>(
          *this, grid, block, shmem,
          m_policy.space().impl_internal_space_instance(),
          false);  // copy to device and execute

      if (!m_result_ptr_device_accessible && m_result_ptr) {
        const int size = reducer.value_size();
        DeepCopy<HostSpace, HIPSpace, HIP>(m_policy.space(), m_result_ptr,
                                           m_scratch_space, size);
      }
    } else {
      if (m_result_ptr) {
        reducer.init(m_result_ptr);
      }
    }
  }

  template <class ViewType>
  ParallelReduce(const CombinedFunctorReducerType& arg_functor_reducer,
                 const Policy& arg_policy, const ViewType& arg_result)
      : m_functor_reducer(arg_functor_reducer),
        m_policy(arg_policy),
        m_result_ptr(arg_result.data()),
        m_result_ptr_device_accessible(
            MemorySpaceAccess<HIPSpace,
                              typename ViewType::memory_space>::accessible),
        m_result_ptr_host_accessible(
            MemorySpaceAccess<Kokkos::HostSpace,
                              typename ViewType::memory_space>::accessible) {}
};

}  // namespace Impl
}  // namespace Kokkos

#endif
