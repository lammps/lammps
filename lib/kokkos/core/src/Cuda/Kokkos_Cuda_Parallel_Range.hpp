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

#ifndef KOKKOS_CUDA_PARALLEL_RANGE_HPP
#define KOKKOS_CUDA_PARALLEL_RANGE_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_CUDA)

#include <algorithm>
#include <string>

#include <Kokkos_Parallel.hpp>

#include <Cuda/Kokkos_Cuda_KernelLaunch.hpp>
#include <Cuda/Kokkos_Cuda_ReduceScan.hpp>
#include <Cuda/Kokkos_Cuda_BlockSize_Deduction.hpp>
#include <Kokkos_MinMaxClamp.hpp>

#include <impl/Kokkos_Tools.hpp>
#include <typeinfo>

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::RangePolicy<Traits...>, Kokkos::Cuda> {
 public:
  using Policy = Kokkos::RangePolicy<Traits...>;

 private:
  using Member       = typename Policy::member_type;
  using WorkTag      = typename Policy::work_tag;
  using LaunchBounds = typename Policy::launch_bounds;

  const FunctorType m_functor;
  const Policy m_policy;

  ParallelFor()        = delete;
  ParallelFor& operator=(const ParallelFor&) = delete;

  template <class TagType>
  inline __device__ std::enable_if_t<std::is_void<TagType>::value> exec_range(
      const Member i) const {
    m_functor(i);
  }

  template <class TagType>
  inline __device__ std::enable_if_t<!std::is_void<TagType>::value> exec_range(
      const Member i) const {
    m_functor(TagType(), i);
  }

 public:
  using functor_type = FunctorType;

  Policy const& get_policy() const { return m_policy; }

  inline __device__ void operator()() const {
    const Member work_stride = blockDim.y * gridDim.x;
    const Member work_end    = m_policy.end();

    for (Member iwork =
             m_policy.begin() + threadIdx.y + blockDim.y * blockIdx.x;
         iwork < work_end;
         iwork = iwork < work_end - work_stride ? iwork + work_stride
                                                : work_end) {
      this->template exec_range<WorkTag>(iwork);
    }
  }

  inline void execute() const {
    const typename Policy::index_type nwork = m_policy.end() - m_policy.begin();

    cudaFuncAttributes attr =
        CudaParallelLaunch<ParallelFor,
                           LaunchBounds>::get_cuda_func_attributes();
    const int block_size =
        Kokkos::Impl::cuda_get_opt_block_size<FunctorType, LaunchBounds>(
            m_policy.space().impl_internal_space_instance(), attr, m_functor, 1,
            0, 0);
    KOKKOS_ASSERT(block_size > 0);
    dim3 block(1, block_size, 1);
    dim3 grid(
        std::min(
            typename Policy::index_type((nwork + block.y - 1) / block.y),
            typename Policy::index_type(cuda_internal_maximum_grid_count()[0])),
        1, 1);
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
    if (Kokkos::Impl::CudaInternal::cuda_use_serial_execution()) {
      block = dim3(1, 1, 1);
      grid  = dim3(1, 1, 1);
    }
#endif

    CudaParallelLaunch<ParallelFor, LaunchBounds>(
        *this, grid, block, 0, m_policy.space().impl_internal_space_instance());
  }

  ParallelFor(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

template <class CombinedFunctorReducerType, class... Traits>
class ParallelReduce<CombinedFunctorReducerType, Kokkos::RangePolicy<Traits...>,
                     Kokkos::Cuda> {
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
  // Conditionally set word_size_type to int16_t or int8_t if value_type is
  // smaller than int32_t (Kokkos::Cuda::size_type)
  // word_size_type is used to determine the word count, shared memory buffer
  // size, and global memory buffer size before the reduction is performed.
  // Within the reduction, the word count is recomputed based on word_size_type
  // and when calculating indexes into the shared/global memory buffers for
  // performing the reduction, word_size_type is used again.
  // For scalars > 4 bytes in size, indexing into shared/global memory relies
  // on the block and grid dimensions to ensure that we index at the correct
  // offset rather than at every 4 byte word; such that, when the join is
  // performed, we have the correct data that was copied over in chunks of 4
  // bytes.
  using word_size_type = std::conditional_t<
      sizeof(value_type) < sizeof(Kokkos::Cuda::size_type),
      std::conditional_t<sizeof(value_type) == 2, int16_t, int8_t>,
      Kokkos::Cuda::size_type>;
  using index_type   = typename Policy::index_type;
  using reducer_type = ReducerType;

  // Algorithmic constraints: blockSize is a power of two AND blockDim.y ==
  // blockDim.z == 1

  const CombinedFunctorReducerType m_functor_reducer;
  const Policy m_policy;
  const pointer_type m_result_ptr;
  const bool m_result_ptr_device_accessible;
  const bool m_result_ptr_host_accessible;
  word_size_type* m_scratch_space;
  // m_scratch_flags must be of type Cuda::size_type due to use of atomics
  // for tracking metadata in Kokkos_Cuda_ReduceScan.hpp
  Cuda::size_type* m_scratch_flags;
  word_size_type* m_unified_space;

  // FIXME_CUDA Shall we use the shfl based reduction or not (only use it for
  // static sized types of more than 128bit:
  // sizeof(value_type)>2*sizeof(double)) && ReducerType::static_value_size())
  static constexpr bool UseShflReduction = false;

 public:
  Policy const& get_policy() const { return m_policy; }

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

  __device__ inline void operator()() const {
    const integral_nonzero_constant<word_size_type,
                                    ReducerType::static_value_size() /
                                        sizeof(word_size_type)>
        word_count(m_functor_reducer.get_reducer().value_size() /
                   sizeof(word_size_type));

    {
      reference_type value =
          m_functor_reducer.get_reducer().init(reinterpret_cast<pointer_type>(
              kokkos_impl_cuda_shared_memory<word_size_type>() +
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
    bool zero_length        = m_policy.begin() == m_policy.end();
    bool do_final_reduction = true;
    if (!zero_length)
      do_final_reduction = cuda_single_inter_block_reduce_scan<false>(
          m_functor_reducer.get_reducer(), blockIdx.x, gridDim.x,
          kokkos_impl_cuda_shared_memory<word_size_type>(), m_scratch_space,
          m_scratch_flags);

    if (do_final_reduction) {
      // This is the final block with the final result at the final threads'
      // location

      word_size_type* const shared =
          kokkos_impl_cuda_shared_memory<word_size_type>() +
          (blockDim.y - 1) * word_count.value;
      word_size_type* const global =
          m_result_ptr_device_accessible
              ? reinterpret_cast<word_size_type*>(m_result_ptr)
              : (m_unified_space ? m_unified_space : m_scratch_space);

      if (threadIdx.y == 0) {
        m_functor_reducer.get_reducer().final(
            reinterpret_cast<value_type*>(shared));
      }

      if (CudaTraits::WarpSize < word_count.value) {
        __syncthreads();
      }

      for (unsigned i = threadIdx.y; i < word_count.value; i += blockDim.y) {
        global[i] = shared[i];
      }
    }
  }

  // Determine block size constrained by shared memory:
  inline unsigned local_block_size(const FunctorType& f) {
    unsigned n = CudaTraits::WarpSize * 8;
    int shmem_size =
        cuda_single_inter_block_reduce_scan_shmem<false, WorkTag, value_type>(
            f, n);
    using closure_type =
        Impl::ParallelReduce<CombinedFunctorReducer<FunctorType, ReducerType>,
                             Policy, Kokkos::Cuda>;
    cudaFuncAttributes attr =
        CudaParallelLaunch<closure_type,
                           LaunchBounds>::get_cuda_func_attributes();
    while (
        (n &&
         (m_policy.space().impl_internal_space_instance()->m_maxShmemPerBlock <
          shmem_size)) ||
        (n >
         static_cast<unsigned>(
             Kokkos::Impl::cuda_get_max_block_size<FunctorType, LaunchBounds>(
                 m_policy.space().impl_internal_space_instance(), attr, f, 1,
                 shmem_size, 0)))) {
      n >>= 1;
      shmem_size =
          cuda_single_inter_block_reduce_scan_shmem<false, WorkTag, value_type>(
              f, n);
    }
    return n;
  }

  inline void execute() {
    const index_type nwork     = m_policy.end() - m_policy.begin();
    const bool need_device_set = ReducerType::has_init_member_function() ||
                                 ReducerType::has_final_member_function() ||
                                 !m_result_ptr_host_accessible ||
                                 Policy::is_graph_kernel::value ||
                                 !std::is_same<ReducerType, InvalidType>::value;
    if ((nwork > 0) || need_device_set) {
      const int block_size = local_block_size(m_functor_reducer.get_functor());

      KOKKOS_ASSERT(block_size > 0);

      // TODO: down casting these uses more space than required?
      m_scratch_space = (word_size_type*)cuda_internal_scratch_space(
          m_policy.space(), m_functor_reducer.get_reducer().value_size() *
                                block_size /* block_size == max block_count */);

      // Intentionally do not downcast to word_size_type since we use Cuda
      // atomics in Kokkos_Cuda_ReduceScan.hpp
      m_scratch_flags = cuda_internal_scratch_flags(m_policy.space(),
                                                    sizeof(Cuda::size_type));
      m_unified_space =
          reinterpret_cast<word_size_type*>(cuda_internal_scratch_unified(
              m_policy.space(), m_functor_reducer.get_reducer().value_size()));

      // REQUIRED ( 1 , N , 1 )
      dim3 block(1, block_size, 1);
      // Required grid.x <= block.y
      dim3 grid(std::min(int(block.y), int((nwork + block.y - 1) / block.y)), 1,
                1);

      // TODO @graph We need to effectively insert this in to the graph
      const int shmem =
          UseShflReduction
              ? 0
              : cuda_single_inter_block_reduce_scan_shmem<false, WorkTag,
                                                          value_type>(
                    m_functor_reducer.get_functor(), block.y);

      if ((nwork == 0)
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
          || Kokkos::Impl::CudaInternal::cuda_use_serial_execution()
#endif
      ) {
        block = dim3(1, 1, 1);
        grid  = dim3(1, 1, 1);
      }

      CudaParallelLaunch<ParallelReduce, LaunchBounds>(
          *this, grid, block, shmem,
          m_policy.space()
              .impl_internal_space_instance());  // copy to device and execute

      if (!m_result_ptr_device_accessible) {
        if (m_result_ptr) {
          if (m_unified_space) {
            m_policy.space().fence(
                "Kokkos::Impl::ParallelReduce<Cuda, RangePolicy>::execute: "
                "Result "
                "Not Device Accessible");
            const int count = m_functor_reducer.get_reducer().value_count();
            for (int i = 0; i < count; ++i) {
              m_result_ptr[i] = pointer_type(m_unified_space)[i];
            }
          } else {
            const int size = m_functor_reducer.get_reducer().value_size();
            DeepCopy<HostSpace, CudaSpace, Cuda>(m_policy.space(), m_result_ptr,
                                                 m_scratch_space, size);
          }
        }
      }
    } else {
      if (m_result_ptr) {
        // TODO @graph We need to effectively insert this in to the graph
        m_functor_reducer.get_reducer().init(m_result_ptr);
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
            MemorySpaceAccess<Kokkos::CudaSpace,
                              typename ViewType::memory_space>::accessible),
        m_result_ptr_host_accessible(
            MemorySpaceAccess<Kokkos::HostSpace,
                              typename ViewType::memory_space>::accessible),
        m_scratch_space(nullptr),
        m_scratch_flags(nullptr),
        m_unified_space(nullptr) {
    check_reduced_view_shmem_size<WorkTag, value_type>(
        m_policy, m_functor_reducer.get_functor());
  }
};

template <class FunctorType, class... Traits>
class ParallelScan<FunctorType, Kokkos::RangePolicy<Traits...>, Kokkos::Cuda> {
 public:
  using Policy = Kokkos::RangePolicy<Traits...>;

 private:
  using Member       = typename Policy::member_type;
  using WorkTag      = typename Policy::work_tag;
  using WorkRange    = typename Policy::WorkRange;
  using LaunchBounds = typename Policy::launch_bounds;

  using Analysis = Kokkos::Impl::FunctorAnalysis<FunctorPatternInterface::SCAN,
                                                 Policy, FunctorType, void>;

 public:
  using pointer_type   = typename Analysis::pointer_type;
  using reference_type = typename Analysis::reference_type;
  using value_type     = typename Analysis::value_type;
  using functor_type   = FunctorType;
  using size_type      = Cuda::size_type;
  // Conditionally set word_size_type to int16_t or int8_t if value_type is
  // smaller than int32_t (Kokkos::Cuda::size_type)
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

 private:
  // Algorithmic constraints:
  //  (a) blockDim.y is a power of two
  //  (b) blockDim.y == blockDim.z == 1
  //  (c) gridDim.x  <= blockDim.y * blockDim.y
  //  (d) gridDim.y  == gridDim.z == 1

  const CombinedFunctorReducer<FunctorType, typename Analysis::Reducer>
      m_functor_reducer;
  const Policy m_policy;
  word_size_type* m_scratch_space;
  size_type* m_scratch_flags;
  size_type m_final;
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
  bool m_run_serial;
#endif

  template <class TagType>
  __device__ inline std::enable_if_t<std::is_void<TagType>::value> exec_range(
      const Member& i, reference_type update, const bool final_result) const {
    m_functor_reducer.get_functor()(i, update, final_result);
  }

  template <class TagType>
  __device__ inline std::enable_if_t<!std::is_void<TagType>::value> exec_range(
      const Member& i, reference_type update, const bool final_result) const {
    m_functor_reducer.get_functor()(TagType(), i, update, final_result);
  }

  //----------------------------------------

  __device__ inline void initial() const {
    const typename Analysis::Reducer& final_reducer =
        m_functor_reducer.get_reducer();

    const integral_nonzero_constant<word_size_type, Analysis::StaticValueSize /
                                                        sizeof(word_size_type)>
        word_count(Analysis::value_size(m_functor_reducer.get_functor()) /
                   sizeof(word_size_type));

    word_size_type* const shared_value =
        kokkos_impl_cuda_shared_memory<word_size_type>() +
        word_count.value * threadIdx.y;

    final_reducer.init(reinterpret_cast<pointer_type>(shared_value));

    // Number of blocks is bounded so that the reduction can be limited to two
    // passes. Each thread block is given an approximately equal amount of work
    // to perform. Accumulate the values for this block. The accumulation
    // ordering does not match the final pass, but is arithmetically equivalent.

    const WorkRange range(m_policy, blockIdx.x, gridDim.x);

    for (Member iwork = range.begin() + threadIdx.y, iwork_end = range.end();
         iwork < iwork_end; iwork += blockDim.y) {
      this->template exec_range<WorkTag>(
          iwork,
          final_reducer.reference(reinterpret_cast<pointer_type>(shared_value)),
          false);
    }

    // Reduce and scan, writing out scan of blocks' totals and block-groups'
    // totals. Blocks' scan values are written to 'blockIdx.x' location.
    // Block-groups' scan values are at: i = ( j * blockDim.y - 1 ) for i <
    // gridDim.x
    cuda_single_inter_block_reduce_scan<true>(
        final_reducer, blockIdx.x, gridDim.x,
        kokkos_impl_cuda_shared_memory<word_size_type>(), m_scratch_space,
        m_scratch_flags);
  }

  //----------------------------------------

  __device__ inline void final() const {
    const typename Analysis::Reducer& final_reducer =
        m_functor_reducer.get_reducer();

    const integral_nonzero_constant<word_size_type, Analysis::StaticValueSize /
                                                        sizeof(word_size_type)>
        word_count(Analysis::value_size(m_functor_reducer.get_functor()) /
                   sizeof(word_size_type));

    // Use shared memory as an exclusive scan: { 0 , value[0] , value[1] ,
    // value[2] , ... }
    word_size_type* const shared_data =
        kokkos_impl_cuda_shared_memory<word_size_type>();
    word_size_type* const shared_prefix =
        shared_data + word_count.value * threadIdx.y;
    word_size_type* const shared_accum =
        shared_data + word_count.value * (blockDim.y + 1);

    // Starting value for this thread block is the previous block's total.
    if (blockIdx.x) {
      word_size_type* const block_total =
          m_scratch_space + word_count.value * (blockIdx.x - 1);
      for (unsigned i = threadIdx.y; i < word_count.value; ++i) {
        shared_accum[i] = block_total[i];
      }
    } else if (0 == threadIdx.y) {
      final_reducer.init(reinterpret_cast<pointer_type>(shared_accum));
    }

    const WorkRange range(m_policy, blockIdx.x, gridDim.x);

    for (typename Policy::member_type iwork_base = range.begin();
         iwork_base < range.end(); iwork_base += blockDim.y) {
      unsigned MASK                            = __activemask();
      const typename Policy::member_type iwork = iwork_base + threadIdx.y;

      __syncthreads();  // Don't overwrite previous iteration values until they
                        // are used

      final_reducer.init(
          reinterpret_cast<pointer_type>(shared_prefix + word_count.value));

      // Copy previous block's accumulation total into thread[0] prefix and
      // inclusive scan value of this block
      for (unsigned i = threadIdx.y; i < word_count.value; ++i) {
        shared_data[i + word_count.value] = shared_data[i] = shared_accum[i];
      }
      __syncwarp(MASK);
      if (CudaTraits::WarpSize < word_count.value) {
        __syncthreads();
      }  // Protect against large scan values.

      // Call functor to accumulate inclusive scan value for this work item
      if (iwork < range.end()) {
        this->template exec_range<WorkTag>(
            iwork,
            final_reducer.reference(reinterpret_cast<pointer_type>(
                shared_prefix + word_count.value)),
            false);
      }

      // Scan block values into locations shared_data[1..blockDim.y]
      cuda_intra_block_reduce_scan<true>(
          final_reducer,
          typename Analysis::pointer_type(shared_data + word_count.value));

      {
        word_size_type* const block_total =
            shared_data + word_count.value * blockDim.y;
        for (unsigned i = threadIdx.y; i < word_count.value; ++i) {
          shared_accum[i] = block_total[i];
        }
      }

      // Call functor with exclusive scan value
      if (iwork < range.end()) {
        this->template exec_range<WorkTag>(
            iwork,
            final_reducer.reference(
                reinterpret_cast<pointer_type>(shared_prefix)),
            true);
      }
    }
  }

 public:
  Policy const& get_policy() const { return m_policy; }

  //----------------------------------------

  __device__ inline void operator()() const {
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
    if (m_run_serial) {
      typename Analysis::value_type value;
      ValueInit::init(m_functor, (void*)&value);
      const WorkRange range(m_policy, blockIdx.x, gridDim.x);

      for (typename Policy::member_type iwork_base = range.begin();
           iwork_base < range.end(); iwork_base++) {
        this->template exec_range<WorkTag>(iwork_base, value, true);
      }
    } else {
#endif
      if (!m_final) {
        initial();
      } else {
        final();
      }
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
    }
#endif
  }

  // Determine block size constrained by shared memory:
  inline unsigned local_block_size(const FunctorType& f) {
    // blockDim.y must be power of two = 128 (4 warps) or 256 (8 warps) or 512
    // (16 warps) gridDim.x <= blockDim.y * blockDim.y
    //
    // 4 warps was 10% faster than 8 warps and 20% faster than 16 warps in unit
    // testing

    unsigned n = CudaTraits::WarpSize * 4;
    while (n &&
           unsigned(m_policy.space()
                        .impl_internal_space_instance()
                        ->m_maxShmemPerBlock) <
               cuda_single_inter_block_reduce_scan_shmem<true, WorkTag,
                                                         value_type>(f, n)) {
      n >>= 1;
    }
    return n;
  }

  inline void execute() {
    const auto nwork = m_policy.end() - m_policy.begin();
    if (nwork) {
      constexpr int GridMaxComputeCapability_2x = 0x0ffff;

      const int block_size = local_block_size(m_functor_reducer.get_functor());
      KOKKOS_ASSERT(block_size > 0);

      const int grid_max =
          (block_size * block_size) < GridMaxComputeCapability_2x
              ? (block_size * block_size)
              : GridMaxComputeCapability_2x;

      // At most 'max_grid' blocks:
      const int max_grid =
          std::min(int(grid_max), int((nwork + block_size - 1) / block_size));

      // How much work per block:
      const int work_per_block = (nwork + max_grid - 1) / max_grid;

      // How many block are really needed for this much work:
      const int grid_x = (nwork + work_per_block - 1) / work_per_block;

      m_scratch_space =
          reinterpret_cast<word_size_type*>(cuda_internal_scratch_space(
              m_policy.space(),
              Analysis::value_size(m_functor_reducer.get_functor()) * grid_x));
      m_scratch_flags =
          cuda_internal_scratch_flags(m_policy.space(), sizeof(size_type) * 1);

      dim3 grid(grid_x, 1, 1);
      dim3 block(1, block_size, 1);  // REQUIRED DIMENSIONS ( 1 , N , 1 )
      const int shmem = Analysis::value_size(m_functor_reducer.get_functor()) *
                        (block_size + 2);

#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
      if (m_run_serial) {
        block = dim3(1, 1, 1);
        grid  = dim3(1, 1, 1);
      } else {
#endif
        m_final = false;
        CudaParallelLaunch<ParallelScan, LaunchBounds>(
            *this, grid, block, shmem,
            m_policy.space()
                .impl_internal_space_instance());  // copy to device and execute
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
      }
#endif
      m_final = true;
      CudaParallelLaunch<ParallelScan, LaunchBounds>(
          *this, grid, block, shmem,
          m_policy.space()
              .impl_internal_space_instance());  // copy to device and execute
    }
  }

  ParallelScan(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_functor_reducer(arg_functor, typename Analysis::Reducer{arg_functor}),
        m_policy(arg_policy),
        m_scratch_space(nullptr),
        m_scratch_flags(nullptr),
        m_final(false)
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
        ,
        m_run_serial(Kokkos::Impl::CudaInternal::cuda_use_serial_execution())
#endif
  {
  }
};

//----------------------------------------------------------------------------
template <class FunctorType, class ReturnType, class... Traits>
class ParallelScanWithTotal<FunctorType, Kokkos::RangePolicy<Traits...>,
                            ReturnType, Kokkos::Cuda> {
 public:
  using Policy = Kokkos::RangePolicy<Traits...>;

 private:
  using Member       = typename Policy::member_type;
  using WorkTag      = typename Policy::work_tag;
  using WorkRange    = typename Policy::WorkRange;
  using LaunchBounds = typename Policy::launch_bounds;

  using Analysis =
      Kokkos::Impl::FunctorAnalysis<FunctorPatternInterface::SCAN, Policy,
                                    FunctorType, ReturnType>;

 public:
  using value_type     = typename Analysis::value_type;
  using pointer_type   = typename Analysis::pointer_type;
  using reference_type = typename Analysis::reference_type;
  using functor_type   = FunctorType;
  using size_type      = Cuda::size_type;
  // Conditionally set word_size_type to int16_t or int8_t if value_type is
  // smaller than int32_t (Kokkos::Cuda::size_type)
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

 private:
  // Algorithmic constraints:
  //  (a) blockDim.y is a power of two
  //  (b) blockDim.y == blockDim.z == 1
  //  (c) gridDim.x  <= blockDim.y * blockDim.y
  //  (d) gridDim.y  == gridDim.z == 1

  const CombinedFunctorReducer<FunctorType, typename Analysis::Reducer>
      m_functor_reducer;
  const Policy m_policy;
  word_size_type* m_scratch_space;
  size_type* m_scratch_flags;
  size_type m_final;
  const pointer_type m_result_ptr;
  const bool m_result_ptr_device_accessible;

#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
  bool m_run_serial;
#endif

  template <class TagType>
  __device__ inline std::enable_if_t<std::is_void<TagType>::value> exec_range(
      const Member& i, reference_type update, const bool final_result) const {
    m_functor_reducer.get_functor()(i, update, final_result);
  }

  template <class TagType>
  __device__ inline std::enable_if_t<!std::is_void<TagType>::value> exec_range(
      const Member& i, reference_type update, const bool final_result) const {
    m_functor_reducer.get_functor()(TagType(), i, update, final_result);
  }

  //----------------------------------------

  __device__ inline void initial() const {
    const typename Analysis::Reducer& final_reducer =
        m_functor_reducer.get_reducer();

    const integral_nonzero_constant<word_size_type, Analysis::StaticValueSize /
                                                        sizeof(word_size_type)>
        word_count(Analysis::value_size(m_functor_reducer.get_functor()) /
                   sizeof(word_size_type));

    word_size_type* const shared_value =
        kokkos_impl_cuda_shared_memory<word_size_type>() +
        word_count.value * threadIdx.y;

    final_reducer.init(reinterpret_cast<pointer_type>(shared_value));

    // Number of blocks is bounded so that the reduction can be limited to two
    // passes. Each thread block is given an approximately equal amount of work
    // to perform. Accumulate the values for this block. The accumulation
    // ordering does not match the final pass, but is arithmetically equivalent.

    const WorkRange range(m_policy, blockIdx.x, gridDim.x);

    for (Member iwork = range.begin() + threadIdx.y, iwork_end = range.end();
         iwork < iwork_end; iwork += blockDim.y) {
      this->template exec_range<WorkTag>(
          iwork,
          final_reducer.reference(reinterpret_cast<pointer_type>(shared_value)),
          false);
    }

    // Reduce and scan, writing out scan of blocks' totals and block-groups'
    // totals. Blocks' scan values are written to 'blockIdx.x' location.
    // Block-groups' scan values are at: i = ( j * blockDim.y - 1 ) for i <
    // gridDim.x
    cuda_single_inter_block_reduce_scan<true>(
        final_reducer, blockIdx.x, gridDim.x,
        kokkos_impl_cuda_shared_memory<word_size_type>(), m_scratch_space,
        m_scratch_flags);
  }

  //----------------------------------------

  __device__ inline void final() const {
    const typename Analysis::Reducer& final_reducer =
        m_functor_reducer.get_reducer();

    const integral_nonzero_constant<word_size_type, Analysis::StaticValueSize /
                                                        sizeof(word_size_type)>
        word_count(final_reducer.value_size() / sizeof(word_size_type));

    // Use shared memory as an exclusive scan: { 0 , value[0] , value[1] ,
    // value[2] , ... }
    word_size_type* const shared_data =
        kokkos_impl_cuda_shared_memory<word_size_type>();
    word_size_type* const shared_prefix =
        shared_data + word_count.value * threadIdx.y;
    word_size_type* const shared_accum =
        shared_data + word_count.value * (blockDim.y + 1);

    // Starting value for this thread block is the previous block's total.
    if (blockIdx.x) {
      word_size_type* const block_total =
          m_scratch_space + word_count.value * (blockIdx.x - 1);
      for (unsigned i = threadIdx.y; i < word_count.value; ++i) {
        shared_accum[i] = block_total[i];
      }
    } else if (0 == threadIdx.y) {
      final_reducer.init(reinterpret_cast<pointer_type>(shared_accum));
    }

    const WorkRange range(m_policy, blockIdx.x, gridDim.x);

    for (typename Policy::member_type iwork_base = range.begin();
         iwork_base < range.end(); iwork_base += blockDim.y) {
      unsigned MASK = __activemask();

      const typename Policy::member_type iwork = iwork_base + threadIdx.y;

      __syncthreads();  // Don't overwrite previous iteration values until they
                        // are used

      final_reducer.init(
          reinterpret_cast<pointer_type>(shared_prefix + word_count.value));

      // Copy previous block's accumulation total into thread[0] prefix and
      // inclusive scan value of this block
      for (unsigned i = threadIdx.y; i < word_count.value; ++i) {
        shared_data[i + word_count.value] = shared_data[i] = shared_accum[i];
      }

      __syncwarp(MASK);
      if (CudaTraits::WarpSize < word_count.value) {
        __syncthreads();
      }  // Protect against large scan values.

      // Call functor to accumulate inclusive scan value for this work item
      if (iwork < range.end()) {
        this->template exec_range<WorkTag>(
            iwork,
            final_reducer.reference(reinterpret_cast<pointer_type>(
                shared_prefix + word_count.value)),
            false);
      }

      // Scan block values into locations shared_data[1..blockDim.y]
      cuda_intra_block_reduce_scan<true>(
          final_reducer,
          typename Analysis::pointer_type(shared_data + word_count.value));

      {
        word_size_type* const block_total =
            shared_data + word_count.value * blockDim.y;
        for (unsigned i = threadIdx.y; i < word_count.value; ++i) {
          shared_accum[i] = block_total[i];
        }
      }

      // Call functor with exclusive scan value
      if (iwork < range.end()) {
        this->template exec_range<WorkTag>(
            iwork,
            final_reducer.reference(
                reinterpret_cast<pointer_type>(shared_prefix)),
            true);
      }
      if (iwork + 1 == m_policy.end() && m_policy.end() == range.end() &&
          m_result_ptr_device_accessible)
        *m_result_ptr = *reinterpret_cast<pointer_type>(shared_prefix);
    }
  }

 public:
  Policy const& get_policy() const { return m_policy; }

  //----------------------------------------

  __device__ inline void operator()() const {
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
    if (m_run_serial) {
      typename Analysis::value_type value;
      ValueInit::init(m_functor, (void*)&value);
      const WorkRange range(m_policy, blockIdx.x, gridDim.x);

      for (typename Policy::member_type iwork_base = range.begin();
           iwork_base < range.end(); iwork_base++) {
        this->template exec_range<WorkTag>(iwork_base, value, true);
      }
      *((typename Analysis::value_type*)m_scratch_space) = value;
    } else {
#endif
      if (!m_final) {
        initial();
      } else {
        final();
      }
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
    }
#endif
  }

  // Determine block size constrained by shared memory:
  inline unsigned local_block_size(const FunctorType& f) {
    // blockDim.y must be power of two = 128 (4 warps) or 256 (8 warps) or 512
    // (16 warps) gridDim.x <= blockDim.y * blockDim.y
    //
    // 4 warps was 10% faster than 8 warps and 20% faster than 16 warps in unit
    // testing

    unsigned n = CudaTraits::WarpSize * 4;
    while (n &&
           unsigned(m_policy.space()
                        .impl_internal_space_instance()
                        ->m_maxShmemPerBlock) <
               cuda_single_inter_block_reduce_scan_shmem<true, WorkTag,
                                                         value_type>(f, n)) {
      n >>= 1;
    }
    return n;
  }

  inline void execute() {
    const auto nwork = m_policy.end() - m_policy.begin();
    if (nwork) {
      enum { GridMaxComputeCapability_2x = 0x0ffff };

      const int block_size = local_block_size(m_functor_reducer.get_functor());
      KOKKOS_ASSERT(block_size > 0);

      const int grid_max =
          (block_size * block_size) < GridMaxComputeCapability_2x
              ? (block_size * block_size)
              : GridMaxComputeCapability_2x;

      // At most 'max_grid' blocks:
      const int max_grid =
          std::min(int(grid_max), int((nwork + block_size - 1) / block_size));

      // How much work per block:
      const int work_per_block = (nwork + max_grid - 1) / max_grid;

      // How many block are really needed for this much work:
      const int grid_x = (nwork + work_per_block - 1) / work_per_block;

      const typename Analysis::Reducer& final_reducer =
          m_functor_reducer.get_reducer();
      m_scratch_space =
          reinterpret_cast<word_size_type*>(cuda_internal_scratch_space(
              m_policy.space(), final_reducer.value_size() * grid_x));
      m_scratch_flags =
          cuda_internal_scratch_flags(m_policy.space(), sizeof(size_type) * 1);

      dim3 grid(grid_x, 1, 1);
      dim3 block(1, block_size, 1);  // REQUIRED DIMENSIONS ( 1 , N , 1 )
      const int shmem = final_reducer.value_size() * (block_size + 2);

#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
      if (m_run_serial) {
        block = dim3(1, 1, 1);
        grid  = dim3(1, 1, 1);
      } else
#endif
      {
        m_final = false;
        CudaParallelLaunch<ParallelScanWithTotal, LaunchBounds>(
            *this, grid, block, shmem,
            m_policy.space()
                .impl_internal_space_instance());  // copy to device and execute
      }
      m_final = true;
      CudaParallelLaunch<ParallelScanWithTotal, LaunchBounds>(
          *this, grid, block, shmem,
          m_policy.space()
              .impl_internal_space_instance());  // copy to device and execute

      const int size = final_reducer.value_size();
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
      if (m_run_serial)
        DeepCopy<HostSpace, CudaSpace, Cuda>(m_policy.space(), &m_returnvalue,
                                             m_scratch_space, size);
      else
#endif
      {
        if (!m_result_ptr_device_accessible)
          DeepCopy<HostSpace, CudaSpace, Cuda>(
              m_policy.space(), m_result_ptr,
              m_scratch_space + (grid_x - 1) * size / sizeof(word_size_type),
              size);
      }
    }
  }

  template <class ViewType>
  ParallelScanWithTotal(const FunctorType& arg_functor,
                        const Policy& arg_policy,
                        const ViewType& arg_result_view)
      : m_functor_reducer(arg_functor, typename Analysis::Reducer{arg_functor}),
        m_policy(arg_policy),
        m_scratch_space(nullptr),
        m_scratch_flags(nullptr),
        m_final(false),
        m_result_ptr(arg_result_view.data()),
        m_result_ptr_device_accessible(
            MemorySpaceAccess<Kokkos::CudaSpace,
                              typename ViewType::memory_space>::accessible)
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
        ,
        m_run_serial(Kokkos::Impl::CudaInternal::cuda_use_serial_execution())
#endif
  {
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif
#endif
