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

#ifndef KOKKOS_CUDA_REDUCESCAN_HPP
#define KOKKOS_CUDA_REDUCESCAN_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_CUDA)

#include <utility>

#include <Kokkos_Parallel.hpp>
#include <impl/Kokkos_Error.hpp>
#include <Cuda/Kokkos_Cuda_Vectorization.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
/*
 *  Algorithmic constraints:
 *   (a) threads with same threadIdx.y have same value
 *   (b) blockDim.x == power of two
 *   (c) blockDim.z == 1
 */

template <class ValueType, class ReducerType>
__device__ inline void cuda_intra_warp_reduction(
    ValueType& result, const ReducerType& reducer,
    const uint32_t max_active_thread = blockDim.y) {
  unsigned int shift = 1;

  // Reduce over values from threads with different threadIdx.y
  while (blockDim.x * shift < 32) {
    const ValueType tmp = shfl_down(result, blockDim.x * shift, 32u);
    // Only join if upper thread is active (this allows non power of two for
    // blockDim.y
    if (threadIdx.y + shift < max_active_thread) reducer.join(&result, &tmp);
    shift *= 2;
  }

  result = shfl(result, 0, 32);
}

template <class ValueType, class ReducerType>
__device__ inline void cuda_inter_warp_reduction(
    ValueType& value, const ReducerType& reducer,
    const int max_active_thread = blockDim.y) {
#define STEP_WIDTH 4
  // Depending on the ValueType _shared__ memory must be aligned up to 8byte
  // boundaries The reason not to use ValueType directly is that for types with
  // constructors it could lead to race conditions
  alignas(alignof(ValueType) > alignof(double) ? alignof(ValueType)
                                               : alignof(double))
      __shared__ double sh_result[(sizeof(ValueType) + 7) / 8 * STEP_WIDTH];
  ValueType* result = (ValueType*)&sh_result;
  const int step    = 32 / blockDim.x;
  int shift         = STEP_WIDTH;
  const int id      = threadIdx.y % step == 0 ? threadIdx.y / step : 65000;
  if (id < STEP_WIDTH) {
    result[id] = value;
  }
  __syncthreads();
  while (shift <= max_active_thread / step) {
    if (shift <= id && shift + STEP_WIDTH > id && threadIdx.x == 0) {
      reducer.join(&result[id % STEP_WIDTH], &value);
    }
    __syncthreads();
    shift += STEP_WIDTH;
  }

  value = result[0];
  for (int i = 1; (i * step < max_active_thread) && i < STEP_WIDTH; i++)
    reducer.join(&value, &result[i]);
  __syncthreads();
}

template <class ValueType, class ReducerType>
__device__ inline void cuda_intra_block_reduction(
    ValueType& value, const ReducerType& reducer,
    const int max_active_thread = blockDim.y) {
  cuda_intra_warp_reduction(value, reducer, max_active_thread);
  cuda_inter_warp_reduction(value, reducer, max_active_thread);
}

template <class FunctorType>
__device__ bool cuda_inter_block_reduction(
    typename FunctorType::reference_type value,
    typename FunctorType::reference_type neutral, const FunctorType& reducer,
    Cuda::size_type* const m_scratch_space,
    typename FunctorType::pointer_type const /*result*/,
    Cuda::size_type* const m_scratch_flags,
    const int max_active_thread = blockDim.y) {
  using pointer_type = typename FunctorType::pointer_type;
  using value_type   = typename FunctorType::value_type;

  // Do the intra-block reduction with shfl operations and static shared memory
  cuda_intra_block_reduction(value, reducer, max_active_thread);

  const int id = threadIdx.y * blockDim.x + threadIdx.x;

  // One thread in the block writes block result to global scratch_memory
  if (id == 0) {
    pointer_type global = ((pointer_type)m_scratch_space) + blockIdx.x;
    *global             = value;
  }

  // One warp of last block performs inter block reduction through loading the
  // block values from global scratch_memory
  bool last_block = false;
  __threadfence();
  __syncthreads();
  if (id < 32) {
    Cuda::size_type count;

    // Figure out whether this is the last block
    if (id == 0) count = Kokkos::atomic_fetch_add(m_scratch_flags, 1);
    count = Kokkos::shfl(count, 0, 32);

    // Last block does the inter block reduction
    if (count == gridDim.x - 1) {
      // set flag back to zero
      if (id == 0) *m_scratch_flags = 0;
      last_block = true;
      value      = neutral;

      pointer_type const volatile global = (pointer_type)m_scratch_space;

      // Reduce all global values with splitting work over threads in one warp
      const int step_size =
          blockDim.x * blockDim.y < 32 ? blockDim.x * blockDim.y : 32;
      for (int i = id; i < (int)gridDim.x; i += step_size) {
        value_type tmp = global[i];
        reducer.join(&value, &tmp);
      }

      // Perform shfl reductions within the warp only join if contribution is
      // valid (allows gridDim.x non power of two and <32)
      if (int(blockDim.x * blockDim.y) > 1) {
        value_type tmp = Kokkos::shfl_down(value, 1, 32);
        if (id + 1 < int(gridDim.x)) reducer.join(&value, &tmp);
      }
      unsigned int mask = __activemask();
      __syncwarp(mask);
      if (int(blockDim.x * blockDim.y) > 2) {
        value_type tmp = Kokkos::shfl_down(value, 2, 32);
        if (id + 2 < int(gridDim.x)) reducer.join(&value, &tmp);
      }
      __syncwarp(mask);
      if (int(blockDim.x * blockDim.y) > 4) {
        value_type tmp = Kokkos::shfl_down(value, 4, 32);
        if (id + 4 < int(gridDim.x)) reducer.join(&value, &tmp);
      }
      __syncwarp(mask);
      if (int(blockDim.x * blockDim.y) > 8) {
        value_type tmp = Kokkos::shfl_down(value, 8, 32);
        if (id + 8 < int(gridDim.x)) reducer.join(&value, &tmp);
      }
      __syncwarp(mask);
      if (int(blockDim.x * blockDim.y) > 16) {
        value_type tmp = Kokkos::shfl_down(value, 16, 32);
        if (id + 16 < int(gridDim.x)) reducer.join(&value, &tmp);
      }
      __syncwarp(mask);
    }
  }
  // The last block has in its thread=0 the global reduction value through
  // "value"
  return last_block;
}

template <class FunctorType, bool DoScan, bool UseShfl>
struct CudaReductionsFunctor;

template <class FunctorType>
struct CudaReductionsFunctor<FunctorType, false, true> {
  using pointer_type = typename FunctorType::pointer_type;
  using Scalar       = typename FunctorType::value_type;

  __device__ static inline void scalar_intra_warp_reduction(
      const FunctorType& functor,
      Scalar value,            // Contribution
      const bool skip_vector,  // Skip threads if Kokkos vector lanes are not
                               // part of the reduction
      const int width,         // How much of the warp participates
      Scalar& result) {
    unsigned mask =
        width == 32
            ? 0xffffffff
            : ((1 << width) - 1)
                  << ((threadIdx.y * blockDim.x + threadIdx.x) / width) * width;
    for (int delta = skip_vector ? blockDim.x : 1; delta < width; delta *= 2) {
      Scalar tmp = Kokkos::shfl_down(value, delta, width, mask);
      functor.join(&value, &tmp);
    }

    Impl::in_place_shfl(result, value, 0, width, mask);
  }

  __device__ static inline void scalar_intra_block_reduction(
      const FunctorType& functor, Scalar value, const bool skip,
      Scalar* my_global_team_buffer_element, const int shared_elements,
      Scalar* shared_team_buffer_element) {
    const int warp_id = (threadIdx.y * blockDim.x) / 32;
    Scalar* const my_shared_team_buffer_element =
        shared_team_buffer_element + warp_id % shared_elements;

    // Warp Level Reduction, ignoring Kokkos vector entries
    scalar_intra_warp_reduction(functor, value, skip, 32, value);

    if (warp_id < shared_elements) {
      *my_shared_team_buffer_element = value;
    }
    // Wait for every warp to be done before using one warp to do final cross
    // warp reduction
    __syncthreads();

    const int num_warps = blockDim.x * blockDim.y / 32;
    for (int w = shared_elements; w < num_warps; w += shared_elements) {
      if (warp_id >= w && warp_id < w + shared_elements) {
        if ((threadIdx.y * blockDim.x + threadIdx.x) % 32 == 0)
          functor.join(my_shared_team_buffer_element, &value);
      }
      __syncthreads();
    }

    if (warp_id == 0) {
      functor.init(&value);
      for (unsigned int i = threadIdx.y * blockDim.x + threadIdx.x;
           i < blockDim.y * blockDim.x / 32; i += 32)
        functor.join(&value, &shared_team_buffer_element[i]);
      scalar_intra_warp_reduction(functor, value, false, 32,
                                  *my_global_team_buffer_element);
    }
  }

  __device__ static inline bool scalar_inter_block_reduction(
      const FunctorType& functor, const Cuda::size_type /*block_id*/,
      const Cuda::size_type block_count, Cuda::size_type* const shared_data,
      Cuda::size_type* const global_data, Cuda::size_type* const global_flags) {
    Scalar* const global_team_buffer_element = ((Scalar*)global_data);
    Scalar* const my_global_team_buffer_element =
        global_team_buffer_element + blockIdx.x;
    Scalar* shared_team_buffer_elements = ((Scalar*)shared_data);
    Scalar value        = shared_team_buffer_elements[threadIdx.y];
    int shared_elements = blockDim.x * blockDim.y / 32;
    int global_elements = block_count;
    __syncthreads();

    scalar_intra_block_reduction(functor, value, true,
                                 my_global_team_buffer_element, shared_elements,
                                 shared_team_buffer_elements);
    __threadfence();
    __syncthreads();
    unsigned int num_teams_done = 0;
    // The cast in the atomic call is necessary to find matching call with
    // MSVC/NVCC
    if (threadIdx.x + threadIdx.y == 0) {
      num_teams_done =
          Kokkos::atomic_fetch_add(global_flags, static_cast<unsigned int>(1)) +
          1;
    }
    bool is_last_block = false;
    if (__syncthreads_or(num_teams_done == gridDim.x)) {
      is_last_block = true;
      *global_flags = 0;
      functor.init(&value);
      for (int i = threadIdx.y * blockDim.x + threadIdx.x; i < global_elements;
           i += blockDim.x * blockDim.y) {
        functor.join(&value, &global_team_buffer_element[i]);
      }
      scalar_intra_block_reduction(
          functor, value, false, shared_team_buffer_elements + (blockDim.y - 1),
          shared_elements, shared_team_buffer_elements);
    }
    return is_last_block;
  }
};

template <class FunctorType>
struct CudaReductionsFunctor<FunctorType, false, false> {
  using pointer_type = typename FunctorType::pointer_type;
  using Scalar       = typename FunctorType::value_type;

  __device__ static inline void scalar_intra_warp_reduction(
      const FunctorType& functor,
      Scalar* value,           // Contribution
      const bool skip_vector,  // Skip threads if Kokkos vector lanes are not
                               // part of the reduction
      const int width)         // How much of the warp participates
  {
    unsigned mask =
        width == 32
            ? 0xffffffff
            : ((1 << width) - 1)
                  << ((threadIdx.y * blockDim.x + threadIdx.x) / width) * width;
    const int lane_id = (threadIdx.y * blockDim.x + threadIdx.x) % 32;

    __syncwarp(mask);

    for (int delta = skip_vector ? blockDim.x : 1; delta < width; delta *= 2) {
      if (lane_id + delta < 32) {
        functor.join(value, value + delta);
      }
      __syncwarp(mask);
    }
    if (lane_id != 0) {
      *value = *(value - lane_id);
    }
  }

  __device__ static inline void scalar_intra_block_reduction(
      const FunctorType& functor, Scalar value, const bool skip, Scalar* result,
      const int /*shared_elements*/, Scalar* shared_team_buffer_element) {
    const int warp_id = (threadIdx.y * blockDim.x) / 32;
    Scalar* const my_shared_team_buffer_element =
        shared_team_buffer_element + threadIdx.y * blockDim.x + threadIdx.x;
    *my_shared_team_buffer_element = value;
    // Warp Level Reduction, ignoring Kokkos vector entries
    scalar_intra_warp_reduction(functor, my_shared_team_buffer_element, skip,
                                32);
    // Wait for every warp to be done before using one warp to do final cross
    // warp reduction
    __syncthreads();

    if (warp_id == 0) {
      const unsigned int delta = (threadIdx.y * blockDim.x + threadIdx.x) * 32;
      if (delta < blockDim.x * blockDim.y)
        *my_shared_team_buffer_element = shared_team_buffer_element[delta];
      __syncwarp(0xffffffff);
      scalar_intra_warp_reduction(functor, my_shared_team_buffer_element, false,
                                  blockDim.x * blockDim.y / 32);
      if (threadIdx.x + threadIdx.y == 0) *result = *shared_team_buffer_element;
    }
  }

  template <class SizeType = Cuda::size_type>
  __device__ static inline bool scalar_inter_block_reduction(
      const FunctorType& functor, const Cuda::size_type /*block_id*/,
      const Cuda::size_type block_count, SizeType* const shared_data,
      SizeType* const global_data, Cuda::size_type* const global_flags) {
    Scalar* const global_team_buffer_element = ((Scalar*)global_data);
    Scalar* const my_global_team_buffer_element =
        global_team_buffer_element + blockIdx.x;
    Scalar* shared_team_buffer_elements = ((Scalar*)shared_data);
    Scalar value        = shared_team_buffer_elements[threadIdx.y];
    int shared_elements = blockDim.x * blockDim.y / 32;
    int global_elements = block_count;
    __syncthreads();

    scalar_intra_block_reduction(functor, value, true,
                                 my_global_team_buffer_element, shared_elements,
                                 shared_team_buffer_elements);
    __threadfence();
    __syncthreads();

    unsigned int num_teams_done = 0;
    // The cast in the atomic call is necessary to find matching call with
    // MSVC/NVCC
    if (threadIdx.x + threadIdx.y == 0) {
      num_teams_done =
          Kokkos::atomic_fetch_add(global_flags, static_cast<unsigned int>(1)) +
          1;
    }
    bool is_last_block = false;
    if (__syncthreads_or(num_teams_done == gridDim.x)) {
      is_last_block = true;
      *global_flags = 0;
      functor.init(&value);
      for (int i = threadIdx.y * blockDim.x + threadIdx.x; i < global_elements;
           i += blockDim.x * blockDim.y) {
        functor.join(&value, &global_team_buffer_element[i]);
      }
      scalar_intra_block_reduction(
          functor, value, false, shared_team_buffer_elements + (blockDim.y - 1),
          shared_elements, shared_team_buffer_elements);
    }
    return is_last_block;
  }
};
//----------------------------------------------------------------------------
// See section B.17 of Cuda C Programming Guide Version 3.2
// for discussion of
//   __launch_bounds__(maxThreadsPerBlock,minBlocksPerMultiprocessor)
// function qualifier which could be used to improve performance.
//----------------------------------------------------------------------------
/*
 *  Algorithmic constraints:
 *   (a) blockDim.y <= 1024
 *   (b) blockDim.x == blockDim.z == 1
 */

template <bool DoScan, class FunctorType>
__device__ void cuda_intra_block_reduce_scan(
    const FunctorType& functor,
    const typename FunctorType::pointer_type base_data) {
  using pointer_type = typename FunctorType::pointer_type;

  const unsigned value_count = functor.length();
  const unsigned not_less_power_of_two =
      (1 << (Impl::int_log2(blockDim.y - 1) + 1));
  const unsigned BlockSizeMask = not_less_power_of_two - 1;
  // There is at most one warp that is neither completely full or empty.
  // For that warp, we shift all indices logically to the end and ignore join
  // operations with unassigned indices in the warp when performing the intra
  // warp reduction/scan.
  const bool is_full_warp = (((threadIdx.y >> CudaTraits::WarpIndexShift) + 1)
                             << CudaTraits::WarpIndexShift) <= blockDim.y;

  const unsigned mapped_idx =
      threadIdx.y + (is_full_warp ? 0
                                  : (not_less_power_of_two - blockDim.y) &
                                        (CudaTraits::WarpSize - 1));
  const pointer_type tdata_intra = base_data + value_count * threadIdx.y;
  const pointer_type warp_start =
      base_data + value_count * ((threadIdx.y >> CudaTraits::WarpIndexShift)
                                 << CudaTraits::WarpIndexShift);

  auto block_reduce_step = [&functor, value_count](
                               int const R, pointer_type const TD, int const S,
                               pointer_type memory_start, int index_shift) {
    const auto join_ptr = TD - (value_count << S) + value_count * index_shift;
    if (((R + 1) & ((1 << (S + 1)) - 1)) == 0 && join_ptr >= memory_start) {
      functor.join(TD, join_ptr);
    }
  };

  auto block_scan_step = [&functor, value_count](
                             int const R, pointer_type const TD, int const S,
                             pointer_type memory_start, int index_shift) {
    const auto N        = (1 << (S + 1));
    const auto join_ptr = TD - (value_count << S) + value_count * index_shift;
    if (R >= N && ((R + 1) & (N - 1)) == (N >> 1) && join_ptr >= memory_start) {
      functor.join(TD, join_ptr);
    }
  };

  {  // Intra-warp reduction:
    __syncwarp(0xffffffff);
    block_reduce_step(mapped_idx, tdata_intra, 0, warp_start, 0);
    __syncwarp(0xffffffff);
    block_reduce_step(mapped_idx, tdata_intra, 1, warp_start, 0);
    __syncwarp(0xffffffff);
    block_reduce_step(mapped_idx, tdata_intra, 2, warp_start, 0);
    __syncwarp(0xffffffff);
    block_reduce_step(mapped_idx, tdata_intra, 3, warp_start, 0);
    __syncwarp(0xffffffff);
    block_reduce_step(mapped_idx, tdata_intra, 4, warp_start, 0);
    __syncwarp(0xffffffff);
  }

  __syncthreads();  // Wait for all warps to reduce

  // Inter-warp reduce-scan by a single warp to avoid extra synchronizations.
  {
    // There is at most one warp where the memory address to be used is not
    // (CudaTraits::WarpSize - 1) away from the warp start adress. For the
    // following reduction, we shift all indices logically to the end of the
    // next power-of-two to the number of warps.
    const unsigned n_active_warps =
        ((blockDim.y - 1) >> CudaTraits::WarpIndexShift) + 1;
    const unsigned inner_mask =
        __ballot_sync(0xffffffff, (threadIdx.y < n_active_warps));
    if (threadIdx.y < n_active_warps) {
      const bool is_full_warp_inter =
          threadIdx.y < (blockDim.y >> CudaTraits::WarpIndexShift);
      const pointer_type tdata_inter =
          base_data +
          value_count * (is_full_warp_inter
                             ? (threadIdx.y << CudaTraits::WarpIndexShift) +
                                   (CudaTraits::WarpSize - 1)
                             : blockDim.y - 1);
      const unsigned index_shift =
          is_full_warp_inter
              ? 0
              : blockDim.y - (threadIdx.y << CudaTraits::WarpIndexShift);
      const int rtid_inter = (threadIdx.y << CudaTraits::WarpIndexShift) +
                             (CudaTraits::WarpSize - 1) - index_shift;

      if ((1 << 5) < BlockSizeMask) {
        __syncwarp(inner_mask);
        block_reduce_step(rtid_inter, tdata_inter, 5, base_data, index_shift);
      }
      if ((1 << 6) < BlockSizeMask) {
        __syncwarp(inner_mask);
        block_reduce_step(rtid_inter, tdata_inter, 6, base_data, index_shift);
      }
      if ((1 << 7) < BlockSizeMask) {
        __syncwarp(inner_mask);
        block_reduce_step(rtid_inter, tdata_inter, 7, base_data, index_shift);
      }
      if ((1 << 8) < BlockSizeMask) {
        __syncwarp(inner_mask);
        block_reduce_step(rtid_inter, tdata_inter, 8, base_data, index_shift);
      }
      if ((1 << 9) < BlockSizeMask) {
        __syncwarp(inner_mask);
        block_reduce_step(rtid_inter, tdata_inter, 9, base_data, index_shift);
      }

      if (DoScan) {
        __syncwarp(inner_mask);
        block_scan_step(rtid_inter, tdata_inter, 8, base_data, index_shift);
        __syncwarp(inner_mask);
        block_scan_step(rtid_inter, tdata_inter, 7, base_data, index_shift);
        __syncwarp(inner_mask);
        block_scan_step(rtid_inter, tdata_inter, 6, base_data, index_shift);
        __syncwarp(inner_mask);
        block_scan_step(rtid_inter, tdata_inter, 5, base_data, index_shift);
      }
    }
  }

  __syncthreads();  // Wait for inter-warp reduce-scan to complete

  if (DoScan) {
    block_scan_step(mapped_idx, tdata_intra, 4, warp_start, 0);
    __threadfence_block();
    __syncwarp(0xffffffff);
    block_scan_step(mapped_idx, tdata_intra, 3, warp_start, 0);
    __threadfence_block();
    __syncwarp(0xffffffff);
    block_scan_step(mapped_idx, tdata_intra, 2, warp_start, 0);
    __threadfence_block();
    __syncwarp(0xffffffff);
    block_scan_step(mapped_idx, tdata_intra, 1, warp_start, 0);
    __threadfence_block();
    __syncwarp(0xffffffff);
    block_scan_step(mapped_idx, tdata_intra, 0, warp_start, 0);
    __threadfence_block();
    __syncwarp(0xffffffff);
    // Update with total from previous warps
    if (mapped_idx >= CudaTraits::WarpSize &&
        (mapped_idx & (CudaTraits::WarpSize - 1)) != (CudaTraits::WarpSize - 1))
      functor.join(tdata_intra, warp_start - value_count);
    __syncwarp(0xffffffff);
  }
}

//----------------------------------------------------------------------------
/**\brief  Input value-per-thread starting at 'shared_data'.
 *         Reduction value at last thread's location.
 *
 *  If 'DoScan' then write blocks' scan values and block-groups' scan values.
 *
 *  Global reduce result is in the last threads' 'shared_data' location.
 */

template <bool DoScan, class FunctorType, class SizeType = Cuda::size_type>
__device__ bool cuda_single_inter_block_reduce_scan2(
    const FunctorType& functor, const Cuda::size_type block_id,
    const Cuda::size_type block_count, SizeType* const shared_data,
    SizeType* const global_data, Cuda::size_type* const global_flags) {
  using size_type    = SizeType;
  using value_type   = typename FunctorType::value_type;
  using pointer_type = typename FunctorType::pointer_type;

  // '__ffs' = position of the least significant bit set to 1.
  // 'blockDim.y' is guaranteed to be a power of two so this
  // is the integral shift value that can replace an integral divide.
  const unsigned BlockSizeShift = __ffs(blockDim.y) - 1;
  const unsigned BlockSizeMask  = blockDim.y - 1;

  // Must have power of two thread count
  if (BlockSizeMask & blockDim.y) {
    Kokkos::abort(
        "Cuda::cuda_single_inter_block_reduce_scan requires power-of-two "
        "blockDim");
  }

  const integral_nonzero_constant<
      size_type, std::is_pointer<typename FunctorType::reference_type>::value
                     ? 0
                     : sizeof(value_type) / sizeof(size_type)>
      word_count((sizeof(value_type) * functor.length()) / sizeof(size_type));

  // Reduce the accumulation for the entire block.
  cuda_intra_block_reduce_scan<false>(functor, pointer_type(shared_data));
  {
    // Write accumulation total to global scratch space.
    // Accumulation total is the last thread's data.
    size_type* const shared = shared_data + word_count.value * BlockSizeMask;
    size_type* const global = global_data + word_count.value * block_id;

    for (int i = int(threadIdx.y); i < int(word_count.value);
         i += int(blockDim.y)) {
      global[i] = shared[i];
    }
  }
  __threadfence();

  // Contributing blocks note that their contribution has been completed via an
  // atomic-increment flag If this block is not the last block to contribute to
  // this group then the block is done.
  const bool is_last_block = !__syncthreads_or(
      threadIdx.y
          ? 0
          : (1 + atomicInc(global_flags, block_count - 1) < block_count));

  if (is_last_block) {
    const size_type b =
        (long(block_count) * long(threadIdx.y)) >> BlockSizeShift;
    const size_type e =
        (long(block_count) * long(threadIdx.y + 1)) >> BlockSizeShift;

    {
      void* const shared_ptr = shared_data + word_count.value * threadIdx.y;
      /* reference_type shared_value = */ functor.init(
          static_cast<pointer_type>(shared_ptr));

      for (size_type i = b; i < e; ++i) {
        functor.join(
            static_cast<pointer_type>(shared_ptr),
            reinterpret_cast<pointer_type>(global_data + word_count.value * i));
      }
    }

    cuda_intra_block_reduce_scan<DoScan>(functor, pointer_type(shared_data));

    if (DoScan) {
      pointer_type const shared_value = reinterpret_cast<pointer_type>(
          shared_data +
          word_count.value * (threadIdx.y ? threadIdx.y - 1 : blockDim.y));

      if (!threadIdx.y) {
        functor.init(shared_value);
      }

      // Join previous inclusive scan value to each member
      for (size_type i = b; i < e; ++i) {
        size_type* const global_value = global_data + word_count.value * i;
        functor.join(shared_value,
                     reinterpret_cast<pointer_type>(global_value));
        functor.copy(reinterpret_cast<pointer_type>(global_value),
                     reinterpret_cast<pointer_type>(shared_value));
      }
    }
  }

  return is_last_block;
}

template <bool DoScan, class FunctorType, class SizeType = Cuda::size_type>
__device__ bool cuda_single_inter_block_reduce_scan(
    const FunctorType& functor, const Cuda::size_type block_id,
    const Cuda::size_type block_count, SizeType* const shared_data,
    SizeType* const global_data, Cuda::size_type* const global_flags) {
  if (!DoScan && !std::is_pointer<typename FunctorType::reference_type>::value)
    return Kokkos::Impl::CudaReductionsFunctor<
        FunctorType, false, (sizeof(typename FunctorType::value_type) > 16)>::
        scalar_inter_block_reduction(functor, block_id, block_count,
                                     shared_data, global_data, global_flags);
  else
    return cuda_single_inter_block_reduce_scan2<DoScan>(
        functor, block_id, block_count, shared_data, global_data, global_flags);
}

// Size in bytes required for inter block reduce or scan
template <bool DoScan, class FunctorType, class ArgTag>
inline std::enable_if_t<DoScan, unsigned>
cuda_single_inter_block_reduce_scan_shmem(const FunctorType& functor,
                                          const unsigned BlockSize) {
  using Analysis =
      Impl::FunctorAnalysis<Impl::FunctorPatternInterface::SCAN,
                            RangePolicy<Cuda, ArgTag>, FunctorType>;

  return (BlockSize + 2) * Analysis::value_size(functor);
}

template <bool DoScan, class FunctorType, class ArgTag>
inline std::enable_if_t<!DoScan, unsigned>
cuda_single_inter_block_reduce_scan_shmem(const FunctorType& functor,
                                          const unsigned BlockSize) {
  using Analysis =
      Impl::FunctorAnalysis<Impl::FunctorPatternInterface::REDUCE,
                            RangePolicy<Cuda, ArgTag>, FunctorType>;

  return (BlockSize + 2) * Analysis::value_size(functor);
}

template <typename WorkTag, typename Policy, typename FunctorType>
inline void check_reduced_view_shmem_size(const Policy& policy,
                                          const FunctorType& functor) {
  size_t minBlockSize = CudaTraits::WarpSize * 1;
  unsigned reqShmemSize =
      cuda_single_inter_block_reduce_scan_shmem<false, FunctorType, WorkTag>(
          functor, minBlockSize);
  size_t maxShmemPerBlock =
      policy.space().impl_internal_space_instance()->m_maxShmemPerBlock;

  if (reqShmemSize > maxShmemPerBlock) {
    Kokkos::Impl::throw_runtime_exception(
        "Kokkos::Impl::ParallelReduce< Cuda > requested too much L0 scratch "
        "memory");
  }
}

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined(KOKKOS_ENABLE_CUDA) */
#endif /* KOKKOS_CUDA_REDUCESCAN_HPP */
