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

#ifndef KOKKOS_HIP_SHUFFLE_REDUCE_HPP
#define KOKKOS_HIP_SHUFFLE_REDUCE_HPP

#include <Kokkos_Macros.hpp>

#if defined(__HIPCC__)

#include <HIP/Kokkos_HIP_Vectorization.hpp>

#include <climits>

namespace Kokkos {
namespace Impl {

/* Algorithmic constraints:
 *   (a) threads with the same threadIdx.x have same value
 *   (b) blockDim.x == power of two
 *   (x) blockDim.z == 1
 */
template <typename ValueType, typename ReducerType>
__device__ inline void hip_intra_warp_shuffle_reduction(
    ValueType& result, ReducerType const& reducer,
    uint32_t const max_active_thread = blockDim.y) {
  unsigned int shift = 1;

  // Reduce over values from threads with different threadIdx.y
  unsigned int constexpr warp_size = HIPTraits::WarpSize;
  while (blockDim.x * shift < warp_size) {
    ValueType const tmp = shfl_down(result, blockDim.x * shift, warp_size);
    // Only join if upper thread is active (this allows non power of two for
    // blockDim.y)
    if (threadIdx.y + shift < max_active_thread) {
      reducer.join(&result, &tmp);
    }
    shift *= 2;
  }

  // Broadcast the result to all the threads in the warp
  result = shfl(result, 0, warp_size);
}

template <typename ValueType, typename ReducerType>
__device__ inline void hip_inter_warp_shuffle_reduction(
    ValueType& value, const ReducerType& reducer,
    const int max_active_thread = blockDim.y) {
  unsigned int constexpr warp_size = HIPTraits::WarpSize;
  int constexpr step_width         = 8;
  // Depending on the ValueType __shared__ memory must be aligned up to 8 byte
  // boundaries. The reason not to use ValueType directly is that for types with
  // constructors it could lead to race conditions.
  __shared__ double sh_result[(sizeof(ValueType) + 7) / 8 * step_width];
  ValueType* result = reinterpret_cast<ValueType*>(&sh_result);
  int const step    = warp_size / blockDim.x;
  int shift         = step_width;
  // Skip the code below if  threadIdx.y % step != 0
  int const id = threadIdx.y % step == 0 ? threadIdx.y / step : INT_MAX;
  if (id < step_width) {
    result[id] = value;
  }
  __syncthreads();
  while (shift <= max_active_thread / step) {
    if (shift <= id && shift + step_width > id && threadIdx.x == 0) {
      reducer.join(&result[id % step_width], &value);
    }
    __syncthreads();
    shift += step_width;
  }

  value = result[0];
  for (int i = 1; (i * step < max_active_thread) && (i < step_width); ++i)
    reducer.join(&value, &result[i]);
  __syncthreads();
}

template <typename ValueType, typename ReducerType>
__device__ inline void hip_intra_block_shuffle_reduction(
    ValueType& value, ReducerType const& reducer,
    int const max_active_thread = blockDim.y) {
  hip_intra_warp_shuffle_reduction(value, reducer, max_active_thread);
  hip_inter_warp_shuffle_reduction(value, reducer, max_active_thread);
}

template <class FunctorType>
__device__ inline bool hip_inter_block_shuffle_reduction(
    typename FunctorType::reference_type value,
    typename FunctorType::reference_type neutral, FunctorType const& reducer,
    HIP::size_type* const m_scratch_space,
    typename FunctorType::pointer_type const /*result*/,
    HIP::size_type* const m_scratch_flags,
    int const max_active_thread = blockDim.y) {
  using pointer_type = typename FunctorType::pointer_type;
  using value_type   = typename FunctorType::value_type;

  // Do the intra-block reduction with shfl operations for the intra warp
  // reduction and static shared memory for the inter warp reduction
  hip_intra_block_shuffle_reduction(value, reducer, max_active_thread);

  int const id = threadIdx.y * blockDim.x + threadIdx.x;

  // One thread in the block writes block result to global scratch_memory
  if (id == 0) {
    pointer_type global =
        reinterpret_cast<pointer_type>(m_scratch_space) + blockIdx.x;
    *global = value;
  }

  // One warp of last block performs inter block reduction through loading the
  // block values from global scratch_memory
  bool last_block = false;
  __syncthreads();
  int constexpr warp_size = HIPTraits::WarpSize;
  if (id < warp_size) {
    HIP::size_type count;

    // Figure out whether this is the last block
    if (id == 0) count = Kokkos::atomic_fetch_add(m_scratch_flags, 1);
    count = shfl(count, 0, warp_size);

    // Last block does the inter block reduction
    if (count == gridDim.x - 1) {
      // set flag back to zero
      if (id == 0) *m_scratch_flags = 0;
      last_block = true;
      value      = neutral;

      pointer_type const global =
          reinterpret_cast<pointer_type>(m_scratch_space);

      // Reduce all global values with splitting work over threads in one warp
      const int step_size = blockDim.x * blockDim.y < warp_size
                                ? blockDim.x * blockDim.y
                                : warp_size;
      for (int i = id; i < static_cast<int>(gridDim.x); i += step_size) {
        value_type tmp = global[i];
        reducer.join(&value, &tmp);
      }

      // Perform shfl reductions within the warp only join if contribution is
      // valid (allows gridDim.x non power of two and <warp_size)
      for (unsigned int i = 1; i < warp_size; i *= 2) {
        if ((blockDim.x * blockDim.y) > i) {
          value_type tmp = shfl_down(value, i, warp_size);
          if (id + i < gridDim.x) reducer.join(&value, &tmp);
        }
      }
    }
  }
  // The last block has in its thread=0 the global reduction value through
  // "value"
  return last_block;
}
}  // namespace Impl
}  // namespace Kokkos

#endif

#endif
