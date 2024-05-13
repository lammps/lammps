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

#ifndef KOKKOS_HIP_REDUCESCAN_HPP
#define KOKKOS_HIP_REDUCESCAN_HPP

#include <Kokkos_Macros.hpp>

#if defined(__HIPCC__)

#include <HIP/Kokkos_HIP.hpp>
#include <HIP/Kokkos_HIP_Vectorization.hpp>

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
// Reduction-only implementation
//----------------------------------------------------------------------------

template <class FunctorType, bool UseShfl>
struct HIPReductionsFunctor;

template <typename FunctorType>
struct HIPReductionsFunctor<FunctorType, true> {
  using pointer_type = typename FunctorType::pointer_type;
  using Scalar       = typename FunctorType::value_type;

  __device__ static inline void scalar_intra_warp_reduction(
      FunctorType const& functor,
      Scalar value,            // Contribution
      bool const skip_vector,  // Skip threads if Kokkos vector lanes are not
                               // part of the reduction
      int const width,         // How much of the warp participates
      Scalar& result) {
    for (int delta = skip_vector ? blockDim.x : 1; delta < width; delta *= 2) {
      Scalar tmp = shfl_down(value, delta, width);
      functor.join(&value, &tmp);
    }

    in_place_shfl(result, value, 0, width);
  }

  __device__ static inline void scalar_intra_block_reduction(
      FunctorType const& functor, Scalar value, bool const skip,
      Scalar* my_global_team_buffer_element, int const shared_elements,
      Scalar* shared_team_buffer_element) {
    constexpr unsigned int warp_size = HIPTraits::WarpSize;
    int const warp_id                = (threadIdx.y * blockDim.x) / warp_size;
    Scalar* const my_shared_team_buffer_element =
        shared_team_buffer_element + warp_id % shared_elements;

    // Warp Level Reduction, ignoring Kokkos vector entries
    scalar_intra_warp_reduction(functor, value, skip, warp_size, value);

    if (warp_id < shared_elements) {
      *my_shared_team_buffer_element = value;
    }
    // Wait for every warp to be done before using one warp to do the final
    // cross warp reduction
    __syncthreads();

    int const num_warps = blockDim.x * blockDim.y / warp_size;
    for (int w = shared_elements; w < num_warps; w += shared_elements) {
      if (warp_id >= w && warp_id < w + shared_elements) {
        if ((threadIdx.y * blockDim.x + threadIdx.x) % warp_size == 0)
          functor.join(my_shared_team_buffer_element, &value);
      }
      __syncthreads();
    }

    if (warp_id == 0) {
      functor.init(&value);
      for (unsigned int i = threadIdx.y * blockDim.x + threadIdx.x;
           i < blockDim.y * blockDim.x / warp_size; i += warp_size) {
        functor.join(&value, &shared_team_buffer_element[i]);
      }
      scalar_intra_warp_reduction(functor, value, false, warp_size,
                                  *my_global_team_buffer_element);
      __threadfence();
    }
  }

  __device__ static inline bool scalar_inter_block_reduction(
      FunctorType const& functor, HIP::size_type const block_count,
      HIP::size_type* const shared_data, HIP::size_type* const global_data,
      HIP::size_type* const global_flags) {
    Scalar* const global_team_buffer_element =
        reinterpret_cast<Scalar*>(global_data);
    Scalar* const my_global_team_buffer_element =
        global_team_buffer_element + blockIdx.x;
    Scalar* shared_team_buffer_elements =
        reinterpret_cast<Scalar*>(shared_data);
    Scalar value                     = shared_team_buffer_elements[threadIdx.y];
    constexpr unsigned int warp_size = Impl::HIPTraits::WarpSize;
    int shared_elements              = blockDim.x * blockDim.y / warp_size;
    int global_elements              = block_count;
    __syncthreads();

    scalar_intra_block_reduction(functor, value, true,
                                 my_global_team_buffer_element, shared_elements,
                                 shared_team_buffer_elements);
    __syncthreads();

    // Use the last block that is done to do the do the reduction across the
    // block
    unsigned int num_teams_done = 0;
    if (threadIdx.x + threadIdx.y == 0) {
      num_teams_done = Kokkos::atomic_fetch_add(global_flags, 1) + 1;
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
          functor, value, false, shared_team_buffer_elements + blockDim.y - 1,
          shared_elements, shared_team_buffer_elements);
    }

    return is_last_block;
  }
};

template <typename FunctorType>
struct HIPReductionsFunctor<FunctorType, false> {
  using pointer_type = typename FunctorType::pointer_type;
  using Scalar       = typename FunctorType::value_type;

  __device__ static inline void scalar_intra_warp_reduction(
      FunctorType const& functor,
      Scalar* value,           // Contribution
      bool const skip_vector,  // Skip threads if Kokkos vector lanes are not
                               // part of the reduction
      int const width)         // How much of the warp participates
  {
    int const lane_id =
        (threadIdx.y * blockDim.x + threadIdx.x) % HIPTraits::WarpSize;
    for (int delta = skip_vector ? blockDim.x : 1; delta < width; delta *= 2) {
      if (lane_id + delta < HIPTraits::WarpSize &&
          (lane_id % (delta * 2) == 0)) {
        functor.join(value, value + delta);
      }
    }
    *value = *(value - lane_id);
  }

  __device__ static inline void scalar_intra_block_reduction(
      FunctorType const& functor, Scalar value, bool const skip, Scalar* result,
      int const /*shared_elements*/, Scalar* shared_team_buffer_element) {
    int const warp_id = (threadIdx.y * blockDim.x) / HIPTraits::WarpSize;
    Scalar* const my_shared_team_buffer_element =
        shared_team_buffer_element + threadIdx.y * blockDim.x + threadIdx.x;
    *my_shared_team_buffer_element = value;
    // Warp Level Reduction, ignoring Kokkos vector entries
    scalar_intra_warp_reduction(functor, my_shared_team_buffer_element, skip,
                                HIPTraits::WarpSize);
    // Wait for every warp to be done before using one warp to do final cross
    // warp reduction
    __syncthreads();

    if (warp_id == 0) {
      const unsigned int delta =
          (threadIdx.y * blockDim.x + threadIdx.x) * HIPTraits::WarpSize;
      if (delta < blockDim.x * blockDim.y)
        *my_shared_team_buffer_element = shared_team_buffer_element[delta];
      scalar_intra_warp_reduction(
          functor, my_shared_team_buffer_element, false,
          blockDim.x * blockDim.y / HIPTraits::WarpSize);
      if (threadIdx.x + threadIdx.y == 0) {
        *result = *shared_team_buffer_element;
        if (skip) __threadfence();
      }
    }
  }

  template <typename SizeType>
  __device__ static inline bool scalar_inter_block_reduction(
      FunctorType const& functor, HIP::size_type const block_count,
      SizeType* const shared_data, SizeType* const global_data,
      HIP::size_type* const global_flags) {
    Scalar* const global_team_buffer_element =
        reinterpret_cast<Scalar*>(global_data);
    Scalar* const my_global_team_buffer_element =
        global_team_buffer_element + blockIdx.x;
    Scalar* shared_team_buffer_elements =
        reinterpret_cast<Scalar*>(shared_data);
    Scalar value        = shared_team_buffer_elements[threadIdx.y];
    int shared_elements = (blockDim.x * blockDim.y) / HIPTraits::WarpSize;
    int global_elements = block_count;
    __syncthreads();

    // Do the scalar reduction inside each block
    scalar_intra_block_reduction(functor, value, true,
                                 my_global_team_buffer_element, shared_elements,
                                 shared_team_buffer_elements);
    __syncthreads();

    // Use the last block that is done to do the do the reduction across the
    // block
    unsigned int num_teams_done = 0;
    if (threadIdx.x + threadIdx.y == 0) {
      num_teams_done = Kokkos::atomic_fetch_add(global_flags, 1) + 1;
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
// Fused reduction and scan implementation
//----------------------------------------------------------------------------
/*
 *  Algorithmic constraints:
 *   (a) blockDim.y <= 1024
 *   (b) blockDim.x == blockDim.z == 1
 */

template <bool DoScan, class FunctorType>
__device__ void hip_intra_block_reduce_scan(
    FunctorType const& functor,
    typename FunctorType::pointer_type const base_data) {
  using pointer_type = typename FunctorType::pointer_type;

  const unsigned value_count = functor.length();
  const unsigned not_less_power_of_two =
      (1 << (Impl::int_log2(blockDim.y - 1) + 1));
  const unsigned BlockSizeMask = not_less_power_of_two - 1;
  // There is at most one warp that is neither completely full or empty.
  // For that warp, we shift all indices logically to the end and ignore join
  // operations with unassigned indices in the warp when performing the intra
  // warp reduction/scan.
  const bool is_full_warp = (((threadIdx.y >> HIPTraits::WarpIndexShift) + 1)
                             << HIPTraits::WarpIndexShift) <= blockDim.y;

  auto block_reduce_step = [&functor, value_count](
                               int const R, pointer_type const TD, int const S,
                               pointer_type memory_start, int index_shift) {
    const auto join_ptr = TD - (value_count << S) + value_count * index_shift;
    if (R > ((1 << S) - 1) && join_ptr >= memory_start) {
      functor.join(TD, join_ptr);
    }
  };

  // Intra-warp reduction:
  int bit_shift = 0;
  {
    const unsigned mapped_idx =
        threadIdx.y + (is_full_warp ? 0
                                    : (not_less_power_of_two - blockDim.y) &
                                          (HIPTraits::WarpSize - 1));
    const pointer_type tdata_intra = base_data + value_count * threadIdx.y;
    const pointer_type warp_start =
        base_data + value_count * ((threadIdx.y >> HIPTraits::WarpIndexShift)
                                   << HIPTraits::WarpIndexShift);
    for (; (1 << bit_shift) < HIPTraits::WarpSize; ++bit_shift) {
      block_reduce_step(mapped_idx, tdata_intra, bit_shift, warp_start, 0);
    }
  }

  __syncthreads();  // Wait for all warps to reduce

  // Inter-warp reduce-scan by a single warp to avoid extra synchronizations
  {
    // There is at most one warp where the memory address to be used is not
    // (HIPTraits::WarpSize - 1) away from the warp start adress. For the
    // following reduction, we shift all indices logically to the end of the
    // next power-of-two to the number of warps.
    const unsigned n_active_warps =
        ((blockDim.y - 1) >> HIPTraits::WarpIndexShift) + 1;
    if (threadIdx.y < n_active_warps) {
      const bool is_full_warp_inter =
          threadIdx.y < (blockDim.y >> HIPTraits::WarpIndexShift);
      pointer_type const tdata_inter =
          base_data +
          value_count * (is_full_warp_inter
                             ? (threadIdx.y << HIPTraits::WarpIndexShift) +
                                   (HIPTraits::WarpSize - 1)
                             : blockDim.y - 1);
      const unsigned index_shift =
          is_full_warp_inter
              ? 0
              : blockDim.y - (threadIdx.y << HIPTraits::WarpIndexShift);
      const int rtid_inter = (threadIdx.y << HIPTraits::WarpIndexShift) +
                             (HIPTraits::WarpSize - 1) - index_shift;

      for (; (1 << bit_shift) < BlockSizeMask; ++bit_shift) {
        block_reduce_step(rtid_inter, tdata_inter, bit_shift, base_data,
                          index_shift);
      }
    }
  }

  __syncthreads();  // Wait for inter-warp reduce-scan to complete

  if (DoScan) {
    // Update all the values for the respective warps (except for the last one)
    // by adding from the last value of the previous warp.
    const unsigned int WarpMask = HIPTraits::WarpSize - 1;
    const int is_last_thread_in_warp =
        is_full_warp ? ((threadIdx.y & WarpMask) == HIPTraits::WarpSize - 1)
                     : (threadIdx.y == blockDim.y - 1);
    if (threadIdx.y >= HIPTraits::WarpSize && !is_last_thread_in_warp) {
      const int offset_to_previous_warp_total = (threadIdx.y & (~WarpMask)) - 1;
      functor.join(base_data + value_count * threadIdx.y,
                   base_data + value_count * offset_to_previous_warp_total);
    }
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

template <bool DoScan, typename FunctorType, typename SizeType>
__device__ bool hip_single_inter_block_reduce_scan_impl(
    FunctorType const& functor, HIP::size_type const block_id,
    HIP::size_type const block_count, SizeType* const shared_data,
    SizeType* const global_data, HIP::size_type* const global_flags) {
  using size_type    = SizeType;
  using value_type   = typename FunctorType::value_type;
  using pointer_type = typename FunctorType::pointer_type;

  // '__ffs' = position of the least significant bit set to 1.
  // 'blockDim.y' is guaranteed to be a power of two so this
  // is the integral shift value that can replace an integral divide.
  unsigned int const BlockSizeShift = __ffs(blockDim.y) - 1;
  unsigned int const BlockSizeMask  = blockDim.y - 1;

  // Must have power of two thread count
  if (BlockSizeMask & blockDim.y) {
    Kokkos::abort(
        "HIP::hip_single_inter_block_reduce_scan requires power-of-two "
        "blockDim");
  }

  const integral_nonzero_constant<
      size_type, std::is_pointer<typename FunctorType::reference_type>::value
                     ? 0
                     : sizeof(value_type) / sizeof(size_type)>
      word_count((sizeof(value_type) * functor.length()) / sizeof(size_type));

  // Reduce the accumulation for the entire block.
  hip_intra_block_reduce_scan<false>(functor, pointer_type(shared_data));

  {
    // Write accumulation total to global scratch space.
    // Accumulation total is the last thread's data.
    size_type* const shared = shared_data + word_count.value * BlockSizeMask;
    size_type* const global = global_data + word_count.value * block_id;

    for (size_t i = threadIdx.y; i < word_count.value; i += blockDim.y) {
      global[i] = shared[i];
    }
    __threadfence();
  }

  // Contributing blocks note that their contribution has been completed via an
  // atomic-increment flag If this block is not the last block to contribute to
  // this group then the block is done.
  const bool is_last_block = !__syncthreads_or(
      threadIdx.y
          ? 0
          : (1 + atomicInc(global_flags, block_count - 1) < block_count));
  if (is_last_block) {
    size_type const b = (static_cast<long long int>(block_count) *
                         static_cast<long long int>(threadIdx.y)) >>
                        BlockSizeShift;
    size_type const e = (static_cast<long long int>(block_count) *
                         static_cast<long long int>(threadIdx.y + 1)) >>
                        BlockSizeShift;

    {
      pointer_type const shared_data_thread = reinterpret_cast<pointer_type>(
          shared_data + word_count.value * threadIdx.y);
      /* reference_type shared_value = */ functor.init(shared_data_thread);

      for (size_type i = b; i < e; ++i) {
        functor.join(
            shared_data_thread,
            reinterpret_cast<pointer_type>(global_data + word_count.value * i));
      }
    }

    hip_intra_block_reduce_scan<DoScan>(functor, pointer_type(shared_data));

    if (DoScan) {
      pointer_type const shared_value = reinterpret_cast<pointer_type>(
          shared_data +
          word_count.value * (threadIdx.y ? threadIdx.y - 1 : blockDim.y));

      if (!threadIdx.y) {
        functor.init(shared_value);
      }

      // Join previous inclusive scan value to each member
      for (size_type i = b; i < e; ++i) {
        pointer_type const global_value =
            reinterpret_cast<pointer_type>(global_data + word_count.value * i);
        functor.join(shared_value, global_value);
        functor.copy(global_value, shared_value);
      }
    }
  }

  return is_last_block;
}

template <bool DoScan, typename FunctorType, typename SizeType>
__device__ bool hip_single_inter_block_reduce_scan(
    FunctorType const& functor, HIP::size_type const block_id,
    HIP::size_type const block_count, SizeType* const shared_data,
    SizeType* const global_data, HIP::size_type* const global_flags) {
  // If we are doing a reduction and we don't do an array reduction, we use the
  // reduction-only path. Otherwise, we use the common path between reduction
  // and scan.
  if (!DoScan && !std::is_pointer<typename FunctorType::reference_type>::value)
    // FIXME_HIP_PERFORMANCE I don't know where 16 comes from. This inequality
    // determines if we use shared memory (false) or shuffle (true)
    return Kokkos::Impl::HIPReductionsFunctor<
        FunctorType, (sizeof(typename FunctorType::value_type) >
                      16)>::scalar_inter_block_reduction(functor, block_count,
                                                         shared_data,
                                                         global_data,
                                                         global_flags);
  else {
    return hip_single_inter_block_reduce_scan_impl<DoScan>(
        functor, block_id, block_count, shared_data, global_data, global_flags);
  }
}

// Size in bytes required for inter block reduce or scan
template <bool DoScan, class ArgTag, class ValueType, class FunctorType>
inline std::enable_if_t<DoScan, unsigned>
hip_single_inter_block_reduce_scan_shmem(const FunctorType& functor,
                                         const unsigned BlockSize) {
  using Analysis =
      Impl::FunctorAnalysis<Impl::FunctorPatternInterface::SCAN,
                            RangePolicy<HIP, ArgTag>, FunctorType, ValueType>;

  return (BlockSize + 2) * Analysis::value_size(functor);
}

template <bool DoScan, class ArgTag, class ValueType, class FunctorType>
inline std::enable_if_t<!DoScan, unsigned>
hip_single_inter_block_reduce_scan_shmem(const FunctorType& functor,
                                         const unsigned BlockSize) {
  using Analysis =
      Impl::FunctorAnalysis<Impl::FunctorPatternInterface::REDUCE,
                            RangePolicy<HIP, ArgTag>, FunctorType, ValueType>;

  return (BlockSize + 2) * Analysis::value_size(functor);
}

}  // namespace Impl
}  // namespace Kokkos

#endif

#endif
