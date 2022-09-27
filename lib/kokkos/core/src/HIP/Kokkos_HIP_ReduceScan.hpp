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

#ifndef KOKKOS_HIP_REDUCESCAN_HPP
#define KOKKOS_HIP_REDUCESCAN_HPP

#include <Kokkos_Macros.hpp>

#if defined(__HIPCC__)

#include <HIP/Kokkos_HIP_Vectorization.hpp>

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
// Reduction-only implementation
//----------------------------------------------------------------------------

template <class FunctorType, class ArgTag, bool UseShfl>
struct HIPReductionsFunctor;

template <typename FunctorType, typename ArgTag>
struct HIPReductionsFunctor<FunctorType, ArgTag, true> {
  using ValueTraits  = FunctorValueTraits<FunctorType, ArgTag>;
  using ValueJoin    = FunctorValueJoin<FunctorType, ArgTag>;
  using ValueInit    = FunctorValueInit<FunctorType, ArgTag>;
  using ValueOps     = FunctorValueOps<FunctorType, ArgTag>;
  using pointer_type = typename ValueTraits::pointer_type;
  using Scalar       = typename ValueTraits::value_type;

  __device__ static inline void scalar_intra_warp_reduction(
      FunctorType const& functor,
      Scalar value,            // Contribution
      bool const skip_vector,  // Skip threads if Kokkos vector lanes are not
                               // part of the reduction
      int const width,         // How much of the warp participates
      Scalar& result) {
    for (int delta = skip_vector ? blockDim.x : 1; delta < width; delta *= 2) {
      Scalar tmp = Kokkos::Experimental::shfl_down(value, delta, width);
      ValueJoin::join(functor, &value, &tmp);
    }

    Experimental::Impl::in_place_shfl(result, value, 0, width);
  }

  __device__ static inline void scalar_intra_block_reduction(
      FunctorType const& functor, Scalar value, bool const skip,
      Scalar* my_global_team_buffer_element, int const shared_elements,
      Scalar* shared_team_buffer_element) {
    unsigned int constexpr warp_size =
        Kokkos::Experimental::Impl::HIPTraits::WarpSize;
    int const warp_id = (threadIdx.y * blockDim.x) / warp_size;
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
          ValueJoin::join(functor, my_shared_team_buffer_element, &value);
      }
      __syncthreads();
    }

    if (warp_id == 0) {
      ValueInit::init(functor, &value);
      for (unsigned int i = threadIdx.y * blockDim.x + threadIdx.x;
           i < blockDim.y * blockDim.x / warp_size; i += warp_size) {
        ValueJoin::join(functor, &value, &shared_team_buffer_element[i]);
      }
      scalar_intra_warp_reduction(functor, value, false, warp_size,
                                  *my_global_team_buffer_element);
    }
  }

  __device__ static inline bool scalar_inter_block_reduction(
      FunctorType const& functor,
      ::Kokkos::Experimental::HIP::size_type const block_count,
      ::Kokkos::Experimental::HIP::size_type* const shared_data,
      ::Kokkos::Experimental::HIP::size_type* const global_data,
      ::Kokkos::Experimental::HIP::size_type* const global_flags) {
    Scalar* const global_team_buffer_element =
        reinterpret_cast<Scalar*>(global_data);
    Scalar* const my_global_team_buffer_element =
        global_team_buffer_element + blockIdx.x;
    Scalar* shared_team_buffer_elements =
        reinterpret_cast<Scalar*>(shared_data);
    Scalar value = shared_team_buffer_elements[threadIdx.y];
    unsigned int constexpr warp_size =
        Kokkos::Experimental::Impl::HIPTraits::WarpSize;
    int shared_elements = blockDim.x * blockDim.y / warp_size;
    int global_elements = block_count;
    __syncthreads();

    scalar_intra_block_reduction(functor, value, true,
                                 my_global_team_buffer_element, shared_elements,
                                 shared_team_buffer_elements);
    __syncthreads();

    // Use the last block that is done to do the do the reduction across the
    // block
    __shared__ unsigned int num_teams_done;
    if (threadIdx.x + threadIdx.y == 0) {
      num_teams_done = Kokkos::atomic_fetch_add(global_flags, 1) + 1;
    }
    bool is_last_block = false;
    // FIXME_HIP HIP does not support syncthreads_or. That's why we need to make
    // num_teams_done __shared__
    // if (__syncthreads_or(num_teams_done == gridDim.x)) {*/
    __syncthreads();
    if (num_teams_done == gridDim.x) {
      is_last_block = true;
      *global_flags = 0;
      ValueInit::init(functor, &value);
      for (int i = threadIdx.y * blockDim.x + threadIdx.x; i < global_elements;
           i += blockDim.x * blockDim.y) {
        ValueJoin::join(functor, &value, &global_team_buffer_element[i]);
      }
      scalar_intra_block_reduction(
          functor, value, false, shared_team_buffer_elements + blockDim.y - 1,
          shared_elements, shared_team_buffer_elements);
    }

    return is_last_block;
  }
};

template <typename FunctorType, typename ArgTag>
struct HIPReductionsFunctor<FunctorType, ArgTag, false> {
  using ValueTraits  = FunctorValueTraits<FunctorType, ArgTag>;
  using ValueJoin    = FunctorValueJoin<FunctorType, ArgTag>;
  using ValueInit    = FunctorValueInit<FunctorType, ArgTag>;
  using ValueOps     = FunctorValueOps<FunctorType, ArgTag>;
  using pointer_type = typename ValueTraits::pointer_type;
  using Scalar       = typename ValueTraits::value_type;

  __device__ static inline void scalar_intra_warp_reduction(
      FunctorType const& functor,
      Scalar* value,           // Contribution
      bool const skip_vector,  // Skip threads if Kokkos vector lanes are not
                               // part of the reduction
      int const width)         // How much of the warp participates
  {
    int const lane_id = (threadIdx.y * blockDim.x + threadIdx.x) %
                        ::Kokkos::Experimental::Impl::HIPTraits::WarpSize;
    for (int delta = skip_vector ? blockDim.x : 1; delta < width; delta *= 2) {
      if (lane_id + delta < ::Kokkos::Experimental::Impl::HIPTraits::WarpSize) {
        ValueJoin::join(functor, value, value + delta);
      }
    }
    *value = *(value - lane_id);
  }

  __device__ static inline void scalar_intra_block_reduction(
      FunctorType const& functor, Scalar value, bool const skip, Scalar* result,
      int const /*shared_elements*/, Scalar* shared_team_buffer_element) {
    int const warp_id = (threadIdx.y * blockDim.x) /
                        ::Kokkos::Experimental::Impl::HIPTraits::WarpSize;
    Scalar* const my_shared_team_buffer_element =
        shared_team_buffer_element + threadIdx.y * blockDim.x + threadIdx.x;
    *my_shared_team_buffer_element = value;
    // Warp Level Reduction, ignoring Kokkos vector entries
    scalar_intra_warp_reduction(
        functor, my_shared_team_buffer_element, skip,
        ::Kokkos::Experimental::Impl::HIPTraits::WarpSize);
    // Wait for every warp to be done before using one warp to do final cross
    // warp reduction
    __syncthreads();

    if (warp_id == 0) {
      const unsigned int delta =
          (threadIdx.y * blockDim.x + threadIdx.x) *
          ::Kokkos::Experimental::Impl::HIPTraits::WarpSize;
      if (delta < blockDim.x * blockDim.y)
        *my_shared_team_buffer_element = shared_team_buffer_element[delta];
      scalar_intra_warp_reduction(
          functor, my_shared_team_buffer_element, false,
          blockDim.x * blockDim.y /
              ::Kokkos::Experimental::Impl::HIPTraits::WarpSize);
      if (threadIdx.x + threadIdx.y == 0) *result = *shared_team_buffer_element;
    }
  }

  __device__ static inline bool scalar_inter_block_reduction(
      FunctorType const& functor,
      ::Kokkos::Experimental::HIP::size_type const block_count,
      ::Kokkos::Experimental::HIP::size_type* const shared_data,
      ::Kokkos::Experimental::HIP::size_type* const global_data,
      ::Kokkos::Experimental::HIP::size_type* const global_flags) {
    Scalar* const global_team_buffer_element =
        reinterpret_cast<Scalar*>(global_data);
    Scalar* const my_global_team_buffer_element =
        global_team_buffer_element + blockIdx.x;
    Scalar* shared_team_buffer_elements =
        reinterpret_cast<Scalar*>(shared_data);
    Scalar value        = shared_team_buffer_elements[threadIdx.y];
    int shared_elements = (blockDim.x * blockDim.y) /
                          ::Kokkos::Experimental::Impl::HIPTraits::WarpSize;
    int global_elements = block_count;
    __syncthreads();

    // Do the scalar reduction inside each block
    scalar_intra_block_reduction(functor, value, true,
                                 my_global_team_buffer_element, shared_elements,
                                 shared_team_buffer_elements);
    __syncthreads();

    // Use the last block that is done to do the do the reduction across the
    // block
    __shared__ unsigned int num_teams_done;
    if (threadIdx.x + threadIdx.y == 0) {
      num_teams_done = Kokkos::atomic_fetch_add(global_flags, 1) + 1;
    }
    bool is_last_block = false;
    // FIXME_HIP HIP does not support syncthreads_or. That's why we need to make
    // num_teams_done __shared__
    // if (__syncthreads_or(num_teams_done == gridDim.x)) {*/
    __syncthreads();
    if (num_teams_done == gridDim.x) {
      is_last_block = true;
      *global_flags = 0;
      ValueInit::init(functor, &value);
      for (int i = threadIdx.y * blockDim.x + threadIdx.x; i < global_elements;
           i += blockDim.x * blockDim.y) {
        ValueJoin::join(functor, &value, &global_team_buffer_element[i]);
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
 *   (a) blockDim.y is a power of two
 *   (b) blockDim.y <= 1024
 *   (c) blockDim.x == blockDim.z == 1
 */

template <bool DoScan, class FunctorType, class ArgTag>
__device__ void hip_intra_block_reduce_scan(
    FunctorType const& functor,
    typename FunctorValueTraits<FunctorType, ArgTag>::pointer_type const
        base_data) {
  using ValueTraits = FunctorValueTraits<FunctorType, ArgTag>;
  using ValueJoin   = FunctorValueJoin<FunctorType, ArgTag>;

  using pointer_type = typename ValueTraits::pointer_type;

  unsigned int const value_count   = ValueTraits::value_count(functor);
  unsigned int const BlockSizeMask = blockDim.y - 1;
  int const WarpMask = Experimental::Impl::HIPTraits::WarpSize - 1;

  // Must have power of two thread count
  if ((blockDim.y - 1) & blockDim.y) {
    Kokkos::abort(
        "HIP::hip_intra_block_reduce_scan requires power-of-two "
        "blockDim.y\n");
  }

  auto block_reduce_step =
      [&functor, value_count](int const R, pointer_type const TD, int const S) {
        if (R > ((1 << S) - 1)) {
          ValueJoin::join(functor, TD, (TD - (value_count << S)));
        }
      };

  {  // Intra-warp reduction:
    const unsigned rtid_intra      = threadIdx.y & WarpMask;
    const pointer_type tdata_intra = base_data + value_count * threadIdx.y;

    block_reduce_step(rtid_intra, tdata_intra, 0);
    block_reduce_step(rtid_intra, tdata_intra, 1);
    block_reduce_step(rtid_intra, tdata_intra, 2);
    block_reduce_step(rtid_intra, tdata_intra, 3);
    block_reduce_step(rtid_intra, tdata_intra, 4);
    block_reduce_step(rtid_intra, tdata_intra, 5);
  }

  __syncthreads();  // Wait for all warps to reduce

  {  // Inter-warp reduce-scan by a single warp to avoid extra synchronizations
    unsigned int const rtid_inter =
        ((threadIdx.y + 1) << Experimental::Impl::HIPTraits::WarpIndexShift) -
        1;

    if (rtid_inter < blockDim.y) {
      pointer_type const tdata_inter = base_data + value_count * rtid_inter;

      if ((1 << 6) < BlockSizeMask) {
        block_reduce_step(rtid_inter, tdata_inter, 6);
      }
      if ((1 << 7) < BlockSizeMask) {
        block_reduce_step(rtid_inter, tdata_inter, 7);
      }
      if ((1 << 8) < BlockSizeMask) {
        block_reduce_step(rtid_inter, tdata_inter, 8);
      }
      if ((1 << 9) < BlockSizeMask) {
        block_reduce_step(rtid_inter, tdata_inter, 9);
      }
      if ((1 << 10) < BlockSizeMask) {
        block_reduce_step(rtid_inter, tdata_inter, 10);
      }
    }
  }

  __syncthreads();  // Wait for inter-warp reduce-scan to complete

  if (DoScan) {
    // Update all the values for the respective warps (except for the last one)
    // by adding from the last value of the previous warp.
    if (threadIdx.y >= Experimental::Impl::HIPTraits::WarpSize &&
        (threadIdx.y & WarpMask) !=
            Experimental::Impl::HIPTraits::WarpSize - 1) {
      const int offset_to_previous_warp_total = (threadIdx.y & (~WarpMask)) - 1;
      ValueJoin::join(functor, base_data + value_count * threadIdx.y,
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

template <bool DoScan, class FunctorType, class ArgTag>
__device__ bool hip_single_inter_block_reduce_scan_impl(
    FunctorType const& functor,
    ::Kokkos::Experimental::HIP::size_type const block_id,
    ::Kokkos::Experimental::HIP::size_type const block_count,
    ::Kokkos::Experimental::HIP::size_type* const shared_data,
    ::Kokkos::Experimental::HIP::size_type* const global_data,
    ::Kokkos::Experimental::HIP::size_type* const global_flags) {
  using size_type   = ::Kokkos::Experimental::HIP::size_type;
  using ValueTraits = FunctorValueTraits<FunctorType, ArgTag>;
  using ValueJoin   = FunctorValueJoin<FunctorType, ArgTag>;
  using ValueInit   = FunctorValueInit<FunctorType, ArgTag>;
  using ValueOps    = FunctorValueOps<FunctorType, ArgTag>;

  using pointer_type = typename ValueTraits::pointer_type;

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

  integral_nonzero_constant<size_type, ValueTraits::StaticValueSize /
                                           sizeof(size_type)> const
      word_count(ValueTraits::value_size(functor) / sizeof(size_type));

  // Reduce the accumulation for the entire block.
  hip_intra_block_reduce_scan<false, FunctorType, ArgTag>(
      functor, pointer_type(shared_data));

  {
    // Write accumulation total to global scratch space.
    // Accumulation total is the last thread's data.
    size_type* const shared = shared_data + word_count.value * BlockSizeMask;
    size_type* const global = global_data + word_count.value * block_id;

    for (size_t i = threadIdx.y; i < word_count.value; i += blockDim.y) {
      global[i] = shared[i];
    }
  }

  // Contributing blocks note that their contribution has been completed via an
  // atomic-increment flag If this block is not the last block to contribute to
  // this group then the block is done.
  // FIXME_HIP __syncthreads_or is not supported by HIP yet.
  // const bool is_last_block = !__syncthreads_or(
  //    threadIdx.y
  //        ? 0
  //        : (1 + atomicInc(global_flags, block_count - 1) < block_count));
  __shared__ int n_done;
  n_done = 0;
  __syncthreads();
  if (threadIdx.y == 0) {
    n_done = 1 + atomicInc(global_flags, block_count - 1);
  }
  __syncthreads();
  bool const is_last_block = (n_done == static_cast<int>(block_count));

  if (is_last_block) {
    size_type const b = (static_cast<long long int>(block_count) *
                         static_cast<long long int>(threadIdx.y)) >>
                        BlockSizeShift;
    size_type const e = (static_cast<long long int>(block_count) *
                         static_cast<long long int>(threadIdx.y + 1)) >>
                        BlockSizeShift;

    {
      void* const shared_ptr = shared_data + word_count.value * threadIdx.y;
      /* reference_type shared_value = */ ValueInit::init(functor, shared_ptr);

      for (size_type i = b; i < e; ++i) {
        ValueJoin::join(functor, shared_ptr,
                        global_data + word_count.value * i);
      }
    }

    hip_intra_block_reduce_scan<DoScan, FunctorType, ArgTag>(
        functor, pointer_type(shared_data));

    if (DoScan) {
      size_type* const shared_value =
          shared_data +
          word_count.value * (threadIdx.y ? threadIdx.y - 1 : blockDim.y);

      if (!threadIdx.y) {
        ValueInit::init(functor, shared_value);
      }

      // Join previous inclusive scan value to each member
      for (size_type i = b; i < e; ++i) {
        size_type* const global_value = global_data + word_count.value * i;
        ValueJoin::join(functor, shared_value, global_value);
        ValueOps::copy(functor, global_value, shared_value);
      }
    }
  }

  return is_last_block;
}

template <bool DoScan, typename FunctorType, typename ArgTag>
__device__ bool hip_single_inter_block_reduce_scan(
    FunctorType const& functor,
    ::Kokkos::Experimental::HIP::size_type const block_id,
    ::Kokkos::Experimental::HIP::size_type const block_count,
    ::Kokkos::Experimental::HIP::size_type* const shared_data,
    ::Kokkos::Experimental::HIP::size_type* const global_data,
    ::Kokkos::Experimental::HIP::size_type* const global_flags) {
  using ValueTraits = FunctorValueTraits<FunctorType, ArgTag>;
  // If we are doing a reduction and StaticValueSize is true, we use the
  // reduction-only path. Otherwise, we use the common path between reduction
  // and scan.
  if (!DoScan && static_cast<bool>(ValueTraits::StaticValueSize))
    // FIXME_HIP_PERFORMANCE I don't know where 16 comes from. This inequality
    // determines if we use shared memory (false) or shuffle (true)
    return Kokkos::Impl::HIPReductionsFunctor<
        FunctorType, ArgTag, (ValueTraits::StaticValueSize > 16)>::
        scalar_inter_block_reduction(functor, block_count, shared_data,
                                     global_data, global_flags);
  else {
    return hip_single_inter_block_reduce_scan_impl<DoScan, FunctorType, ArgTag>(
        functor, block_id, block_count, shared_data, global_data, global_flags);
  }
}

// Size in bytes required for inter block reduce or scan
template <bool DoScan, class FunctorType, class ArgTag>
inline unsigned hip_single_inter_block_reduce_scan_shmem(
    const FunctorType& functor, const unsigned BlockSize) {
  return (BlockSize + 2) *
         Impl::FunctorValueTraits<FunctorType, ArgTag>::value_size(functor);
}

}  // namespace Impl
}  // namespace Kokkos

#endif

#endif
