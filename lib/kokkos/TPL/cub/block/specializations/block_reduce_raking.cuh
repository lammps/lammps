/******************************************************************************
 * Copyright (c) 2011, Duane Merrill.  All rights reserved.
 * Copyright (c) 2011-2013, NVIDIA CORPORATION.  All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the NVIDIA CORPORATION nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ******************************************************************************/

/**
 * \file
 * cub::BlockReduceRaking provides raking-based methods of parallel reduction across a CUDA threadblock
 */

#pragma once

#include "../../block/block_raking_layout.cuh"
#include "../../warp/warp_reduce.cuh"
#include "../../thread/thread_reduce.cuh"
#include "../../util_namespace.cuh"

/// Optional outer namespace(s)
CUB_NS_PREFIX

/// CUB namespace
namespace cub {


/**
 * \brief BlockReduceRaking provides raking-based methods of parallel reduction across a CUDA threadblock
 */
template <
    typename    T,              ///< Data type being reduced
    int         BLOCK_THREADS>  ///< The thread block size in threads
struct BlockReduceRaking
{
    /// Layout type for padded threadblock raking grid
    typedef BlockRakingLayout<T, BLOCK_THREADS, 1> BlockRakingLayout;

    ///  WarpReduce utility type
    typedef typename WarpReduce<T, 1, BlockRakingLayout::RAKING_THREADS>::InternalWarpReduce WarpReduce;

    /// Constants
    enum
    {
        /// Number of raking threads
        RAKING_THREADS = BlockRakingLayout::RAKING_THREADS,

        /// Number of raking elements per warp synchronous raking thread
        SEGMENT_LENGTH = BlockRakingLayout::SEGMENT_LENGTH,

        /// Cooperative work can be entirely warp synchronous
        WARP_SYNCHRONOUS = (RAKING_THREADS == BLOCK_THREADS),

        /// Whether or not warp-synchronous reduction should be unguarded (i.e., the warp-reduction elements is a power of two
        WARP_SYNCHRONOUS_UNGUARDED = ((RAKING_THREADS & (RAKING_THREADS - 1)) == 0),

        /// Whether or not accesses into smem are unguarded
        RAKING_UNGUARDED = BlockRakingLayout::UNGUARDED,

    };


    /// Shared memory storage layout type
    struct _TempStorage
    {
        typename WarpReduce::TempStorage            warp_storage;        ///< Storage for warp-synchronous reduction
        typename BlockRakingLayout::TempStorage     raking_grid;         ///< Padded threadblock raking grid
    };


    /// Alias wrapper allowing storage to be unioned
    struct TempStorage : Uninitialized<_TempStorage> {};


    // Thread fields
    _TempStorage &temp_storage;
    int linear_tid;


    /// Constructor
    __device__ __forceinline__ BlockReduceRaking(
        TempStorage &temp_storage,
        int linear_tid)
    :
        temp_storage(temp_storage.Alias()),
        linear_tid(linear_tid)
    {}


    /// Computes a threadblock-wide reduction using addition (+) as the reduction operator. The first num_valid threads each contribute one reduction partial.  The return value is only valid for thread<sub>0</sub>.
    template <bool FULL_TILE>
    __device__ __forceinline__ T Sum(
        T                   partial,            ///< [in] Calling thread's input partial reductions
        int                 num_valid)          ///< [in] Number of valid elements (may be less than BLOCK_THREADS)
    {
        cub::Sum reduction_op;

        if (WARP_SYNCHRONOUS)
        {
            // Short-circuit directly to warp synchronous reduction (unguarded if active threads is a power-of-two)
            partial = WarpReduce(temp_storage.warp_storage, 0, linear_tid).template Sum<FULL_TILE, SEGMENT_LENGTH>(
                partial,
                num_valid);
        }
        else
        {
            // Place partial into shared memory grid.
            *BlockRakingLayout::PlacementPtr(temp_storage.raking_grid, linear_tid) = partial;

            __syncthreads();

            // Reduce parallelism to one warp
            if (linear_tid < RAKING_THREADS)
            {
                // Raking reduction in grid
                T *raking_segment = BlockRakingLayout::RakingPtr(temp_storage.raking_grid, linear_tid);
                partial = raking_segment[0];

                #pragma unroll
                for (int ITEM = 1; ITEM < SEGMENT_LENGTH; ITEM++)
                {
                    // Update partial if addend is in range
                    if ((FULL_TILE && RAKING_UNGUARDED) || ((linear_tid * SEGMENT_LENGTH) + ITEM < num_valid))
                    {
                        partial = reduction_op(partial, raking_segment[ITEM]);
                    }
                }

                partial = WarpReduce(temp_storage.warp_storage, 0, linear_tid).template Sum<FULL_TILE && RAKING_UNGUARDED, SEGMENT_LENGTH>(
                    partial,
                    num_valid);
            }
        }

        return partial;
    }


    /// Computes a threadblock-wide reduction using the specified reduction operator. The first num_valid threads each contribute one reduction partial.  The return value is only valid for thread<sub>0</sub>.
    template <
        bool                FULL_TILE,
        typename            ReductionOp>
    __device__ __forceinline__ T Reduce(
        T                   partial,            ///< [in] Calling thread's input partial reductions
        int                 num_valid,          ///< [in] Number of valid elements (may be less than BLOCK_THREADS)
        ReductionOp         reduction_op)       ///< [in] Binary reduction operator
    {
        if (WARP_SYNCHRONOUS)
        {
            // Short-circuit directly to warp synchronous reduction (unguarded if active threads is a power-of-two)
            partial = WarpReduce(temp_storage.warp_storage, 0, linear_tid).template Reduce<FULL_TILE, SEGMENT_LENGTH>(
                partial,
                num_valid,
                reduction_op);
        }
        else
        {
            // Place partial into shared memory grid.
            *BlockRakingLayout::PlacementPtr(temp_storage.raking_grid, linear_tid) = partial;

            __syncthreads();

            // Reduce parallelism to one warp
            if (linear_tid < RAKING_THREADS)
            {
                // Raking reduction in grid
                T *raking_segment = BlockRakingLayout::RakingPtr(temp_storage.raking_grid, linear_tid);
                partial = raking_segment[0];

                #pragma unroll
                for (int ITEM = 1; ITEM < SEGMENT_LENGTH; ITEM++)
                {
                    // Update partial if addend is in range
                    if ((FULL_TILE && RAKING_UNGUARDED) || ((linear_tid * SEGMENT_LENGTH) + ITEM < num_valid))
                    {
                        partial = reduction_op(partial, raking_segment[ITEM]);
                    }
                }

                partial = WarpReduce(temp_storage.warp_storage, 0, linear_tid).template Reduce<FULL_TILE && RAKING_UNGUARDED, SEGMENT_LENGTH>(
                    partial,
                    num_valid,
                    reduction_op);
            }
        }

        return partial;
    }

};

}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)

