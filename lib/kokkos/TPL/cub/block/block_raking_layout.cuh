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
 * cub::BlockRakingLayout provides a conflict-free shared memory layout abstraction for warp-raking across thread block data.
 */


#pragma once

#include "../util_macro.cuh"
#include "../util_arch.cuh"
#include "../util_namespace.cuh"

/// Optional outer namespace(s)
CUB_NS_PREFIX

/// CUB namespace
namespace cub {

/**
 * \brief BlockRakingLayout provides a conflict-free shared memory layout abstraction for raking across thread block data.    ![](raking.png)
 * \ingroup BlockModule
 *
 * \par Overview
 * This type facilitates a shared memory usage pattern where a block of CUDA
 * threads places elements into shared memory and then reduces the active
 * parallelism to one "raking" warp of threads for serially aggregating consecutive
 * sequences of shared items.  Padding is inserted to eliminate bank conflicts
 * (for most data types).
 *
 * \tparam T                    The data type to be exchanged.
 * \tparam BLOCK_THREADS        The thread block size in threads.
 * \tparam BLOCK_STRIPS         When strip-mining, the number of threadblock-strips per tile
 */
template <
    typename    T,
    int         BLOCK_THREADS,
    int         BLOCK_STRIPS = 1>
struct BlockRakingLayout
{
    //---------------------------------------------------------------------
    // Constants and typedefs
    //---------------------------------------------------------------------

    enum
    {
        /// The total number of elements that need to be cooperatively reduced
        SHARED_ELEMENTS =
            BLOCK_THREADS * BLOCK_STRIPS,

        /// Maximum number of warp-synchronous raking threads
        MAX_RAKING_THREADS =
            CUB_MIN(BLOCK_THREADS, PtxArchProps::WARP_THREADS),

        /// Number of raking elements per warp-synchronous raking thread (rounded up)
        SEGMENT_LENGTH =
            (SHARED_ELEMENTS + MAX_RAKING_THREADS - 1) / MAX_RAKING_THREADS,

        /// Never use a raking thread that will have no valid data (e.g., when BLOCK_THREADS is 62 and SEGMENT_LENGTH is 2, we should only use 31 raking threads)
        RAKING_THREADS =
            (SHARED_ELEMENTS + SEGMENT_LENGTH - 1) / SEGMENT_LENGTH,

        /// Pad each segment length with one element if it evenly divides the number of banks
        SEGMENT_PADDING =
            (PtxArchProps::SMEM_BANKS % SEGMENT_LENGTH == 0) ? 1 : 0,

        /// Total number of elements in the raking grid
        GRID_ELEMENTS =
            RAKING_THREADS * (SEGMENT_LENGTH + SEGMENT_PADDING),

        /// Whether or not we need bounds checking during raking (the number of reduction elements is not a multiple of the warp size)
        UNGUARDED =
            (SHARED_ELEMENTS % RAKING_THREADS == 0),
    };


    /**
     * \brief Shared memory storage type
     */
    typedef T TempStorage[BlockRakingLayout::GRID_ELEMENTS];


    /**
     * \brief Returns the location for the calling thread to place data into the grid
     */
    static __device__ __forceinline__ T* PlacementPtr(
        TempStorage &temp_storage,
        int linear_tid,
        int block_strip = 0)
    {
        // Offset for partial
        unsigned int offset = (block_strip * BLOCK_THREADS) + linear_tid;

        // Add in one padding element for every segment
        if (SEGMENT_PADDING > 0)
        {
            offset += offset / SEGMENT_LENGTH;
        }

        // Incorporating a block of padding partials every shared memory segment
        return temp_storage + offset;
    }


    /**
     * \brief Returns the location for the calling thread to begin sequential raking
     */
    static __device__ __forceinline__ T* RakingPtr(
        TempStorage &temp_storage,
        int linear_tid)
    {
        return temp_storage + (linear_tid * (SEGMENT_LENGTH + SEGMENT_PADDING));
    }
};

}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)

