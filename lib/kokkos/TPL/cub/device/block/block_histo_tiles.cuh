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
 * cub::BlockHistogramTiles implements a stateful abstraction of CUDA thread blocks for histogramming multiple tiles as part of device-wide histogram.
 */

#pragma once

#include <iterator>

#include "specializations/block_histo_tiles_gatomic.cuh"
#include "specializations/block_histo_tiles_satomic.cuh"
#include "specializations/block_histo_tiles_sort.cuh"
#include "../../util_type.cuh"
#include "../../grid/grid_mapping.cuh"
#include "../../grid/grid_even_share.cuh"
#include "../../grid/grid_queue.cuh"
#include "../../util_namespace.cuh"

/// Optional outer namespace(s)
CUB_NS_PREFIX

/// CUB namespace
namespace cub {


/******************************************************************************
 * Algorithmic variants
 ******************************************************************************/


/**
 * \brief BlockHistogramTilesAlgorithm enumerates alternative algorithms for BlockHistogramTiles.
 */
enum BlockHistogramTilesAlgorithm
{

    /**
     * \par Overview
     * A two-kernel approach in which:
     * -# Thread blocks in the first kernel aggregate their own privatized
     *    histograms using block-wide sorting (see BlockHistogramAlgorithm::BLOCK_HISTO_SORT).
     * -# A single thread block in the second kernel reduces them into the output histogram(s).
     *
     * \par Performance Considerations
     * Delivers consistent throughput regardless of sample bin distribution.
     *
     * However, because histograms are privatized in shared memory, a large
     * number of bins (e.g., thousands) may adversely affect occupancy and
     * performance (or even the ability to launch).
     */
    GRID_HISTO_SORT,


    /**
     * \par Overview
     * A two-kernel approach in which:
     * -# Thread blocks in the first kernel aggregate their own privatized
     *    histograms using shared-memory \p atomicAdd().
     * -# A single thread block in the second kernel reduces them into the
     *    output histogram(s).
     *
     * \par Performance Considerations
     * Performance is strongly tied to the hardware implementation of atomic
     * addition, and may be significantly degraded for non uniformly-random
     * input distributions where many concurrent updates are likely to be
     * made to the same bin counter.
     *
     * However, because histograms are privatized in shared memory, a large
     * number of bins (e.g., thousands) may adversely affect occupancy and
     * performance (or even the ability to launch).
     */
    GRID_HISTO_SHARED_ATOMIC,


    /**
     * \par Overview
     * A single-kernel approach in which thread blocks update the output histogram(s) directly
     * using global-memory \p atomicAdd().
     *
     * \par Performance Considerations
     * Performance is strongly tied to the hardware implementation of atomic
     * addition, and may be significantly degraded for non uniformly-random
     * input distributions where many concurrent updates are likely to be
     * made to the same bin counter.
     *
     * Performance is not significantly impacted when computing histograms having large
     * numbers of bins (e.g., thousands).
     */
    GRID_HISTO_GLOBAL_ATOMIC,

};


/******************************************************************************
 * Tuning policy
 ******************************************************************************/

/**
 * Tuning policy for BlockHistogramTiles
 */
template <
    int                             _BLOCK_THREADS,
    int                             _ITEMS_PER_THREAD,
    BlockHistogramTilesAlgorithm    _GRID_ALGORITHM,
    GridMappingStrategy             _GRID_MAPPING,
    int                             _SM_OCCUPANCY>
struct BlockHistogramTilesPolicy
{
    enum
    {
        BLOCK_THREADS       = _BLOCK_THREADS,
        ITEMS_PER_THREAD    = _ITEMS_PER_THREAD,
        SM_OCCUPANCY        = _SM_OCCUPANCY,
    };

    static const BlockHistogramTilesAlgorithm   GRID_ALGORITHM      = _GRID_ALGORITHM;
    static const GridMappingStrategy            GRID_MAPPING        = _GRID_MAPPING;
};



/******************************************************************************
 * Thread block abstractions
 ******************************************************************************/


/**
 * Implements a stateful abstraction of CUDA thread blocks for histogramming multiple tiles as part of device-wide histogram using global atomics
 */
template <
    typename    BlockHistogramTilesPolicy,          ///< Tuning policy
    int         BINS,                           ///< Number of histogram bins per channel
    int         CHANNELS,                       ///< Number of channels interleaved in the input data (may be greater than the number of active channels being histogrammed)
    int         ACTIVE_CHANNELS,                ///< Number of channels actively being histogrammed
    typename    InputIteratorRA,                ///< The input iterator type (may be a simple pointer type).  Must have a value type that can be cast as an integer in the range [0..BINS-1]
    typename    HistoCounter,                   ///< Integral type for counting sample occurrences per histogram bin
    typename    SizeT>                          ///< Integer type for offsets
struct BlockHistogramTiles
{
    //---------------------------------------------------------------------
    // Types and constants
    //---------------------------------------------------------------------

    // Histogram grid algorithm
    static const BlockHistogramTilesAlgorithm GRID_ALGORITHM = BlockHistogramTilesPolicy::GRID_ALGORITHM;

    // Alternative internal implementation types
    typedef BlockHistogramTilesSort<            BlockHistogramTilesPolicy, BINS, CHANNELS, ACTIVE_CHANNELS, InputIteratorRA, HistoCounter, SizeT>   BlockHistogramTilesSortT;
    typedef BlockHistogramTilesSharedAtomic<    BlockHistogramTilesPolicy, BINS, CHANNELS, ACTIVE_CHANNELS, InputIteratorRA, HistoCounter, SizeT>   BlockHistogramTilesSharedAtomicT;
    typedef BlockHistogramTilesGlobalAtomic<    BlockHistogramTilesPolicy, BINS, CHANNELS, ACTIVE_CHANNELS, InputIteratorRA, HistoCounter, SizeT>   BlockHistogramTilesGlobalAtomicT;

    // Internal block sweep histogram type
    typedef typename If<(GRID_ALGORITHM == GRID_HISTO_SORT),
        BlockHistogramTilesSortT,
        typename If<(GRID_ALGORITHM == GRID_HISTO_SHARED_ATOMIC),
            BlockHistogramTilesSharedAtomicT,
            BlockHistogramTilesGlobalAtomicT>::Type>::Type InternalBlockDelegate;

    enum
    {
        TILE_ITEMS = InternalBlockDelegate::TILE_ITEMS,
    };


    // Temporary storage type
    typedef typename InternalBlockDelegate::TempStorage TempStorage;

    //---------------------------------------------------------------------
    // Per-thread fields
    //---------------------------------------------------------------------

    // Internal block delegate
    InternalBlockDelegate internal_delegate;


    //---------------------------------------------------------------------
    // Interface
    //---------------------------------------------------------------------

    /**
     * Constructor
     */
    __device__ __forceinline__ BlockHistogramTiles(
        TempStorage         &temp_storage,                                  ///< Reference to temp_storage
        InputIteratorRA     d_in,                                           ///< Input data to reduce
        HistoCounter*       (&d_out_histograms)[ACTIVE_CHANNELS])           ///< Reference to output histograms
    :
        internal_delegate(temp_storage, d_in, d_out_histograms)
    {}


    /**
     * \brief Reduce a consecutive segment of input tiles
     */
    __device__ __forceinline__ void ConsumeTiles(
        SizeT   block_offset,                       ///< [in] Threadblock begin offset (inclusive)
        SizeT   block_oob)                          ///< [in] Threadblock end offset (exclusive)
    {
        // Consume subsequent full tiles of input
        while (block_offset + TILE_ITEMS <= block_oob)
        {
            internal_delegate.ConsumeTile<true>(block_offset);
            block_offset += TILE_ITEMS;
        }

        // Consume a partially-full tile
        if (block_offset < block_oob)
        {
            int valid_items = block_oob - block_offset;
            internal_delegate.ConsumeTile<false>(block_offset, valid_items);
        }

        // Aggregate output
        internal_delegate.AggregateOutput();
    }


    /**
     * Reduce a consecutive segment of input tiles
     */
    __device__ __forceinline__ void ConsumeTiles(
        SizeT                               num_items,          ///< [in] Total number of global input items
        GridEvenShare<SizeT>                &even_share,        ///< [in] GridEvenShare descriptor
        GridQueue<SizeT>                    &queue,             ///< [in,out] GridQueue descriptor
        Int2Type<GRID_MAPPING_EVEN_SHARE>   is_even_share)      ///< [in] Marker type indicating this is an even-share mapping
    {
        even_share.BlockInit();
        ConsumeTiles(even_share.block_offset, even_share.block_oob);
    }


    /**
     * Dequeue and reduce tiles of items as part of a inter-block scan
     */
    __device__ __forceinline__ void ConsumeTiles(
        int                 num_items,          ///< Total number of input items
        GridQueue<SizeT>    queue)              ///< Queue descriptor for assigning tiles of work to thread blocks
    {
        // Shared block offset
        __shared__ SizeT shared_block_offset;

        // We give each thread block at least one tile of input.
        SizeT block_offset      = blockIdx.x * TILE_ITEMS;
        SizeT even_share_base   = gridDim.x * TILE_ITEMS;

        // Process full tiles of input
        while (block_offset + TILE_ITEMS <= num_items)
        {
            internal_delegate.ConsumeTile<true>(block_offset);

            // Dequeue up to TILE_ITEMS
            if (threadIdx.x == 0)
                shared_block_offset = queue.Drain(TILE_ITEMS) + even_share_base;

            __syncthreads();

            block_offset = shared_block_offset;

            __syncthreads();
        }

        // Consume a partially-full tile
        if (block_offset < num_items)
        {
            int valid_items = num_items - block_offset;
            internal_delegate.ConsumeTile<false>(block_offset, valid_items);
        }

        // Aggregate output
        internal_delegate.AggregateOutput();
    }


    /**
     * Dequeue and reduce tiles of items as part of a inter-block scan
     */
    __device__ __forceinline__ void ConsumeTiles(
        SizeT                               num_items,          ///< [in] Total number of global input items
        GridEvenShare<SizeT>                &even_share,        ///< [in] GridEvenShare descriptor
        GridQueue<SizeT>                    &queue,             ///< [in,out] GridQueue descriptor
        Int2Type<GRID_MAPPING_DYNAMIC>      is_dynamic)         ///< [in] Marker type indicating this is a dynamic mapping
    {
        ConsumeTiles(num_items, queue);
    }


};




}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)

