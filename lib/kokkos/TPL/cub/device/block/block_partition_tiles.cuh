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
 * cub::BlockPartitionTiles implements a stateful abstraction of CUDA thread blocks for participating in device-wide list partitioning.
 */

#pragma once

#include <iterator>

#include "scan_tiles_types.cuh"
#include "../../thread/thread_operators.cuh"
#include "../../block/block_load.cuh"
#include "../../block/block_store.cuh"
#include "../../block/block_scan.cuh"
#include "../../grid/grid_queue.cuh"
#include "../../util_vector.cuh"
#include "../../util_namespace.cuh"

/// Optional outer namespace(s)
CUB_NS_PREFIX

/// CUB namespace
namespace cub {


/******************************************************************************
 * Tuning policy types
 ******************************************************************************/

/**
 * Tuning policy for BlockPartitionTiles
 */
template <
    int                         _PARTITIONS,
    int                         _BLOCK_THREADS,
    int                         _ITEMS_PER_THREAD,
    PtxLoadModifier             _LOAD_MODIFIER,
    BlockScanAlgorithm          _SCAN_ALGORITHM>
struct BlockPartitionTilesPolicy
{
    enum
    {
        PARTITIONS              = _PARTITIONS,
        BLOCK_THREADS           = _BLOCK_THREADS,
        ITEMS_PER_THREAD        = _ITEMS_PER_THREAD,
    };

    static const PtxLoadModifier        LOAD_MODIFIER       = _LOAD_MODIFIER;
    static const BlockScanAlgorithm     SCAN_ALGORITHM      = _SCAN_ALGORITHM;
};



/**
 * Tuple type for scanning partition membership flags
 */
template <
    typename    SizeT,
    int         PARTITIONS>
struct PartitionScanTuple;


/**
 * Tuple type for scanning partition membership flags (specialized for 1 output partition)
 */
template <typename SizeT>
struct PartitionScanTuple<SizeT, 1> : VectorHelper<SizeT, 1>::Type
{
    __device__ __forceinline__ PartitionScanTuple operator+(const PartitionScanTuple &other)
    {
        PartitionScanTuple retval;
        retval.x = x + other.x;
        return retval;
    }

    template <typename PredicateOp, typename T>
    __device__ __forceinline__ void SetFlags(PredicateOp pred_op, T val)
    {
        this->x = pred_op(val);
    }

    template <typename PredicateOp, typename T, typename OutputIteratorRA, SizeT num_items>
    __device__ __forceinline__ void Scatter(PredicateOp pred_op, T val, OutputIteratorRA d_out, SizeT num_items)
    {
        if (pred_op(val))
            d_out[this->x - 1] = val;
    }

};


/**
 * Tuple type for scanning partition membership flags (specialized for 2 output partitions)
 */
template <typename SizeT>
struct PartitionScanTuple<SizeT, 2> : VectorHelper<SizeT, 2>::Type
{
    __device__ __forceinline__ PartitionScanTuple operator+(const PartitionScanTuple &other)
    {
        PartitionScanTuple retval;
        retval.x = x + other.x;
        retval.y = y + other.y;
        return retval;
    }

    template <typename PredicateOp, typename T>
    __device__ __forceinline__ void SetFlags(PredicateOp pred_op, T val)
    {
        bool pred = pred_op(val);
        this->x = pred;
        this->y = !pred;
    }

    template <typename PredicateOp, typename T, typename OutputIteratorRA, SizeT num_items>
    __device__ __forceinline__ void Scatter(PredicateOp pred_op, T val, OutputIteratorRA d_out, SizeT num_items)
    {
        SizeT scatter_offset = (pred_op(val)) ?
            this->x - 1 :
            num_items - this->y;

        d_out[scatter_offset] = val;
    }
};




/******************************************************************************
 * Thread block abstractions
 ******************************************************************************/

/**
 * \brief BlockPartitionTiles implements a stateful abstraction of CUDA thread blocks for participating in device-wide list partitioning.
 *
 * Implements a single-pass "domino" strategy with adaptive prefix lookback.
 */
template <
    typename BlockPartitionTilesPolicy, ///< Tuning policy
    typename InputIteratorRA,           ///< Input iterator type
    typename OutputIteratorRA,          ///< Output iterator type
    typename PredicateOp,               ///< Partition predicate functor type
    typename SizeT>                     ///< Offset integer type
struct BlockPartitionTiles
{
    //---------------------------------------------------------------------
    // Types and constants
    //---------------------------------------------------------------------

    // Constants
    enum
    {
        PARTITIONS          = BlockPartitionTilesPolicy::PARTITIONS,
        BLOCK_THREADS       = BlockPartitionTilesPolicy::BLOCK_THREADS,
        ITEMS_PER_THREAD    = BlockPartitionTilesPolicy::ITEMS_PER_THREAD,
        TILE_ITEMS          = BLOCK_THREADS * ITEMS_PER_THREAD,
    };

    // Load modifier
    static const PtxLoadModifier LOAD_MODIFIER = BlockPartitionTilesPolicy::LOAD_MODIFIER;

    // Data type of input iterator
    typedef typename std::iterator_traits<InputIteratorRA>::value_type T;

    // Tuple type for scanning partition membership flags
    typedef PartitionScanTuple<SizeT, PARTITIONS> PartitionScanTuple;

    // Tile status descriptor type
    typedef ScanTileDescriptor<PartitionScanTuple> ScanTileDescriptorT;

    // Block scan type for scanning membership flag scan_tuples
    typedef BlockScan<
        PartitionScanTuple,
        BlockPartitionTilesPolicy::BLOCK_THREADS,
        BlockPartitionTilesPolicy::SCAN_ALGORITHM> BlockScanT;

    // Callback type for obtaining inter-tile prefix during block scan
    typedef DeviceScanBlockPrefixOp<PartitionScanTuple, Sum> InterblockPrefixOp;

    // Shared memory type for this threadblock
    struct TempStorage
    {
        typename InterblockPrefixOp::TempStorage    prefix;         // Smem needed for cooperative prefix callback
        typename BlockScanT::TempStorage            scan;           // Smem needed for tile scanning
        SizeT                                       tile_idx;       // Shared tile index
    };


    //---------------------------------------------------------------------
    // Per-thread fields
    //---------------------------------------------------------------------

    TempStorage                 &temp_storage;      ///< Reference to temp_storage
    InputIteratorRA             d_in;               ///< Input data
    OutputIteratorRA            d_out;              ///< Output data
    ScanTileDescriptorT         *d_tile_status;     ///< Global list of tile status
    PredicateOp                 pred_op;            ///< Unary predicate operator indicating membership in the first partition
    SizeT                       num_items;          ///< Total number of input items


    //---------------------------------------------------------------------
    // Constructor
    //---------------------------------------------------------------------

    // Constructor
    __device__ __forceinline__
    BlockPartitionTiles(
        TempStorage                 &temp_storage,      ///< Reference to temp_storage
        InputIteratorRA             d_in,               ///< Input data
        OutputIteratorRA            d_out,              ///< Output data
        ScanTileDescriptorT         *d_tile_status,     ///< Global list of tile status
        PredicateOp                 pred_op,            ///< Unary predicate operator indicating membership in the first partition
        SizeT                       num_items)          ///< Total number of input items
    :
        temp_storage(temp_storage.Alias()),
        d_in(d_in),
        d_out(d_out),
        d_tile_status(d_tile_status),
        pred_op(pred_op),
        num_items(num_items)
    {}


    //---------------------------------------------------------------------
    // Domino scan
    //---------------------------------------------------------------------

    /**
     * Process a tile of input
     */
    template <bool FULL_TILE>
    __device__ __forceinline__ void ConsumeTile(
        int                 tile_idx,           ///< Tile index
        SizeT               block_offset,       ///< Tile offset
        PartitionScanTuple  &partition_ends)    ///< Running total
    {
        T                   items[ITEMS_PER_THREAD];
        PartitionScanTuple  scan_tuples[ITEMS_PER_THREAD];

        // Load items
        int valid_items = num_items - block_offset;
        if (FULL_TILE)
            LoadStriped<LOAD_MODIFIER, BLOCK_THREADS>(threadIdx.x, d_in + block_offset, items);
        else
            LoadStriped<LOAD_MODIFIER, BLOCK_THREADS>(threadIdx.x, d_in + block_offset, items, valid_items);

        // Prevent hoisting
//        __syncthreads();
//        __threadfence_block();

        // Set partition membership flags in scan scan_tuples
        #pragma unroll
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ++ITEM)
        {
            scan_tuples[ITEM].SetFlags(pred_op, items[ITEM]);
        }

        // Perform inclusive scan over scan scan_tuples
        PartitionScanTuple block_aggregate;
        if (tile_idx == 0)
        {
            BlockScanT(temp_storage.scan).InclusiveScan(scan_tuples, scan_tuples, Sum(), block_aggregate);
            partition_ends = block_aggregate;

            // Update tile status if there are successor tiles
            if (FULL_TILE && (threadIdx.x == 0))
                ScanTileDescriptorT::SetPrefix(d_tile_status, block_aggregate);
        }
        else
        {
            InterblockPrefixOp prefix_op(d_tile_status, temp_storage.prefix, Sum(), tile_idx);
            BlockScanT(temp_storage.scan).InclusiveScan(scan_tuples, scan_tuples, Sum(), block_aggregate, prefix_op);
            partition_ends = prefix_op.inclusive_prefix;
        }

        // Scatter items
        #pragma unroll
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ++ITEM)
        {
            // Scatter if not out-of-bounds
            if (FULL_TILE || (threadIdx.x + (ITEM * BLOCK_THREADS) < valid_items))
            {
                scan_tuples[ITEM].Scatter(pred_op, items[ITEM], d_out, num_items);
            }
        }
    }


    /**
     * Dequeue and scan tiles of items as part of a domino scan
     */
    __device__ __forceinline__ void ConsumeTiles(
        GridQueue<int>      queue,              ///< [in] Queue descriptor for assigning tiles of work to thread blocks
        SizeT               num_tiles,          ///< [in] Total number of input tiles
        PartitionScanTuple  &partition_ends,    ///< [out] Running partition end offsets
        bool                &is_last_tile)      ///< [out] Whether or not this block handled the last tile (i.e., partition_ends is valid for the entire input)
    {
#if CUB_PTX_ARCH < 200

        // No concurrent kernels allowed and blocks are launched in increasing order, so just assign one tile per block (up to 65K blocks)
        int     tile_idx        = blockIdx.x;
        SizeT   block_offset    = SizeT(TILE_ITEMS) * tile_idx;

        if (block_offset + TILE_ITEMS <= num_items)
        {
            ConsumeTile<true>(tile_idx, block_offset, partition_ends);
        }
        else if (block_offset < num_items)
        {
            ConsumeTile<false>(tile_idx, block_offset, partition_ends);
        }
        is_last_tile = (tile_idx == num_tiles - 1);

#else

        // Get first tile
        if (threadIdx.x == 0)
            temp_storage.tile_idx = queue.Drain(1);

        __syncthreads();

        int tile_idx = temp_storage.tile_idx;
        SizeT block_offset = SizeT(TILE_ITEMS) * tile_idx;

        while (block_offset + TILE_ITEMS <= num_items)
        {
            // Consume full tile
            ConsumeTile<true>(tile_idx, block_offset, partition_ends);
            is_last_tile = (tile_idx == num_tiles - 1);

            // Get next tile
            if (threadIdx.x == 0)
                temp_storage.tile_idx = queue.Drain(1);

            __syncthreads();

            tile_idx = temp_storage.tile_idx;
            block_offset = SizeT(TILE_ITEMS) * tile_idx;
        }

        // Consume a partially-full tile
        if (block_offset < num_items)
        {
            ConsumeTile<false>(tile_idx, block_offset, partition_ends);
            is_last_tile = (tile_idx == num_tiles - 1);
        }
#endif
    }
};


}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)

