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
 * cub::BlockReduceByKeyiles implements a stateful abstraction of CUDA thread blocks for participating in device-wide reduce-value-by-key.
 */

#pragma once

#include <iterator>

#include "scan_tiles_types.cuh"
#include "../../block/block_load.cuh"
#include "../../block/block_discontinuity.cuh"
#include "../../block/block_scan.cuh"
#include "../../util_namespace.cuh"

/// Optional outer namespace(s)
CUB_NS_PREFIX

/// CUB namespace
namespace cub {


/******************************************************************************
 * Utility data types
 ******************************************************************************/

/// Scan tuple data type for reduce-value-by-key
template <typename Value, typename SizeT>
struct ReduceByKeyuple
{
    Value   value;      // Initially set as value, contains segment aggregate after prefix scan
    SizeT   flag;       // Initially set as a tail flag, contains scatter offset after prefix scan
};


/// Binary reduce-by-key scan operator
template <typename ReductionOp>
struct ReduceByKeyScanOp
{
    /// Reduction functor
    ReductionOp reduction_op;

    /// Constructor
    ReduceByKeyScanOp(ReductionOp reduction_op) : reduction_op(reduction_op)
    {}

    /// Binary scan operator
    template <typename ReduceByKeyuple>
    __device__ __forceinline__ ReduceByKeyuple operator()(
        const ReduceByKeyuple &first,
        const ReduceByKeyuple &second)
    {
        ReduceByKeyuple retval;
        retval.val = (second.flag) ? second.val : reduction_op(first.val, second.val);
        retval.flag = first.flag + second.flag;
        return retval;
    }
};



/******************************************************************************
 * Tuning policy types
 ******************************************************************************/

/**
 * Tuning policy for BlockReduceByKeyiles
 */
template <
    int                         _BLOCK_THREADS,
    int                         _ITEMS_PER_THREAD,
    BlockLoadAlgorithm          _LOAD_ALGORITHM,
    bool                        _LOAD_WARP_TIME_SLICING,
    PtxLoadModifier             _LOAD_MODIFIER,
    BlockScanAlgorithm          _SCAN_ALGORITHM>
struct BlockReduceByKeyilesPolicy
{
    enum
    {
        BLOCK_THREADS           = _BLOCK_THREADS,
        ITEMS_PER_THREAD        = _ITEMS_PER_THREAD,
        LOAD_WARP_TIME_SLICING  = _LOAD_WARP_TIME_SLICING,
    };

    static const BlockLoadAlgorithm     LOAD_ALGORITHM      = _LOAD_ALGORITHM;
    static const PtxLoadModifier        LOAD_MODIFIER       = _LOAD_MODIFIER;
    static const BlockScanAlgorithm     SCAN_ALGORITHM      = _SCAN_ALGORITHM;
};


/******************************************************************************
 * Thread block abstractions
 ******************************************************************************/

/**
 * \brief BlockReduceByKeyiles implements a stateful abstraction of CUDA thread blocks for participating in device-wide prefix scan.
 */
template <
    typename BlockReduceByKeyilesPolicy,   ///< Tuning policy
    typename KeyInputIteratorRA,            ///< Random-access input iterator type for keys
    typename KeyOutputIteratorRA,           ///< Random-access output iterator type for keys
    typename ValueInputIteratorRA,          ///< Random-access input iterator type for values
    typename ValueOutputIteratorRA,         ///< Random-access output iterator type for values
    typename ReductionOp,                   ///< Reduction functor type
    typename SizeT>                         ///< Offset integer type
struct BlockReduceByKeyiles
{
    //---------------------------------------------------------------------
    // Types and constants
    //---------------------------------------------------------------------

    // Data types of input iterators
    typedef typename std::iterator_traits<KeyInputIteratorRA>::value_type   Key;    // Key data type
    typedef typename std::iterator_traits<ValueInputIteratorRA>::value_type Value;  // Value data type

    // Constants
    enum
    {
        BLOCK_THREADS       = BlockReduceByKeyilesPolicy::BLOCK_THREADS,
        ITEMS_PER_THREAD    = BlockReduceByKeyilesPolicy::ITEMS_PER_THREAD,
        TILE_ITEMS          = BLOCK_THREADS * ITEMS_PER_THREAD,
        STATUS_PADDING      = PtxArchProps::WARP_THREADS,
    };

    // Block load type for keys
    typedef BlockLoad<
        KeyInputIteratorRA,
        BlockReduceByKeyilesPolicy::BLOCK_THREADS,
        BlockReduceByKeyilesPolicy::ITEMS_PER_THREAD,
        BlockReduceByKeyilesPolicy::LOAD_ALGORITHM,
        BlockReduceByKeyilesPolicy::LOAD_MODIFIER,
        BlockReduceByKeyilesPolicy::LOAD_WARP_TIME_SLICING>    BlockLoadKeys;

    // Block load type for values
    typedef BlockLoad<
        ValueInputIteratorRA,
        BlockReduceByKeyilesPolicy::BLOCK_THREADS,
        BlockReduceByKeyilesPolicy::ITEMS_PER_THREAD,
        BlockReduceByKeyilesPolicy::LOAD_ALGORITHM,
        BlockReduceByKeyilesPolicy::LOAD_MODIFIER,
        BlockReduceByKeyilesPolicy::LOAD_WARP_TIME_SLICING>    BlockLoadValues;

    // Block discontinuity type for setting tail flags
    typedef BlockDiscontinuity<Key, BLOCK_THREADS>              BlockDiscontinuityKeys;

    // Scan tuple type
    typedef ReduceByKeyuple<Value, SizeT>                      ScanTuple;

    // Tile status descriptor type
    typedef ScanTileDescriptor<ScanTuple>                 ScanTileDescriptorT;

    // Block scan functor type
    typedef ReduceByKeyScanOp<ReductionOp>                      ScanOp;

    // Block scan prefix callback type
    typedef DeviceScanBlockPrefixOp<ScanTuple, ScanOp>          PrefixCallback;

    // Block scan type
    typedef BlockScan<
        ScanTuple,
        BlockReduceByKeyilesPolicy::BLOCK_THREADS,
        BlockReduceByKeyilesPolicy::SCAN_ALGORITHM>            BlockScanT;

    /// Shared memory type for this threadblock
    struct _TempStorage
    {
        union
        {
            typename BlockLoadKeys::TempStorage         load_keys;      // Smem needed for loading tiles of keys
            typename BlockLoadValues::TempStorage       load_values;    // Smem needed for loading tiles of values
            struct
            {
                typename BlockScanT::TempStorage        scan;           // Smem needed for tile scanning
                typename PrefixCallback::TempStorage    prefix;         // Smem needed for cooperative prefix callback
            };
        };

        typename BlockDiscontinuityKeys::TempStorage    flagging;       // Smem needed for tile scanning
        SizeT                                           tile_idx;       // Shared tile index
    };

    /// Alias wrapper allowing storage to be unioned
    struct TempStorage : Uninitialized<_TempStorage> {};



    //---------------------------------------------------------------------
    // Per-thread fields
    //---------------------------------------------------------------------

    _TempStorage                &temp_storage;      ///< Reference to temp_storage
    KeyInputIteratorRA          d_keys_in;          ///< Key input data
    KeyOutputIteratorRA         d_keys_out;         ///< Key output data
    ValueInputIteratorRA        d_values_in;        ///< Value input data
    ValueOutputIteratorRA       d_values_out;       ///< Value output data
    ScanTileDescriptorT         *d_tile_status;     ///< Global list of tile status
    ScanOp                      scan_op;            ///< Binary scan operator
    int                         num_tiles;          ///< Total number of input tiles for the entire problem
    SizeT                       num_items;          ///< Total number of scan items for the entire problem


    //---------------------------------------------------------------------
    // Interface
    //---------------------------------------------------------------------

    // Constructor
    __device__ __forceinline__
    BlockReduceByKeyiles(
        TempStorage                 &temp_storage,      ///< Reference to temp_storage
        KeyInputIteratorRA          d_keys_in,          ///< Key input data
        KeyOutputIteratorRA         d_keys_out,         ///< Key output data
        ValueInputIteratorRA        d_values_in,        ///< Value input data
        ValueOutputIteratorRA       d_values_out,       ///< Value output data
        ScanTileDescriptorT       *d_tile_status,     ///< Global list of tile status
        ReductionOp                 reduction_op,       ///< Binary scan operator
        int                         num_tiles,          ///< Total number of input tiles for the entire problem
        SizeT                       num_items)          ///< Total number of scan items for the entire problem
    :
        temp_storage(temp_storage.Alias()),
        d_keys_in(d_keys_in),
        d_keys_out(d_keys_out),
        d_values_in(d_values_in),
        d_values_out(d_values_out),
        d_tile_status(d_tile_status),
        scan_op(reduction_op),
        num_tiles(num_tiles),
        num_items(num_items)
    {}


    /**
     * Process a tile of input
     */
    template <bool FULL_TILE>
    __device__ __forceinline__ void ConsumeTile(
        int     tile_idx,                   ///< Tile index
        SizeT   block_offset,               ///< Tile offset
        int     valid_items = TILE_ITEMS)   ///< Number of valid items in the tile
    {
        Key         keys[ITEMS_PER_THREAD];
        Value       values[ITEMS_PER_THREAD];
        int         tail_flags[ITEMS_PER_THREAD];
        ScanTuple   scan_tuples[ITEMS_PER_THREAD];

        // Load keys
        if (FULL_TILE)
            BlockLoadKeys(temp_storage.load_keys).Load(d_keys_in + block_offset, keys);
        else
            BlockLoadKeys(temp_storage.load_keys).Load(d_keys_in + block_offset, keys, valid_items);

        // Set tail flags
        if (tile_idx == num_tiles - 1)
        {
            // Last tile
            BlockDiscontinuityKeys(temp_storage.flagging).FlagTails(tail_flags, keys, Equality());
        }
        else
        {
            // Preceding tiles require the first element of the next tile
            Key tile_suffix_item;
            if (threadIdx.x == 0)
                tile_suffix_item = d_keys_in[block_offset + TILE_ITEMS];

            BlockDiscontinuityKeys(temp_storage.flagging).FlagTails(tail_flags, keys, Equality(), tile_suffix_item);
        }

        __syncthreads();

        // Load values
        if (FULL_TILE)
            BlockLoadValues(temp_storage.load_values).Load(d_values_in + block_offset, values);
        else
            BlockLoadValues(temp_storage.load_values).Load(d_values_in + block_offset, values, valid_items);

        // Assemble scan tuples
        #pragma unroll
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ++ITEM)
        {
            scan_tuples[ITEM].value     = values[ITEM];
            scan_tuples[ITEM].flag      = tail_flags[ITEM];
        }

        __syncthreads();

        // Perform inclusive prefix scan
        ScanTuple block_aggregate;
        if (tile_idx == 0)
        {
            // Without prefix callback
            BlockScanT(temp_storage.scan).InclusiveScan(scan_tuples, scan_tuples, scan_op, block_aggregate);

            // Update tile status
            if (threadIdx.x == 0)
                ScanTileDescriptorT::SetPrefix(d_tile_status, block_aggregate);
        }
        else
        {
            // With prefix callback
            PrefixCallback prefix_op(d_tile_status, temp_storage.prefix, scan_op, tile_idx);
            BlockScanT(temp_storage.scan).InclusiveScan(scan_tuples, scan_tuples, scan_op, block_aggregate, prefix_op);
        }

        // Scatter flagged keys and values to output
        #pragma unroll
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ++ITEM)
        {
            int tile_item = (threadIdx.x * ITEMS_PER_THREAD) + ITEM;

            // Set the head flag on the last item in a partially-full tile
            if (!FULL_TILE && (tile_item == valid_items - 1))
                tail_flags[ITEM] = 1;

            // Decrement scatter offset
            scan_tuples[ITEM].flag--;

            // Scatter key and aggregate value if flagged and in range
            if ((FULL_TILE || (tile_item < valid_items)) && (tail_flags[ITEM]))
            {
                d_keys_out[scan_tuples[ITEM].flag]      = keys[ITEM];
                d_values_out[scan_tuples[ITEM].flag]    = scan_tuples[ITEM].value;
            }
        }
    }



    /**
     * Dequeue and scan tiles of elements
     */
    __device__ __forceinline__ void ProcessTiles(GridQueue<int> queue)          ///< Queue descriptor for assigning tiles of work to thread blocks
    {
        // We give each thread block at least one tile of input
        int tile_idx = blockIdx.x;

        // Consume full tiles of input
        SizeT block_offset = SizeT(TILE_ITEMS) * tile_idx;
        while (block_offset + TILE_ITEMS <= num_items)
        {
            ConsumeTile<true>(tile_idx, block_offset);

            // Get next tile
#if CUB_PTX_ARCH < 200
            // No concurrent kernels allowed, so just stripe tiles
            tile_idx += gridDim.x;
#else
            // Concurrent kernels are allowed, so we must only use active blocks to dequeue tile indices
            if (threadIdx.x == 0)
                temp_storage.tile_idx = queue.Drain(1) + gridDim.x;

            __syncthreads();

            tile_idx = temp_storage.tile_idx;
#endif
            block_offset = SizeT(TILE_ITEMS) * tile_idx;
        }

        // Consume a partially-full tile
        if (block_offset < num_items)
        {
            // Consume a partially-full tile
            int valid_items = num_items - block_offset;
            ConsumeTile<false>(tile_idx, block_offset, valid_items);
        }
    }

};


}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)

