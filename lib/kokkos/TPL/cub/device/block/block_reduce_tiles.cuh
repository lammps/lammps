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
 * cub::BlockReduceTiles implements a stateful abstraction of CUDA thread blocks for participating in device-wide reduction.
 */

#pragma once

#include <iterator>

#include "../../block/block_load.cuh"
#include "../../block/block_reduce.cuh"
#include "../../grid/grid_mapping.cuh"
#include "../../grid/grid_queue.cuh"
#include "../../grid/grid_even_share.cuh"
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
 * Tuning policy for BlockReduceTiles
 */
template <
    int                     _BLOCK_THREADS,         ///< Threads per thread block
    int                     _ITEMS_PER_THREAD,      ///< Items per thread per tile of input
    int                     _VECTOR_LOAD_LENGTH,    ///< Number of items per vectorized load
    BlockReduceAlgorithm    _BLOCK_ALGORITHM,       ///< Cooperative block-wide reduction algorithm to use
    PtxLoadModifier         _LOAD_MODIFIER,         ///< PTX load modifier
    GridMappingStrategy     _GRID_MAPPING>          ///< How to map tiles of input onto thread blocks
struct BlockReduceTilesPolicy
{
    enum
    {
        BLOCK_THREADS       = _BLOCK_THREADS,
        ITEMS_PER_THREAD    = _ITEMS_PER_THREAD,
        VECTOR_LOAD_LENGTH  = _VECTOR_LOAD_LENGTH,
    };

    static const BlockReduceAlgorithm  BLOCK_ALGORITHM      = _BLOCK_ALGORITHM;
    static const GridMappingStrategy   GRID_MAPPING         = _GRID_MAPPING;
    static const PtxLoadModifier       LOAD_MODIFIER        = _LOAD_MODIFIER;
};



/******************************************************************************
 * Thread block abstractions
 ******************************************************************************/

/**
 * \brief BlockReduceTiles implements a stateful abstraction of CUDA thread blocks for participating in device-wide reduction.
 *
 * Each thread reduces only the values it loads. If \p FIRST_TILE, this
 * partial reduction is stored into \p thread_aggregate.  Otherwise it is
 * accumulated into \p thread_aggregate.
 */
template <
    typename BlockReduceTilesPolicy,
    typename InputIteratorRA,
    typename SizeT,
    typename ReductionOp>
struct BlockReduceTiles
{

    //---------------------------------------------------------------------
    // Types and constants
    //---------------------------------------------------------------------

    typedef typename std::iterator_traits<InputIteratorRA>::value_type  T;              // Type of input iterator
    typedef VectorHelper<T, BlockReduceTilesPolicy::VECTOR_LOAD_LENGTH> VecHelper;      // Helper type for vectorizing loads of T
    typedef typename VecHelper::Type                                    VectorT;        // Vector of T

    // Constants
    enum
    {
        BLOCK_THREADS       = BlockReduceTilesPolicy::BLOCK_THREADS,
        ITEMS_PER_THREAD    = BlockReduceTilesPolicy::ITEMS_PER_THREAD,
        TILE_ITEMS          = BLOCK_THREADS * ITEMS_PER_THREAD,
        VECTOR_LOAD_LENGTH  = BlockReduceTilesPolicy::VECTOR_LOAD_LENGTH,

        // Can vectorize according to the policy if the input iterator is a native pointer to a built-in primitive
        CAN_VECTORIZE       = (BlockReduceTilesPolicy::VECTOR_LOAD_LENGTH > 1) &&
                                (IsPointer<InputIteratorRA>::VALUE) &&
                                (VecHelper::BUILT_IN),

    };

    static const PtxLoadModifier      LOAD_MODIFIER   = BlockReduceTilesPolicy::LOAD_MODIFIER;
    static const BlockReduceAlgorithm BLOCK_ALGORITHM = BlockReduceTilesPolicy::BLOCK_ALGORITHM;

    // Parameterized BlockReduce primitive
    typedef BlockReduce<T, BLOCK_THREADS, BlockReduceTilesPolicy::BLOCK_ALGORITHM> BlockReduceT;

    /// Shared memory type required by this thread block
    typedef typename BlockReduceT::TempStorage _TempStorage;

    /// Alias wrapper allowing storage to be unioned
    struct TempStorage : Uninitialized<_TempStorage> {};


    //---------------------------------------------------------------------
    // Per-thread fields
    //---------------------------------------------------------------------

    T                       thread_aggregate;   ///< Each thread's partial reduction
    _TempStorage&           temp_storage;       ///< Reference to temp_storage
    InputIteratorRA         d_in;               ///< Input data to reduce
    ReductionOp             reduction_op;       ///< Binary reduction operator
    int                     first_tile_size;    ///< Size of first tile consumed
    bool                    input_aligned;      ///< Whether or not input is vector-aligned


    //---------------------------------------------------------------------
    // Interface
    //---------------------------------------------------------------------

    /**
     * Constructor
     */
    __device__ __forceinline__ BlockReduceTiles(
        TempStorage&            temp_storage,       ///< Reference to temp_storage
        InputIteratorRA         d_in,               ///< Input data to reduce
        ReductionOp             reduction_op)       ///< Binary reduction operator
    :
        temp_storage(temp_storage.Alias()),
        d_in(d_in),
        reduction_op(reduction_op),
        first_tile_size(0),
        input_aligned(CAN_VECTORIZE && ((size_t(d_in) & (sizeof(VectorT) - 1)) == 0))
    {}


    /**
     * Process a single tile of input
     */
    template <bool FULL_TILE>
    __device__ __forceinline__ void ConsumeTile(
        SizeT   block_offset,                   ///< The offset the tile to consume
        int     valid_items = TILE_ITEMS)       ///< The number of valid items in the tile
    {
        if (FULL_TILE)
        {
            T stripe_partial;

            // Load full tile
            if (input_aligned)
            {
                // Alias items as an array of VectorT and load it in striped fashion
                enum { WORDS =  ITEMS_PER_THREAD / VECTOR_LOAD_LENGTH };

                VectorT vec_items[WORDS];

                // Load striped into vec items
                VectorT* alias_ptr = reinterpret_cast<VectorT*>(d_in + block_offset + (threadIdx.x * VECTOR_LOAD_LENGTH));

                #pragma unroll
                for (int i = 0; i < WORDS; ++i)
                    vec_items[i] = alias_ptr[BLOCK_THREADS * i];

                // Reduce items within each thread stripe
                stripe_partial = ThreadReduce<ITEMS_PER_THREAD>(
                    reinterpret_cast<T*>(vec_items),
                    reduction_op);
            }
            else
            {
                T items[ITEMS_PER_THREAD];

                // Load items in striped fashion
                LoadStriped<LOAD_MODIFIER, BLOCK_THREADS>(threadIdx.x, d_in + block_offset, items);

                // Reduce items within each thread stripe
                stripe_partial = ThreadReduce(items, reduction_op);
            }

            // Update running thread aggregate
            thread_aggregate = (first_tile_size) ?
                reduction_op(thread_aggregate, stripe_partial) :       // Update
                stripe_partial;                                        // Assign
        }
        else
        {

            // Partial tile
            int thread_offset = threadIdx.x;

            if (!first_tile_size && (thread_offset < valid_items))
            {
                // Assign thread_aggregate
                thread_aggregate = ThreadLoad<LOAD_MODIFIER>(d_in + block_offset + thread_offset);
                thread_offset += BLOCK_THREADS;
            }

            while (thread_offset < valid_items)
            {
                // Update thread aggregate
                T item = ThreadLoad<LOAD_MODIFIER>(d_in + block_offset + thread_offset);
                thread_aggregate = reduction_op(thread_aggregate, item);
                thread_offset += BLOCK_THREADS;
            }
        }

        // Set first tile size if necessary
        if (!first_tile_size)
            first_tile_size = valid_items;
    }


    //---------------------------------------------------------------------
    // Consume a contiguous segment of tiles
    //---------------------------------------------------------------------

    /**
     * \brief Reduce a contiguous segment of input tiles
     */
    __device__ __forceinline__ void ConsumeTiles(
        SizeT   block_offset,                       ///< [in] Threadblock begin offset (inclusive)
        SizeT   block_oob,                          ///< [in] Threadblock end offset (exclusive)
        T       &block_aggregate)                   ///< [out] Running total
    {
        // Consume subsequent full tiles of input
        while (block_offset + TILE_ITEMS <= block_oob)
        {
            ConsumeTile<true>(block_offset);
            block_offset += TILE_ITEMS;
        }

        // Consume a partially-full tile
        if (block_offset < block_oob)
        {
            int valid_items = block_oob - block_offset;
            ConsumeTile<false>(block_offset, valid_items);
        }

        // Compute block-wide reduction
        block_aggregate = (first_tile_size < TILE_ITEMS) ?
            BlockReduceT(temp_storage).Reduce(thread_aggregate, reduction_op, first_tile_size) :
            BlockReduceT(temp_storage).Reduce(thread_aggregate, reduction_op);
    }


    /**
     * Reduce a contiguous segment of input tiles
     */
    __device__ __forceinline__ void ConsumeTiles(
        SizeT                               num_items,          ///< [in] Total number of global input items
        GridEvenShare<SizeT>                &even_share,        ///< [in] GridEvenShare descriptor
        GridQueue<SizeT>                    &queue,             ///< [in,out] GridQueue descriptor
        T                                   &block_aggregate,   ///< [out] Running total
        Int2Type<GRID_MAPPING_EVEN_SHARE>   is_even_share)      ///< [in] Marker type indicating this is an even-share mapping
    {
        // Initialize even-share descriptor for this thread block
        even_share.BlockInit();

        // Consume input tiles
        ConsumeTiles(even_share.block_offset, even_share.block_oob, block_aggregate);
    }


    //---------------------------------------------------------------------
    // Dynamically consume tiles
    //---------------------------------------------------------------------

    /**
     * Dequeue and reduce tiles of items as part of a inter-block scan
     */
    __device__ __forceinline__ void ConsumeTiles(
        int                 num_items,          ///< Total number of input items
        GridQueue<SizeT>    queue,              ///< Queue descriptor for assigning tiles of work to thread blocks
        T                   &block_aggregate)   ///< [out] Running total
    {
        // Shared dequeue offset
        __shared__ SizeT dequeue_offset;

        // We give each thread block at least one tile of input.
        SizeT block_offset = blockIdx.x * TILE_ITEMS;
        SizeT even_share_base = gridDim.x * TILE_ITEMS;

        if (block_offset + TILE_ITEMS <= num_items)
        {
            // Consume full tile of input
            ConsumeTile<true>(block_offset);

            // Dequeue more tiles
            while (true)
            {
                 // Dequeue a tile of items
                if (threadIdx.x == 0)
                    dequeue_offset = queue.Drain(TILE_ITEMS) + even_share_base;

                __syncthreads();

                // Grab tile offset and check if we're done with full tiles
                block_offset = dequeue_offset;

                __syncthreads();

                if (block_offset + TILE_ITEMS > num_items)
                    break;

                // Consume a full tile
                ConsumeTile<true>(block_offset);
            }
        }

        if (block_offset < num_items)
        {
            int valid_items = num_items - block_offset;
            ConsumeTile<false>(block_offset, valid_items);
        }

        // Compute block-wide reduction
        block_aggregate = (first_tile_size < TILE_ITEMS) ?
            BlockReduceT(temp_storage).Reduce(thread_aggregate, reduction_op, first_tile_size) :
            BlockReduceT(temp_storage).Reduce(thread_aggregate, reduction_op);
    }


    /**
     * Dequeue and reduce tiles of items as part of a inter-block scan
     */
    __device__ __forceinline__ void ConsumeTiles(
        SizeT                               num_items,          ///< [in] Total number of global input items
        GridEvenShare<SizeT>                &even_share,        ///< [in] GridEvenShare descriptor
        GridQueue<SizeT>                    &queue,             ///< [in,out] GridQueue descriptor
        T                                   &block_aggregate,   ///< [out] Running total
        Int2Type<GRID_MAPPING_DYNAMIC>      is_dynamic)         ///< [in] Marker type indicating this is a dynamic mapping
    {
        ConsumeTiles(num_items, queue, block_aggregate);
    }

};


}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)

