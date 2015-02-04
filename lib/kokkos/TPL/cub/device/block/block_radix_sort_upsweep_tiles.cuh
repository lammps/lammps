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
 * BlockRadixSortUpsweepTiles implements a stateful abstraction of CUDA thread blocks for participating in device-wide radix sort upsweep.
 */

#pragma once

#include "../../thread/thread_reduce.cuh"
#include "../../thread/thread_load.cuh"
#include "../../block/block_load.cuh"
#include "../../util_type.cuh"
#include "../../util_namespace.cuh"

/// Optional outer namespace(s)
CUB_NS_PREFIX

/// CUB namespace
namespace cub {

/******************************************************************************
 * Tuning policy types
 ******************************************************************************/

/**
 * Tuning policy for BlockRadixSortUpsweepTiles
 */
template <
    int                 _BLOCK_THREADS,     ///< The number of threads per CTA
    int                 _ITEMS_PER_THREAD,  ///< The number of items to load per thread per tile
    PtxLoadModifier     _LOAD_MODIFIER,     ///< Load cache-modifier
    int                 _RADIX_BITS>        ///< The number of radix bits, i.e., log2(bins)
struct BlockRadixSortUpsweepTilesPolicy
{
    enum
    {
        BLOCK_THREADS       = _BLOCK_THREADS,
        ITEMS_PER_THREAD    = _ITEMS_PER_THREAD,
        RADIX_BITS          = _RADIX_BITS,
        TILE_ITEMS          = BLOCK_THREADS * ITEMS_PER_THREAD,
    };

    static const PtxLoadModifier LOAD_MODIFIER = _LOAD_MODIFIER;

    typedef BlockRadixSortUpsweepTilesPolicy<
        BLOCK_THREADS,
        ITEMS_PER_THREAD,
        LOAD_MODIFIER,
        CUB_MAX(1, RADIX_BITS - 1)> AltPolicy;
};


/******************************************************************************
 * Thread block abstractions
 ******************************************************************************/

/**
 * \brief BlockRadixSortUpsweepTiles implements a stateful abstraction of CUDA thread blocks for participating in device-wide radix sort upsweep.
 *
 * Computes radix digit histograms over a range of input tiles.
 */
template <
    typename BlockRadixSortUpsweepTilesPolicy,
    typename Key,
    typename SizeT>
struct BlockRadixSortUpsweepTiles
{

    //---------------------------------------------------------------------
    // Type definitions and constants
    //---------------------------------------------------------------------

    typedef typename Traits<Key>::UnsignedBits UnsignedBits;

    // Integer type for digit counters (to be packed into words of PackedCounters)
    typedef unsigned char DigitCounter;

    // Integer type for packing DigitCounters into columns of shared memory banks
    typedef unsigned int PackedCounter;

    static const PtxLoadModifier LOAD_MODIFIER = BlockRadixSortUpsweepTilesPolicy::LOAD_MODIFIER;

    enum
    {
        RADIX_BITS              = BlockRadixSortUpsweepTilesPolicy::RADIX_BITS,
        BLOCK_THREADS           = BlockRadixSortUpsweepTilesPolicy::BLOCK_THREADS,
        KEYS_PER_THREAD         = BlockRadixSortUpsweepTilesPolicy::ITEMS_PER_THREAD,

        RADIX_DIGITS            = 1 << RADIX_BITS,

        LOG_WARP_THREADS        = PtxArchProps::LOG_WARP_THREADS,
        WARP_THREADS            = 1 << LOG_WARP_THREADS,
        WARPS                   = (BLOCK_THREADS + WARP_THREADS - 1) / WARP_THREADS,

        TILE_ITEMS              = BLOCK_THREADS * KEYS_PER_THREAD,

        BYTES_PER_COUNTER       = sizeof(DigitCounter),
        LOG_BYTES_PER_COUNTER   = Log2<BYTES_PER_COUNTER>::VALUE,

        PACKING_RATIO           = sizeof(PackedCounter) / sizeof(DigitCounter),
        LOG_PACKING_RATIO       = Log2<PACKING_RATIO>::VALUE,

        LOG_COUNTER_LANES       = CUB_MAX(0, RADIX_BITS - LOG_PACKING_RATIO),
        COUNTER_LANES           = 1 << LOG_COUNTER_LANES,

        // To prevent counter overflow, we must periodically unpack and aggregate the
        // digit counters back into registers.  Each counter lane is assigned to a
        // warp for aggregation.

        LANES_PER_WARP          = CUB_MAX(1, (COUNTER_LANES + WARPS - 1) / WARPS),

        // Unroll tiles in batches without risk of counter overflow
        UNROLL_COUNT            = CUB_MIN(64, 255 / KEYS_PER_THREAD),
        UNROLLED_ELEMENTS       = UNROLL_COUNT * TILE_ITEMS,
    };



    /**
     * Shared memory storage layout
     */
    struct _TempStorage
    {
        union
        {
            DigitCounter    digit_counters[COUNTER_LANES][BLOCK_THREADS][PACKING_RATIO];
            PackedCounter   packed_counters[COUNTER_LANES][BLOCK_THREADS];
            SizeT           digit_partials[RADIX_DIGITS][WARP_THREADS + 1];
        };
    };


    /// Alias wrapper allowing storage to be unioned
    struct TempStorage : Uninitialized<_TempStorage> {};


    //---------------------------------------------------------------------
    // Thread fields (aggregate state bundle)
    //---------------------------------------------------------------------

    // Shared storage for this CTA
    _TempStorage    &temp_storage;

    // Thread-local counters for periodically aggregating composite-counter lanes
    SizeT           local_counts[LANES_PER_WARP][PACKING_RATIO];

    // Input and output device pointers
    UnsignedBits    *d_keys_in;

    // The least-significant bit position of the current digit to extract
    int             current_bit;



    //---------------------------------------------------------------------
    // Helper structure for templated iteration
    //---------------------------------------------------------------------

    // Iterate
    template <int COUNT, int MAX>
    struct Iterate
    {
        enum {
            HALF = (MAX / 2),
        };

        // BucketKeys
        static __device__ __forceinline__ void BucketKeys(
            BlockRadixSortUpsweepTiles &cta,
            UnsignedBits keys[KEYS_PER_THREAD])
        {
            cta.Bucket(keys[COUNT]);

            // Next
            Iterate<COUNT + 1, MAX>::BucketKeys(cta, keys);
        }

        // ProcessTiles
        static __device__ __forceinline__ void ProcessTiles(BlockRadixSortUpsweepTiles &cta, SizeT block_offset)
        {
            // Next
            Iterate<1, HALF>::ProcessTiles(cta, block_offset);
            Iterate<1, MAX - HALF>::ProcessTiles(cta, block_offset + (HALF * TILE_ITEMS));
        }
    };

    // Terminate
    template <int MAX>
    struct Iterate<MAX, MAX>
    {
        // BucketKeys
        static __device__ __forceinline__ void BucketKeys(BlockRadixSortUpsweepTiles &cta, UnsignedBits keys[KEYS_PER_THREAD]) {}

        // ProcessTiles
        static __device__ __forceinline__ void ProcessTiles(BlockRadixSortUpsweepTiles &cta, SizeT block_offset)
        {
            cta.ProcessFullTile(block_offset);
        }
    };


    //---------------------------------------------------------------------
    // Utility methods
    //---------------------------------------------------------------------

    /**
     * Decode a key and increment corresponding smem digit counter
     */
    __device__ __forceinline__ void Bucket(UnsignedBits key)
    {
        // Perform transform op
        UnsignedBits converted_key = Traits<Key>::TwiddleIn(key);

        // Add in sub-counter offset
        UnsignedBits sub_counter = BFE(converted_key, current_bit, LOG_PACKING_RATIO);

        // Add in row offset
        UnsignedBits row_offset = BFE(converted_key, current_bit + LOG_PACKING_RATIO, LOG_COUNTER_LANES);

        // Increment counter
        temp_storage.digit_counters[row_offset][threadIdx.x][sub_counter]++;

    }


    /**
     * Reset composite counters
     */
    __device__ __forceinline__ void ResetDigitCounters()
    {
        #pragma unroll
        for (int LANE = 0; LANE < COUNTER_LANES; LANE++)
        {
            temp_storage.packed_counters[LANE][threadIdx.x] = 0;
        }
    }


    /**
     * Reset the unpacked counters in each thread
     */
    __device__ __forceinline__ void ResetUnpackedCounters()
    {
        #pragma unroll
        for (int LANE = 0; LANE < LANES_PER_WARP; LANE++)
        {
            #pragma unroll
            for (int UNPACKED_COUNTER = 0; UNPACKED_COUNTER < PACKING_RATIO; UNPACKED_COUNTER++)
            {
                local_counts[LANE][UNPACKED_COUNTER] = 0;
            }
        }
    }


    /**
     * Extracts and aggregates the digit counters for each counter lane
     * owned by this warp
     */
    __device__ __forceinline__ void UnpackDigitCounts()
    {
        unsigned int warp_id = threadIdx.x >> LOG_WARP_THREADS;
        unsigned int warp_tid = threadIdx.x & (WARP_THREADS - 1);

        #pragma unroll
        for (int LANE = 0; LANE < LANES_PER_WARP; LANE++)
        {
            const int counter_lane = (LANE * WARPS) + warp_id;
            if (counter_lane < COUNTER_LANES)
            {
                #pragma unroll
                for (int PACKED_COUNTER = 0; PACKED_COUNTER < BLOCK_THREADS; PACKED_COUNTER += WARP_THREADS)
                {
                    #pragma unroll
                    for (int UNPACKED_COUNTER = 0; UNPACKED_COUNTER < PACKING_RATIO; UNPACKED_COUNTER++)
                    {
                        SizeT counter = temp_storage.digit_counters[counter_lane][warp_tid + PACKED_COUNTER][UNPACKED_COUNTER];
                        local_counts[LANE][UNPACKED_COUNTER] += counter;
                    }
                }
            }
        }
    }


    /**
     * Places unpacked counters into smem for final digit reduction
     */
    __device__ __forceinline__ void ReduceUnpackedCounts(SizeT &bin_count)
    {
        unsigned int warp_id = threadIdx.x >> LOG_WARP_THREADS;
        unsigned int warp_tid = threadIdx.x & (WARP_THREADS - 1);

        // Place unpacked digit counters in shared memory
        #pragma unroll
        for (int LANE = 0; LANE < LANES_PER_WARP; LANE++)
        {
            int counter_lane = (LANE * WARPS) + warp_id;
            if (counter_lane < COUNTER_LANES)
            {
                int digit_row = counter_lane << LOG_PACKING_RATIO;

                #pragma unroll
                for (int UNPACKED_COUNTER = 0; UNPACKED_COUNTER < PACKING_RATIO; UNPACKED_COUNTER++)
                {
                    temp_storage.digit_partials[digit_row + UNPACKED_COUNTER][warp_tid] =
                        local_counts[LANE][UNPACKED_COUNTER];
                }
            }
        }

        __syncthreads();

        // Rake-reduce bin_count reductions
        if (threadIdx.x < RADIX_DIGITS)
        {
            bin_count = ThreadReduce<WARP_THREADS>(
                temp_storage.digit_partials[threadIdx.x],
                Sum());
        }
    }


    /**
     * Processes a single, full tile
     */
    __device__ __forceinline__ void ProcessFullTile(SizeT block_offset)
    {
        // Tile of keys
        UnsignedBits keys[KEYS_PER_THREAD];

        LoadStriped<LOAD_MODIFIER, BLOCK_THREADS>(threadIdx.x, d_keys_in + block_offset, keys);

        // Prevent hoisting
//        __threadfence_block();
//        __syncthreads();

        // Bucket tile of keys
        Iterate<0, KEYS_PER_THREAD>::BucketKeys(*this, keys);
    }


    /**
     * Processes a single load (may have some threads masked off)
     */
    __device__ __forceinline__ void ProcessPartialTile(
        SizeT block_offset,
        const SizeT &block_oob)
    {
        // Process partial tile if necessary using single loads
        block_offset += threadIdx.x;
        while (block_offset < block_oob)
        {
            // Load and bucket key
            UnsignedBits key = ThreadLoad<LOAD_MODIFIER>(d_keys_in + block_offset);
            Bucket(key);
            block_offset += BLOCK_THREADS;
        }
    }


    //---------------------------------------------------------------------
    // Interface
    //---------------------------------------------------------------------

    /**
     * Constructor
     */
    __device__ __forceinline__ BlockRadixSortUpsweepTiles(
        TempStorage &temp_storage,
        Key         *d_keys_in,
        int         current_bit)
    :
        temp_storage(temp_storage.Alias()),
        d_keys_in(reinterpret_cast<UnsignedBits*>(d_keys_in)),
        current_bit(current_bit)
    {}


    /**
     * Compute radix digit histograms from a segment of input tiles.
     */
    __device__ __forceinline__ void ProcessTiles(
        SizeT           block_offset,
        const SizeT     &block_oob,
        SizeT           &bin_count)                ///< [out] The digit count for tid'th bin (output param, valid in the first RADIX_DIGITS threads)
    {
        // Reset digit counters in smem and unpacked counters in registers
        ResetDigitCounters();
        ResetUnpackedCounters();

        // Unroll batches of full tiles
        while (block_offset + UNROLLED_ELEMENTS <= block_oob)
        {
            Iterate<0, UNROLL_COUNT>::ProcessTiles(*this, block_offset);
            block_offset += UNROLLED_ELEMENTS;

            __syncthreads();

            // Aggregate back into local_count registers to prevent overflow
            UnpackDigitCounts();

            __syncthreads();

            // Reset composite counters in lanes
            ResetDigitCounters();
        }

        // Unroll single full tiles
        while (block_offset + TILE_ITEMS <= block_oob)
        {
            ProcessFullTile(block_offset);
            block_offset += TILE_ITEMS;
        }

        // Process partial tile if necessary
        ProcessPartialTile(
            block_offset,
            block_oob);

        __syncthreads();

        // Aggregate back into local_count registers
        UnpackDigitCounts();

        __syncthreads();

        // Final raking reduction of counts by bin
        ReduceUnpackedCounts(bin_count);
    }

};


}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)

