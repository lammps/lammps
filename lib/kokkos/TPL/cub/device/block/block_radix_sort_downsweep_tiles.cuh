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
 * BlockRadixSortDownsweepTiles implements a stateful abstraction of CUDA thread blocks for participating in device-wide radix sort downsweep.
 */


#pragma once

#include "../../thread/thread_load.cuh"
#include "../../block/block_load.cuh"
#include "../../block/block_store.cuh"
#include "../../block/block_radix_rank.cuh"
#include "../../block/block_exchange.cuh"
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
 * Types of scattering strategies
 */
enum RadixSortScatterAlgorithm
{
    RADIX_SORT_SCATTER_DIRECT,      ///< Scatter directly from registers to global bins
    RADIX_SORT_SCATTER_TWO_PHASE,   ///< First scatter from registers into shared memory bins, then into global bins
};


/**
 * Tuning policy for BlockRadixSortDownsweepTiles
 */
template <
    int                         _BLOCK_THREADS,             ///< The number of threads per CTA
    int                         _ITEMS_PER_THREAD,          ///< The number of consecutive downsweep keys to process per thread
    BlockLoadAlgorithm          _LOAD_ALGORITHM,            ///< The BlockLoad algorithm to use
    PtxLoadModifier             _LOAD_MODIFIER,             ///< The PTX cache-modifier to use for loads
    bool                        _EXCHANGE_TIME_SLICING,     ///< Whether or not to time-slice key/value exchanges through shared memory to lower shared memory pressure
    bool                        _MEMOIZE_OUTER_SCAN,        ///< Whether or not to buffer outer raking scan partials to incur fewer shared memory reads at the expense of higher register pressure.  See BlockScanAlgorithm::BLOCK_SCAN_RAKING_MEMOIZE for more details.
    BlockScanAlgorithm          _INNER_SCAN_ALGORITHM,      ///< The cub::BlockScanAlgorithm algorithm to use
    RadixSortScatterAlgorithm   _SCATTER_ALGORITHM,         ///< The scattering strategy to use
    cudaSharedMemConfig         _SMEM_CONFIG,               ///< Shared memory bank mode (default: \p cudaSharedMemBankSizeFourByte)
    int                         _RADIX_BITS>                ///< The number of radix bits, i.e., log2(bins)
struct BlockRadixSortDownsweepTilesPolicy
{
    enum
    {
        BLOCK_THREADS           = _BLOCK_THREADS,
        ITEMS_PER_THREAD        = _ITEMS_PER_THREAD,
        EXCHANGE_TIME_SLICING   = _EXCHANGE_TIME_SLICING,
        RADIX_BITS              = _RADIX_BITS,
        MEMOIZE_OUTER_SCAN      = _MEMOIZE_OUTER_SCAN,
        TILE_ITEMS              = BLOCK_THREADS * ITEMS_PER_THREAD,
    };

    static const BlockLoadAlgorithm         LOAD_ALGORITHM          = _LOAD_ALGORITHM;
    static const PtxLoadModifier            LOAD_MODIFIER           = _LOAD_MODIFIER;
    static const BlockScanAlgorithm         INNER_SCAN_ALGORITHM    = _INNER_SCAN_ALGORITHM;
    static const RadixSortScatterAlgorithm  SCATTER_ALGORITHM       = _SCATTER_ALGORITHM;
    static const cudaSharedMemConfig        SMEM_CONFIG             = _SMEM_CONFIG;

    typedef BlockRadixSortDownsweepTilesPolicy<
        BLOCK_THREADS,
        ITEMS_PER_THREAD,
        LOAD_ALGORITHM,
        LOAD_MODIFIER,
        EXCHANGE_TIME_SLICING,
        MEMOIZE_OUTER_SCAN,
        INNER_SCAN_ALGORITHM,
        SCATTER_ALGORITHM,
        SMEM_CONFIG,
        CUB_MAX(1, RADIX_BITS - 1)> AltPolicy;
};


/******************************************************************************
 * Thread block abstractions
 ******************************************************************************/

/**
 * CTA-wide "downsweep" abstraction for distributing keys from
 * a range of input tiles.
 */
template <
    typename BlockRadixSortDownsweepTilesPolicy,
    typename Key,
    typename Value,
    typename SizeT>
struct BlockRadixSortDownsweepTiles
{
    //---------------------------------------------------------------------
    // Type definitions and constants
    //---------------------------------------------------------------------

    // Appropriate unsigned-bits representation of Key
    typedef typename Traits<Key>::UnsignedBits UnsignedBits;

    static const UnsignedBits MIN_KEY = Traits<Key>::MIN_KEY;
    static const UnsignedBits MAX_KEY = Traits<Key>::MAX_KEY;

    static const BlockLoadAlgorithm         LOAD_ALGORITHM          = BlockRadixSortDownsweepTilesPolicy::LOAD_ALGORITHM;
    static const PtxLoadModifier            LOAD_MODIFIER           = BlockRadixSortDownsweepTilesPolicy::LOAD_MODIFIER;
    static const BlockScanAlgorithm         INNER_SCAN_ALGORITHM    = BlockRadixSortDownsweepTilesPolicy::INNER_SCAN_ALGORITHM;
    static const RadixSortScatterAlgorithm  SCATTER_ALGORITHM       = BlockRadixSortDownsweepTilesPolicy::SCATTER_ALGORITHM;
    static const cudaSharedMemConfig        SMEM_CONFIG             = BlockRadixSortDownsweepTilesPolicy::SMEM_CONFIG;

    enum
    {
        BLOCK_THREADS           = BlockRadixSortDownsweepTilesPolicy::BLOCK_THREADS,
        ITEMS_PER_THREAD        = BlockRadixSortDownsweepTilesPolicy::ITEMS_PER_THREAD,
        EXCHANGE_TIME_SLICING   = BlockRadixSortDownsweepTilesPolicy::EXCHANGE_TIME_SLICING,
        RADIX_BITS              = BlockRadixSortDownsweepTilesPolicy::RADIX_BITS,
        MEMOIZE_OUTER_SCAN      = BlockRadixSortDownsweepTilesPolicy::MEMOIZE_OUTER_SCAN,
        TILE_ITEMS              = BLOCK_THREADS * ITEMS_PER_THREAD,

        RADIX_DIGITS            = 1 << RADIX_BITS,
        KEYS_ONLY               = Equals<Value, NullType>::VALUE,

        WARP_THREADS            = PtxArchProps::LOG_WARP_THREADS,
        WARPS                   = (BLOCK_THREADS + WARP_THREADS - 1) / WARP_THREADS,

        BYTES_PER_SIZET         = sizeof(SizeT),
        LOG_BYTES_PER_SIZET     = Log2<BYTES_PER_SIZET>::VALUE,

        LOG_SMEM_BANKS          = PtxArchProps::LOG_SMEM_BANKS,
        SMEM_BANKS              = 1 << LOG_SMEM_BANKS,

        DIGITS_PER_SCATTER_PASS = BLOCK_THREADS / SMEM_BANKS,
        SCATTER_PASSES          = RADIX_DIGITS / DIGITS_PER_SCATTER_PASS,

        LOG_STORE_TXN_THREADS   = LOG_SMEM_BANKS,
        STORE_TXN_THREADS       = 1 << LOG_STORE_TXN_THREADS,
    };

    // BlockRadixRank type
    typedef BlockRadixRank<
        BLOCK_THREADS,
        RADIX_BITS,
        MEMOIZE_OUTER_SCAN,
        INNER_SCAN_ALGORITHM,
        SMEM_CONFIG> BlockRadixRank;

    // BlockLoad type (keys)
    typedef BlockLoad<
        UnsignedBits*,
        BLOCK_THREADS,
        ITEMS_PER_THREAD,
        LOAD_ALGORITHM,
        LOAD_MODIFIER,
        EXCHANGE_TIME_SLICING> BlockLoadKeys;

    // BlockLoad type (values)
    typedef BlockLoad<
        Value*,
        BLOCK_THREADS,
        ITEMS_PER_THREAD,
        LOAD_ALGORITHM,
        LOAD_MODIFIER,
        EXCHANGE_TIME_SLICING> BlockLoadValues;

    // BlockExchange type (keys)
    typedef BlockExchange<
        UnsignedBits,
        BLOCK_THREADS,
        ITEMS_PER_THREAD,
        EXCHANGE_TIME_SLICING> BlockExchangeKeys;

    // BlockExchange type (values)
    typedef BlockExchange<
        Value,
        BLOCK_THREADS,
        ITEMS_PER_THREAD,
        EXCHANGE_TIME_SLICING> BlockExchangeValues;


    /**
     * Shared memory storage layout
     */
    struct _TempStorage
    {
        SizeT   relative_bin_offsets[RADIX_DIGITS + 1];
        bool    short_circuit;

        union
        {
            typename BlockRadixRank::TempStorage        ranking;
            typename BlockLoadKeys::TempStorage         load_keys;
            typename BlockLoadValues::TempStorage       load_values;
            typename BlockExchangeKeys::TempStorage     exchange_keys;
            typename BlockExchangeValues::TempStorage   exchange_values;
        };
    };


    /// Alias wrapper allowing storage to be unioned
    struct TempStorage : Uninitialized<_TempStorage> {};


    //---------------------------------------------------------------------
    // Thread fields
    //---------------------------------------------------------------------

    // Shared storage for this CTA
    _TempStorage    &temp_storage;

    // Input and output device pointers
    UnsignedBits    *d_keys_in;
    UnsignedBits    *d_keys_out;
    Value           *d_values_in;
    Value           *d_values_out;

    // The global scatter base offset for each digit (valid in the first RADIX_DIGITS threads)
    SizeT           bin_offset;

    // The least-significant bit position of the current digit to extract
    int             current_bit;

    // Whether to short-ciruit
    bool            short_circuit;



    //---------------------------------------------------------------------
    // Utility methods
    //---------------------------------------------------------------------

    /**
     * Decodes given keys to lookup digit offsets in shared memory
     */
    __device__ __forceinline__ void DecodeRelativeBinOffsets(
        UnsignedBits    (&twiddled_keys)[ITEMS_PER_THREAD],
        SizeT           (&relative_bin_offsets)[ITEMS_PER_THREAD])
    {
        #pragma unroll
        for (int KEY = 0; KEY < ITEMS_PER_THREAD; KEY++)
        {
            UnsignedBits digit = BFE(twiddled_keys[KEY], current_bit, RADIX_BITS);

            // Lookup base digit offset from shared memory
            relative_bin_offsets[KEY] = temp_storage.relative_bin_offsets[digit];
        }
    }


    /**
     * Scatter ranked items to global memory
     */
    template <bool FULL_TILE, typename T>
    __device__ __forceinline__ void ScatterItems(
        T       (&items)[ITEMS_PER_THREAD],
        int     (&local_ranks)[ITEMS_PER_THREAD],
        SizeT   (&relative_bin_offsets)[ITEMS_PER_THREAD],
        T       *d_out,
        SizeT   valid_items)
    {
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ++ITEM)
        {
            // Scatter if not out-of-bounds
            if (FULL_TILE || (local_ranks[ITEM] < valid_items))
            {
                d_out[relative_bin_offsets[ITEM] + local_ranks[ITEM]] = items[ITEM];
            }
        }
    }


    /**
     * Scatter ranked keys directly to global memory
     */
    template <bool FULL_TILE>
    __device__ __forceinline__ void ScatterKeys(
        UnsignedBits                            (&twiddled_keys)[ITEMS_PER_THREAD],
        SizeT                                   (&relative_bin_offsets)[ITEMS_PER_THREAD],
        int                                     (&ranks)[ITEMS_PER_THREAD],
        SizeT                                   valid_items,
        Int2Type<RADIX_SORT_SCATTER_DIRECT>     scatter_algorithm)
    {
        // Compute scatter offsets
        DecodeRelativeBinOffsets(twiddled_keys, relative_bin_offsets);

        // Untwiddle keys before outputting
        UnsignedBits keys[ITEMS_PER_THREAD];

        #pragma unroll
        for (int KEY = 0; KEY < ITEMS_PER_THREAD; KEY++)
        {
            keys[KEY] = Traits<Key>::TwiddleOut(twiddled_keys[KEY]);
        }

        // Scatter to global
        ScatterItems<FULL_TILE>(keys, ranks, relative_bin_offsets, d_keys_out, valid_items);
    }


    /**
     * Scatter ranked keys through shared memory, then to global memory
     */
    template <bool FULL_TILE>
    __device__ __forceinline__ void ScatterKeys(
        UnsignedBits                            (&twiddled_keys)[ITEMS_PER_THREAD],
        SizeT                                   (&relative_bin_offsets)[ITEMS_PER_THREAD],
        int                                     (&ranks)[ITEMS_PER_THREAD],
        SizeT                                   valid_items,
        Int2Type<RADIX_SORT_SCATTER_TWO_PHASE>  scatter_algorithm)
    {
        // Exchange keys through shared memory
        BlockExchangeKeys(temp_storage.exchange_keys).ScatterToStriped(twiddled_keys, ranks);

        // Compute striped local ranks
        int local_ranks[ITEMS_PER_THREAD];
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ++ITEM)
        {
            local_ranks[ITEM] = threadIdx.x + (ITEM * BLOCK_THREADS);
        }

        // Scatter directly
        ScatterKeys<FULL_TILE>(
            twiddled_keys,
            relative_bin_offsets,
            local_ranks,
            valid_items,
            Int2Type<RADIX_SORT_SCATTER_DIRECT>());
    }


    /**
     * Scatter ranked values directly to global memory
     */
    template <bool FULL_TILE>
    __device__ __forceinline__ void ScatterValues(
        Value                                   (&values)[ITEMS_PER_THREAD],
        SizeT                                   (&relative_bin_offsets)[ITEMS_PER_THREAD],
        int                                     (&ranks)[ITEMS_PER_THREAD],
        SizeT                                   valid_items,
        Int2Type<RADIX_SORT_SCATTER_DIRECT>     scatter_algorithm)
    {
        // Scatter to global
        ScatterItems<FULL_TILE>(values, ranks, relative_bin_offsets, d_values_out, valid_items);
    }


    /**
     * Scatter ranked values through shared memory, then to global memory
     */
    template <bool FULL_TILE>
    __device__ __forceinline__ void ScatterValues(
        Value                                   (&values)[ITEMS_PER_THREAD],
        SizeT                                   (&relative_bin_offsets)[ITEMS_PER_THREAD],
        int                                     (&ranks)[ITEMS_PER_THREAD],
        SizeT                                   valid_items,
        Int2Type<RADIX_SORT_SCATTER_TWO_PHASE>  scatter_algorithm)
    {
        __syncthreads();

        // Exchange keys through shared memory
        BlockExchangeValues(temp_storage.exchange_values).ScatterToStriped(values, ranks);

        // Compute striped local ranks
        int local_ranks[ITEMS_PER_THREAD];
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ++ITEM)
        {
            local_ranks[ITEM] = threadIdx.x + (ITEM * BLOCK_THREADS);
        }

        // Scatter directly
        ScatterValues<FULL_TILE>(
            values,
            relative_bin_offsets,
            local_ranks,
            valid_items,
            Int2Type<RADIX_SORT_SCATTER_DIRECT>());
    }


    /**
     * Load a tile of items (specialized for full tile)
     */
    template <typename BlockLoadT, typename T>
    __device__ __forceinline__ void LoadItems(
        BlockLoadT      &block_loader, 
        T               (&items)[ITEMS_PER_THREAD],
        T               *d_in, 
        SizeT           valid_items, 
        Int2Type<true>  is_full_tile)
    {
        block_loader.Load(d_in, items);
    }


    /**
     * Load a tile of items (specialized for partial tile)
     */
    template <typename BlockLoadT, typename T>
    __device__ __forceinline__ void LoadItems(
        BlockLoadT      &block_loader, 
        T               (&items)[ITEMS_PER_THREAD],
        T               *d_in, 
        SizeT           valid_items, 
        Int2Type<false> is_full_tile)
    {
        block_loader.Load(d_in, items, valid_items);
    }


    /**
     * Truck along associated values
     */
    template <bool FULL_TILE, typename _Value>
    __device__ __forceinline__ void GatherScatterValues(
        _Value      (&values)[ITEMS_PER_THREAD],
        SizeT       (&relative_bin_offsets)[ITEMS_PER_THREAD],
        int         (&ranks)[ITEMS_PER_THREAD],
        SizeT       block_offset,
        SizeT       valid_items)
    {
        BlockLoadValues loader(temp_storage.load_values);
        LoadItems(
            loader,
            values,
            d_values_in + block_offset,
            valid_items,
            Int2Type<FULL_TILE>());

        ScatterValues<FULL_TILE>(
            values,
            relative_bin_offsets,
            ranks,
            valid_items,
            Int2Type<SCATTER_ALGORITHM>());
    }


    /**
     * Truck along associated values (specialized for key-only sorting)
     */
    template <bool FULL_TILE>
    __device__ __forceinline__ void GatherScatterValues(
        NullType    (&values)[ITEMS_PER_THREAD],
        SizeT       (&relative_bin_offsets)[ITEMS_PER_THREAD],
        int         (&ranks)[ITEMS_PER_THREAD],
        SizeT       block_offset,
        SizeT       valid_items)
    {}


    /**
     * Process tile
     */
    template <bool FULL_TILE>
    __device__ __forceinline__ void ProcessTile(
        SizeT block_offset,
        const SizeT &valid_items = TILE_ITEMS)
    {
        // Per-thread tile data
        UnsignedBits    keys[ITEMS_PER_THREAD];                     // Keys
        UnsignedBits    twiddled_keys[ITEMS_PER_THREAD];            // Twiddled keys
        int             ranks[ITEMS_PER_THREAD];                    // For each key, the local rank within the CTA
        SizeT           relative_bin_offsets[ITEMS_PER_THREAD];     // For each key, the global scatter base offset of the corresponding digit

        if (LOAD_ALGORITHM != BLOCK_LOAD_DIRECT) __syncthreads();

        // Assign max-key to all keys
        #pragma unroll
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ++ITEM)
        {
            keys[ITEM] = MAX_KEY;
        }

        // Load tile of keys
        BlockLoadKeys loader(temp_storage.load_keys);
        LoadItems(
            loader,
            keys,
            d_keys_in + block_offset,
            valid_items, 
            Int2Type<FULL_TILE>());

        __syncthreads();

        // Twiddle key bits if necessary
        #pragma unroll
        for (int KEY = 0; KEY < ITEMS_PER_THREAD; KEY++)
        {
            twiddled_keys[KEY] = Traits<Key>::TwiddleIn(keys[KEY]);
        }

        // Rank the twiddled keys
        int inclusive_digit_prefix;
        BlockRadixRank(temp_storage.ranking).RankKeys(
            twiddled_keys,
            ranks,
            current_bit,
            inclusive_digit_prefix);

        // Update global scatter base offsets for each digit
        if ((BLOCK_THREADS == RADIX_DIGITS) || (threadIdx.x < RADIX_DIGITS))
        {
            int exclusive_digit_prefix;

            // Get exclusive digit prefix from inclusive prefix
#if CUB_PTX_ARCH >= 300
            exclusive_digit_prefix = ShuffleUp(inclusive_digit_prefix, 1);
            if (threadIdx.x == 0)
                exclusive_digit_prefix = 0;
#else
            volatile int* exchange = reinterpret_cast<int *>(temp_storage.relative_bin_offsets);
            exchange[threadIdx.x] = 0;
            exchange[threadIdx.x + 1] = inclusive_digit_prefix;
            exclusive_digit_prefix = exchange[threadIdx.x];
#endif

            bin_offset -= exclusive_digit_prefix;
            temp_storage.relative_bin_offsets[threadIdx.x] = bin_offset;
            bin_offset += inclusive_digit_prefix;
        }

        __syncthreads();

        // Scatter keys
        ScatterKeys<FULL_TILE>(twiddled_keys, relative_bin_offsets, ranks, valid_items, Int2Type<SCATTER_ALGORITHM>());

        // Gather/scatter values
        Value values[ITEMS_PER_THREAD];
        GatherScatterValues<FULL_TILE>(values, relative_bin_offsets, ranks, block_offset, valid_items);
    }


    /**
     * Copy tiles within the range of input
     */
    template <typename T>
    __device__ __forceinline__ void Copy(
        T       *d_in,
        T       *d_out,
        SizeT   block_offset,
        SizeT   block_oob)
    {
        // Simply copy the input
        while (block_offset + TILE_ITEMS <= block_oob)
        {
            T items[ITEMS_PER_THREAD];

            LoadStriped<LOAD_DEFAULT, BLOCK_THREADS>(threadIdx.x, d_in + block_offset, items);
            __syncthreads();
            StoreStriped<STORE_DEFAULT, BLOCK_THREADS>(threadIdx.x, d_out + block_offset, items);

            block_offset += TILE_ITEMS;
        }

        // Clean up last partial tile with guarded-I/O
        if (block_offset < block_oob)
        {
            SizeT valid_items = block_oob - block_offset;

            T items[ITEMS_PER_THREAD];

            LoadStriped<LOAD_DEFAULT, BLOCK_THREADS>(threadIdx.x, d_in + block_offset, items, valid_items);
            __syncthreads();
            StoreStriped<STORE_DEFAULT, BLOCK_THREADS>(threadIdx.x, d_out + block_offset, items, valid_items);
        }
    }


    /**
     * Copy tiles within the range of input (specialized for NullType)
     */
    __device__ __forceinline__ void Copy(
        NullType    *d_in,
        NullType    *d_out,
        SizeT       block_offset,
        SizeT       block_oob)
    {}


    //---------------------------------------------------------------------
    // Interface
    //---------------------------------------------------------------------

    /**
     * Constructor
     */
    __device__ __forceinline__ BlockRadixSortDownsweepTiles(
        TempStorage &temp_storage,
        SizeT       bin_offset,
        Key         *d_keys_in,
        Key         *d_keys_out,
        Value       *d_values_in,
        Value       *d_values_out,
        int         current_bit)
    :
        temp_storage(temp_storage.Alias()),
        bin_offset(bin_offset),
        d_keys_in(reinterpret_cast<UnsignedBits*>(d_keys_in)),
        d_keys_out(reinterpret_cast<UnsignedBits*>(d_keys_out)),
        d_values_in(d_values_in),
        d_values_out(d_values_out),
        current_bit(current_bit),
        short_circuit(false)
    {}


    /**
     * Constructor
     */
    __device__ __forceinline__ BlockRadixSortDownsweepTiles(
        TempStorage &temp_storage,
        SizeT       num_items,
        SizeT       *d_spine,
        Key         *d_keys_in,
        Key         *d_keys_out,
        Value       *d_values_in,
        Value       *d_values_out,
        int         current_bit)
    :
        temp_storage(temp_storage.Alias()),
        d_keys_in(reinterpret_cast<UnsignedBits*>(d_keys_in)),
        d_keys_out(reinterpret_cast<UnsignedBits*>(d_keys_out)),
        d_values_in(d_values_in),
        d_values_out(d_values_out),
        current_bit(current_bit)
    {
        // Load digit bin offsets (each of the first RADIX_DIGITS threads will load an offset for that digit)
        if (threadIdx.x < RADIX_DIGITS)
        {
            // Short circuit if the first block's histogram has only bin counts of only zeros or problem-size
            SizeT first_block_bin_offset = d_spine[gridDim.x * threadIdx.x];
            int predicate = ((first_block_bin_offset == 0) || (first_block_bin_offset == num_items));
            this->temp_storage.short_circuit = WarpAll(predicate);

            // Load my block's bin offset for my bin
            bin_offset = d_spine[(gridDim.x * threadIdx.x) + blockIdx.x];
        }

        __syncthreads();

        short_circuit = this->temp_storage.short_circuit;
    }


    /**
     * Distribute keys from a segment of input tiles.
     */
    __device__ __forceinline__ void ProcessTiles(
        SizeT           block_offset,
        const SizeT     &block_oob)
    {
        if (short_circuit)
        {
            // Copy keys
            Copy(d_keys_in, d_keys_out, block_offset, block_oob);

            // Copy values
            Copy(d_values_in, d_values_out, block_offset, block_oob);
        }
        else
        {
            // Process full tiles of tile_items
            while (block_offset + TILE_ITEMS <= block_oob)
            {
                ProcessTile<true>(block_offset);
                block_offset += TILE_ITEMS;
            }

            // Clean up last partial tile with guarded-I/O
            if (block_offset < block_oob)
            {
                ProcessTile<false>(block_offset, block_oob - block_offset);
            }
        }
    }
};



}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)

