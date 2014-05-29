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
 * cub::BlockRadixRank provides operations for ranking unsigned integer types within a CUDA threadblock
 */

#pragma once

#include "../util_arch.cuh"
#include "../util_type.cuh"
#include "../thread/thread_reduce.cuh"
#include "../thread/thread_scan.cuh"
#include "../block/block_scan.cuh"
#include "../util_namespace.cuh"


/// Optional outer namespace(s)
CUB_NS_PREFIX

/// CUB namespace
namespace cub {

/**
 * \brief BlockRadixRank provides operations for ranking unsigned integer types within a CUDA threadblock.
 * \ingroup BlockModule
 *
 * \par Overview
 * Blah...
 *
 * \tparam BLOCK_THREADS        The thread block size in threads
 * \tparam RADIX_BITS           <b>[optional]</b> The number of radix bits per digit place (default: 5 bits)
 * \tparam MEMOIZE_OUTER_SCAN   <b>[optional]</b> Whether or not to buffer outer raking scan partials to incur fewer shared memory reads at the expense of higher register pressure (default: true for architectures SM35 and newer, false otherwise).  See BlockScanAlgorithm::BLOCK_SCAN_RAKING_MEMOIZE for more details.
 * \tparam INNER_SCAN_ALGORITHM <b>[optional]</b> The cub::BlockScanAlgorithm algorithm to use (default: cub::BLOCK_SCAN_WARP_SCANS)
 * \tparam SMEM_CONFIG          <b>[optional]</b> Shared memory bank mode (default: \p cudaSharedMemBankSizeFourByte)
 *
 * \par Usage Considerations
 * - Keys must be in a form suitable for radix ranking (i.e., unsigned bits).
 * - Assumes a [<em>blocked arrangement</em>](index.html#sec5sec4) of elements across threads
 * - \smemreuse{BlockRadixRank::TempStorage}
 *
 * \par Performance Considerations
 *
 * \par Algorithm
 * These parallel radix ranking variants have <em>O</em>(<em>n</em>) work complexity and are implemented in XXX phases:
 * -# blah
 * -# blah
 *
 * \par Examples
 * \par
 * - <b>Example 1:</b> Simple radix rank of 32-bit integer keys
 *      \code
 *      #include <cub/cub.cuh>
 *
 *      template <int BLOCK_THREADS>
 *      __global__ void ExampleKernel(...)
 *      {
 *
 *      \endcode
 */
template <
    int                     BLOCK_THREADS,
    int                     RADIX_BITS,
    bool                    MEMOIZE_OUTER_SCAN      = (CUB_PTX_ARCH >= 350) ? true : false,
    BlockScanAlgorithm      INNER_SCAN_ALGORITHM    = BLOCK_SCAN_WARP_SCANS,
    cudaSharedMemConfig     SMEM_CONFIG             = cudaSharedMemBankSizeFourByte>
class BlockRadixRank
{
private:

    /******************************************************************************
     * Type definitions and constants
     ******************************************************************************/

    // Integer type for digit counters (to be packed into words of type PackedCounters)
    typedef unsigned short DigitCounter;

    // Integer type for packing DigitCounters into columns of shared memory banks
    typedef typename If<(SMEM_CONFIG == cudaSharedMemBankSizeEightByte),
        unsigned long long,
        unsigned int>::Type PackedCounter;

    enum
    {
        RADIX_DIGITS                 = 1 << RADIX_BITS,

        LOG_WARP_THREADS             = PtxArchProps::LOG_WARP_THREADS,
        WARP_THREADS                 = 1 << LOG_WARP_THREADS,
        WARPS                        = (BLOCK_THREADS + WARP_THREADS - 1) / WARP_THREADS,

        BYTES_PER_COUNTER            = sizeof(DigitCounter),
        LOG_BYTES_PER_COUNTER        = Log2<BYTES_PER_COUNTER>::VALUE,

        PACKING_RATIO                = sizeof(PackedCounter) / sizeof(DigitCounter),
        LOG_PACKING_RATIO            = Log2<PACKING_RATIO>::VALUE,

        LOG_COUNTER_LANES            = CUB_MAX((RADIX_BITS - LOG_PACKING_RATIO), 0),                // Always at least one lane
        COUNTER_LANES                = 1 << LOG_COUNTER_LANES,

        // The number of packed counters per thread (plus one for padding)
        RAKING_SEGMENT               = COUNTER_LANES + 1,

        LOG_SMEM_BANKS               = PtxArchProps::LOG_SMEM_BANKS,
        SMEM_BANKS                   = 1 << LOG_SMEM_BANKS,
    };


    /// BlockScan type
    typedef BlockScan<PackedCounter, BLOCK_THREADS, INNER_SCAN_ALGORITHM> BlockScan;


    /// Shared memory storage layout type for BlockRadixRank
    struct _TempStorage
    {
        // Storage for scanning local ranks
        typename BlockScan::TempStorage block_scan;

        union
        {
            DigitCounter            digit_counters[COUNTER_LANES + 1][BLOCK_THREADS][PACKING_RATIO];
            PackedCounter           raking_grid[BLOCK_THREADS][RAKING_SEGMENT];
        };
    };


    /******************************************************************************
     * Thread fields
     ******************************************************************************/

    /// Shared storage reference
    _TempStorage &temp_storage;

    /// Linear thread-id
    int linear_tid;

    /// Copy of raking segment, promoted to registers
    PackedCounter cached_segment[RAKING_SEGMENT];


    /******************************************************************************
     * Templated iteration
     ******************************************************************************/

    // General template iteration
    template <int COUNT, int MAX>
    struct Iterate
    {
        /**
         * Decode keys.  Decodes the radix digit from the current digit place
         * and increments the thread's corresponding counter in shared
         * memory for that digit.
         *
         * Saves both (1) the prior value of that counter (the key's
         * thread-local exclusive prefix sum for that digit), and (2) the shared
         * memory offset of the counter (for later use).
         */
        template <typename UnsignedBits, int KEYS_PER_THREAD>
        static __device__ __forceinline__ void DecodeKeys(
            BlockRadixRank  &cta,                                   // BlockRadixRank instance
            UnsignedBits    (&keys)[KEYS_PER_THREAD],               // Key to decode
            DigitCounter    (&thread_prefixes)[KEYS_PER_THREAD],    // Prefix counter value (out parameter)
            DigitCounter*   (&digit_counters)[KEYS_PER_THREAD],     // Counter smem offset (out parameter)
            int             current_bit)                            // The least-significant bit position of the current digit to extract
        {
            // Add in sub-counter offset
            UnsignedBits sub_counter = BFE(keys[COUNT], current_bit + LOG_COUNTER_LANES, LOG_PACKING_RATIO);

            // Add in row offset
            UnsignedBits row_offset = BFE(keys[COUNT], current_bit, LOG_COUNTER_LANES);

            // Pointer to smem digit counter
            digit_counters[COUNT] = &cta.temp_storage.digit_counters[row_offset][cta.linear_tid][sub_counter];

            // Load thread-exclusive prefix
            thread_prefixes[COUNT] = *digit_counters[COUNT];

            // Store inclusive prefix
            *digit_counters[COUNT] = thread_prefixes[COUNT] + 1;

            // Iterate next key
            Iterate<COUNT + 1, MAX>::DecodeKeys(cta, keys, thread_prefixes, digit_counters, current_bit);
        }


        // Termination
        template <int KEYS_PER_THREAD>
        static __device__ __forceinline__ void UpdateRanks(
            int             (&ranks)[KEYS_PER_THREAD],              // Local ranks (out parameter)
            DigitCounter    (&thread_prefixes)[KEYS_PER_THREAD],    // Prefix counter value
            DigitCounter*   (&digit_counters)[KEYS_PER_THREAD])     // Counter smem offset
        {
            // Add in threadblock exclusive prefix
            ranks[COUNT] = thread_prefixes[COUNT] + *digit_counters[COUNT];

            // Iterate next key
            Iterate<COUNT + 1, MAX>::UpdateRanks(ranks, thread_prefixes, digit_counters);
        }
    };


    // Termination
    template <int MAX>
    struct Iterate<MAX, MAX>
    {
        // DecodeKeys
        template <typename UnsignedBits, int KEYS_PER_THREAD>
        static __device__ __forceinline__ void DecodeKeys(
            BlockRadixRank  &cta,
            UnsignedBits    (&keys)[KEYS_PER_THREAD],
            DigitCounter    (&thread_prefixes)[KEYS_PER_THREAD],
            DigitCounter*   (&digit_counters)[KEYS_PER_THREAD],
            int             current_bit) {}


        // UpdateRanks
        template <int KEYS_PER_THREAD>
        static __device__ __forceinline__ void UpdateRanks(
            int             (&ranks)[KEYS_PER_THREAD],
            DigitCounter    (&thread_prefixes)[KEYS_PER_THREAD],
            DigitCounter    *(&digit_counters)[KEYS_PER_THREAD]) {}
    };


    /******************************************************************************
     * Utility methods
     ******************************************************************************/

    /**
     * Internal storage allocator
     */
    __device__ __forceinline__ _TempStorage& PrivateStorage()
    {
        __shared__ _TempStorage private_storage;
        return private_storage;
    }


    /**
     * Performs upsweep raking reduction, returning the aggregate
     */
    __device__ __forceinline__ PackedCounter Upsweep()
    {
        PackedCounter *smem_raking_ptr = temp_storage.raking_grid[linear_tid];
        PackedCounter *raking_ptr;

        if (MEMOIZE_OUTER_SCAN)
        {
            // Copy data into registers
            #pragma unroll
            for (int i = 0; i < RAKING_SEGMENT; i++)
            {
                cached_segment[i] = smem_raking_ptr[i];
            }
            raking_ptr = cached_segment;
        }
        else
        {
            raking_ptr = smem_raking_ptr;
        }

        return ThreadReduce<RAKING_SEGMENT>(raking_ptr, Sum());
    }


    /// Performs exclusive downsweep raking scan
    __device__ __forceinline__ void ExclusiveDownsweep(
        PackedCounter raking_partial)
    {
        PackedCounter *smem_raking_ptr = temp_storage.raking_grid[linear_tid];

        PackedCounter *raking_ptr = (MEMOIZE_OUTER_SCAN) ?
            cached_segment :
            smem_raking_ptr;

        // Exclusive raking downsweep scan
        ThreadScanExclusive<RAKING_SEGMENT>(raking_ptr, raking_ptr, Sum(), raking_partial);

        if (MEMOIZE_OUTER_SCAN)
        {
            // Copy data back to smem
            #pragma unroll
            for (int i = 0; i < RAKING_SEGMENT; i++)
            {
                smem_raking_ptr[i] = cached_segment[i];
            }
        }
    }


    /**
     * Reset shared memory digit counters
     */
    __device__ __forceinline__ void ResetCounters()
    {
        // Reset shared memory digit counters
        #pragma unroll
        for (int LANE = 0; LANE < COUNTER_LANES + 1; LANE++)
        {
            *((PackedCounter*) temp_storage.digit_counters[LANE][linear_tid]) = 0;
        }
    }


    /**
     * Scan shared memory digit counters.
     */
    __device__ __forceinline__ void ScanCounters()
    {
        // Upsweep scan
        PackedCounter raking_partial = Upsweep();

        // Compute inclusive sum
        PackedCounter inclusive_partial;
        PackedCounter packed_aggregate;
        BlockScan(temp_storage.block_scan, linear_tid).InclusiveSum(raking_partial, inclusive_partial, packed_aggregate);

        // Propagate totals in packed fields
        #pragma unroll
        for (int PACKED = 1; PACKED < PACKING_RATIO; PACKED++)
        {
            inclusive_partial += packed_aggregate << (sizeof(DigitCounter) * 8 * PACKED);
        }

        // Downsweep scan with exclusive partial
        PackedCounter exclusive_partial = inclusive_partial - raking_partial;
        ExclusiveDownsweep(exclusive_partial);
    }

public:

    /// \smemstorage{BlockScan}
    struct TempStorage : Uninitialized<_TempStorage> {};


    /******************************************************************//**
     * \name Collective constructors
     *********************************************************************/
    //@{

    /**
     * \brief Collective constructor for 1D thread blocks using a private static allocation of shared memory as temporary storage.  Threads are identified using <tt>threadIdx.x</tt>.
     */
    __device__ __forceinline__ BlockRadixRank()
    :
        temp_storage(PrivateStorage()),
        linear_tid(threadIdx.x)
    {}


    /**
     * \brief Collective constructor for 1D thread blocks using the specified memory allocation as temporary storage.  Threads are identified using <tt>threadIdx.x</tt>.
     */
    __device__ __forceinline__ BlockRadixRank(
        TempStorage &temp_storage)             ///< [in] Reference to memory allocation having layout type TempStorage
    :
        temp_storage(temp_storage.Alias()),
        linear_tid(threadIdx.x)
    {}


    /**
     * \brief Collective constructor using a private static allocation of shared memory as temporary storage.  Each thread is identified using the supplied linear thread identifier
     */
    __device__ __forceinline__ BlockRadixRank(
        int linear_tid)                        ///< [in] A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    :
        temp_storage(PrivateStorage()),
        linear_tid(linear_tid)
    {}


    /**
     * \brief Collective constructor using the specified memory allocation as temporary storage.  Each thread is identified using the supplied linear thread identifier.
     */
    __device__ __forceinline__ BlockRadixRank(
        TempStorage &temp_storage,             ///< [in] Reference to memory allocation having layout type TempStorage
        int linear_tid)                        ///< [in] <b>[optional]</b> A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    :
        temp_storage(temp_storage.Alias()),
        linear_tid(linear_tid)
    {}



    //@}  end member group
    /******************************************************************//**
     * \name Raking
     *********************************************************************/
    //@{

    /**
     * \brief Rank keys.
     */
    template <
        typename        UnsignedBits,
        int             KEYS_PER_THREAD>
    __device__ __forceinline__ void RankKeys(
        UnsignedBits    (&keys)[KEYS_PER_THREAD],           ///< [in] Keys for this tile
        int             (&ranks)[KEYS_PER_THREAD],          ///< [out] For each key, the local rank within the tile
        int             current_bit)                        ///< [in] The least-significant bit position of the current digit to extract
    {
        DigitCounter    thread_prefixes[KEYS_PER_THREAD];   // For each key, the count of previous keys in this tile having the same digit
        DigitCounter*   digit_counters[KEYS_PER_THREAD];    // For each key, the byte-offset of its corresponding digit counter in smem

        // Reset shared memory digit counters
        ResetCounters();

        // Decode keys and update digit counters
        Iterate<0, KEYS_PER_THREAD>::DecodeKeys(*this, keys, thread_prefixes, digit_counters, current_bit);

        __syncthreads();

        // Scan shared memory counters
        ScanCounters();

        __syncthreads();

        // Extract the local ranks of each key
        Iterate<0, KEYS_PER_THREAD>::UpdateRanks(ranks, thread_prefixes, digit_counters);
    }


    /**
     * \brief Rank keys.  For the lower \p RADIX_DIGITS threads, digit counts for each digit are provided for the corresponding thread.
     */
    template <
        typename        UnsignedBits,
        int             KEYS_PER_THREAD>
    __device__ __forceinline__ void RankKeys(
        UnsignedBits    (&keys)[KEYS_PER_THREAD],           ///< [in] Keys for this tile
        int             (&ranks)[KEYS_PER_THREAD],          ///< [out] For each key, the local rank within the tile (out parameter)
        int             current_bit,                        ///< [in] The least-significant bit position of the current digit to extract
        int             &inclusive_digit_prefix)            ///< [out] The incluisve prefix sum for the digit threadIdx.x
    {
        // Rank keys
        RankKeys(keys, ranks, current_bit);

        // Get the inclusive and exclusive digit totals corresponding to the calling thread.
        if ((BLOCK_THREADS == RADIX_DIGITS) || (linear_tid < RADIX_DIGITS))
        {
            // Obtain ex/inclusive digit counts.  (Unfortunately these all reside in the
            // first counter column, resulting in unavoidable bank conflicts.)
            int counter_lane = (linear_tid & (COUNTER_LANES - 1));
            int sub_counter = linear_tid >> (LOG_COUNTER_LANES);
            inclusive_digit_prefix = temp_storage.digit_counters[counter_lane + 1][0][sub_counter];
        }
    }
};

}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)


