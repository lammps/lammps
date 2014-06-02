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
 * The cub::BlockExchange class provides [<em>collective</em>](index.html#sec0) methods for rearranging data partitioned across a CUDA thread block.
 */

#pragma once

#include "../util_arch.cuh"
#include "../util_macro.cuh"
#include "../util_type.cuh"
#include "../util_namespace.cuh"

/// Optional outer namespace(s)
CUB_NS_PREFIX

/// CUB namespace
namespace cub {

/**
 * \brief The BlockExchange class provides [<em>collective</em>](index.html#sec0) methods for rearranging data partitioned across a CUDA thread block. ![](transpose_logo.png)
 * \ingroup BlockModule
 *
 * \par Overview
 * It is commonplace for blocks of threads to rearrange data items between
 * threads.  For example, the global memory subsystem prefers access patterns
 * where data items are "striped" across threads (where consecutive threads access consecutive items),
 * yet most block-wide operations prefer a "blocked" partitioning of items across threads
 * (where consecutive items belong to a single thread).
 *
 * \par
 * BlockExchange supports the following types of data exchanges:
 * - Transposing between [<em>blocked</em>](index.html#sec5sec4) and [<em>striped</em>](index.html#sec5sec4) arrangements
 * - Transposing between [<em>blocked</em>](index.html#sec5sec4) and [<em>warp-striped</em>](index.html#sec5sec4) arrangements
 * - Scattering ranked items to a [<em>blocked arrangement</em>](index.html#sec5sec4)
 * - Scattering ranked items to a [<em>striped arrangement</em>](index.html#sec5sec4)
 *
 * \tparam T                    The data type to be exchanged.
 * \tparam BLOCK_THREADS        The thread block size in threads.
 * \tparam ITEMS_PER_THREAD     The number of items partitioned onto each thread.
 * \tparam WARP_TIME_SLICING    <b>[optional]</b> When \p true, only use enough shared memory for a single warp's worth of tile data, time-slicing the block-wide exchange over multiple synchronized rounds.  Yields a smaller memory footprint at the expense of decreased parallelism.  (Default: false)
 *
 * \par A Simple Example
 * \blockcollective{BlockExchange}
 * \par
 * The code snippet below illustrates the conversion from a "blocked" to a "striped" arrangement
 * of 512 integer items partitioned across 128 threads where each thread owns 4 items.
 * \par
 * \code
 * #include <cub/cub.cuh>
 *
 * __global__ void ExampleKernel(int *d_data, ...)
 * {
 *     // Specialize BlockExchange for 128 threads owning 4 integer items each
 *     typedef cub::BlockExchange<int, 128, 4> BlockExchange;
 *
 *     // Allocate shared memory for BlockExchange
 *     __shared__ typename BlockExchange::TempStorage temp_storage;
 *
 *     // Load a tile of data striped across threads
 *     int thread_data[4];
 *     cub::LoadStriped<LOAD_DEFAULT, 128>(threadIdx.x, d_data, thread_data);
 *
 *     // Collectively exchange data into a blocked arrangement across threads
 *     BlockExchange(temp_storage).StripedToBlocked(thread_data);
 *
 * \endcode
 * \par
 * Suppose the set of striped input \p thread_data across the block of threads is
 * <tt>{ [0,128,256,384], [1,129,257,385], ..., [127,255,383,511] }</tt>.
 * The corresponding output \p thread_data in those threads will be
 * <tt>{ [0,1,2,3], [4,5,6,7], [8,9,10,11], ..., [508,509,510,511] }</tt>.
 *
 * \par Performance Considerations
 * - Proper device-specific padding ensures zero bank conflicts for most types.
 *
 */
template <
    typename        T,
    int             BLOCK_THREADS,
    int             ITEMS_PER_THREAD,
    bool            WARP_TIME_SLICING = false>
class BlockExchange
{
private:

    /******************************************************************************
     * Constants
     ******************************************************************************/

    enum
    {
        LOG_WARP_THREADS            = PtxArchProps::LOG_WARP_THREADS,
        WARP_THREADS                = 1 << LOG_WARP_THREADS,
        WARPS                       = (BLOCK_THREADS + PtxArchProps::WARP_THREADS - 1) / PtxArchProps::WARP_THREADS,

        LOG_SMEM_BANKS              = PtxArchProps::LOG_SMEM_BANKS,
        SMEM_BANKS                  = 1 << LOG_SMEM_BANKS,

        TILE_ITEMS                  = BLOCK_THREADS * ITEMS_PER_THREAD,

        TIME_SLICES                 = (WARP_TIME_SLICING) ? WARPS : 1,

        TIME_SLICED_THREADS         = (WARP_TIME_SLICING) ? CUB_MIN(BLOCK_THREADS, WARP_THREADS) : BLOCK_THREADS,
        TIME_SLICED_ITEMS           = TIME_SLICED_THREADS * ITEMS_PER_THREAD,

        WARP_TIME_SLICED_THREADS    = CUB_MIN(BLOCK_THREADS, WARP_THREADS),
        WARP_TIME_SLICED_ITEMS      = WARP_TIME_SLICED_THREADS * ITEMS_PER_THREAD,

        // Insert padding if the number of items per thread is a power of two
        INSERT_PADDING              = ((ITEMS_PER_THREAD & (ITEMS_PER_THREAD - 1)) == 0),
        PADDING_ITEMS               = (INSERT_PADDING) ? (TIME_SLICED_ITEMS >> LOG_SMEM_BANKS) : 0,
    };

    /******************************************************************************
     * Type definitions
     ******************************************************************************/

    /// Shared memory storage layout type
    typedef T _TempStorage[TIME_SLICED_ITEMS + PADDING_ITEMS];

public:

    /// \smemstorage{BlockExchange}
    struct TempStorage : Uninitialized<_TempStorage> {};

private:


    /******************************************************************************
     * Thread fields
     ******************************************************************************/

    /// Shared storage reference
    _TempStorage &temp_storage;

    /// Linear thread-id
    int linear_tid;
    int warp_lane;
    int warp_id;
    int warp_offset;


    /******************************************************************************
     * Utility methods
     ******************************************************************************/

    /// Internal storage allocator
    __device__ __forceinline__ _TempStorage& PrivateStorage()
    {
        __shared__ _TempStorage private_storage;
        return private_storage;
    }


    /**
     * Transposes data items from <em>blocked</em> arrangement to <em>striped</em> arrangement.  Specialized for no timeslicing.
     */
    __device__ __forceinline__ void BlockedToStriped(
        T               items[ITEMS_PER_THREAD],    ///< [in-out] Items to exchange, converting between <em>blocked</em> and <em>striped</em> arrangements.
        Int2Type<false> time_slicing)
    {
        #pragma unroll
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
        {
            int item_offset = (linear_tid * ITEMS_PER_THREAD) + ITEM;
            if (INSERT_PADDING) item_offset += item_offset >> LOG_SMEM_BANKS;
            temp_storage[item_offset] = items[ITEM];
        }

        __syncthreads();

        #pragma unroll
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
        {
            int item_offset = int(ITEM * BLOCK_THREADS) + linear_tid;
            if (INSERT_PADDING) item_offset += item_offset >> LOG_SMEM_BANKS;
            items[ITEM] = temp_storage[item_offset];
        }
    }


    /**
     * Transposes data items from <em>blocked</em> arrangement to <em>striped</em> arrangement.  Specialized for warp-timeslicing.
     */
    __device__ __forceinline__ void BlockedToStriped(
        T               items[ITEMS_PER_THREAD],    ///< [in-out] Items to exchange, converting between <em>blocked</em> and <em>striped</em> arrangements.
        Int2Type<true>  time_slicing)
    {
        T temp_items[ITEMS_PER_THREAD];

        #pragma unroll
        for (int SLICE = 0; SLICE < TIME_SLICES; SLICE++)
        {
            const int SLICE_OFFSET  = SLICE * TIME_SLICED_ITEMS;
            const int SLICE_OOB     = SLICE_OFFSET + TIME_SLICED_ITEMS;

            __syncthreads();

            if (warp_id == SLICE)
            {
                #pragma unroll
                for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
                {
                    int item_offset = (warp_lane * ITEMS_PER_THREAD) + ITEM;
                    if (INSERT_PADDING) item_offset += item_offset >> LOG_SMEM_BANKS;
                    temp_storage[item_offset] = items[ITEM];
                }
            }

            __syncthreads();

            #pragma unroll
            for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
            {
                // Read a strip of items
                const int STRIP_OFFSET  = ITEM * BLOCK_THREADS;
                const int STRIP_OOB     = STRIP_OFFSET + BLOCK_THREADS;

                if ((SLICE_OFFSET < STRIP_OOB) && (SLICE_OOB > STRIP_OFFSET))
                {
                    int item_offset = STRIP_OFFSET + linear_tid - SLICE_OFFSET;
                    if ((item_offset >= 0) && (item_offset < TIME_SLICED_ITEMS))
                    {
                        if (INSERT_PADDING) item_offset += item_offset >> LOG_SMEM_BANKS;
                        temp_items[ITEM] = temp_storage[item_offset];
                    }
                }
            }
        }

        // Copy
        #pragma unroll
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
        {
            items[ITEM] = temp_items[ITEM];
        }
    }


    /**
     * Transposes data items from <em>blocked</em> arrangement to <em>warp-striped</em> arrangement. Specialized for no timeslicing
     */
    __device__ __forceinline__ void BlockedToWarpStriped(
        T               items[ITEMS_PER_THREAD],   ///< [in-out] Items to exchange, converting between <em>blocked</em> and <em>warp-striped</em> arrangements.
        Int2Type<false> time_slicing)
    {
        #pragma unroll
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
        {
            int item_offset = warp_offset + ITEM + (warp_lane * ITEMS_PER_THREAD);
            if (INSERT_PADDING) item_offset += item_offset >> LOG_SMEM_BANKS;
            temp_storage[item_offset] = items[ITEM];
        }

        #pragma unroll
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
        {
            int item_offset = warp_offset + (ITEM * WARP_TIME_SLICED_THREADS) + warp_lane;
            if (INSERT_PADDING) item_offset += item_offset >> LOG_SMEM_BANKS;
            items[ITEM] = temp_storage[item_offset];
        }
    }

    /**
     * Transposes data items from <em>blocked</em> arrangement to <em>warp-striped</em> arrangement. Specialized for warp-timeslicing
     */
    __device__ __forceinline__ void BlockedToWarpStriped(
        T               items[ITEMS_PER_THREAD],   ///< [in-out] Items to exchange, converting between <em>blocked</em> and <em>warp-striped</em> arrangements.
        Int2Type<true>  time_slicing)
    {
        #pragma unroll
        for (int SLICE = 0; SLICE < TIME_SLICES; ++SLICE)
        {
            __syncthreads();

            if (warp_id == SLICE)
            {
                #pragma unroll
                for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
                {
                    int item_offset = ITEM + (warp_lane * ITEMS_PER_THREAD);
                    if (INSERT_PADDING) item_offset += item_offset >> LOG_SMEM_BANKS;
                    temp_storage[item_offset] = items[ITEM];
                }

                #pragma unroll
                for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
                {
                    int item_offset = (ITEM * WARP_TIME_SLICED_THREADS) + warp_lane;
                    if (INSERT_PADDING) item_offset += item_offset >> LOG_SMEM_BANKS;
                    items[ITEM] = temp_storage[item_offset];
                }
            }
        }
    }


    /**
     * Transposes data items from <em>striped</em> arrangement to <em>blocked</em> arrangement.  Specialized for no timeslicing.
     */
    __device__ __forceinline__ void StripedToBlocked(
        T               items[ITEMS_PER_THREAD],   ///< [in-out] Items to exchange, converting between <em>striped</em> and <em>blocked</em> arrangements.
        Int2Type<false> time_slicing)
    {
        #pragma unroll
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
        {
            int item_offset = int(ITEM * BLOCK_THREADS) + linear_tid;
            if (INSERT_PADDING) item_offset += item_offset >> LOG_SMEM_BANKS;
            temp_storage[item_offset] = items[ITEM];
        }

        __syncthreads();

        // No timeslicing
        #pragma unroll
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
        {
            int item_offset = (linear_tid * ITEMS_PER_THREAD) + ITEM;
            if (INSERT_PADDING) item_offset += item_offset >> LOG_SMEM_BANKS;
            items[ITEM] = temp_storage[item_offset];
        }
    }


    /**
     * Transposes data items from <em>striped</em> arrangement to <em>blocked</em> arrangement.  Specialized for warp-timeslicing.
     */
    __device__ __forceinline__ void StripedToBlocked(
        T               items[ITEMS_PER_THREAD],   ///< [in-out] Items to exchange, converting between <em>striped</em> and <em>blocked</em> arrangements.
        Int2Type<true>  time_slicing)
    {
        // Warp time-slicing
        T temp_items[ITEMS_PER_THREAD];

        #pragma unroll
        for (int SLICE = 0; SLICE < TIME_SLICES; SLICE++)
        {
            const int SLICE_OFFSET  = SLICE * TIME_SLICED_ITEMS;
            const int SLICE_OOB     = SLICE_OFFSET + TIME_SLICED_ITEMS;

            __syncthreads();

            #pragma unroll
            for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
            {
                // Write a strip of items
                const int STRIP_OFFSET  = ITEM * BLOCK_THREADS;
                const int STRIP_OOB     = STRIP_OFFSET + BLOCK_THREADS;

                if ((SLICE_OFFSET < STRIP_OOB) && (SLICE_OOB > STRIP_OFFSET))
                {
                    int item_offset = STRIP_OFFSET + linear_tid - SLICE_OFFSET;
                    if ((item_offset >= 0) && (item_offset < TIME_SLICED_ITEMS))
                    {
                        if (INSERT_PADDING) item_offset += item_offset >> LOG_SMEM_BANKS;
                        temp_storage[item_offset] = items[ITEM];
                    }
                }
            }

            __syncthreads();

            if (warp_id == SLICE)
            {
                #pragma unroll
                for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
                {
                    int item_offset = (warp_lane * ITEMS_PER_THREAD) + ITEM;
                    if (INSERT_PADDING) item_offset += item_offset >> LOG_SMEM_BANKS;
                    temp_items[ITEM] = temp_storage[item_offset];
                }
            }
        }

        // Copy
        #pragma unroll
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
        {
            items[ITEM] = temp_items[ITEM];
        }
    }


    /**
     * Transposes data items from <em>warp-striped</em> arrangement to <em>blocked</em> arrangement.  Specialized for no timeslicing
     */
    __device__ __forceinline__ void WarpStripedToBlocked(
        T               items[ITEMS_PER_THREAD],   ///< [in-out] Items to exchange, converting between <em>warp-striped</em> and <em>blocked</em> arrangements.
        Int2Type<false> time_slicing)
    {
        #pragma unroll
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
        {
            int item_offset = warp_offset + (ITEM * WARP_TIME_SLICED_THREADS) + warp_lane;
            if (INSERT_PADDING) item_offset += item_offset >> LOG_SMEM_BANKS;
            temp_storage[item_offset] = items[ITEM];
        }

        #pragma unroll
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
        {
            int item_offset = warp_offset + ITEM + (warp_lane * ITEMS_PER_THREAD);
            if (INSERT_PADDING) item_offset += item_offset >> LOG_SMEM_BANKS;
            items[ITEM] = temp_storage[item_offset];
        }
    }


    /**
     * Transposes data items from <em>warp-striped</em> arrangement to <em>blocked</em> arrangement.  Specialized for warp-timeslicing
     */
    __device__ __forceinline__ void WarpStripedToBlocked(
        T               items[ITEMS_PER_THREAD],   ///< [in-out] Items to exchange, converting between <em>warp-striped</em> and <em>blocked</em> arrangements.
        Int2Type<true>  time_slicing)
    {
        #pragma unroll
        for (int SLICE = 0; SLICE < TIME_SLICES; ++SLICE)
        {
            __syncthreads();

            if (warp_id == SLICE)
            {
                #pragma unroll
                for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
                {
                    int item_offset = (ITEM * WARP_TIME_SLICED_THREADS) + warp_lane;
                    if (INSERT_PADDING) item_offset += item_offset >> LOG_SMEM_BANKS;
                    temp_storage[item_offset] = items[ITEM];
                }

                #pragma unroll
                for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
                {
                    int item_offset = ITEM + (warp_lane * ITEMS_PER_THREAD);
                    if (INSERT_PADDING) item_offset += item_offset >> LOG_SMEM_BANKS;
                    items[ITEM] = temp_storage[item_offset];
                }
            }
        }
    }


    /**
     * Exchanges data items annotated by rank into <em>blocked</em> arrangement.  Specialized for no timeslicing.
     */
    __device__ __forceinline__ void ScatterToBlocked(
        T               items[ITEMS_PER_THREAD],    ///< [in-out] Items to exchange
        int             ranks[ITEMS_PER_THREAD],    ///< [in] Corresponding scatter ranks
        Int2Type<false> time_slicing)
    {
        #pragma unroll
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
        {
            int item_offset = ranks[ITEM];
            if (INSERT_PADDING) item_offset = SHR_ADD(item_offset, LOG_SMEM_BANKS, item_offset);
            temp_storage[item_offset] = items[ITEM];
        }

        __syncthreads();

        #pragma unroll
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
        {
            int item_offset = (linear_tid * ITEMS_PER_THREAD) + ITEM;
            if (INSERT_PADDING) item_offset = SHR_ADD(item_offset, LOG_SMEM_BANKS, item_offset);
            items[ITEM] = temp_storage[item_offset];
        }
    }

    /**
     * Exchanges data items annotated by rank into <em>blocked</em> arrangement.  Specialized for warp-timeslicing.
     */
    __device__ __forceinline__ void ScatterToBlocked(
        T               items[ITEMS_PER_THREAD],    ///< [in-out] Items to exchange
        int             ranks[ITEMS_PER_THREAD],    ///< [in] Corresponding scatter ranks
        Int2Type<true>  time_slicing)
    {
        T temp_items[ITEMS_PER_THREAD];

        #pragma unroll
        for (int SLICE = 0; SLICE < TIME_SLICES; SLICE++)
        {
            __syncthreads();

            const int SLICE_OFFSET = TIME_SLICED_ITEMS * SLICE;

            #pragma unroll
            for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
            {
                int item_offset = ranks[ITEM] - SLICE_OFFSET;
                if ((item_offset >= 0) && (item_offset < WARP_TIME_SLICED_ITEMS))
                {
                    if (INSERT_PADDING) item_offset = SHR_ADD(item_offset, LOG_SMEM_BANKS, item_offset);
                    temp_storage[item_offset] = items[ITEM];
                }
            }

            __syncthreads();

            if (warp_id == SLICE)
            {
                #pragma unroll
                for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
                {
                    int item_offset = (warp_lane * ITEMS_PER_THREAD) + ITEM;
                    if (INSERT_PADDING) item_offset = SHR_ADD(item_offset, LOG_SMEM_BANKS, item_offset);
                    temp_items[ITEM] = temp_storage[item_offset];
                }
            }
        }

        // Copy
        #pragma unroll
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
        {
            items[ITEM] = temp_items[ITEM];
        }
    }


    /**
     * Exchanges data items annotated by rank into <em>striped</em> arrangement.  Specialized for no timeslicing.
     */
    __device__ __forceinline__ void ScatterToStriped(
        T               items[ITEMS_PER_THREAD],    ///< [in-out] Items to exchange
        int             ranks[ITEMS_PER_THREAD],    ///< [in] Corresponding scatter ranks
        Int2Type<false> time_slicing)
    {
        #pragma unroll
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
        {
            int item_offset = ranks[ITEM];
            if (INSERT_PADDING) item_offset = SHR_ADD(item_offset, LOG_SMEM_BANKS, item_offset);
            temp_storage[item_offset] = items[ITEM];
        }

        __syncthreads();

        #pragma unroll
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
        {
            int item_offset = int(ITEM * BLOCK_THREADS) + linear_tid;
            if (INSERT_PADDING) item_offset = SHR_ADD(item_offset, LOG_SMEM_BANKS, item_offset);
            items[ITEM] = temp_storage[item_offset];
        }
    }


    /**
     * Exchanges data items annotated by rank into <em>striped</em> arrangement.  Specialized for warp-timeslicing.
     */
    __device__ __forceinline__ void ScatterToStriped(
        T               items[ITEMS_PER_THREAD],    ///< [in-out] Items to exchange
        int             ranks[ITEMS_PER_THREAD],    ///< [in] Corresponding scatter ranks
        Int2Type<true> time_slicing)
    {
        T temp_items[ITEMS_PER_THREAD];

        #pragma unroll
        for (int SLICE = 0; SLICE < TIME_SLICES; SLICE++)
        {
            const int SLICE_OFFSET  = SLICE * TIME_SLICED_ITEMS;
            const int SLICE_OOB     = SLICE_OFFSET + TIME_SLICED_ITEMS;

            __syncthreads();

            #pragma unroll
            for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
            {
                int item_offset = ranks[ITEM] - SLICE_OFFSET;
                if ((item_offset >= 0) && (item_offset < WARP_TIME_SLICED_ITEMS))
                {
                    if (INSERT_PADDING) item_offset = SHR_ADD(item_offset, LOG_SMEM_BANKS, item_offset);
                    temp_storage[item_offset] = items[ITEM];
                }
            }

            __syncthreads();

            #pragma unroll
            for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
            {
                // Read a strip of items
                const int STRIP_OFFSET  = ITEM * BLOCK_THREADS;
                const int STRIP_OOB     = STRIP_OFFSET + BLOCK_THREADS;

                if ((SLICE_OFFSET < STRIP_OOB) && (SLICE_OOB > STRIP_OFFSET))
                {
                    int item_offset = STRIP_OFFSET + linear_tid - SLICE_OFFSET;
                    if ((item_offset >= 0) && (item_offset < TIME_SLICED_ITEMS))
                    {
                        if (INSERT_PADDING) item_offset += item_offset >> LOG_SMEM_BANKS;
                        temp_items[ITEM] = temp_storage[item_offset];
                    }
                }
            }
        }

        // Copy
        #pragma unroll
        for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
        {
            items[ITEM] = temp_items[ITEM];
        }
    }


public:

    /******************************************************************//**
     * \name Collective constructors
     *********************************************************************/
    //@{

    /**
     * \brief Collective constructor for 1D thread blocks using a private static allocation of shared memory as temporary storage.  Threads are identified using <tt>threadIdx.x</tt>.
     */
    __device__ __forceinline__ BlockExchange()
    :
        temp_storage(PrivateStorage()),
        linear_tid(threadIdx.x),
        warp_lane(linear_tid & (WARP_THREADS - 1)),
        warp_id(linear_tid >> LOG_WARP_THREADS),
        warp_offset(warp_id * WARP_TIME_SLICED_ITEMS)
    {}


    /**
     * \brief Collective constructor for 1D thread blocks using the specified memory allocation as temporary storage.  Threads are identified using <tt>threadIdx.x</tt>.
     */
    __device__ __forceinline__ BlockExchange(
        TempStorage &temp_storage)             ///< [in] Reference to memory allocation having layout type TempStorage
    :
        temp_storage(temp_storage.Alias()),
        linear_tid(threadIdx.x),
        warp_lane(linear_tid & (WARP_THREADS - 1)),
        warp_id(linear_tid >> LOG_WARP_THREADS),
        warp_offset(warp_id * WARP_TIME_SLICED_ITEMS)
    {}


    /**
     * \brief Collective constructor using a private static allocation of shared memory as temporary storage.  Each thread is identified using the supplied linear thread identifier
     */
    __device__ __forceinline__ BlockExchange(
        int linear_tid)                        ///< [in] A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    :
        temp_storage(PrivateStorage()),
        linear_tid(linear_tid),
        warp_lane(linear_tid & (WARP_THREADS - 1)),
        warp_id(linear_tid >> LOG_WARP_THREADS),
        warp_offset(warp_id * WARP_TIME_SLICED_ITEMS)
    {}


    /**
     * \brief Collective constructor using the specified memory allocation as temporary storage.  Each thread is identified using the supplied linear thread identifier.
     */
    __device__ __forceinline__ BlockExchange(
        TempStorage &temp_storage,              ///< [in] Reference to memory allocation having layout type TempStorage
        int         linear_tid)                 ///< [in] <b>[optional]</b> A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    :
        temp_storage(temp_storage.Alias()),
        linear_tid(linear_tid),
        warp_lane(linear_tid & (WARP_THREADS - 1)),
        warp_id(linear_tid >> LOG_WARP_THREADS),
        warp_offset(warp_id * WARP_TIME_SLICED_ITEMS)
    {}


    //@}  end member group
    /******************************************************************//**
     * \name Structured exchanges
     *********************************************************************/
    //@{

    /**
     * \brief Transposes data items from <em>striped</em> arrangement to <em>blocked</em> arrangement.
     *
     * \smemreuse
     *
     * The code snippet below illustrates the conversion from a "striped" to a "blocked" arrangement
     * of 512 integer items partitioned across 128 threads where each thread owns 4 items.
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(int *d_data, ...)
     * {
     *     // Specialize BlockExchange for 128 threads owning 4 integer items each
     *     typedef cub::BlockExchange<int, 128, 4> BlockExchange;
     *
     *     // Allocate shared memory for BlockExchange
     *     __shared__ typename BlockExchange::TempStorage temp_storage;
     *
     *     // Load a tile of ordered data into a striped arrangement across block threads
     *     int thread_data[4];
     *     cub::LoadStriped<LOAD_DEFAULT, 128>(threadIdx.x, d_data, thread_data);
     *
     *     // Collectively exchange data into a blocked arrangement across threads
     *     BlockExchange(temp_storage).StripedToBlocked(thread_data);
     *
     * \endcode
     * \par
     * Suppose the set of striped input \p thread_data across the block of threads is
     * <tt>{ [0,128,256,384], [1,129,257,385], ..., [127,255,383,511] }</tt> after loading from global memory.
     * The corresponding output \p thread_data in those threads will be
     * <tt>{ [0,1,2,3], [4,5,6,7], [8,9,10,11], ..., [508,509,510,511] }</tt>.
     *
     */
    __device__ __forceinline__ void StripedToBlocked(
        T                items[ITEMS_PER_THREAD])   ///< [in-out] Items to exchange, converting between <em>striped</em> and <em>blocked</em> arrangements.
    {
        StripedToBlocked(items, Int2Type<WARP_TIME_SLICING>());
    }

    /**
     * \brief Transposes data items from <em>blocked</em> arrangement to <em>striped</em> arrangement.
     *
     * \smemreuse
     *
     * The code snippet below illustrates the conversion from a "blocked" to a "striped" arrangement
     * of 512 integer items partitioned across 128 threads where each thread owns 4 items.
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(int *d_data, ...)
     * {
     *     // Specialize BlockExchange for 128 threads owning 4 integer items each
     *     typedef cub::BlockExchange<int, 128, 4> BlockExchange;
     *
     *     // Allocate shared memory for BlockExchange
     *     __shared__ typename BlockExchange::TempStorage temp_storage;
     *
     *     // Obtain a segment of consecutive items that are blocked across threads
     *     int thread_data[4];
     *     ...
     *
     *     // Collectively exchange data into a striped arrangement across threads
     *     BlockExchange(temp_storage).BlockedToStriped(thread_data);
     *
     *     // Store data striped across block threads into an ordered tile
     *     cub::StoreStriped<STORE_DEFAULT, 128>(threadIdx.x, d_data, thread_data);
     *
     * \endcode
     * \par
     * Suppose the set of blocked input \p thread_data across the block of threads is
     * <tt>{ [0,1,2,3], [4,5,6,7], [8,9,10,11], ..., [508,509,510,511] }</tt>.
     * The corresponding output \p thread_data in those threads will be
     * <tt>{ [0,128,256,384], [1,129,257,385], ..., [127,255,383,511] }</tt> in
     * preparation for storing to global memory.
     *
     */
    __device__ __forceinline__ void BlockedToStriped(
        T               items[ITEMS_PER_THREAD])    ///< [in-out] Items to exchange, converting between <em>blocked</em> and <em>striped</em> arrangements.
    {
        BlockedToStriped(items, Int2Type<WARP_TIME_SLICING>());
    }


    /**
     * \brief Transposes data items from <em>warp-striped</em> arrangement to <em>blocked</em> arrangement.
     *
     * \smemreuse
     *
     * The code snippet below illustrates the conversion from a "warp-striped" to a "blocked" arrangement
     * of 512 integer items partitioned across 128 threads where each thread owns 4 items.
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(int *d_data, ...)
     * {
     *     // Specialize BlockExchange for 128 threads owning 4 integer items each
     *     typedef cub::BlockExchange<int, 128, 4> BlockExchange;
     *
     *     // Allocate shared memory for BlockExchange
     *     __shared__ typename BlockExchange::TempStorage temp_storage;
     *
     *     // Load a tile of ordered data into a warp-striped arrangement across warp threads
     *     int thread_data[4];
     *     cub::LoadSWarptriped<LOAD_DEFAULT>(threadIdx.x, d_data, thread_data);
     *
     *     // Collectively exchange data into a blocked arrangement across threads
     *     BlockExchange(temp_storage).WarpStripedToBlocked(thread_data);
     *
     * \endcode
     * \par
     * Suppose the set of warp-striped input \p thread_data across the block of threads is
     * <tt>{ [0,32,64,96], [1,33,65,97], [2,34,66,98], ..., [415,447,479,511] }</tt>
     * after loading from global memory.  (The first 128 items are striped across
     * the first warp of 32 threads, the second 128 items are striped across the second warp, etc.)
     * The corresponding output \p thread_data in those threads will be
     * <tt>{ [0,1,2,3], [4,5,6,7], [8,9,10,11], ..., [508,509,510,511] }</tt>.
     *
     */
    __device__ __forceinline__ void WarpStripedToBlocked(
        T                items[ITEMS_PER_THREAD])   ///< [in-out] Items to exchange, converting between <em>warp-striped</em> and <em>blocked</em> arrangements.
    {
        WarpStripedToBlocked(items, Int2Type<WARP_TIME_SLICING>());
    }

    /**
     * \brief Transposes data items from <em>blocked</em> arrangement to <em>warp-striped</em> arrangement.
     *
     * \smemreuse
     *
     * The code snippet below illustrates the conversion from a "blocked" to a "warp-striped" arrangement
     * of 512 integer items partitioned across 128 threads where each thread owns 4 items.
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(int *d_data, ...)
     * {
     *     // Specialize BlockExchange for 128 threads owning 4 integer items each
     *     typedef cub::BlockExchange<int, 128, 4> BlockExchange;
     *
     *     // Allocate shared memory for BlockExchange
     *     __shared__ typename BlockExchange::TempStorage temp_storage;
     *
     *     // Obtain a segment of consecutive items that are blocked across threads
     *     int thread_data[4];
     *     ...
     *
     *     // Collectively exchange data into a warp-striped arrangement across threads
     *     BlockExchange(temp_storage).BlockedToWarpStriped(thread_data);
     *
     *     // Store data striped across warp threads into an ordered tile
     *     cub::StoreStriped<STORE_DEFAULT, 128>(threadIdx.x, d_data, thread_data);
     *
     * \endcode
     * \par
     * Suppose the set of blocked input \p thread_data across the block of threads is
     * <tt>{ [0,1,2,3], [4,5,6,7], [8,9,10,11], ..., [508,509,510,511] }</tt>.
     * The corresponding output \p thread_data in those threads will be
     * <tt>{ [0,32,64,96], [1,33,65,97], [2,34,66,98], ..., [415,447,479,511] }</tt>
     * in preparation for storing to global memory. (The first 128 items are striped across
     * the first warp of 32 threads, the second 128 items are striped across the second warp, etc.)
     *
     */
    __device__ __forceinline__ void BlockedToWarpStriped(
        T                items[ITEMS_PER_THREAD])   ///< [in-out] Items to exchange, converting between <em>blocked</em> and <em>warp-striped</em> arrangements.
    {
        BlockedToWarpStriped(items, Int2Type<WARP_TIME_SLICING>());
    }


    //@}  end member group
    /******************************************************************//**
     * \name Scatter exchanges
     *********************************************************************/
    //@{


    /**
     * \brief Exchanges data items annotated by rank into <em>blocked</em> arrangement.
     *
     * \smemreuse
     */
    __device__ __forceinline__ void ScatterToBlocked(
        T               items[ITEMS_PER_THREAD],    ///< [in-out] Items to exchange
        int             ranks[ITEMS_PER_THREAD])    ///< [in] Corresponding scatter ranks
    {
        ScatterToBlocked(items, ranks, Int2Type<WARP_TIME_SLICING>());
    }


    /**
     * \brief Exchanges data items annotated by rank into <em>striped</em> arrangement.
     *
     * \smemreuse
     */
    __device__ __forceinline__ void ScatterToStriped(
        T               items[ITEMS_PER_THREAD],    ///< [in-out] Items to exchange
        int             ranks[ITEMS_PER_THREAD])    ///< [in] Corresponding scatter ranks
    {
        ScatterToStriped(items, ranks, Int2Type<WARP_TIME_SLICING>());
    }

    //@}  end member group


};

}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)

