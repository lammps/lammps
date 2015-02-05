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
 * Operations for writing linear segments of data from the CUDA thread block
 */

#pragma once

#include <iterator>

#include "../util_namespace.cuh"
#include "../util_macro.cuh"
#include "../util_type.cuh"
#include "../util_vector.cuh"
#include "../thread/thread_store.cuh"
#include "block_exchange.cuh"

/// Optional outer namespace(s)
CUB_NS_PREFIX

/// CUB namespace
namespace cub {

/**
 * \addtogroup IoModule
 * @{
 */


/******************************************************************//**
 * \name Blocked I/O
 *********************************************************************/
//@{

/**
 * \brief Store a blocked arrangement of items across a thread block into a linear segment of items using the specified cache modifier.
 *
 * \blocked
 *
 * \tparam MODIFIER             cub::PtxStoreModifier cache modifier.
 * \tparam T                    <b>[inferred]</b> The data type to store.
 * \tparam ITEMS_PER_THREAD     <b>[inferred]</b> The number of consecutive items partitioned onto each thread.
 * \tparam OutputIteratorRA     <b>[inferred]</b> The random-access iterator type for output (may be a simple pointer type).
 */
template <
    PtxStoreModifier    MODIFIER,
    typename            T,
    int                 ITEMS_PER_THREAD,
    typename            OutputIteratorRA>
__device__ __forceinline__ void StoreBlocked(
    int                 linear_tid,                 ///< [in] A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    OutputIteratorRA    block_itr,                  ///< [in] The thread block's base output iterator for storing to
    T                   (&items)[ITEMS_PER_THREAD]) ///< [in] Data to store
{
    // Store directly in thread-blocked order
    #pragma unroll
    for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
    {
        ThreadStore<MODIFIER>(block_itr + (linear_tid * ITEMS_PER_THREAD) + ITEM, items[ITEM]);
    }
}


/**
 * \brief Store a blocked arrangement of items across a thread block into a linear segment of items using the specified cache modifier, guarded by range
 *
 * \blocked
 *
 * \tparam MODIFIER             cub::PtxStoreModifier cache modifier.
 * \tparam T                    <b>[inferred]</b> The data type to store.
 * \tparam ITEMS_PER_THREAD     <b>[inferred]</b> The number of consecutive items partitioned onto each thread.
 * \tparam OutputIteratorRA     <b>[inferred]</b> The random-access iterator type for output (may be a simple pointer type).
 */
template <
    PtxStoreModifier    MODIFIER,
    typename            T,
    int                 ITEMS_PER_THREAD,
    typename            OutputIteratorRA>
__device__ __forceinline__ void StoreBlocked(
    int                 linear_tid,                 ///< [in] A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    OutputIteratorRA    block_itr,                  ///< [in] The thread block's base output iterator for storing to
    T                   (&items)[ITEMS_PER_THREAD], ///< [in] Data to store
    int                 valid_items)                ///< [in] Number of valid items to write
{
    // Store directly in thread-blocked order
    #pragma unroll
    for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
    {
        if (ITEM + (linear_tid * ITEMS_PER_THREAD) < valid_items)
        {
            ThreadStore<MODIFIER>(block_itr + (linear_tid * ITEMS_PER_THREAD) + ITEM, items[ITEM]);
        }
    }
}



//@}  end member group
/******************************************************************//**
 * \name Striped I/O
 *********************************************************************/
//@{


/**
 * \brief Store a striped arrangement of data across the thread block into a linear segment of items using the specified cache modifier.
 *
 * \striped
 *
 * \tparam MODIFIER             cub::PtxStoreModifier cache modifier.
 * \tparam BLOCK_THREADS        The thread block size in threads
 * \tparam T                    <b>[inferred]</b> The data type to store.
 * \tparam ITEMS_PER_THREAD     <b>[inferred]</b> The number of consecutive items partitioned onto each thread.
 * \tparam OutputIteratorRA     <b>[inferred]</b> The random-access iterator type for output (may be a simple pointer type).
 */
template <
    PtxStoreModifier    MODIFIER,
    int                 BLOCK_THREADS,
    typename            T,
    int                 ITEMS_PER_THREAD,
    typename            OutputIteratorRA>
__device__ __forceinline__ void StoreStriped(
    int                 linear_tid,                 ///< [in] A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    OutputIteratorRA    block_itr,                  ///< [in] The thread block's base output iterator for storing to
    T                   (&items)[ITEMS_PER_THREAD]) ///< [in] Data to store
{
    // Store directly in striped order
    #pragma unroll
    for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
    {
        ThreadStore<MODIFIER>(block_itr + (ITEM * BLOCK_THREADS) + linear_tid, items[ITEM]);
    }
}


/**
 * \brief Store a striped arrangement of data across the thread block into a linear segment of items using the specified cache modifier, guarded by range
 *
 * \striped
 *
 * \tparam MODIFIER             cub::PtxStoreModifier cache modifier.
 * \tparam BLOCK_THREADS        The thread block size in threads
 * \tparam T                    <b>[inferred]</b> The data type to store.
 * \tparam ITEMS_PER_THREAD     <b>[inferred]</b> The number of consecutive items partitioned onto each thread.
 * \tparam OutputIteratorRA     <b>[inferred]</b> The random-access iterator type for output (may be a simple pointer type).
 */
template <
    PtxStoreModifier    MODIFIER,
    int                 BLOCK_THREADS,
    typename            T,
    int                 ITEMS_PER_THREAD,
    typename            OutputIteratorRA>
__device__ __forceinline__ void StoreStriped(
    int                 linear_tid,                 ///< [in] A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    OutputIteratorRA    block_itr,                  ///< [in] The thread block's base output iterator for storing to
    T                   (&items)[ITEMS_PER_THREAD], ///< [in] Data to store
    int                 valid_items)                ///< [in] Number of valid items to write
{
    // Store directly in striped order
    #pragma unroll
    for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
    {
        if ((ITEM * BLOCK_THREADS) + linear_tid < valid_items)
        {
            ThreadStore<MODIFIER>(block_itr + (ITEM * BLOCK_THREADS) + linear_tid, items[ITEM]);
        }
    }
}



//@}  end member group
/******************************************************************//**
 * \name Warp-striped I/O
 *********************************************************************/
//@{


/**
 * \brief Store a warp-striped arrangement of data across the thread block into a linear segment of items using the specified cache modifier.
 *
 * \warpstriped
 *
 * \par Usage Considerations
 * The number of threads in the thread block must be a multiple of the architecture's warp size.
 *
 * \tparam MODIFIER             cub::PtxStoreModifier cache modifier.
 * \tparam T                    <b>[inferred]</b> The data type to store.
 * \tparam ITEMS_PER_THREAD     <b>[inferred]</b> The number of consecutive items partitioned onto each thread.
 * \tparam OutputIteratorRA     <b>[inferred]</b> The random-access iterator type for output (may be a simple pointer type).
 */
template <
    PtxStoreModifier    MODIFIER,
    typename            T,
    int                 ITEMS_PER_THREAD,
    typename            OutputIteratorRA>
__device__ __forceinline__ void StoreWarpStriped(
    int                 linear_tid,                 ///< [in] A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    OutputIteratorRA    block_itr,                  ///< [in] The thread block's base output iterator for storing to
    T                   (&items)[ITEMS_PER_THREAD]) ///< [out] Data to load
{
    int tid         = linear_tid & (PtxArchProps::WARP_THREADS - 1);
    int wid         = linear_tid >> PtxArchProps::LOG_WARP_THREADS;
    int warp_offset = wid * PtxArchProps::WARP_THREADS * ITEMS_PER_THREAD;

    // Store directly in warp-striped order
    #pragma unroll
    for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
    {
        ThreadStore<MODIFIER>(block_itr + warp_offset + tid + (ITEM * PtxArchProps::WARP_THREADS), items[ITEM]);
    }
}


/**
 * \brief Store a warp-striped arrangement of data across the thread block into a linear segment of items using the specified cache modifier, guarded by range
 *
 * \warpstriped
 *
 * \par Usage Considerations
 * The number of threads in the thread block must be a multiple of the architecture's warp size.
 *
 * \tparam MODIFIER             cub::PtxStoreModifier cache modifier.
 * \tparam T                    <b>[inferred]</b> The data type to store.
 * \tparam ITEMS_PER_THREAD     <b>[inferred]</b> The number of consecutive items partitioned onto each thread.
 * \tparam OutputIteratorRA     <b>[inferred]</b> The random-access iterator type for output (may be a simple pointer type).
 */
template <
    PtxStoreModifier    MODIFIER,
    typename            T,
    int                 ITEMS_PER_THREAD,
    typename            OutputIteratorRA>
__device__ __forceinline__ void StoreWarpStriped(
    int                 linear_tid,                 ///< [in] A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    OutputIteratorRA    block_itr,                  ///< [in] The thread block's base output iterator for storing to
    T                   (&items)[ITEMS_PER_THREAD], ///< [in] Data to store
    int                 valid_items)                ///< [in] Number of valid items to write
{
    int tid         = linear_tid & (PtxArchProps::WARP_THREADS - 1);
    int wid         = linear_tid >> PtxArchProps::LOG_WARP_THREADS;
    int warp_offset = wid * PtxArchProps::WARP_THREADS * ITEMS_PER_THREAD;

    // Store directly in warp-striped order
    #pragma unroll
    for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
    {
        if (warp_offset + tid + (ITEM * PtxArchProps::WARP_THREADS) < valid_items)
        {
            ThreadStore<MODIFIER>(block_itr + warp_offset + tid + (ITEM * PtxArchProps::WARP_THREADS), items[ITEM]);
        }
    }
}



//@}  end member group
/******************************************************************//**
 * \name Blocked, vectorized I/O
 *********************************************************************/
//@{

/**
 * \brief Store a blocked arrangement of items across a thread block into a linear segment of items using the specified cache modifier.
 *
 * \blocked
 *
 * The output offset (\p block_ptr + \p block_offset) must be quad-item aligned,
 * which is the default starting offset returned by \p cudaMalloc()
 *
 * \par
 * The following conditions will prevent vectorization and storing will fall back to cub::BLOCK_STORE_DIRECT:
 *   - \p ITEMS_PER_THREAD is odd
 *   - The data type \p T is not a built-in primitive or CUDA vector type (e.g., \p short, \p int2, \p double, \p float2, etc.)
 *
 * \tparam MODIFIER             cub::PtxStoreModifier cache modifier.
 * \tparam T                    <b>[inferred]</b> The data type to store.
 * \tparam ITEMS_PER_THREAD     <b>[inferred]</b> The number of consecutive items partitioned onto each thread.
 *
 */
template <
    PtxStoreModifier    MODIFIER,
    typename            T,
    int                 ITEMS_PER_THREAD>
__device__ __forceinline__ void StoreBlockedVectorized(
    int                 linear_tid,                 ///< [in] A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    T                   *block_ptr,                 ///< [in] Input pointer for storing from
    T                   (&items)[ITEMS_PER_THREAD]) ///< [in] Data to store
{
    enum
    {
        // Maximum CUDA vector size is 4 elements
        MAX_VEC_SIZE = CUB_MIN(4, ITEMS_PER_THREAD),

        // Vector size must be a power of two and an even divisor of the items per thread
        VEC_SIZE = ((((MAX_VEC_SIZE - 1) & MAX_VEC_SIZE) == 0) && ((ITEMS_PER_THREAD % MAX_VEC_SIZE) == 0)) ?
            MAX_VEC_SIZE :
            1,

        VECTORS_PER_THREAD = ITEMS_PER_THREAD / VEC_SIZE,
    };

    // Vector type
    typedef typename VectorHelper<T, VEC_SIZE>::Type Vector;

    // Alias global pointer
    Vector *block_ptr_vectors = reinterpret_cast<Vector *>(block_ptr);

    // Alias pointers (use "raw" array here which should get optimized away to prevent conservative PTXAS lmem spilling)
    Vector raw_vector[VECTORS_PER_THREAD];
    T *raw_items = reinterpret_cast<T*>(raw_vector);

    // Copy
    #pragma unroll
    for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
    {
        raw_items[ITEM] = items[ITEM];
    }

    // Direct-store using vector types
    StoreBlocked<MODIFIER>(linear_tid, block_ptr_vectors, raw_vector);
}


//@}  end member group


/** @} */       // end group IoModule


//-----------------------------------------------------------------------------
// Generic BlockStore abstraction
//-----------------------------------------------------------------------------

/**
 * \brief cub::BlockStoreAlgorithm enumerates alternative algorithms for cub::BlockStore to write a blocked arrangement of items across a CUDA thread block to a linear segment of memory.
 */
enum BlockStoreAlgorithm
{
    /**
     * \par Overview
     *
     * A [<em>blocked arrangement</em>](index.html#sec5sec4) of data is written
     * directly to memory.  The thread block writes items in a parallel "raking" fashion:
     * thread<sub><em>i</em></sub> writes the <em>i</em><sup>th</sup> segment of consecutive elements.
     *
     * \par Performance Considerations
     * - The utilization of memory transactions (coalescing) decreases as the
     *   access stride between threads increases (i.e., the number items per thread).
     */
    BLOCK_STORE_DIRECT,

    /**
     * \par Overview
     *
     * A [<em>blocked arrangement</em>](index.html#sec5sec4) of data is written directly
     * to memory using CUDA's built-in vectorized stores as a coalescing optimization.
     * The thread block writes items in a parallel "raking" fashion: thread<sub><em>i</em></sub> uses vector stores to
     * write the <em>i</em><sup>th</sup> segment of consecutive elements.
     *
     * For example, <tt>st.global.v4.s32</tt> instructions will be generated when \p T = \p int and \p ITEMS_PER_THREAD > 4.
     *
     * \par Performance Considerations
     * - The utilization of memory transactions (coalescing) remains high until the the
     *   access stride between threads (i.e., the number items per thread) exceeds the
     *   maximum vector store width (typically 4 items or 64B, whichever is lower).
     * - The following conditions will prevent vectorization and writing will fall back to cub::BLOCK_STORE_DIRECT:
     *   - \p ITEMS_PER_THREAD is odd
     *   - The \p OutputIteratorRA is not a simple pointer type
     *   - The block output offset is not quadword-aligned
     *   - The data type \p T is not a built-in primitive or CUDA vector type (e.g., \p short, \p int2, \p double, \p float2, etc.)
     */
    BLOCK_STORE_VECTORIZE,

    /**
     * \par Overview
     * A [<em>blocked arrangement</em>](index.html#sec5sec4) is locally
     * transposed into a [<em>striped arrangement</em>](index.html#sec5sec4)
     * which is then written to memory.  More specifically, cub::BlockExchange
     * used to locally reorder the items into a
     * [<em>striped arrangement</em>](index.html#sec5sec4), after which the
     * thread block writes items in a parallel "strip-mining" fashion: consecutive
     * items owned by thread<sub><em>i</em></sub> are written to memory with
     * stride \p BLOCK_THREADS between them.
     *
     * \par Performance Considerations
     * - The utilization of memory transactions (coalescing) remains high regardless
     *   of items written per thread.
     * - The local reordering incurs slightly longer latencies and throughput than the
     *   direct cub::BLOCK_STORE_DIRECT and cub::BLOCK_STORE_VECTORIZE alternatives.
     */
    BLOCK_STORE_TRANSPOSE,

    /**
     * \par Overview
     * A [<em>blocked arrangement</em>](index.html#sec5sec4) is locally
     * transposed into a [<em>warp-striped arrangement</em>](index.html#sec5sec4)
     * which is then written to memory.  More specifically, cub::BlockExchange used
     * to locally reorder the items into a
     * [<em>warp-striped arrangement</em>](index.html#sec5sec4), after which
     * each warp writes its own contiguous segment in a parallel "strip-mining" fashion:
     * consecutive items owned by lane<sub><em>i</em></sub> are written to memory
     * with stride \p WARP_THREADS between them.
     *
     * \par Performance Considerations
     * - The utilization of memory transactions (coalescing) remains high regardless
     *   of items written per thread.
     * - The local reordering incurs slightly longer latencies and throughput than the
     *   direct cub::BLOCK_STORE_DIRECT and cub::BLOCK_STORE_VECTORIZE alternatives.
     */
    BLOCK_STORE_WARP_TRANSPOSE,
};



/**
 * \addtogroup BlockModule
 * @{
 */


/**
 * \brief The BlockStore class provides [<em>collective</em>](index.html#sec0) data movement methods for writing a [<em>blocked arrangement</em>](index.html#sec5sec4) of items partitioned across a CUDA thread block to a linear segment of memory.  ![](block_store_logo.png)
 *
 * \par Overview
 * The BlockStore class provides a single data movement abstraction that can be specialized
 * to implement different cub::BlockStoreAlgorithm strategies.  This facilitates different
 * performance policies for different architectures, data types, granularity sizes, etc.
 *
 * \par Optionally, BlockStore can be specialized by different data movement strategies:
 *   -# <b>cub::BLOCK_STORE_DIRECT</b>.  A [<em>blocked arrangement</em>](index.html#sec5sec4) of data is written
 *      directly to memory. [More...](\ref cub::BlockStoreAlgorithm)
 *   -# <b>cub::BLOCK_STORE_VECTORIZE</b>.  A [<em>blocked arrangement</em>](index.html#sec5sec4)
 *      of data is written directly to memory using CUDA's built-in vectorized stores as a
 *      coalescing optimization.  [More...](\ref cub::BlockStoreAlgorithm)
 *   -# <b>cub::BLOCK_STORE_TRANSPOSE</b>.  A [<em>blocked arrangement</em>](index.html#sec5sec4)
 *      is locally transposed into a [<em>striped arrangement</em>](index.html#sec5sec4) which is
 *      then written to memory.  [More...](\ref cub::BlockStoreAlgorithm)
 *   -# <b>cub::BLOCK_STORE_WARP_TRANSPOSE</b>.  A [<em>blocked arrangement</em>](index.html#sec5sec4)
 *      is locally transposed into a [<em>warp-striped arrangement</em>](index.html#sec5sec4) which is
 *      then written to memory.  [More...](\ref cub::BlockStoreAlgorithm)
 *
 * \tparam OutputIteratorRA     The input iterator type (may be a simple pointer type).
 * \tparam BLOCK_THREADS        The thread block size in threads.
 * \tparam ITEMS_PER_THREAD     The number of consecutive items partitioned onto each thread.
 * \tparam ALGORITHM            <b>[optional]</b> cub::BlockStoreAlgorithm tuning policy enumeration.  default: cub::BLOCK_STORE_DIRECT.
 * \tparam MODIFIER             <b>[optional]</b> cub::PtxStoreModifier cache modifier.  default: cub::STORE_DEFAULT.
 * \tparam WARP_TIME_SLICING    <b>[optional]</b> For transposition-based cub::BlockStoreAlgorithm parameterizations that utilize shared memory: When \p true, only use enough shared memory for a single warp's worth of data, time-slicing the block-wide exchange over multiple synchronized rounds (default: false)
 *
 * \par A Simple Example
 * \blockcollective{BlockStore}
 * \par
 * The code snippet below illustrates the storing of a "blocked" arrangement
 * of 512 integers across 128 threads (where each thread owns 4 consecutive items)
 * into a linear segment of memory.  The store is specialized for \p BLOCK_STORE_WARP_TRANSPOSE,
 * meaning items are locally reordered among threads so that memory references will be
 * efficiently coalesced using a warp-striped access pattern.
 * \par
 * \code
 * #include <cub/cub.cuh>
 *
 * __global__ void ExampleKernel(int *d_data, ...)
 * {
 *     // Specialize BlockStore for 128 threads owning 4 integer items each
 *     typedef cub::BlockStore<int*, 128, 4, BLOCK_STORE_WARP_TRANSPOSE> BlockStore;
 *
 *     // Allocate shared memory for BlockStore
 *     __shared__ typename BlockStore::TempStorage temp_storage;
 *
 *     // Obtain a segment of consecutive items that are blocked across threads
 *     int thread_data[4];
 *     ...
 *
 *     // Store items to linear memory
 *     int thread_data[4];
 *     BlockStore(temp_storage).Store(d_data, thread_data);
 *
 * \endcode
 * \par
 * Suppose the set of \p thread_data across the block of threads is
 * <tt>{ [0,1,2,3], [4,5,6,7], ..., [508,509,510,511] }</tt>.
 * The output \p d_data will be <tt>0, 1, 2, 3, 4, 5, ...</tt>.
 *
 */
template <
    typename                OutputIteratorRA,
    int                     BLOCK_THREADS,
    int                     ITEMS_PER_THREAD,
    BlockStoreAlgorithm     ALGORITHM           = BLOCK_STORE_DIRECT,
    PtxStoreModifier        MODIFIER            = STORE_DEFAULT,
    bool                    WARP_TIME_SLICING   = false>
class BlockStore
{
private:
    /******************************************************************************
     * Constants and typed definitions
     ******************************************************************************/

    // Data type of input iterator
    typedef typename std::iterator_traits<OutputIteratorRA>::value_type T;


    /******************************************************************************
     * Algorithmic variants
     ******************************************************************************/

    /// Store helper
    template <BlockStoreAlgorithm _POLICY, int DUMMY = 0>
    struct StoreInternal;


    /**
     * BLOCK_STORE_DIRECT specialization of store helper
     */
    template <int DUMMY>
    struct StoreInternal<BLOCK_STORE_DIRECT, DUMMY>
    {
        /// Shared memory storage layout type
        typedef NullType TempStorage;

        /// Linear thread-id
        int linear_tid;

        /// Constructor
        __device__ __forceinline__ StoreInternal(
            TempStorage &temp_storage,
            int linear_tid)
        :
            linear_tid(linear_tid)
        {}

        /// Store items into a linear segment of memory
        __device__ __forceinline__ void Store(
            OutputIteratorRA    block_itr,                  ///< [in] The thread block's base output iterator for storing to
            T                   (&items)[ITEMS_PER_THREAD]) ///< [in] Data to store
        {
            StoreBlocked<MODIFIER>(linear_tid, block_itr, items);
        }

        /// Store items into a linear segment of memory, guarded by range
        __device__ __forceinline__ void Store(
            OutputIteratorRA    block_itr,                  ///< [in] The thread block's base output iterator for storing to
            T                   (&items)[ITEMS_PER_THREAD], ///< [in] Data to store
            int                 valid_items)                ///< [in] Number of valid items to write
        {
            StoreBlocked<MODIFIER>(linear_tid, block_itr, items, valid_items);
        }
    };


    /**
     * BLOCK_STORE_VECTORIZE specialization of store helper
     */
    template <int DUMMY>
    struct StoreInternal<BLOCK_STORE_VECTORIZE, DUMMY>
    {
        /// Shared memory storage layout type
        typedef NullType TempStorage;

        /// Linear thread-id
        int linear_tid;

        /// Constructor
        __device__ __forceinline__ StoreInternal(
            TempStorage &temp_storage,
            int linear_tid)
        :
            linear_tid(linear_tid)
        {}

        /// Store items into a linear segment of memory, specialized for native pointer types (attempts vectorization)
        __device__ __forceinline__ void Store(
            T                   *block_ptr,                 ///< [in] The thread block's base output iterator for storing to
            T                   (&items)[ITEMS_PER_THREAD]) ///< [in] Data to store
        {
            StoreBlockedVectorized<MODIFIER>(linear_tid, block_ptr, items);
        }

        /// Store items into a linear segment of memory, specialized for opaque input iterators (skips vectorization)
        template <typename _OutputIteratorRA>
        __device__ __forceinline__ void Store(
            _OutputIteratorRA   block_itr,                  ///< [in] The thread block's base output iterator for storing to
            T                   (&items)[ITEMS_PER_THREAD]) ///< [in] Data to store
        {
            StoreBlocked<MODIFIER>(linear_tid, block_itr, items);
        }

        /// Store items into a linear segment of memory, guarded by range
        __device__ __forceinline__ void Store(
            OutputIteratorRA    block_itr,                  ///< [in] The thread block's base output iterator for storing to
            T                   (&items)[ITEMS_PER_THREAD], ///< [in] Data to store
            int                 valid_items)                ///< [in] Number of valid items to write
        {
            StoreBlocked<MODIFIER>(linear_tid, block_itr, items, valid_items);
        }
    };


    /**
     * BLOCK_STORE_TRANSPOSE specialization of store helper
     */
    template <int DUMMY>
    struct StoreInternal<BLOCK_STORE_TRANSPOSE, DUMMY>
    {
        // BlockExchange utility type for keys
        typedef BlockExchange<T, BLOCK_THREADS, ITEMS_PER_THREAD, WARP_TIME_SLICING> BlockExchange;

        /// Shared memory storage layout type
        typedef typename BlockExchange::TempStorage _TempStorage;

        /// Alias wrapper allowing storage to be unioned
        struct TempStorage : Uninitialized<_TempStorage> {};

        /// Thread reference to shared storage
        _TempStorage &temp_storage;

        /// Linear thread-id
        int linear_tid;

        /// Constructor
        __device__ __forceinline__ StoreInternal(
            TempStorage &temp_storage,
            int linear_tid)
        :
            temp_storage(temp_storage.Alias()),
            linear_tid(linear_tid)
        {}

        /// Store items into a linear segment of memory
        __device__ __forceinline__ void Store(
            OutputIteratorRA    block_itr,                  ///< [in] The thread block's base output iterator for storing to
            T                   (&items)[ITEMS_PER_THREAD]) ///< [in] Data to store
        {
            BlockExchange(temp_storage).BlockedToStriped(items);
            StoreStriped<MODIFIER, BLOCK_THREADS>(linear_tid, block_itr, items);
        }

        /// Store items into a linear segment of memory, guarded by range
        __device__ __forceinline__ void Store(
            OutputIteratorRA    block_itr,                  ///< [in] The thread block's base output iterator for storing to
            T                   (&items)[ITEMS_PER_THREAD], ///< [in] Data to store
            int                 valid_items)                ///< [in] Number of valid items to write
        {
            BlockExchange(temp_storage).BlockedToStriped(items);
            StoreStriped<MODIFIER, BLOCK_THREADS>(linear_tid, block_itr, items, valid_items);
        }
    };


    /**
     * BLOCK_STORE_WARP_TRANSPOSE specialization of store helper
     */
    template <int DUMMY>
    struct StoreInternal<BLOCK_STORE_WARP_TRANSPOSE, DUMMY>
    {
        enum
        {
            WARP_THREADS = PtxArchProps::WARP_THREADS
        };

        // Assert BLOCK_THREADS must be a multiple of WARP_THREADS
        CUB_STATIC_ASSERT((BLOCK_THREADS % WARP_THREADS == 0), "BLOCK_THREADS must be a multiple of WARP_THREADS");

        // BlockExchange utility type for keys
        typedef BlockExchange<T, BLOCK_THREADS, ITEMS_PER_THREAD, WARP_TIME_SLICING> BlockExchange;

        /// Shared memory storage layout type
        typedef typename BlockExchange::TempStorage _TempStorage;

        /// Alias wrapper allowing storage to be unioned
        struct TempStorage : Uninitialized<_TempStorage> {};

        /// Thread reference to shared storage
        _TempStorage &temp_storage;

        /// Linear thread-id
        int linear_tid;

        /// Constructor
        __device__ __forceinline__ StoreInternal(
            TempStorage &temp_storage,
            int linear_tid)
        :
            temp_storage(temp_storage.Alias()),
            linear_tid(linear_tid)
        {}

        /// Store items into a linear segment of memory
        __device__ __forceinline__ void Store(
            OutputIteratorRA    block_itr,                  ///< [in] The thread block's base output iterator for storing to
            T                   (&items)[ITEMS_PER_THREAD]) ///< [in] Data to store
        {
            BlockExchange(temp_storage).BlockedToWarpStriped(items);
            StoreWarpStriped<MODIFIER>(linear_tid, block_itr, items);
        }

        /// Store items into a linear segment of memory, guarded by range
        __device__ __forceinline__ void Store(
            OutputIteratorRA    block_itr,                  ///< [in] The thread block's base output iterator for storing to
            T                   (&items)[ITEMS_PER_THREAD], ///< [in] Data to store
            int                 valid_items)                ///< [in] Number of valid items to write
        {
            BlockExchange(temp_storage).BlockedToWarpStriped(items);
            StoreWarpStriped<MODIFIER>(linear_tid, block_itr, items, valid_items);
        }
    };

    /******************************************************************************
     * Type definitions
     ******************************************************************************/

    /// Internal load implementation to use
    typedef StoreInternal<ALGORITHM> InternalStore;


    /// Shared memory storage layout type
    typedef typename InternalStore::TempStorage _TempStorage;


    /******************************************************************************
     * Utility methods
     ******************************************************************************/

    /// Internal storage allocator
    __device__ __forceinline__ _TempStorage& PrivateStorage()
    {
        __shared__ _TempStorage private_storage;
        return private_storage;
    }


    /******************************************************************************
     * Thread fields
     ******************************************************************************/

    /// Thread reference to shared storage
    _TempStorage &temp_storage;

    /// Linear thread-id
    int linear_tid;

public:


    /// \smemstorage{BlockStore}
    struct TempStorage : Uninitialized<_TempStorage> {};


    /******************************************************************//**
     * \name Collective constructors
     *********************************************************************/
    //@{

    /**
     * \brief Collective constructor for 1D thread blocks using a private static allocation of shared memory as temporary storage.  Threads are identified using <tt>threadIdx.x</tt>.
     */
    __device__ __forceinline__ BlockStore()
    :
        temp_storage(PrivateStorage()),
        linear_tid(threadIdx.x)
    {}


    /**
     * \brief Collective constructor for 1D thread blocks using the specified memory allocation as temporary storage.  Threads are identified using <tt>threadIdx.x</tt>.
     */
    __device__ __forceinline__ BlockStore(
        TempStorage &temp_storage)             ///< [in] Reference to memory allocation having layout type TempStorage
    :
        temp_storage(temp_storage.Alias()),
        linear_tid(threadIdx.x)
    {}


    /**
     * \brief Collective constructor using a private static allocation of shared memory as temporary storage.  Each thread is identified using the supplied linear thread identifier
     */
    __device__ __forceinline__ BlockStore(
        int linear_tid)                        ///< [in] A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    :
        temp_storage(PrivateStorage()),
        linear_tid(linear_tid)
    {}


    /**
     * \brief Collective constructor using the specified memory allocation as temporary storage.  Each thread is identified using the supplied linear thread identifier.
     */
    __device__ __forceinline__ BlockStore(
        TempStorage &temp_storage,             ///< [in] Reference to memory allocation having layout type TempStorage
        int linear_tid)                        ///< [in] <b>[optional]</b> A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    :
        temp_storage(temp_storage.Alias()),
        linear_tid(linear_tid)
    {}


    //@}  end member group
    /******************************************************************//**
     * \name Data movement
     *********************************************************************/
    //@{


    /**
     * \brief Store items into a linear segment of memory.
     *
     * \blocked
     *
     * The code snippet below illustrates the storing of a "blocked" arrangement
     * of 512 integers across 128 threads (where each thread owns 4 consecutive items)
     * into a linear segment of memory.  The store is specialized for \p BLOCK_STORE_WARP_TRANSPOSE,
     * meaning items are locally reordered among threads so that memory references will be
     * efficiently coalesced using a warp-striped access pattern.
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(int *d_data, ...)
     * {
     *     // Specialize BlockStore for 128 threads owning 4 integer items each
     *     typedef cub::BlockStore<int*, 128, 4, BLOCK_STORE_WARP_TRANSPOSE> BlockStore;
     *
     *     // Allocate shared memory for BlockStore
     *     __shared__ typename BlockStore::TempStorage temp_storage;
     *
     *     // Obtain a segment of consecutive items that are blocked across threads
     *     int thread_data[4];
     *     ...
     *
     *     // Store items to linear memory
     *     int thread_data[4];
     *     BlockStore(temp_storage).Store(d_data, thread_data);
     *
     * \endcode
     * \par
     * Suppose the set of \p thread_data across the block of threads is
     * <tt>{ [0,1,2,3], [4,5,6,7], ..., [508,509,510,511] }</tt>.
     * The output \p d_data will be <tt>0, 1, 2, 3, 4, 5, ...</tt>.
     *
     */
    __device__ __forceinline__ void Store(
        OutputIteratorRA    block_itr,                  ///< [in] The thread block's base output iterator for storing to
        T                   (&items)[ITEMS_PER_THREAD]) ///< [in] Data to store
    {
        InternalStore(temp_storage, linear_tid).Store(block_itr, items);
    }

    /**
     * \brief Store items into a linear segment of memory, guarded by range.
     *
     * \blocked
     *
     * The code snippet below illustrates the guarded storing of a "blocked" arrangement
     * of 512 integers across 128 threads (where each thread owns 4 consecutive items)
     * into a linear segment of memory.  The store is specialized for \p BLOCK_STORE_WARP_TRANSPOSE,
     * meaning items are locally reordered among threads so that memory references will be
     * efficiently coalesced using a warp-striped access pattern.
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(int *d_data, int valid_items, ...)
     * {
     *     // Specialize BlockStore for 128 threads owning 4 integer items each
     *     typedef cub::BlockStore<int*, 128, 4, BLOCK_STORE_WARP_TRANSPOSE> BlockStore;
     *
     *     // Allocate shared memory for BlockStore
     *     __shared__ typename BlockStore::TempStorage temp_storage;
     *
     *     // Obtain a segment of consecutive items that are blocked across threads
     *     int thread_data[4];
     *     ...
     *
     *     // Store items to linear memory
     *     int thread_data[4];
     *     BlockStore(temp_storage).Store(d_data, thread_data, valid_items);
     *
     * \endcode
     * \par
     * Suppose the set of \p thread_data across the block of threads is
     * <tt>{ [0,1,2,3], [4,5,6,7], ..., [508,509,510,511] }</tt> and \p valid_items is \p 5.
     * The output \p d_data will be <tt>0, 1, 2, 3, 4, ?, ?, ?, ...</tt>, with
     * only the first two threads being unmasked to store portions of valid data.
     *
     */
    __device__ __forceinline__ void Store(
        OutputIteratorRA    block_itr,                  ///< [in] The thread block's base output iterator for storing to
        T                   (&items)[ITEMS_PER_THREAD], ///< [in] Data to store
        int                 valid_items)                ///< [in] Number of valid items to write
    {
        InternalStore(temp_storage, linear_tid).Store(block_itr, items, valid_items);
    }
};

/** @} */       // end group BlockModule

}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)

