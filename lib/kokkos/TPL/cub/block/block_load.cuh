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
 * Operations for reading linear tiles of data into the CUDA thread block.
 */

#pragma once

#include <iterator>

#include "../util_namespace.cuh"
#include "../util_macro.cuh"
#include "../util_type.cuh"
#include "../util_vector.cuh"
#include "../thread/thread_load.cuh"
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
 * \brief Load a linear segment of items into a blocked arrangement across the thread block using the specified cache modifier.
 *
 * \blocked
 *
 * \tparam MODIFIER             cub::PtxLoadModifier cache modifier.
 * \tparam T                    <b>[inferred]</b> The data type to load.
 * \tparam ITEMS_PER_THREAD     <b>[inferred]</b> The number of consecutive items partitioned onto each thread.
 * \tparam InputIteratorRA      <b>[inferred]</b> The random-access iterator type for input (may be a simple pointer type).
 */
template <
    PtxLoadModifier MODIFIER,
    typename        T,
    int             ITEMS_PER_THREAD,
    typename        InputIteratorRA>
__device__ __forceinline__ void LoadBlocked(
    int             linear_tid,                 ///< [in] A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    InputIteratorRA block_itr,                  ///< [in] The thread block's base input iterator for loading from
    T               (&items)[ITEMS_PER_THREAD]) ///< [out] Data to load
{
    // Load directly in thread-blocked order
    #pragma unroll
    for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
    {
        items[ITEM] = ThreadLoad<MODIFIER>(block_itr + (linear_tid * ITEMS_PER_THREAD) + ITEM);
    }
}


/**
 * \brief Load a linear segment of items into a blocked arrangement across the thread block using the specified cache modifier, guarded by range.
 *
 * \blocked
 *
 * \tparam MODIFIER             cub::PtxLoadModifier cache modifier.
 * \tparam T                    <b>[inferred]</b> The data type to load.
 * \tparam ITEMS_PER_THREAD     <b>[inferred]</b> The number of consecutive items partitioned onto each thread.
 * \tparam InputIteratorRA      <b>[inferred]</b> The random-access iterator type for input (may be a simple pointer type).
 */
template <
    PtxLoadModifier MODIFIER,
    typename        T,
    int             ITEMS_PER_THREAD,
    typename        InputIteratorRA>
__device__ __forceinline__ void LoadBlocked(
    int             linear_tid,                 ///< [in] A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    InputIteratorRA block_itr,                  ///< [in] The thread block's base input iterator for loading from
    T               (&items)[ITEMS_PER_THREAD], ///< [out] Data to load
    int             valid_items)                ///< [in] Number of valid items to load
{
    int bounds = valid_items - (linear_tid * ITEMS_PER_THREAD);

    #pragma unroll
    for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
    {
        if (ITEM < bounds)
        {
            items[ITEM] = ThreadLoad<MODIFIER>(block_itr + (linear_tid * ITEMS_PER_THREAD) + ITEM);
        }
    }
}


/**
 * \brief Load a linear segment of items into a blocked arrangement across the thread block using the specified cache modifier, guarded by range, with a fall-back assignment of out-of-bound elements..
 *
 * \blocked
 *
 * \tparam MODIFIER             cub::PtxLoadModifier cache modifier.
 * \tparam T                    <b>[inferred]</b> The data type to load.
 * \tparam ITEMS_PER_THREAD     <b>[inferred]</b> The number of consecutive items partitioned onto each thread.
 * \tparam InputIteratorRA      <b>[inferred]</b> The random-access iterator type for input (may be a simple pointer type).
 */
template <
    PtxLoadModifier MODIFIER,
    typename        T,
    int             ITEMS_PER_THREAD,
    typename        InputIteratorRA>
__device__ __forceinline__ void LoadBlocked(
    int             linear_tid,                 ///< [in] A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    InputIteratorRA block_itr,                  ///< [in] The thread block's base input iterator for loading from
    T               (&items)[ITEMS_PER_THREAD], ///< [out] Data to load
    int             valid_items,                ///< [in] Number of valid items to load
    T               oob_default)                ///< [in] Default value to assign out-of-bound items
{
    int bounds = valid_items - (linear_tid * ITEMS_PER_THREAD);

    #pragma unroll
    for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
    {
        items[ITEM] = (ITEM < bounds) ?
            ThreadLoad<MODIFIER>(block_itr + (linear_tid * ITEMS_PER_THREAD) + ITEM) :
            oob_default;
    }
}



//@}  end member group
/******************************************************************//**
 * \name Striped I/O
 *********************************************************************/
//@{


/**
 * \brief Load a linear segment of items into a striped arrangement across the thread block using the specified cache modifier.
 *
 * \striped
 *
 * \tparam MODIFIER             cub::PtxLoadModifier cache modifier.
 * \tparam BLOCK_THREADS        The thread block size in threads
 * \tparam T                    <b>[inferred]</b> The data type to load.
 * \tparam ITEMS_PER_THREAD     <b>[inferred]</b> The number of consecutive items partitioned onto each thread.
 * \tparam InputIteratorRA      <b>[inferred]</b> The random-access iterator type for input (may be a simple pointer type).
 */
template <
    PtxLoadModifier MODIFIER,
    int             BLOCK_THREADS,
    typename        T,
    int             ITEMS_PER_THREAD,
    typename        InputIteratorRA>
__device__ __forceinline__ void LoadStriped(
    int             linear_tid,                 ///< [in] A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    InputIteratorRA block_itr,                  ///< [in] The thread block's base input iterator for loading from
    T               (&items)[ITEMS_PER_THREAD]) ///< [out] Data to load
{
    #pragma unroll
    for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
    {
        items[ITEM] = ThreadLoad<MODIFIER>(block_itr + (ITEM * BLOCK_THREADS) + linear_tid);
    }
}


/**
 * \brief Load a linear segment of items into a striped arrangement across the thread block using the specified cache modifier, guarded by range
 *
 * \striped
 *
 * \tparam MODIFIER             cub::PtxLoadModifier cache modifier.
 * \tparam BLOCK_THREADS        The thread block size in threads
 * \tparam T                    <b>[inferred]</b> The data type to load.
 * \tparam ITEMS_PER_THREAD     <b>[inferred]</b> The number of consecutive items partitioned onto each thread.
 * \tparam InputIteratorRA      <b>[inferred]</b> The random-access iterator type for input (may be a simple pointer type).
 */
template <
    PtxLoadModifier MODIFIER,
    int             BLOCK_THREADS,
    typename        T,
    int             ITEMS_PER_THREAD,
    typename        InputIteratorRA>
__device__ __forceinline__ void LoadStriped(
    int             linear_tid,                 ///< [in] A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    InputIteratorRA block_itr,                  ///< [in] The thread block's base input iterator for loading from
    T               (&items)[ITEMS_PER_THREAD], ///< [out] Data to load
    int             valid_items)                ///< [in] Number of valid items to load
{
    int bounds = valid_items - linear_tid;

    #pragma unroll
    for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
    {
        if (ITEM * BLOCK_THREADS < bounds)
        {
            items[ITEM] = ThreadLoad<MODIFIER>(block_itr + linear_tid + (ITEM * BLOCK_THREADS));
        }
    }
}


/**
 * \brief Load a linear segment of items into a striped arrangement across the thread block using the specified cache modifier, guarded by range, with a fall-back assignment of out-of-bound elements.
 *
 * \striped
 *
 * \tparam MODIFIER             cub::PtxLoadModifier cache modifier.
 * \tparam BLOCK_THREADS        The thread block size in threads
 * \tparam T                    <b>[inferred]</b> The data type to load.
 * \tparam ITEMS_PER_THREAD     <b>[inferred]</b> The number of consecutive items partitioned onto each thread.
 * \tparam InputIteratorRA      <b>[inferred]</b> The random-access iterator type for input (may be a simple pointer type).
 */
template <
    PtxLoadModifier MODIFIER,
    int             BLOCK_THREADS,
    typename        T,
    int             ITEMS_PER_THREAD,
    typename        InputIteratorRA>
__device__ __forceinline__ void LoadStriped(
    int             linear_tid,                 ///< [in] A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    InputIteratorRA block_itr,                  ///< [in] The thread block's base input iterator for loading from
    T               (&items)[ITEMS_PER_THREAD], ///< [out] Data to load
    int             valid_items,                ///< [in] Number of valid items to load
    T               oob_default)                ///< [in] Default value to assign out-of-bound items
{
    int bounds = valid_items - linear_tid;

    #pragma unroll
    for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
    {
        items[ITEM] = (ITEM * BLOCK_THREADS < bounds) ?
             ThreadLoad<MODIFIER>(block_itr + linear_tid + (ITEM * BLOCK_THREADS)) :
             oob_default;
    }
}



//@}  end member group
/******************************************************************//**
 * \name Warp-striped I/O
 *********************************************************************/
//@{


/**
 * \brief Load a linear segment of items into a warp-striped arrangement across the thread block using the specified cache modifier.
 *
 * \warpstriped
 *
 * \par Usage Considerations
 * The number of threads in the thread block must be a multiple of the architecture's warp size.
 *
 * \tparam MODIFIER             cub::PtxLoadModifier cache modifier.
 * \tparam T                    <b>[inferred]</b> The data type to load.
 * \tparam ITEMS_PER_THREAD     <b>[inferred]</b> The number of consecutive items partitioned onto each thread.
 * \tparam InputIteratorRA      <b>[inferred]</b> The random-access iterator type for input (may be a simple pointer type).
 */
template <
    PtxLoadModifier MODIFIER,
    typename        T,
    int             ITEMS_PER_THREAD,
    typename        InputIteratorRA>
__device__ __forceinline__ void LoadWarpStriped(
    int             linear_tid,                 ///< [in] A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    InputIteratorRA block_itr,                  ///< [in] The thread block's base input iterator for loading from
    T               (&items)[ITEMS_PER_THREAD]) ///< [out] Data to load
{
    int tid         = linear_tid & (PtxArchProps::WARP_THREADS - 1);
    int wid         = linear_tid >> PtxArchProps::LOG_WARP_THREADS;
    int warp_offset = wid * PtxArchProps::WARP_THREADS * ITEMS_PER_THREAD;

    // Load directly in warp-striped order
    #pragma unroll
    for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
    {
        items[ITEM] = ThreadLoad<MODIFIER>(block_itr + warp_offset + tid + (ITEM * PtxArchProps::WARP_THREADS));
    }
}


/**
 * \brief Load a linear segment of items into a warp-striped arrangement across the thread block using the specified cache modifier, guarded by range
 *
 * \warpstriped
 *
 * \par Usage Considerations
 * The number of threads in the thread block must be a multiple of the architecture's warp size.
 *
 * \tparam MODIFIER             cub::PtxLoadModifier cache modifier.
 * \tparam T                    <b>[inferred]</b> The data type to load.
 * \tparam ITEMS_PER_THREAD     <b>[inferred]</b> The number of consecutive items partitioned onto each thread.
 * \tparam InputIteratorRA      <b>[inferred]</b> The random-access iterator type for input (may be a simple pointer type).
 */
template <
    PtxLoadModifier MODIFIER,
    typename        T,
    int             ITEMS_PER_THREAD,
    typename        InputIteratorRA>
__device__ __forceinline__ void LoadWarpStriped(
    int             linear_tid,                 ///< [in] A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    InputIteratorRA block_itr,                  ///< [in] The thread block's base input iterator for loading from
    T               (&items)[ITEMS_PER_THREAD], ///< [out] Data to load
    int             valid_items)               ///< [in] Number of valid items to load
{
    int tid                 = linear_tid & (PtxArchProps::WARP_THREADS - 1);
    int wid                 = linear_tid >> PtxArchProps::LOG_WARP_THREADS;
    int warp_offset         = wid * PtxArchProps::WARP_THREADS * ITEMS_PER_THREAD;
    int bounds              = valid_items - warp_offset - tid;

    // Load directly in warp-striped order
    #pragma unroll
    for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
    {
        if ((ITEM * PtxArchProps::WARP_THREADS) < bounds)
        {
            items[ITEM] = ThreadLoad<MODIFIER>(block_itr + warp_offset + tid + (ITEM * PtxArchProps::WARP_THREADS));
        }
    }
}


/**
 * \brief Load a linear segment of items into a warp-striped arrangement across the thread block using the specified cache modifier, guarded by range, with a fall-back assignment of out-of-bound elements.
 *
 * \warpstriped
 *
 * \par Usage Considerations
 * The number of threads in the thread block must be a multiple of the architecture's warp size.
 *
 * \tparam MODIFIER             cub::PtxLoadModifier cache modifier.
 * \tparam T                    <b>[inferred]</b> The data type to load.
 * \tparam ITEMS_PER_THREAD     <b>[inferred]</b> The number of consecutive items partitioned onto each thread.
 * \tparam InputIteratorRA      <b>[inferred]</b> The random-access iterator type for input (may be a simple pointer type).
 */
template <
    PtxLoadModifier MODIFIER,
    typename        T,
    int             ITEMS_PER_THREAD,
    typename        InputIteratorRA>
__device__ __forceinline__ void LoadWarpStriped(
    int             linear_tid,                 ///< [in] A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    InputIteratorRA block_itr,                  ///< [in] The thread block's base input iterator for loading from
    T               (&items)[ITEMS_PER_THREAD], ///< [out] Data to load
    int             valid_items,               ///< [in] Number of valid items to load
    T               oob_default)                ///< [in] Default value to assign out-of-bound items
{
    int tid         = linear_tid & (PtxArchProps::WARP_THREADS - 1);
    int wid         = linear_tid >> PtxArchProps::LOG_WARP_THREADS;
    int warp_offset = wid * PtxArchProps::WARP_THREADS * ITEMS_PER_THREAD;
    int bounds      = valid_items - warp_offset - tid;

    // Load directly in warp-striped order
    #pragma unroll
    for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
    {
        items[ITEM] = ((ITEM * PtxArchProps::WARP_THREADS) < bounds) ?
            ThreadLoad<MODIFIER>(block_itr + warp_offset + tid + (ITEM * PtxArchProps::WARP_THREADS)) :
            oob_default;
    }
}



//@}  end member group
/******************************************************************//**
 * \name Blocked, vectorized I/O
 *********************************************************************/
//@{

/**
 * \brief Load a linear segment of items into a blocked arrangement across the thread block using the specified cache modifier.
 *
 * \blocked
 *
 * The input offset (\p block_ptr + \p block_offset) must be quad-item aligned
 *
 * The following conditions will prevent vectorization and loading will fall back to cub::BLOCK_LOAD_DIRECT:
 *   - \p ITEMS_PER_THREAD is odd
 *   - The data type \p T is not a built-in primitive or CUDA vector type (e.g., \p short, \p int2, \p double, \p float2, etc.)
 *
 * \tparam MODIFIER             cub::PtxLoadModifier cache modifier.
 * \tparam T                    <b>[inferred]</b> The data type to load.
 * \tparam ITEMS_PER_THREAD     <b>[inferred]</b> The number of consecutive items partitioned onto each thread.
 */
template <
    PtxLoadModifier MODIFIER,
    typename        T,
    int             ITEMS_PER_THREAD>
__device__ __forceinline__ void LoadBlockedVectorized(
    int             linear_tid,                 ///< [in] A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    T               *block_ptr,                 ///< [in] Input pointer for loading from
    T               (&items)[ITEMS_PER_THREAD]) ///< [out] Data to load
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

    // Alias local data (use raw_items array here which should get optimized away to prevent conservative PTXAS lmem spilling)
    T raw_items[ITEMS_PER_THREAD];

    // Direct-load using vector types
    LoadBlocked<MODIFIER>(
        linear_tid,
        reinterpret_cast<Vector *>(block_ptr),
        reinterpret_cast<Vector (&)[VECTORS_PER_THREAD]>(raw_items));

    // Copy
    #pragma unroll
    for (int ITEM = 0; ITEM < ITEMS_PER_THREAD; ITEM++)
    {
        items[ITEM] = raw_items[ITEM];
    }
}


//@}  end member group

/** @} */       // end group IoModule



//-----------------------------------------------------------------------------
// Generic BlockLoad abstraction
//-----------------------------------------------------------------------------

/**
 * \brief cub::BlockLoadAlgorithm enumerates alternative algorithms for cub::BlockLoad to read a linear segment of data from memory into a blocked arrangement across a CUDA thread block.
 */
enum BlockLoadAlgorithm
{
    /**
     * \par Overview
     *
     * A [<em>blocked arrangement</em>](index.html#sec5sec4) of data is read
     * directly from memory.  The thread block reads items in a parallel "raking" fashion: thread<sub><em>i</em></sub>
     * reads the <em>i</em><sup>th</sup> segment of consecutive elements.
     *
     * \par Performance Considerations
     * - The utilization of memory transactions (coalescing) decreases as the
     *   access stride between threads increases (i.e., the number items per thread).
     */
    BLOCK_LOAD_DIRECT,

    /**
     * \par Overview
     *
     * A [<em>blocked arrangement</em>](index.html#sec5sec4) of data is read directly
     * from memory using CUDA's built-in vectorized loads as a coalescing optimization.
     * The thread block reads items in a parallel "raking" fashion: thread<sub><em>i</em></sub> uses vector loads to
     * read the <em>i</em><sup>th</sup> segment of consecutive elements.
     *
     * For example, <tt>ld.global.v4.s32</tt> instructions will be generated when \p T = \p int and \p ITEMS_PER_THREAD > 4.
     *
     * \par Performance Considerations
     * - The utilization of memory transactions (coalescing) remains high until the the
     *   access stride between threads (i.e., the number items per thread) exceeds the
     *   maximum vector load width (typically 4 items or 64B, whichever is lower).
     * - The following conditions will prevent vectorization and loading will fall back to cub::BLOCK_LOAD_DIRECT:
     *   - \p ITEMS_PER_THREAD is odd
     *   - The \p InputIteratorRA is not a simple pointer type
     *   - The block input offset is not quadword-aligned
     *   - The data type \p T is not a built-in primitive or CUDA vector type (e.g., \p short, \p int2, \p double, \p float2, etc.)
     */
    BLOCK_LOAD_VECTORIZE,

    /**
     * \par Overview
     *
     * A [<em>striped arrangement</em>](index.html#sec5sec4) of data is read
     * directly from memory and then is locally transposed into a
     * [<em>blocked arrangement</em>](index.html#sec5sec4). The thread block
     * reads items in a parallel "strip-mining" fashion:
     * thread<sub><em>i</em></sub> reads items having stride \p BLOCK_THREADS
     * between them. cub::BlockExchange is then used to locally reorder the items
     * into a [<em>blocked arrangement</em>](index.html#sec5sec4).
     *
     * \par Performance Considerations
     * - The utilization of memory transactions (coalescing) remains high regardless
     *   of items loaded per thread.
     * - The local reordering incurs slightly longer latencies and throughput than the
     *   direct cub::BLOCK_LOAD_DIRECT and cub::BLOCK_LOAD_VECTORIZE alternatives.
     */
    BLOCK_LOAD_TRANSPOSE,


    /**
     * \par Overview
     *
     * A [<em>warp-striped arrangement</em>](index.html#sec5sec4) of data is read
     * directly from memory and then is locally transposed into a
     * [<em>blocked arrangement</em>](index.html#sec5sec4). Each warp reads its own
     * contiguous segment in a parallel "strip-mining" fashion: lane<sub><em>i</em></sub>
     * reads items having stride \p WARP_THREADS between them. cub::BlockExchange
     * is then used to locally reorder the items into a
     * [<em>blocked arrangement</em>](index.html#sec5sec4).
     *
     * \par Usage Considerations
     * - BLOCK_THREADS must be a multiple of WARP_THREADS
     *
     * \par Performance Considerations
     * - The utilization of memory transactions (coalescing) remains high regardless
     *   of items loaded per thread.
     * - The local reordering incurs slightly longer latencies and throughput than the
     *   direct cub::BLOCK_LOAD_DIRECT and cub::BLOCK_LOAD_VECTORIZE alternatives.
     */
    BLOCK_LOAD_WARP_TRANSPOSE,
};


/**
 * \brief The BlockLoad class provides [<em>collective</em>](index.html#sec0) data movement methods for loading a linear segment of items from memory into a [<em>blocked arrangement</em>](index.html#sec5sec4) across a CUDA thread block.  ![](block_load_logo.png)
 * \ingroup BlockModule
 *
 * \par Overview
 * The BlockLoad class provides a single data movement abstraction that can be specialized
 * to implement different cub::BlockLoadAlgorithm strategies.  This facilitates different
 * performance policies for different architectures, data types, granularity sizes, etc.
 *
 * \par
 * Optionally, BlockLoad can be specialized by different data movement strategies:
 *   -# <b>cub::BLOCK_LOAD_DIRECT</b>.  A [<em>blocked arrangement</em>](index.html#sec5sec4)
 *      of data is read directly from memory.  [More...](\ref cub::BlockLoadAlgorithm)
 *   -# <b>cub::BLOCK_LOAD_VECTORIZE</b>.  A [<em>blocked arrangement</em>](index.html#sec5sec4)
 *      of data is read directly from memory using CUDA's built-in vectorized loads as a
 *      coalescing optimization.    [More...](\ref cub::BlockLoadAlgorithm)
 *   -# <b>cub::BLOCK_LOAD_TRANSPOSE</b>.  A [<em>striped arrangement</em>](index.html#sec5sec4)
 *      of data is read directly from memory and is then locally transposed into a
 *      [<em>blocked arrangement</em>](index.html#sec5sec4).  [More...](\ref cub::BlockLoadAlgorithm)
 *   -# <b>cub::BLOCK_LOAD_WARP_TRANSPOSE</b>.  A [<em>warp-striped arrangement</em>](index.html#sec5sec4)
 *      of data is read directly from memory and is then locally transposed into a
 *      [<em>blocked arrangement</em>](index.html#sec5sec4).  [More...](\ref cub::BlockLoadAlgorithm)
 *
 * \tparam InputIteratorRA      The input iterator type (may be a simple pointer type).
 * \tparam BLOCK_THREADS        The thread block size in threads.
 * \tparam ITEMS_PER_THREAD     The number of consecutive items partitioned onto each thread.
 * \tparam ALGORITHM            <b>[optional]</b> cub::BlockLoadAlgorithm tuning policy.  default: cub::BLOCK_LOAD_DIRECT.
 * \tparam MODIFIER             <b>[optional]</b> cub::PtxLoadModifier cache modifier.  default: cub::LOAD_DEFAULT.
 * \tparam WARP_TIME_SLICING    <b>[optional]</b> For transposition-based cub::BlockLoadAlgorithm parameterizations that utilize shared memory: When \p true, only use enough shared memory for a single warp's worth of data, time-slicing the block-wide exchange over multiple synchronized rounds (default: false)
 *
 * \par A Simple Example
 * \blockcollective{BlockLoad}
 * \par
 * The code snippet below illustrates the loading of a linear
 * segment of 512 integers into a "blocked" arrangement across 128 threads where each
 * thread owns 4 consecutive items.  The load is specialized for \p BLOCK_LOAD_WARP_TRANSPOSE,
 * meaning memory references are efficiently coalesced using a warp-striped access
 * pattern (after which items are locally reordered among threads).
 * \par
 * \code
 * #include <cub/cub.cuh>
 *
 * __global__ void ExampleKernel(int *d_data, ...)
 * {
 *     // Specialize BlockLoad for 128 threads owning 4 integer items each
 *     typedef cub::BlockLoad<int*, 128, 4, BLOCK_LOAD_WARP_TRANSPOSE> BlockLoad;
 *
 *     // Allocate shared memory for BlockLoad
 *     __shared__ typename BlockLoad::TempStorage temp_storage;
 *
 *     // Load a segment of consecutive items that are blocked across threads
 *     int thread_data[4];
 *     BlockLoad(temp_storage).Load(d_data, thread_data);
 *
 * \endcode
 * \par
 * Suppose the input \p d_data is <tt>0, 1, 2, 3, 4, 5, ...</tt>.
 * The set of \p thread_data across the block of threads in those threads will be
 * <tt>{ [0,1,2,3], [4,5,6,7], ..., [508,509,510,511] }</tt>.
 *
 */
template <
    typename            InputIteratorRA,
    int                 BLOCK_THREADS,
    int                 ITEMS_PER_THREAD,
    BlockLoadAlgorithm  ALGORITHM = BLOCK_LOAD_DIRECT,
    PtxLoadModifier     MODIFIER = LOAD_DEFAULT,
    bool                WARP_TIME_SLICING = false>
class BlockLoad
{
private:

    /******************************************************************************
     * Constants and typed definitions
     ******************************************************************************/

    // Data type of input iterator
    typedef typename std::iterator_traits<InputIteratorRA>::value_type T;


    /******************************************************************************
     * Algorithmic variants
     ******************************************************************************/

    /// Load helper
    template <BlockLoadAlgorithm _POLICY, int DUMMY = 0>
    struct LoadInternal;


    /**
     * BLOCK_LOAD_DIRECT specialization of load helper
     */
    template <int DUMMY>
    struct LoadInternal<BLOCK_LOAD_DIRECT, DUMMY>
    {
        /// Shared memory storage layout type
        typedef NullType TempStorage;

        /// Linear thread-id
        int linear_tid;

        /// Constructor
        __device__ __forceinline__ LoadInternal(
            TempStorage &temp_storage,
            int linear_tid)
        :
            linear_tid(linear_tid)
        {}

        /// Load a linear segment of items from memory
        __device__ __forceinline__ void Load(
            InputIteratorRA block_itr,                      ///< [in] The thread block's base input iterator for loading from
            T               (&items)[ITEMS_PER_THREAD])     ///< [out] Data to load
        {
            LoadBlocked<MODIFIER>(linear_tid, block_itr, items);
        }

        /// Load a linear segment of items from memory, guarded by range
        __device__ __forceinline__ void Load(
            InputIteratorRA block_itr,                      ///< [in] The thread block's base input iterator for loading from
            T               (&items)[ITEMS_PER_THREAD],     ///< [out] Data to load
            int             valid_items)                    ///< [in] Number of valid items to load
        {
            LoadBlocked<MODIFIER>(linear_tid, block_itr, items, valid_items);
        }

        /// Load a linear segment of items from memory, guarded by range, with a fall-back assignment of out-of-bound elements
        __device__ __forceinline__ void Load(
            InputIteratorRA block_itr,                      ///< [in] The thread block's base input iterator for loading from
            T               (&items)[ITEMS_PER_THREAD],     ///< [out] Data to load
            int             valid_items,                    ///< [in] Number of valid items to load
            T               oob_default)                    ///< [in] Default value to assign out-of-bound items
        {
            LoadBlocked<MODIFIER>(linear_tid, block_itr, items, valid_items, oob_default);
        }

    };


    /**
     * BLOCK_LOAD_VECTORIZE specialization of load helper
     */
    template <int DUMMY>
    struct LoadInternal<BLOCK_LOAD_VECTORIZE, DUMMY>
    {
        /// Shared memory storage layout type
        typedef NullType TempStorage;

        /// Linear thread-id
        int linear_tid;

        /// Constructor
        __device__ __forceinline__ LoadInternal(
            TempStorage &temp_storage,
            int linear_tid)
        :
            linear_tid(linear_tid)
        {}

        /// Load a linear segment of items from memory, specialized for native pointer types (attempts vectorization)
        __device__ __forceinline__ void Load(
            T               *block_ptr,                     ///< [in] The thread block's base input iterator for loading from
            T               (&items)[ITEMS_PER_THREAD])     ///< [out] Data to load
        {
            LoadBlockedVectorized<MODIFIER>(linear_tid, block_ptr, items);
        }

        /// Load a linear segment of items from memory, specialized for opaque input iterators (skips vectorization)
        template <
            typename T,
            typename _InputIteratorRA>
        __device__ __forceinline__ void Load(
            _InputIteratorRA    block_itr,                  ///< [in] The thread block's base input iterator for loading from
            T                   (&items)[ITEMS_PER_THREAD]) ///< [out] Data to load
        {
            LoadBlocked<MODIFIER>(linear_tid, block_itr, items);
        }

        /// Load a linear segment of items from memory, guarded by range (skips vectorization)
        __device__ __forceinline__ void Load(
            InputIteratorRA block_itr,                      ///< [in] The thread block's base input iterator for loading from
            T               (&items)[ITEMS_PER_THREAD],     ///< [out] Data to load
            int             valid_items)                    ///< [in] Number of valid items to load
        {
            LoadBlocked<MODIFIER>(linear_tid, block_itr, items, valid_items);
        }

        /// Load a linear segment of items from memory, guarded by range, with a fall-back assignment of out-of-bound elements (skips vectorization)
        __device__ __forceinline__ void Load(
            InputIteratorRA block_itr,                      ///< [in] The thread block's base input iterator for loading from
            T               (&items)[ITEMS_PER_THREAD],     ///< [out] Data to load
            int             valid_items,                    ///< [in] Number of valid items to load
            T               oob_default)                    ///< [in] Default value to assign out-of-bound items
        {
            LoadBlocked<MODIFIER>(linear_tid, block_itr, items, valid_items, oob_default);
        }

    };


    /**
     * BLOCK_LOAD_TRANSPOSE specialization of load helper
     */
    template <int DUMMY>
    struct LoadInternal<BLOCK_LOAD_TRANSPOSE, DUMMY>
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
        __device__ __forceinline__ LoadInternal(
            TempStorage &temp_storage,
            int linear_tid)
        :
            temp_storage(temp_storage.Alias()),
            linear_tid(linear_tid)
        {}

        /// Load a linear segment of items from memory
        __device__ __forceinline__ void Load(
            InputIteratorRA block_itr,                      ///< [in] The thread block's base input iterator for loading from
            T               (&items)[ITEMS_PER_THREAD])     ///< [out] Data to load{
        {
            LoadStriped<MODIFIER, BLOCK_THREADS>(linear_tid, block_itr, items);
            BlockExchange(temp_storage, linear_tid).StripedToBlocked(items);
        }

        /// Load a linear segment of items from memory, guarded by range
        __device__ __forceinline__ void Load(
            InputIteratorRA block_itr,                      ///< [in] The thread block's base input iterator for loading from
            T               (&items)[ITEMS_PER_THREAD],     ///< [out] Data to load
            int             valid_items)                    ///< [in] Number of valid items to load
        {
            LoadStriped<MODIFIER, BLOCK_THREADS>(linear_tid, block_itr, items, valid_items);
            BlockExchange(temp_storage, linear_tid).StripedToBlocked(items);
        }

        /// Load a linear segment of items from memory, guarded by range, with a fall-back assignment of out-of-bound elements
        __device__ __forceinline__ void Load(
            InputIteratorRA block_itr,                      ///< [in] The thread block's base input iterator for loading from
            T               (&items)[ITEMS_PER_THREAD],     ///< [out] Data to load
            int             valid_items,                    ///< [in] Number of valid items to load
            T               oob_default)                    ///< [in] Default value to assign out-of-bound items
        {
            LoadStriped<MODIFIER, BLOCK_THREADS>(linear_tid, block_itr, items, valid_items, oob_default);
            BlockExchange(temp_storage, linear_tid).StripedToBlocked(items);
        }

    };


    /**
     * BLOCK_LOAD_WARP_TRANSPOSE specialization of load helper
     */
    template <int DUMMY>
    struct LoadInternal<BLOCK_LOAD_WARP_TRANSPOSE, DUMMY>
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
        __device__ __forceinline__ LoadInternal(
            TempStorage &temp_storage,
            int linear_tid)
        :
            temp_storage(temp_storage.Alias()),
            linear_tid(linear_tid)
        {}

        /// Load a linear segment of items from memory
        __device__ __forceinline__ void Load(
            InputIteratorRA block_itr,                      ///< [in] The thread block's base input iterator for loading from
            T               (&items)[ITEMS_PER_THREAD])     ///< [out] Data to load{
        {
            LoadWarpStriped<MODIFIER>(linear_tid, block_itr, items);
            BlockExchange(temp_storage, linear_tid).WarpStripedToBlocked(items);
        }

        /// Load a linear segment of items from memory, guarded by range
        __device__ __forceinline__ void Load(
            InputIteratorRA block_itr,                      ///< [in] The thread block's base input iterator for loading from
            T               (&items)[ITEMS_PER_THREAD],     ///< [out] Data to load
            int             valid_items)                    ///< [in] Number of valid items to load
        {
            LoadWarpStriped<MODIFIER>(linear_tid, block_itr, items, valid_items);
            BlockExchange(temp_storage, linear_tid).WarpStripedToBlocked(items);
        }


        /// Load a linear segment of items from memory, guarded by range, with a fall-back assignment of out-of-bound elements
        __device__ __forceinline__ void Load(
            InputIteratorRA block_itr,                      ///< [in] The thread block's base input iterator for loading from
            T               (&items)[ITEMS_PER_THREAD],     ///< [out] Data to load
            int             valid_items,                    ///< [in] Number of valid items to load
            T               oob_default)                    ///< [in] Default value to assign out-of-bound items
        {
            LoadWarpStriped<MODIFIER>(linear_tid, block_itr, items, valid_items, oob_default);
            BlockExchange(temp_storage, linear_tid).WarpStripedToBlocked(items);
        }
    };


    /******************************************************************************
     * Type definitions
     ******************************************************************************/

    /// Internal load implementation to use
    typedef LoadInternal<ALGORITHM> InternalLoad;


    /// Shared memory storage layout type
    typedef typename InternalLoad::TempStorage _TempStorage;


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

    /// \smemstorage{BlockLoad}
    struct TempStorage : Uninitialized<_TempStorage> {};


    /******************************************************************//**
     * \name Collective constructors
     *********************************************************************/
    //@{

    /**
     * \brief Collective constructor for 1D thread blocks using a private static allocation of shared memory as temporary storage.  Threads are identified using <tt>threadIdx.x</tt>.
     */
    __device__ __forceinline__ BlockLoad()
    :
        temp_storage(PrivateStorage()),
        linear_tid(threadIdx.x)
    {}


    /**
     * \brief Collective constructor for 1D thread blocks using the specified memory allocation as temporary storage.  Threads are identified using <tt>threadIdx.x</tt>.
     */
    __device__ __forceinline__ BlockLoad(
        TempStorage &temp_storage)             ///< [in] Reference to memory allocation having layout type TempStorage
    :
        temp_storage(temp_storage.Alias()),
        linear_tid(threadIdx.x)
    {}


    /**
     * \brief Collective constructor using a private static allocation of shared memory as temporary storage.  Each thread is identified using the supplied linear thread identifier
     */
    __device__ __forceinline__ BlockLoad(
        int linear_tid)                        ///< [in] A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    :
        temp_storage(PrivateStorage()),
        linear_tid(linear_tid)
    {}


    /**
     * \brief Collective constructor using the specified memory allocation as temporary storage.  Each thread is identified using the supplied linear thread identifier.
     */
    __device__ __forceinline__ BlockLoad(
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
     * \brief Load a linear segment of items from memory.
     *
     * \blocked
     *
     * The code snippet below illustrates the loading of a linear
     * segment of 512 integers into a "blocked" arrangement across 128 threads where each
     * thread owns 4 consecutive items.  The load is specialized for \p BLOCK_LOAD_WARP_TRANSPOSE,
     * meaning memory references are efficiently coalesced using a warp-striped access
     * pattern (after which items are locally reordered among threads).
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(int *d_data, ...)
     * {
     *     // Specialize BlockLoad for 128 threads owning 4 integer items each
     *     typedef cub::BlockLoad<int*, 128, 4, BLOCK_LOAD_WARP_TRANSPOSE> BlockLoad;
     *
     *     // Allocate shared memory for BlockLoad
     *     __shared__ typename BlockLoad::TempStorage temp_storage;
     *
     *     // Load a segment of consecutive items that are blocked across threads
     *     int thread_data[4];
     *     BlockLoad(temp_storage).Load(d_data, thread_data);
     *
     * \endcode
     * \par
     * Suppose the input \p d_data is <tt>0, 1, 2, 3, 4, 5, ...</tt>.
     * The set of \p thread_data across the block of threads in those threads will be
     * <tt>{ [0,1,2,3], [4,5,6,7], ..., [508,509,510,511] }</tt>.
     *
     */
    __device__ __forceinline__ void Load(
        InputIteratorRA block_itr,                  ///< [in] The thread block's base input iterator for loading from
        T               (&items)[ITEMS_PER_THREAD]) ///< [out] Data to load
    {
        InternalLoad(temp_storage, linear_tid).Load(block_itr, items);
    }


    /**
     * \brief Load a linear segment of items from memory, guarded by range.
     *
     * \blocked
     *
     * The code snippet below illustrates the guarded loading of a linear
     * segment of 512 integers into a "blocked" arrangement across 128 threads where each
     * thread owns 4 consecutive items.  The load is specialized for \p BLOCK_LOAD_WARP_TRANSPOSE,
     * meaning memory references are efficiently coalesced using a warp-striped access
     * pattern (after which items are locally reordered among threads).
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(int *d_data, int valid_items, ...)
     * {
     *     // Specialize BlockLoad for 128 threads owning 4 integer items each
     *     typedef cub::BlockLoad<int*, 128, 4, BLOCK_LOAD_WARP_TRANSPOSE> BlockLoad;
     *
     *     // Allocate shared memory for BlockLoad
     *     __shared__ typename BlockLoad::TempStorage temp_storage;
     *
     *     // Load a segment of consecutive items that are blocked across threads
     *     int thread_data[4];
     *     BlockLoad(temp_storage).Load(d_data, thread_data, valid_items);
     *
     * \endcode
     * \par
     * Suppose the input \p d_data is <tt>0, 1, 2, 3, 4, 5, 6...</tt> and \p valid_items is \p 5.
     * The set of \p thread_data across the block of threads in those threads will be
     * <tt>{ [0,1,2,3], [4,?,?,?], ..., [?,?,?,?] }</tt>, with only the first two threads
     * being unmasked to load portions of valid data (and other items remaining unassigned).
     *
     */
    __device__ __forceinline__ void Load(
        InputIteratorRA block_itr,                  ///< [in] The thread block's base input iterator for loading from
        T               (&items)[ITEMS_PER_THREAD], ///< [out] Data to load
        int             valid_items)                ///< [in] Number of valid items to load
    {
        InternalLoad(temp_storage, linear_tid).Load(block_itr, items, valid_items);
    }


    /**
     * \brief Load a linear segment of items from memory, guarded by range, with a fall-back assignment of out-of-bound elements
     *
     * \blocked
     *
     * The code snippet below illustrates the guarded loading of a linear
     * segment of 512 integers into a "blocked" arrangement across 128 threads where each
     * thread owns 4 consecutive items.  The load is specialized for \p BLOCK_LOAD_WARP_TRANSPOSE,
     * meaning memory references are efficiently coalesced using a warp-striped access
     * pattern (after which items are locally reordered among threads).
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(int *d_data, int valid_items, ...)
     * {
     *     // Specialize BlockLoad for 128 threads owning 4 integer items each
     *     typedef cub::BlockLoad<int*, 128, 4, BLOCK_LOAD_WARP_TRANSPOSE> BlockLoad;
     *
     *     // Allocate shared memory for BlockLoad
     *     __shared__ typename BlockLoad::TempStorage temp_storage;
     *
     *     // Load a segment of consecutive items that are blocked across threads
     *     int thread_data[4];
     *     BlockLoad(temp_storage).Load(d_data, thread_data, valid_items, -1);
     *
     * \endcode
     * \par
     * Suppose the input \p d_data is <tt>0, 1, 2, 3, 4, 5, 6...</tt>,
     * \p valid_items is \p 5, and the out-of-bounds default is \p -1.
     * The set of \p thread_data across the block of threads in those threads will be
     * <tt>{ [0,1,2,3], [4,-1,-1,-1], ..., [-1,-1,-1,-1] }</tt>, with only the first two threads
     * being unmasked to load portions of valid data (and other items are assigned \p -1)
     *
     */
    __device__ __forceinline__ void Load(
        InputIteratorRA block_itr,                  ///< [in] The thread block's base input iterator for loading from
        T               (&items)[ITEMS_PER_THREAD], ///< [out] Data to load
        int             valid_items,                ///< [in] Number of valid items to load
        T               oob_default)                ///< [in] Default value to assign out-of-bound items
    {
        InternalLoad(temp_storage, linear_tid).Load(block_itr, items, valid_items, oob_default);
    }


    //@}  end member group

};


}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)

