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
 * The cub::BlockReduce class provides [<em>collective</em>](index.html#sec0) methods for computing a parallel reduction of items partitioned across a CUDA thread block.
 */

#pragma once

#include "specializations/block_reduce_raking.cuh"
#include "specializations/block_reduce_warp_reductions.cuh"
#include "../util_type.cuh"
#include "../thread/thread_operators.cuh"
#include "../util_namespace.cuh"

/// Optional outer namespace(s)
CUB_NS_PREFIX

/// CUB namespace
namespace cub {



/******************************************************************************
 * Algorithmic variants
 ******************************************************************************/

/**
 * BlockReduceAlgorithm enumerates alternative algorithms for parallel
 * reduction across a CUDA threadblock.
 */
enum BlockReduceAlgorithm
{

    /**
     * \par Overview
     * An efficient "raking" reduction algorithm.  Execution is comprised of
     * three phases:
     * -# Upsweep sequential reduction in registers (if threads contribute more
     *    than one input each).  Each thread then places the partial reduction
     *    of its item(s) into shared memory.
     * -# Upsweep sequential reduction in shared memory.  Threads within a
     *    single warp rake across segments of shared partial reductions.
     * -# A warp-synchronous Kogge-Stone style reduction within the raking warp.
     *
     * \par
     * \image html block_reduce.png
     * <div class="centercaption">\p BLOCK_REDUCE_RAKING data flow for a hypothetical 16-thread threadblock and 4-thread raking warp.</div>
     *
     * \par Performance Considerations
     * - Although this variant may suffer longer turnaround latencies when the
     *   GPU is under-occupied, it can often provide higher overall throughput
     *   across the GPU when suitably occupied.
     */
    BLOCK_REDUCE_RAKING,


    /**
     * \par Overview
     * A quick "tiled warp-reductions" reduction algorithm.  Execution is
     * comprised of four phases:
     * -# Upsweep sequential reduction in registers (if threads contribute more
     *    than one input each).  Each thread then places the partial reduction
     *    of its item(s) into shared memory.
     * -# Compute a shallow, but inefficient warp-synchronous Kogge-Stone style
     *    reduction within each warp.
     * -# A propagation phase where the warp reduction outputs in each warp are
     *    updated with the aggregate from each preceding warp.
     *
     * \par
     * \image html block_scan_warpscans.png
     * <div class="centercaption">\p BLOCK_REDUCE_WARP_REDUCTIONS data flow for a hypothetical 16-thread threadblock and 4-thread raking warp.</div>
     *
     * \par Performance Considerations
     * - Although this variant may suffer lower overall throughput across the
     *   GPU because due to a heavy reliance on inefficient warp-reductions, it
     *   can often provide lower turnaround latencies when the GPU is
     *   under-occupied.
     */
    BLOCK_REDUCE_WARP_REDUCTIONS,
};


/******************************************************************************
 * Block reduce
 ******************************************************************************/

/**
 * \brief The BlockReduce class provides [<em>collective</em>](index.html#sec0) methods for computing a parallel reduction of items partitioned across a CUDA thread block. ![](reduce_logo.png)
 * \ingroup BlockModule
 *
 * \par Overview
 * A <a href="http://en.wikipedia.org/wiki/Reduce_(higher-order_function)"><em>reduction</em></a> (or <em>fold</em>)
 * uses a binary combining operator to compute a single aggregate from a list of input elements.
 *
 * \par
 * Optionally, BlockReduce can be specialized by algorithm to accommodate different latency/throughput workload profiles:
 *   -# <b>cub::BLOCK_REDUCE_RAKING</b>.  An efficient "raking" reduction algorithm. [More...](\ref cub::BlockReduceAlgorithm)
 *   -# <b>cub::BLOCK_REDUCE_WARP_REDUCTIONS</b>.  A quick "tiled warp-reductions" reduction algorithm. [More...](\ref cub::BlockReduceAlgorithm)
 *
 * \tparam T                Data type being reduced
 * \tparam BLOCK_THREADS    The thread block size in threads
 * \tparam ALGORITHM        <b>[optional]</b> cub::BlockReduceAlgorithm enumerator specifying the underlying algorithm to use (default: cub::BLOCK_REDUCE_RAKING)
 *
 * \par Performance Considerations
 * - Very efficient (only one synchronization barrier).
 * - Zero bank conflicts for most types.
 * - Computation is slightly more efficient (i.e., having lower instruction overhead) for:
 *   - Summation (<b><em>vs.</em></b> generic reduction)
 *   - \p BLOCK_THREADS is a multiple of the architecture's warp size
 *   - Every thread has a valid input (i.e., full <b><em>vs.</em></b> partial-tiles)
 * - See cub::BlockReduceAlgorithm for performance details regarding algorithmic alternatives
 *
 * \par A Simple Example
 * \blockcollective{BlockReduce}
 * \par
 * The code snippet below illustrates a sum reduction of 512 integer items that
 * are partitioned in a [<em>blocked arrangement</em>](index.html#sec5sec4) across 128 threads
 * where each thread owns 4 consecutive items.
 * \par
 * \code
 * #include <cub/cub.cuh>
 *
 * __global__ void ExampleKernel(...)
 * {
 *     // Specialize BlockReduce for 128 threads on type int
 *     typedef cub::BlockReduce<int, 128> BlockReduce;
 *
 *     // Allocate shared memory for BlockReduce
 *     __shared__ typename BlockReduce::TempStorage temp_storage;
 *
 *     // Obtain a segment of consecutive items that are blocked across threads
 *     int thread_data[4];
 *     ...
 *
 *     // Compute the block-wide sum for thread0
 *     int aggregate = BlockReduce(temp_storage).Sum(thread_data);
 *
 * \endcode
 *
 */
template <
    typename                T,
    int                     BLOCK_THREADS,
    BlockReduceAlgorithm    ALGORITHM = BLOCK_REDUCE_RAKING>
class BlockReduce
{
private:

    /******************************************************************************
     * Constants and typedefs
     ******************************************************************************/

    /// Internal specialization.
    typedef typename If<(ALGORITHM == BLOCK_REDUCE_WARP_REDUCTIONS),
        BlockReduceWarpReductions<T, BLOCK_THREADS>,
        BlockReduceRaking<T, BLOCK_THREADS> >::Type InternalBlockReduce;

    /// Shared memory storage layout type for BlockReduce
    typedef typename InternalBlockReduce::TempStorage _TempStorage;


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

    /// Shared storage reference
    _TempStorage &temp_storage;

    /// Linear thread-id
    int linear_tid;


public:

    /// \smemstorage{BlockReduce}
    struct TempStorage : Uninitialized<_TempStorage> {};


    /******************************************************************//**
     * \name Collective constructors
     *********************************************************************/
    //@{

    /**
     * \brief Collective constructor for 1D thread blocks using a private static allocation of shared memory as temporary storage.  Threads are identified using <tt>threadIdx.x</tt>.
     */
    __device__ __forceinline__ BlockReduce()
    :
        temp_storage(PrivateStorage()),
        linear_tid(threadIdx.x)
    {}


    /**
     * \brief Collective constructor for 1D thread blocks using the specified memory allocation as temporary storage.  Threads are identified using <tt>threadIdx.x</tt>.
     */
    __device__ __forceinline__ BlockReduce(
        TempStorage &temp_storage)             ///< [in] Reference to memory allocation having layout type TempStorage
    :
        temp_storage(temp_storage.Alias()),
        linear_tid(threadIdx.x)
    {}


    /**
     * \brief Collective constructor using a private static allocation of shared memory as temporary storage.  Each thread is identified using the supplied linear thread identifier
     */
    __device__ __forceinline__ BlockReduce(
        int linear_tid)                        ///< [in] A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    :
        temp_storage(PrivateStorage()),
        linear_tid(linear_tid)
    {}


    /**
     * \brief Collective constructor using the specified memory allocation as temporary storage.  Each thread is identified using the supplied linear thread identifier.
     */
    __device__ __forceinline__ BlockReduce(
        TempStorage &temp_storage,             ///< [in] Reference to memory allocation having layout type TempStorage
        int linear_tid)                        ///< [in] <b>[optional]</b> A suitable 1D thread-identifier for the calling thread (e.g., <tt>(threadIdx.y * blockDim.x) + linear_tid</tt> for 2D thread blocks)
    :
        temp_storage(temp_storage.Alias()),
        linear_tid(linear_tid)
    {}



    //@}  end member group
    /******************************************************************//**
     * \name Generic reductions
     *********************************************************************/
    //@{


    /**
     * \brief Computes a block-wide reduction for thread<sub>0</sub> using the specified binary reduction functor.  Each thread contributes one input element.
     *
     * The return value is undefined in threads other than thread<sub>0</sub>.
     *
     * Supports non-commutative reduction operators.
     *
     * \smemreuse
     *
     * The code snippet below illustrates a max reduction of 128 integer items that
     * are partitioned across 128 threads.
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(...)
     * {
     *     // Specialize BlockReduce for 128 threads on type int
     *     typedef cub::BlockReduce<int, 128> BlockReduce;
     *
     *     // Allocate shared memory for BlockReduce
     *     __shared__ typename BlockReduce::TempStorage temp_storage;
     *
     *     // Each thread obtains an input item
     *     int thread_data;
     *     ...
     *
     *     // Compute the block-wide max for thread0
     *     int aggregate = BlockReduce(temp_storage).Reduce(thread_data, cub::Max());
     *
     * \endcode
     *
     * \tparam ReductionOp          <b>[inferred]</b> Binary reduction operator type having member <tt>T operator()(const T &a, const T &b)</tt>
     */
    template <typename ReductionOp>
    __device__ __forceinline__ T Reduce(
        T               input,                      ///< [in] Calling thread's input
        ReductionOp     reduction_op)               ///< [in] Binary reduction operator
    {
        return InternalBlockReduce(temp_storage, linear_tid).template Reduce<true>(input, BLOCK_THREADS, reduction_op);
    }


    /**
     * \brief Computes a block-wide reduction for thread<sub>0</sub> using the specified binary reduction functor.  Each thread contributes an array of consecutive input elements.
     *
     * The return value is undefined in threads other than thread<sub>0</sub>.
     *
     * Supports non-commutative reduction operators.
     *
     * \blocked
     *
     * \smemreuse
     *
     * The code snippet below illustrates a max reduction of 512 integer items that
     * are partitioned in a [<em>blocked arrangement</em>](index.html#sec5sec4) across 128 threads
     * where each thread owns 4 consecutive items.
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(...)
     * {
     *     // Specialize BlockReduce for 128 threads on type int
     *     typedef cub::BlockReduce<int, 128> BlockReduce;
     *
     *     // Allocate shared memory for BlockReduce
     *     __shared__ typename BlockReduce::TempStorage temp_storage;
     *
     *     // Obtain a segment of consecutive items that are blocked across threads
     *     int thread_data[4];
     *     ...
     *
     *     // Compute the block-wide max for thread0
     *     int aggregate = BlockReduce(temp_storage).Reduce(thread_data, cub::Max());
     *
     * \endcode
     *
     * \tparam ITEMS_PER_THREAD     <b>[inferred]</b> The number of consecutive items partitioned onto each thread.
     * \tparam ReductionOp          <b>[inferred]</b> Binary reduction operator type having member <tt>T operator()(const T &a, const T &b)</tt>
     */
    template <
        int ITEMS_PER_THREAD,
        typename ReductionOp>
    __device__ __forceinline__ T Reduce(
        T               (&inputs)[ITEMS_PER_THREAD],    ///< [in] Calling thread's input segment
        ReductionOp     reduction_op)                   ///< [in] Binary reduction operator
    {
        // Reduce partials
        T partial = ThreadReduce(inputs, reduction_op);
        return Reduce(partial, reduction_op);
    }


    /**
     * \brief Computes a block-wide reduction for thread<sub>0</sub> using the specified binary reduction functor.  The first \p num_valid threads each contribute one input element.
     *
     * The return value is undefined in threads other than thread<sub>0</sub>.
     *
     * Supports non-commutative reduction operators.
     *
     * \blocked
     *
     * \smemreuse
     *
     * The code snippet below illustrates a max reduction of a partially-full tile of integer items that
     * are partitioned across 128 threads.
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(int num_valid, ...)
     * {
     *     // Specialize BlockReduce for 128 threads on type int
     *     typedef cub::BlockReduce<int, 128> BlockReduce;
     *
     *     // Allocate shared memory for BlockReduce
     *     __shared__ typename BlockReduce::TempStorage temp_storage;
     *
     *     // Each thread obtains an input item
     *     int thread_data;
     *     if (threadIdx.x < num_valid) thread_data = ...
     *
     *     // Compute the block-wide max for thread0
     *     int aggregate = BlockReduce(temp_storage).Reduce(thread_data, cub::Max(), num_valid);
     *
     * \endcode
     *
     * \tparam ReductionOp          <b>[inferred]</b> Binary reduction operator type having member <tt>T operator()(const T &a, const T &b)</tt>
     */
    template <typename ReductionOp>
    __device__ __forceinline__ T Reduce(
        T                   input,                  ///< [in] Calling thread's input
        ReductionOp         reduction_op,           ///< [in] Binary reduction operator
        int                 num_valid)              ///< [in] Number of threads containing valid elements (may be less than BLOCK_THREADS)
    {
        // Determine if we scan skip bounds checking
        if (num_valid >= BLOCK_THREADS)
        {
            return InternalBlockReduce(temp_storage, linear_tid).template Reduce<true>(input, num_valid, reduction_op);
        }
        else
        {
            return InternalBlockReduce(temp_storage, linear_tid).template Reduce<false>(input, num_valid, reduction_op);
        }
    }


    //@}  end member group
    /******************************************************************//**
     * \name Summation reductions
     *********************************************************************/
    //@{


    /**
     * \brief Computes a block-wide reduction for thread<sub>0</sub> using addition (+) as the reduction operator.  Each thread contributes one input element.
     *
     * The return value is undefined in threads other than thread<sub>0</sub>.
     *
     * \smemreuse
     *
     * The code snippet below illustrates a sum reduction of 128 integer items that
     * are partitioned across 128 threads.
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(...)
     * {
     *     // Specialize BlockReduce for 128 threads on type int
     *     typedef cub::BlockReduce<int, 128> BlockReduce;
     *
     *     // Allocate shared memory for BlockReduce
     *     __shared__ typename BlockReduce::TempStorage temp_storage;
     *
     *     // Each thread obtains an input item
     *     int thread_data;
     *     ...
     *
     *     // Compute the block-wide sum for thread0
     *     int aggregate = BlockReduce(temp_storage).Sum(thread_data);
     *
     * \endcode
     *
     */
    __device__ __forceinline__ T Sum(
        T   input)                      ///< [in] Calling thread's input
    {
        return InternalBlockReduce(temp_storage, linear_tid).template Sum<true>(input, BLOCK_THREADS);
    }

    /**
     * \brief Computes a block-wide reduction for thread<sub>0</sub> using addition (+) as the reduction operator.  Each thread contributes an array of consecutive input elements.
     *
     * The return value is undefined in threads other than thread<sub>0</sub>.
     *
     * \smemreuse
     *
     * The code snippet below illustrates a sum reduction of 512 integer items that
     * are partitioned in a [<em>blocked arrangement</em>](index.html#sec5sec4) across 128 threads
     * where each thread owns 4 consecutive items.
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(...)
     * {
     *     // Specialize BlockReduce for 128 threads on type int
     *     typedef cub::BlockReduce<int, 128> BlockReduce;
     *
     *     // Allocate shared memory for BlockReduce
     *     __shared__ typename BlockReduce::TempStorage temp_storage;
     *
     *     // Obtain a segment of consecutive items that are blocked across threads
     *     int thread_data[4];
     *     ...
     *
     *     // Compute the block-wide sum for thread0
     *     int aggregate = BlockReduce(temp_storage).Sum(thread_data);
     *
     * \endcode
     *
     * \tparam ITEMS_PER_THREAD     <b>[inferred]</b> The number of consecutive items partitioned onto each thread.
     */
    template <int ITEMS_PER_THREAD>
    __device__ __forceinline__ T Sum(
        T   (&inputs)[ITEMS_PER_THREAD])    ///< [in] Calling thread's input segment
    {
        // Reduce partials
        T partial = ThreadReduce(inputs, cub::Sum());
        return Sum(partial);
    }


    /**
     * \brief Computes a block-wide reduction for thread<sub>0</sub> using addition (+) as the reduction operator.  The first \p num_valid threads each contribute one input element.
     *
     * The return value is undefined in threads other than thread<sub>0</sub>.
     *
     * \smemreuse
     *
     * The code snippet below illustrates a sum reduction of a partially-full tile of integer items that
     * are partitioned across 128 threads.
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(int num_valid, ...)
     * {
     *     // Specialize BlockReduce for 128 threads on type int
     *     typedef cub::BlockReduce<int, 128> BlockReduce;
     *
     *     // Allocate shared memory for BlockReduce
     *     __shared__ typename BlockReduce::TempStorage temp_storage;
     *
     *     // Each thread obtains an input item (up to num_items)
     *     int thread_data;
     *     if (threadIdx.x < num_valid)
     *         thread_data = ...
     *
     *     // Compute the block-wide sum for thread0
     *     int aggregate = BlockReduce(temp_storage).Sum(thread_data, num_valid);
     *
     * \endcode
     *
     */
    __device__ __forceinline__ T Sum(
        T   input,                  ///< [in] Calling thread's input
        int num_valid)              ///< [in] Number of threads containing valid elements (may be less than BLOCK_THREADS)
    {
        // Determine if we scan skip bounds checking
        if (num_valid >= BLOCK_THREADS)
        {
            return InternalBlockReduce(temp_storage, linear_tid).template Sum<true>(input, num_valid);
        }
        else
        {
            return InternalBlockReduce(temp_storage, linear_tid).template Sum<false>(input, num_valid);
        }
    }


    //@}  end member group
};

}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)

