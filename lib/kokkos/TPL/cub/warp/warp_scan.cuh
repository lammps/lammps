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
 * The cub::WarpScan class provides [<em>collective</em>](index.html#sec0) methods for computing a parallel prefix scan of items partitioned across CUDA warp threads.
 */

#pragma once

#include "specializations/warp_scan_shfl.cuh"
#include "specializations/warp_scan_smem.cuh"
#include "../thread/thread_operators.cuh"
#include "../util_arch.cuh"
#include "../util_type.cuh"
#include "../util_namespace.cuh"

/// Optional outer namespace(s)
CUB_NS_PREFIX

/// CUB namespace
namespace cub {

/**
 * \addtogroup WarpModule
 * @{
 */

/**
 * \brief The WarpScan class provides [<em>collective</em>](index.html#sec0) methods for computing a parallel prefix scan of items partitioned across CUDA warp threads.  ![](warp_scan_logo.png)
 *
 * \par Overview
 * Given a list of input elements and a binary reduction operator, a [<em>prefix scan</em>](http://en.wikipedia.org/wiki/Prefix_sum)
 * produces an output list where each element is computed to be the reduction
 * of the elements occurring earlier in the input list.  <em>Prefix sum</em>
 * connotes a prefix scan with the addition operator. The term \em inclusive indicates
 * that the <em>i</em><sup>th</sup> output reduction incorporates the <em>i</em><sup>th</sup> input.
 * The term \em exclusive indicates the <em>i</em><sup>th</sup> input is not incorporated into
 * the <em>i</em><sup>th</sup> output reduction.
 *
 * \tparam T                        The scan input/output element type
 * \tparam LOGICAL_WARPS            <b>[optional]</b> The number of "logical" warps performing concurrent warp scans. Default is 1.
 * \tparam LOGICAL_WARP_THREADS     <b>[optional]</b> The number of threads per "logical" warp (may be less than the number of hardware warp threads).  Default is the warp size associated with the CUDA Compute Capability targeted by the compiler (e.g., 32 threads for SM20).
 *
 * \par Simple Examples
 * \warpcollective{WarpScan}
 * \par
 * The code snippet below illustrates four concurrent warp prefix sums within a block of
 * 128 threads (one per each of the 32-thread warps).
 * \par
 * \code
 * #include <cub/cub.cuh>
 *
 * __global__ void ExampleKernel(...)
 * {
 *     // Specialize WarpScan for 4 warps on type int
 *     typedef cub::WarpScan<int, 4> WarpScan;
 *
 *     // Allocate shared memory for WarpScan
 *     __shared__ typename WarpScan::TempStorage temp_storage;
 *
 *     // Obtain one input item per thread
 *     int thread_data = ...
 *
 *     // Compute warp-wide prefix sums
 *     WarpScan(temp_storage).ExclusiveSum(thread_data, thread_data);
 *
 * \endcode
 * \par
 * Suppose the set of input \p thread_data across the block of threads is <tt>1, 1, 1, 1, ...</tt>.
 * The corresponding output \p thread_data in each of the four warps of threads will be
 * <tt>0, 1, 2, 3, ..., 31</tt>.
 *
 * \par
 * The code snippet below illustrates a single warp prefix sum within a block of
 * 128 threads.
 * \par
 * \code
 * #include <cub/cub.cuh>
 *
 * __global__ void ExampleKernel(...)
 * {
 *     // Specialize WarpScan for one warp on type int
 *     typedef cub::WarpScan<int, 1> WarpScan;
 *
 *     // Allocate shared memory for WarpScan
 *     __shared__ typename WarpScan::TempStorage temp_storage;
 *     ...
 *
 *     // Only the first warp performs a prefix sum
 *     if (threadIdx.x < 32)
 *     {
 *         // Obtain one input item per thread
 *         int thread_data = ...
 *
 *         // Compute warp-wide prefix sums
 *         WarpScan(temp_storage).ExclusiveSum(thread_data, thread_data);
 *
 * \endcode
 * \par
 * Suppose the set of input \p thread_data across the warp of threads is <tt>1, 1, 1, 1, ...</tt>.
 * The corresponding output \p thread_data will be <tt>0, 1, 2, 3, ..., 31</tt>.
 *
 * \par Usage and Performance Considerations
 * - Supports "logical" warps smaller than the physical warp size (e.g., a logical warp of 8 threads)
 * - The number of entrant threads must be an multiple of \p LOGICAL_WARP_THREADS
 * - Warp scans are concurrent if more than one warp is participating
 * - Uses special instructions when applicable (e.g., warp \p SHFL)
 * - Uses synchronization-free communication between warp lanes when applicable
 * - Zero bank conflicts for most types.
 * - Computation is slightly more efficient (i.e., having lower instruction overhead) for:
 *     - Summation (<b><em>vs.</em></b> generic scan)
 *     - The architecture's warp size is a whole multiple of \p LOGICAL_WARP_THREADS
 *
 */
template <
    typename    T,
    int         LOGICAL_WARPS           = 1,
    int         LOGICAL_WARP_THREADS    = PtxArchProps::WARP_THREADS>
class WarpScan
{
private:

    /******************************************************************************
     * Constants and typedefs
     ******************************************************************************/

    enum
    {
        POW_OF_TWO = ((LOGICAL_WARP_THREADS & (LOGICAL_WARP_THREADS - 1)) == 0),
    };

    /// Internal specialization.  Use SHFL-based reduction if (architecture is >= SM30) and ((only one logical warp) or (LOGICAL_WARP_THREADS is a power-of-two))
    typedef typename If<(CUB_PTX_ARCH >= 300) && ((LOGICAL_WARPS == 1) || POW_OF_TWO),
        WarpScanShfl<T, LOGICAL_WARPS, LOGICAL_WARP_THREADS>,
        WarpScanSmem<T, LOGICAL_WARPS, LOGICAL_WARP_THREADS> >::Type InternalWarpScan;

    /// Shared memory storage layout type for WarpScan
    typedef typename InternalWarpScan::TempStorage _TempStorage;


    /******************************************************************************
     * Thread fields
     ******************************************************************************/

    /// Shared storage reference
    _TempStorage &temp_storage;

    /// Warp ID
    int warp_id;

    /// Lane ID
    int lane_id;


    /******************************************************************************
     * Utility methods
     ******************************************************************************/

    /// Internal storage allocator
    __device__ __forceinline__ _TempStorage& PrivateStorage()
    {
        __shared__ TempStorage private_storage;
        return private_storage;
    }


public:

    /// \smemstorage{WarpScan}
    struct TempStorage : Uninitialized<_TempStorage> {};


    /******************************************************************//**
     * \name Collective constructors
     *********************************************************************/
    //@{

    /**
     * \brief Collective constructor for 1D thread blocks using a private static allocation of shared memory as temporary storage.  Logical warp and lane identifiers are constructed from <tt>threadIdx.x</tt>.
     */
    __device__ __forceinline__ WarpScan()
    :
        temp_storage(PrivateStorage()),
        warp_id((LOGICAL_WARPS == 1) ?
            0 :
            threadIdx.x / LOGICAL_WARP_THREADS),
        lane_id(((LOGICAL_WARPS == 1) || (LOGICAL_WARP_THREADS == PtxArchProps::WARP_THREADS)) ?
            LaneId() :
            threadIdx.x % LOGICAL_WARP_THREADS)
    {}


    /**
     * \brief Collective constructor for 1D thread blocks using the specified memory allocation as temporary storage.  Logical warp and lane identifiers are constructed from <tt>threadIdx.x</tt>.
     */
    __device__ __forceinline__ WarpScan(
        TempStorage &temp_storage)             ///< [in] Reference to memory allocation having layout type TempStorage
    :
        temp_storage(temp_storage.Alias()),
        warp_id((LOGICAL_WARPS == 1) ?
            0 :
            threadIdx.x / LOGICAL_WARP_THREADS),
        lane_id(((LOGICAL_WARPS == 1) || (LOGICAL_WARP_THREADS == PtxArchProps::WARP_THREADS)) ?
            LaneId() :
            threadIdx.x % LOGICAL_WARP_THREADS)
    {}


    /**
     * \brief Collective constructor using a private static allocation of shared memory as temporary storage.  Threads are identified using the given warp and lane identifiers.
     */
    __device__ __forceinline__ WarpScan(
        int warp_id,                           ///< [in] A suitable warp membership identifier
        int lane_id)                           ///< [in] A lane identifier within the warp
    :
        temp_storage(PrivateStorage()),
        warp_id(warp_id),
        lane_id(lane_id)
    {}


    /**
     * \brief Collective constructor using the specified memory allocation as temporary storage.  Threads are identified using the given warp and lane identifiers.
     */
    __device__ __forceinline__ WarpScan(
        TempStorage &temp_storage,             ///< [in] Reference to memory allocation having layout type TempStorage
        int warp_id,                           ///< [in] A suitable warp membership identifier
        int lane_id)                           ///< [in] A lane identifier within the warp
    :
        temp_storage(temp_storage.Alias()),
        warp_id(warp_id),
        lane_id(lane_id)
    {}


    //@}  end member group
    /******************************************************************//**
     * \name Inclusive prefix sums
     *********************************************************************/
    //@{


    /**
     * \brief Computes an inclusive prefix sum in each logical warp.
     *
     * \smemreuse
     *
     * The code snippet below illustrates four concurrent warp-wide inclusive prefix sums within a block of
     * 128 threads (one per each of the 32-thread warps).
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(...)
     * {
     *     // Specialize WarpScan for 4 warps on type int
     *     typedef cub::WarpScan<int, 4> WarpScan;
     *
     *     // Allocate shared memory for WarpScan
     *     __shared__ typename WarpScan::TempStorage temp_storage;
     *
     *     // Obtain one input item per thread
     *     int thread_data = ...
     *
     *     // Compute inclusive warp-wide prefix sums
     *     WarpScan(temp_storage).InclusiveSum(thread_data, thread_data);
     *
     * \endcode
     * \par
     * Suppose the set of input \p thread_data across the block of threads is <tt>1, 1, 1, 1, ...</tt>.
     * The corresponding output \p thread_data in each of the four warps of threads will be
     * <tt>1, 2, 3, ..., 32</tt>.
     */
    __device__ __forceinline__ void InclusiveSum(
        T               input,              ///< [in] Calling thread's input item.
        T               &output)            ///< [out] Calling thread's output item.  May be aliased with \p input.
    {
        InternalWarpScan(temp_storage, warp_id, lane_id).InclusiveSum(input, output);
    }


    /**
     * \brief Computes an inclusive prefix sum in each logical warp.  Also provides every thread with the warp-wide \p warp_aggregate of all inputs.
     *
     * The \p warp_aggregate is undefined in threads other than <em>warp-lane</em><sub>0</sub>.
     *
     * \smemreuse
     *
     * The code snippet below illustrates four concurrent warp-wide inclusive prefix sums within a block of
     * 128 threads (one per each of the 32-thread warps).
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(...)
     * {
     *     // Specialize WarpScan for 4 warps on type int
     *     typedef cub::WarpScan<int, 4> WarpScan;
     *
     *     // Allocate shared memory for WarpScan
     *     __shared__ typename WarpScan::TempStorage temp_storage;
     *
     *     // Obtain one input item per thread
     *     int thread_data = ...
     *
     *     // Compute inclusive warp-wide prefix sums
     *     int warp_aggregate;
     *     WarpScan(temp_storage).InclusiveSum(thread_data, thread_data, warp_aggregate);
     *
     * \endcode
     * \par
     * Suppose the set of input \p thread_data across the block of threads is <tt>1, 1, 1, 1, ...</tt>.
     * The corresponding output \p thread_data in each of the four warps of threads will be
     * <tt>1, 2, 3, ..., 32</tt>.  Furthermore, \p warp_aggregate for all threads in all warps will be \p 32.
     */
    __device__ __forceinline__ void InclusiveSum(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        T               &warp_aggregate)    ///< [out] Warp-wide aggregate reduction of input items.
    {
        InternalWarpScan(temp_storage, warp_id, lane_id).InclusiveSum(input, output, warp_aggregate);
    }


    /**
     * \brief Computes an inclusive prefix sum in each logical warp.  Instead of using 0 as the warp-wide prefix, the call-back functor \p warp_prefix_op is invoked to provide the "seed" value that logically prefixes the warp's scan inputs.  Also provides every thread with the warp-wide \p warp_aggregate of all inputs.
     *
     * The \p warp_aggregate is undefined in threads other than <em>warp-lane</em><sub>0</sub>.
     *
     * The \p warp_prefix_op functor must implement a member function <tt>T operator()(T warp_aggregate)</tt>.
     * The functor's input parameter \p warp_aggregate is the same value also returned by the scan operation.
     * The functor will be invoked by the entire warp of threads, however only the return value from
     * <em>lane</em><sub>0</sub> is applied as the threadblock-wide prefix.  Can be stateful.
     *
     * \smemreuse
     *
     * The code snippet below illustrates a single thread block of 32 threads (one warp) that progressively
     * computes an inclusive prefix sum over multiple "tiles" of input using a
     * prefix functor to maintain a running total between block-wide scans.  Each tile consists
     * of 32 integer items that are partitioned across the warp.
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * // A stateful callback functor that maintains a running prefix to be applied
     * // during consecutive scan operations.
     * struct WarpPrefixOp
     * {
     *     // Running prefix
     *     int running_total;
     *
     *     // Constructor
     *     __device__ WarpPrefixOp(int running_total) : running_total(running_total) {}
     *
     *     // Callback operator to be entered by the entire warp. Lane-0 is responsible
     *     // for returning a value for seeding the warp-wide scan.
     *     __device__ int operator()(int warp_aggregate)
     *     {
     *         int old_prefix = running_total;
     *         running_total += warp_aggregate;
     *         return old_prefix;
     *     }
     * };
     *
     * __global__ void ExampleKernel(int *d_data, int num_items, ...)
     * {
     *     // Specialize WarpScan for one warp
     *     typedef cub::WarpScan<int, 1> WarpScan;
     *
     *     // Allocate shared memory for WarpScan
     *     __shared__ typename WarpScan::TempStorage temp_storage;
     *
     *     // Initialize running total
     *     WarpPrefixOp prefix_op(0);
     *
     *     // Have the warp iterate over segments of items
     *     for (int block_offset = 0; block_offset < num_items; block_offset += 32)
     *     {
     *         // Load a segment of consecutive items
     *         int thread_data = d_data[block_offset];
     *
     *         // Collectively compute the warp-wide inclusive prefix sum
     *         int warp_aggregate;
     *         WarpScan(temp_storage).InclusiveSum(
     *             thread_data, thread_data, warp_aggregate, prefix_op);
     *
     *         // Store scanned items to output segment
     *         d_data[block_offset] = thread_data;
     *     }
     * \endcode
     * \par
     * Suppose the input \p d_data is <tt>1, 1, 1, 1, 1, 1, 1, 1, ...</tt>.
     * The corresponding output for the first segment will be <tt>1, 2, 3, ..., 32</tt>.
     * The output for the second segment will be <tt>33, 34, 35, ..., 64</tt>.  Furthermore,
     * the value \p 32 will be stored in \p warp_aggregate for all threads after each scan.
     *
     * \tparam WarpPrefixOp                 <b>[inferred]</b> Call-back functor type having member <tt>T operator()(T warp_aggregate)</tt>
     */
    template <typename WarpPrefixOp>
    __device__ __forceinline__ void InclusiveSum(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        T               &warp_aggregate,    ///< [out] <b>[<em>warp-lane</em><sub>0</sub> only]</b> Warp-wide aggregate reduction of input items, exclusive of the \p warp_prefix_op value
        WarpPrefixOp    &warp_prefix_op)    ///< [in-out] <b>[<em>warp-lane</em><sub>0</sub> only]</b> Call-back functor for specifying a warp-wide prefix to be applied to all inputs.
    {
        // Compute inclusive warp scan
        InclusiveSum(input, output, warp_aggregate);

        // Compute warp-wide prefix from aggregate, then broadcast to other lanes
        T prefix;
        prefix = warp_prefix_op(warp_aggregate);
        prefix = InternalWarpScan(temp_storage, warp_id, lane_id).Broadcast(prefix, 0);

        // Update output
        output = prefix + output;
    }

    //@}  end member group

private:

    /// Computes an exclusive prefix sum in each logical warp.
    __device__ __forceinline__ void ExclusiveSum(T input, T &output, Int2Type<true> is_primitive)
    {
        // Compute exclusive warp scan from inclusive warp scan
        T inclusive;
        InclusiveSum(input, inclusive);
        output = inclusive - input;
    }

    /// Computes an exclusive prefix sum in each logical warp.  Specialized for non-primitive types.
    __device__ __forceinline__ void ExclusiveSum(T input, T &output, Int2Type<false> is_primitive)
    {
        // Delegate to regular scan for non-primitive types (because we won't be able to use subtraction)
        T identity = T();
        ExclusiveScan(input, output, identity, Sum());
    }

    /// Computes an exclusive prefix sum in each logical warp.  Also provides every thread with the warp-wide \p warp_aggregate of all inputs.
    __device__ __forceinline__ void ExclusiveSum(T input, T &output, T &warp_aggregate, Int2Type<true> is_primitive)
    {
        // Compute exclusive warp scan from inclusive warp scan
        T inclusive;
        InclusiveSum(input, inclusive, warp_aggregate);
        output = inclusive - input;
    }

    /// Computes an exclusive prefix sum in each logical warp.  Also provides every thread with the warp-wide \p warp_aggregate of all inputs.  Specialized for non-primitive types.
    __device__ __forceinline__ void ExclusiveSum(T input, T &output, T &warp_aggregate, Int2Type<false> is_primitive)
    {
        // Delegate to regular scan for non-primitive types (because we won't be able to use subtraction)
        T identity = T();
        ExclusiveScan(input, output, identity, Sum(), warp_aggregate);
    }

    /// Computes an exclusive prefix sum in each logical warp.  Instead of using 0 as the warp-wide prefix, the call-back functor \p warp_prefix_op is invoked to provide the "seed" value that logically prefixes the warp's scan inputs.  Also provides every thread with the warp-wide \p warp_aggregate of all inputs.
    template <typename WarpPrefixOp>
    __device__ __forceinline__ void ExclusiveSum(T input, T &output, T &warp_aggregate, WarpPrefixOp &warp_prefix_op, Int2Type<true> is_primitive)
    {
        // Compute exclusive warp scan from inclusive warp scan
        T inclusive;
        InclusiveSum(input, inclusive, warp_aggregate, warp_prefix_op);
        output = inclusive - input;
    }

    /// Computes an exclusive prefix sum in each logical warp.  Instead of using 0 as the warp-wide prefix, the call-back functor \p warp_prefix_op is invoked to provide the "seed" value that logically prefixes the warp's scan inputs.  Also provides every thread with the warp-wide \p warp_aggregate of all inputs.  Specialized for non-primitive types.
    template <typename WarpPrefixOp>
    __device__ __forceinline__ void ExclusiveSum(T input, T &output, T &warp_aggregate, WarpPrefixOp &warp_prefix_op, Int2Type<false> is_primitive)
    {
        // Delegate to regular scan for non-primitive types (because we won't be able to use subtraction)
        T identity = T();
        ExclusiveScan(input, output, identity, Sum(), warp_aggregate, warp_prefix_op);
    }

public:


    /******************************************************************//**
     * \name Exclusive prefix sums
     *********************************************************************/
    //@{


    /**
     * \brief Computes an exclusive prefix sum in each logical warp.
     *
     * This operation assumes the value of obtained by the <tt>T</tt>'s default
     * constructor (or by zero-initialization if no user-defined default
     * constructor exists) is suitable as the identity value "zero" for
     * addition.
     *
     * \smemreuse
     *
     * The code snippet below illustrates four concurrent warp-wide exclusive prefix sums within a block of
     * 128 threads (one per each of the 32-thread warps).
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(...)
     * {
     *     // Specialize WarpScan for 4 warps on type int
     *     typedef cub::WarpScan<int, 4> WarpScan;
     *
     *     // Allocate shared memory for WarpScan
     *     __shared__ typename WarpScan::TempStorage temp_storage;
     *
     *     // Obtain one input item per thread
     *     int thread_data = ...
     *
     *     // Compute exclusive warp-wide prefix sums
     *     WarpScan(temp_storage).ExclusiveSum(thread_data, thread_data);
     *
     * \endcode
     * \par
     * Suppose the set of input \p thread_data across the block of threads is <tt>1, 1, 1, 1, ...</tt>.
     * The corresponding output \p thread_data in each of the four warps of threads will be
     * <tt>0, 1, 2, ..., 31</tt>.
     *
     */
    __device__ __forceinline__ void ExclusiveSum(
        T               input,              ///< [in] Calling thread's input item.
        T               &output)            ///< [out] Calling thread's output item.  May be aliased with \p input.
    {
        ExclusiveSum(input, output, Int2Type<Traits<T>::PRIMITIVE>());
    }


    /**
     * \brief Computes an exclusive prefix sum in each logical warp.  Also provides every thread with the warp-wide \p warp_aggregate of all inputs.
     *
     * This operation assumes the value of obtained by the <tt>T</tt>'s default
     * constructor (or by zero-initialization if no user-defined default
     * constructor exists) is suitable as the identity value "zero" for
     * addition.
     *
     * \smemreuse
     *
     * The code snippet below illustrates four concurrent warp-wide exclusive prefix sums within a block of
     * 128 threads (one per each of the 32-thread warps).
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(...)
     * {
     *     // Specialize WarpScan for 4 warps on type int
     *     typedef cub::WarpScan<int, 4> WarpScan;
     *
     *     // Allocate shared memory for WarpScan
     *     __shared__ typename WarpScan::TempStorage temp_storage;
     *
     *     // Obtain one input item per thread
     *     int thread_data = ...
     *
     *     // Compute exclusive warp-wide prefix sums
     *     int warp_aggregate;
     *     WarpScan(temp_storage).ExclusiveSum(thread_data, thread_data, warp_aggregate);
     *
     * \endcode
     * \par
     * Suppose the set of input \p thread_data across the block of threads is <tt>1, 1, 1, 1, ...</tt>.
     * The corresponding output \p thread_data in each of the four warps of threads will be
     * <tt>0, 1, 2, ..., 31</tt>.  Furthermore, \p warp_aggregate for all threads in all warps will be \p 32.
     */
    __device__ __forceinline__ void ExclusiveSum(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        T               &warp_aggregate)    ///< [out] Warp-wide aggregate reduction of input items.
    {
        ExclusiveSum(input, output, warp_aggregate, Int2Type<Traits<T>::PRIMITIVE>());
    }


    /**
     * \brief Computes an exclusive prefix sum in each logical warp.  Instead of using 0 as the warp-wide prefix, the call-back functor \p warp_prefix_op is invoked to provide the "seed" value that logically prefixes the warp's scan inputs.  Also provides every thread with the warp-wide \p warp_aggregate of all inputs.
     *
     * This operation assumes the value of obtained by the <tt>T</tt>'s default
     * constructor (or by zero-initialization if no user-defined default
     * constructor exists) is suitable as the identity value "zero" for
     * addition.
     *
     * The \p warp_prefix_op functor must implement a member function <tt>T operator()(T warp_aggregate)</tt>.
     * The functor's input parameter \p warp_aggregate is the same value also returned by the scan operation.
     * The functor will be invoked by the entire warp of threads, however only the return value from
     * <em>lane</em><sub>0</sub> is applied as the threadblock-wide prefix.  Can be stateful.
     *
     * \smemreuse
     *
     * The code snippet below illustrates a single thread block of 32 threads (one warp) that progressively
     * computes an exclusive prefix sum over multiple "tiles" of input using a
     * prefix functor to maintain a running total between block-wide scans.  Each tile consists
     * of 32 integer items that are partitioned across the warp.
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * // A stateful callback functor that maintains a running prefix to be applied
     * // during consecutive scan operations.
     * struct WarpPrefixOp
     * {
     *     // Running prefix
     *     int running_total;
     *
     *     // Constructor
     *     __device__ WarpPrefixOp(int running_total) : running_total(running_total) {}
     *
     *     // Callback operator to be entered by the entire warp. Lane-0 is responsible
     *     // for returning a value for seeding the warp-wide scan.
     *     __device__ int operator()(int warp_aggregate)
     *     {
     *         int old_prefix = running_total;
     *         running_total += warp_aggregate;
     *         return old_prefix;
     *     }
     * };
     *
     * __global__ void ExampleKernel(int *d_data, int num_items, ...)
     * {
     *     // Specialize WarpScan for one warp
     *     typedef cub::WarpScan<int, 1> WarpScan;
     *
     *     // Allocate shared memory for WarpScan
     *     __shared__ typename WarpScan::TempStorage temp_storage;
     *
     *     // Initialize running total
     *     WarpPrefixOp prefix_op(0);
     *
     *     // Have the warp iterate over segments of items
     *     for (int block_offset = 0; block_offset < num_items; block_offset += 32)
     *     {
     *         // Load a segment of consecutive items
     *         int thread_data = d_data[block_offset];
     *
     *         // Collectively compute the warp-wide exclusive prefix sum
     *         int warp_aggregate;
     *         WarpScan(temp_storage).ExclusiveSum(
     *             thread_data, thread_data, warp_aggregate, prefix_op);
     *
     *         // Store scanned items to output segment
     *         d_data[block_offset] = thread_data;
     *     }
     * \endcode
     * \par
     * Suppose the input \p d_data is <tt>1, 1, 1, 1, 1, 1, 1, 1, ...</tt>.
     * The corresponding output for the first segment will be <tt>0, 1, 2, ..., 31</tt>.
     * The output for the second segment will be <tt>32, 33, 34, ..., 63</tt>.  Furthermore,
     * the value \p 32 will be stored in \p warp_aggregate for all threads after each scan.
     *
     * \tparam WarpPrefixOp                 <b>[inferred]</b> Call-back functor type having member <tt>T operator()(T warp_aggregate)</tt>
     */
    template <typename WarpPrefixOp>
    __device__ __forceinline__ void ExclusiveSum(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        T               &warp_aggregate,    ///< [out] <b>[<em>warp-lane</em><sub>0</sub> only]</b> Warp-wide aggregate reduction of input items (exclusive of the \p warp_prefix_op value).
        WarpPrefixOp    &warp_prefix_op)    ///< [in-out] <b>[<em>warp-lane</em><sub>0</sub> only]</b> Call-back functor for specifying a warp-wide prefix to be applied to all inputs.
    {
        ExclusiveSum(input, output, warp_aggregate, warp_prefix_op, Int2Type<Traits<T>::PRIMITIVE>());
    }


    //@}  end member group
    /******************************************************************//**
     * \name Inclusive prefix scans
     *********************************************************************/
    //@{

    /**
     * \brief Computes an inclusive prefix sum using the specified binary scan functor in each logical warp.
     *
     * Supports non-commutative scan operators.
     *
     * \smemreuse
     *
     * The code snippet below illustrates four concurrent warp-wide inclusive prefix max scans within a block of
     * 128 threads (one per each of the 32-thread warps).
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(...)
     * {
     *     // Specialize WarpScan for 4 warps on type int
     *     typedef cub::WarpScan<int, 4> WarpScan;
     *
     *     // Allocate shared memory for WarpScan
     *     __shared__ typename WarpScan::TempStorage temp_storage;
     *
     *     // Obtain one input item per thread
     *     int thread_data = ...
     *
     *     // Compute inclusive warp-wide prefix max scans
     *     WarpScan(temp_storage).InclusiveScan(thread_data, thread_data, cub::Max());
     *
     * \endcode
     * \par
     * Suppose the set of input \p thread_data across the block of threads is <tt>0, -1, 2, -3, ..., 126, -127</tt>.
     * The corresponding output \p thread_data in the first warp would be
     * <tt>0, 0, 2, 2, ..., 30, 30</tt>, the output for the second warp would be <tt>32, 32, 34, 34, ..., 62, 62</tt>, etc.
     *
     * \tparam ScanOp     <b>[inferred]</b> Binary scan operator type having member <tt>T operator()(const T &a, const T &b)</tt>
     */
    template <typename ScanOp>
    __device__ __forceinline__ void InclusiveScan(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        ScanOp          scan_op)            ///< [in] Binary scan operator
    {
        InternalWarpScan(temp_storage, warp_id, lane_id).InclusiveScan(input, output, scan_op);
    }


    /**
     * \brief Computes an inclusive prefix sum using the specified binary scan functor in each logical warp.  Also provides every thread with the warp-wide \p warp_aggregate of all inputs.
     *
     * Supports non-commutative scan operators.
     *
     * \smemreuse
     *
     * The code snippet below illustrates four concurrent warp-wide inclusive prefix max scans within a block of
     * 128 threads (one per each of the 32-thread warps).
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(...)
     * {
     *     // Specialize WarpScan for 4 warps on type int
     *     typedef cub::WarpScan<int, 4> WarpScan;
     *
     *     // Allocate shared memory for WarpScan
     *     __shared__ typename WarpScan::TempStorage temp_storage;
     *
     *     // Obtain one input item per thread
     *     int thread_data = ...
     *
     *     // Compute inclusive warp-wide prefix max scans
     *     int warp_aggregate;
     *     WarpScan(temp_storage).InclusiveScan(
     *         thread_data, thread_data, cub::Max(), warp_aggregate);
     *
     * \endcode
     * \par
     * Suppose the set of input \p thread_data across the block of threads is <tt>0, -1, 2, -3, ..., 126, -127</tt>.
     * The corresponding output \p thread_data in the first warp would be
     * <tt>0, 0, 2, 2, ..., 30, 30</tt>, the output for the second warp would be <tt>32, 32, 34, 34, ..., 62, 62</tt>, etc.
     * Furthermore, \p warp_aggregate would be assigned \p 30 for threads in the first warp, \p 62 for threads
     * in the second warp, etc.
     *
     * \tparam ScanOp     <b>[inferred]</b> Binary scan operator type having member <tt>T operator()(const T &a, const T &b)</tt>
     */
    template <typename ScanOp>
    __device__ __forceinline__ void InclusiveScan(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        ScanOp          scan_op,            ///< [in] Binary scan operator
        T               &warp_aggregate)    ///< [out] Warp-wide aggregate reduction of input items.
    {
        InternalWarpScan(temp_storage, warp_id, lane_id).InclusiveScan(input, output, scan_op, warp_aggregate);
    }


    /**
     * \brief Computes an inclusive prefix sum using the specified binary scan functor in each logical warp.  The call-back functor \p warp_prefix_op is invoked to provide the "seed" value that logically prefixes the warp's scan inputs.  Also provides every thread with the warp-wide \p warp_aggregate of all inputs.
     *
     * The \p warp_prefix_op functor must implement a member function <tt>T operator()(T warp_aggregate)</tt>.
     * The functor's input parameter \p warp_aggregate is the same value also returned by the scan operation.
     * The functor will be invoked by the entire warp of threads, however only the return value from
     * <em>lane</em><sub>0</sub> is applied as the threadblock-wide prefix.  Can be stateful.
     *
     * Supports non-commutative scan operators.
     *
     * \smemreuse
     *
     * The code snippet below illustrates a single thread block of 32 threads (one warp) that progressively
     * computes an inclusive prefix max scan over multiple "tiles" of input using a
     * prefix functor to maintain a running total between block-wide scans.  Each tile consists
     * of 32 integer items that are partitioned across the warp.
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * // A stateful callback functor that maintains a running prefix to be applied
     * // during consecutive scan operations.
     * struct WarpPrefixOp
     * {
     *     // Running prefix
     *     int running_total;
     *
     *     // Constructor
     *     __device__ WarpPrefixOp(int running_total) : running_total(running_total) {}
     *
     *     // Callback operator to be entered by the entire warp. Lane-0 is responsible
     *     // for returning a value for seeding the warp-wide scan.
     *     __device__ int operator()(int warp_aggregate)
     *     {
     *         int old_prefix = running_total;
     *         running_total = (warp_aggregate > old_prefix) ? warp_aggregate : old_prefix;
     *         return old_prefix;
     *     }
     * };
     *
     * __global__ void ExampleKernel(int *d_data, int num_items, ...)
     * {
     *     // Specialize WarpScan for one warp
     *     typedef cub::WarpScan<int, 1> WarpScan;
     *
     *     // Allocate shared memory for WarpScan
     *     __shared__ typename WarpScan::TempStorage temp_storage;
     *
     *     // Initialize running total
     *     WarpPrefixOp prefix_op(0);
     *
     *     // Have the warp iterate over segments of items
     *     for (int block_offset = 0; block_offset < num_items; block_offset += 32)
     *     {
     *         // Load a segment of consecutive items
     *         int thread_data = d_data[block_offset];
     *
     *         // Collectively compute the warp-wide inclusive prefix max scan
     *         int warp_aggregate;
     *         WarpScan(temp_storage).InclusiveScan(
     *             thread_data, thread_data, cub::Max(), warp_aggregate, prefix_op);
     *
     *         // Store scanned items to output segment
     *         d_data[block_offset] = thread_data;
     *     }
     * \endcode
     * \par
     * Suppose the input \p d_data is <tt>0, -1, 2, -3, 4, -5, ...</tt>.
     * The corresponding output for the first segment will be <tt>0, 0, 2, 2, ..., 30, 30</tt>.
     * The output for the second segment will be <tt>32, 32, 34, 34, ..., 62, 62</tt>.  Furthermore,
     * \p block_aggregate will be assigned \p 30 in all threads after the first scan, assigned \p 62 after the second
     * scan, etc.
     *
     * \tparam ScanOp                       <b>[inferred]</b> Binary scan operator type having member <tt>T operator()(const T &a, const T &b)</tt>
     * \tparam WarpPrefixOp                 <b>[inferred]</b> Call-back functor type having member <tt>T operator()(T warp_aggregate)</tt>
     */
    template <
        typename ScanOp,
        typename WarpPrefixOp>
    __device__ __forceinline__ void InclusiveScan(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        ScanOp          scan_op,            ///< [in] Binary scan operator
        T               &warp_aggregate,    ///< [out] <b>[<em>warp-lane</em><sub>0</sub> only]</b> Warp-wide aggregate reduction of input items (exclusive of the \p warp_prefix_op value).
        WarpPrefixOp    &warp_prefix_op)    ///< [in-out] <b>[<em>warp-lane</em><sub>0</sub> only]</b> Call-back functor for specifying a warp-wide prefix to be applied to all inputs.
    {
        // Compute inclusive warp scan
        InclusiveScan(input, output, scan_op, warp_aggregate);

        // Compute warp-wide prefix from aggregate, then broadcast to other lanes
        T prefix;
        prefix = warp_prefix_op(warp_aggregate);
        prefix = InternalWarpScan(temp_storage, warp_id, lane_id).Broadcast(prefix, 0);

        // Update output
        output = scan_op(prefix, output);
    }


    //@}  end member group
    /******************************************************************//**
     * \name Exclusive prefix scans
     *********************************************************************/
    //@{

    /**
     * \brief Computes an exclusive prefix scan using the specified binary scan functor in each logical warp.
     *
     * Supports non-commutative scan operators.
     *
     * \smemreuse
     *
     * The code snippet below illustrates four concurrent warp-wide exclusive prefix max scans within a block of
     * 128 threads (one per each of the 32-thread warps).
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(...)
     * {
     *     // Specialize WarpScan for 4 warps on type int
     *     typedef cub::WarpScan<int, 4> WarpScan;
     *
     *     // Allocate shared memory for WarpScan
     *     __shared__ typename WarpScan::TempStorage temp_storage;
     *
     *     // Obtain one input item per thread
     *     int thread_data = ...
     *
     *     // Compute exclusive warp-wide prefix max scans
     *     WarpScan(temp_storage).ExclusiveScan(thread_data, thread_data, INT_MIN, cub::Max());
     *
     * \endcode
     * \par
     * Suppose the set of input \p thread_data across the block of threads is <tt>0, -1, 2, -3, ..., 126, -127</tt>.
     * The corresponding output \p thread_data in the first warp would be
     * <tt>INT_MIN, 0, 0, 2, ..., 28, 30</tt>, the output for the second warp would be <tt>30, 32, 32, 34, ..., 60, 62</tt>, etc.
     *
     * \tparam ScanOp     <b>[inferred]</b> Binary scan operator type having member <tt>T operator()(const T &a, const T &b)</tt>
     */
    template <typename ScanOp>
    __device__ __forceinline__ void ExclusiveScan(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        T               identity,           ///< [in] Identity value
        ScanOp          scan_op)            ///< [in] Binary scan operator
    {
        InternalWarpScan(temp_storage, warp_id, lane_id).ExclusiveScan(input, output, identity, scan_op);
    }


    /**
     * \brief Computes an exclusive prefix scan using the specified binary scan functor in each logical warp.  Also provides every thread with the warp-wide \p warp_aggregate of all inputs.
     *
     * Supports non-commutative scan operators.
     *
     * \smemreuse
     *
     * The code snippet below illustrates four concurrent warp-wide exclusive prefix max scans within a block of
     * 128 threads (one per each of the 32-thread warps).
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(...)
     * {
     *     // Specialize WarpScan for 4 warps on type int
     *     typedef cub::WarpScan<int, 4> WarpScan;
     *
     *     // Allocate shared memory for WarpScan
     *     __shared__ typename WarpScan::TempStorage temp_storage;
     *
     *     // Obtain one input item per thread
     *     int thread_data = ...
     *
     *     // Compute exclusive warp-wide prefix max scans
     *     WarpScan(temp_storage).ExclusiveScan(thread_data, thread_data, INT_MIN, cub::Max());
     *
     * \endcode
     * \par
     * Suppose the set of input \p thread_data across the block of threads is <tt>0, -1, 2, -3, ..., 126, -127</tt>.
     * The corresponding output \p thread_data in the first warp would be
     * <tt>INT_MIN, 0, 0, 2, ..., 28, 30</tt>, the output for the second warp would be <tt>30, 32, 32, 34, ..., 60, 62</tt>, etc.
     * Furthermore, \p warp_aggregate would be assigned \p 30 for threads in the first warp, \p 62 for threads
     * in the second warp, etc.
     *
     * \tparam ScanOp     <b>[inferred]</b> Binary scan operator type having member <tt>T operator()(const T &a, const T &b)</tt>
     */
    template <typename ScanOp>
    __device__ __forceinline__ void ExclusiveScan(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        T               identity,           ///< [in] Identity value
        ScanOp          scan_op,            ///< [in] Binary scan operator
        T               &warp_aggregate)    ///< [out] Warp-wide aggregate reduction of input items.
    {
        InternalWarpScan(temp_storage, warp_id, lane_id).ExclusiveScan(input, output, identity, scan_op, warp_aggregate);
    }


    /**
     * \brief Computes an exclusive prefix scan using the specified binary scan functor in each logical warp.  The call-back functor \p warp_prefix_op is invoked to provide the "seed" value that logically prefixes the warp's scan inputs.  Also provides every thread with the warp-wide \p warp_aggregate of all inputs.
     *
     * The \p warp_prefix_op functor must implement a member function <tt>T operator()(T warp_aggregate)</tt>.
     * The functor's input parameter \p warp_aggregate is the same value also returned by the scan operation.
     * The functor will be invoked by the entire warp of threads, however only the return value from
     * <em>lane</em><sub>0</sub> is applied as the threadblock-wide prefix.  Can be stateful.
     *
     * Supports non-commutative scan operators.
     *
     * \smemreuse
     *
     * The code snippet below illustrates a single thread block of 32 threads (one warp) that progressively
     * computes an exclusive prefix max scan over multiple "tiles" of input using a
     * prefix functor to maintain a running total between block-wide scans.  Each tile consists
     * of 32 integer items that are partitioned across the warp.
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * // A stateful callback functor that maintains a running prefix to be applied
     * // during consecutive scan operations.
     * struct WarpPrefixOp
     * {
     *     // Running prefix
     *     int running_total;
     *
     *     // Constructor
     *     __device__ WarpPrefixOp(int running_total) : running_total(running_total) {}
     *
     *     // Callback operator to be entered by the entire warp. Lane-0 is responsible
     *     // for returning a value for seeding the warp-wide scan.
     *     __device__ int operator()(int warp_aggregate)
     *     {
     *         int old_prefix = running_total;
     *         running_total = (warp_aggregate > old_prefix) ? warp_aggregate : old_prefix;
     *         return old_prefix;
     *     }
     * };
     *
     * __global__ void ExampleKernel(int *d_data, int num_items, ...)
     * {
     *     // Specialize WarpScan for one warp
     *     typedef cub::WarpScan<int, 1> WarpScan;
     *
     *     // Allocate shared memory for WarpScan
     *     __shared__ typename WarpScan::TempStorage temp_storage;
     *
     *     // Initialize running total
     *     WarpPrefixOp prefix_op(INT_MIN);
     *
     *     // Have the warp iterate over segments of items
     *     for (int block_offset = 0; block_offset < num_items; block_offset += 32)
     *     {
     *         // Load a segment of consecutive items
     *         int thread_data = d_data[block_offset];
     *
     *         // Collectively compute the warp-wide exclusive prefix max scan
     *         int warp_aggregate;
     *         WarpScan(temp_storage).ExclusiveScan(
     *             thread_data, thread_data, INT_MIN, cub::Max(), warp_aggregate, prefix_op);
     *
     *         // Store scanned items to output segment
     *         d_data[block_offset] = thread_data;
     *     }
     * \endcode
     * \par
     * Suppose the input \p d_data is <tt>0, -1, 2, -3, 4, -5, ...</tt>.
     * The corresponding output for the first segment will be <tt>INT_MIN, 0, 0, 2, ..., 28, 30</tt>.
     * The output for the second segment will be <tt>30, 32, 32, 34, ..., 60, 62</tt>.  Furthermore,
     * \p block_aggregate will be assigned \p 30 in all threads after the first scan, assigned \p 62 after the second
     * scan, etc.
     *
     * \tparam ScanOp                       <b>[inferred]</b> Binary scan operator type having member <tt>T operator()(const T &a, const T &b)</tt>
     * \tparam WarpPrefixOp                 <b>[inferred]</b> Call-back functor type having member <tt>T operator()(T warp_aggregate)</tt>
     */
    template <
        typename ScanOp,
        typename WarpPrefixOp>
    __device__ __forceinline__ void ExclusiveScan(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        T               identity,           ///< [in] Identity value
        ScanOp          scan_op,            ///< [in] Binary scan operator
        T               &warp_aggregate,    ///< [out] <b>[<em>warp-lane</em><sub>0</sub> only]</b> Warp-wide aggregate reduction of input items (exclusive of the \p warp_prefix_op value).
        WarpPrefixOp    &warp_prefix_op)    ///< [in-out] <b>[<em>warp-lane</em><sub>0</sub> only]</b> Call-back functor for specifying a warp-wide prefix to be applied to all inputs.
    {
        // Exclusive warp scan
        ExclusiveScan(input, output, identity, scan_op, warp_aggregate);

        // Compute warp-wide prefix from aggregate, then broadcast to other lanes
        T prefix = warp_prefix_op(warp_aggregate);
        prefix = InternalWarpScan(temp_storage, warp_id, lane_id).Broadcast(prefix, 0);

        // Update output
        output = (lane_id == 0) ?
            prefix :
            scan_op(prefix, output);
    }


    //@}  end member group
    /******************************************************************//**
     * \name Identityless exclusive prefix scans
     *********************************************************************/
    //@{


    /**
     * \brief Computes an exclusive prefix scan using the specified binary scan functor in each logical warp.  Because no identity value is supplied, the \p output computed for <em>warp-lane</em><sub>0</sub> is undefined.
     *
     * Supports non-commutative scan operators.
     *
     * \smemreuse
     *
     * The code snippet below illustrates four concurrent warp-wide exclusive prefix max scans within a block of
     * 128 threads (one per each of the 32-thread warps).
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(...)
     * {
     *     // Specialize WarpScan for 4 warps on type int
     *     typedef cub::WarpScan<int, 4> WarpScan;
     *
     *     // Allocate shared memory for WarpScan
     *     __shared__ typename WarpScan::TempStorage temp_storage;
     *
     *     // Obtain one input item per thread
     *     int thread_data = ...
     *
     *     // Compute exclusive warp-wide prefix max scans
     *     WarpScan(temp_storage).ExclusiveScan(thread_data, thread_data, cub::Max());
     *
     * \endcode
     * \par
     * Suppose the set of input \p thread_data across the block of threads is <tt>0, -1, 2, -3, ..., 126, -127</tt>.
     * The corresponding output \p thread_data in the first warp would be
     * <tt>?, 0, 0, 2, ..., 28, 30</tt>, the output for the second warp would be <tt>?, 32, 32, 34, ..., 60, 62</tt>, etc.
     * (The output \p thread_data in each warp lane0 is undefined.)
     *
     * \tparam ScanOp     <b>[inferred]</b> Binary scan operator type having member <tt>T operator()(const T &a, const T &b)</tt>
     */
    template <typename ScanOp>
    __device__ __forceinline__ void ExclusiveScan(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        ScanOp          scan_op)            ///< [in] Binary scan operator
    {
        InternalWarpScan(temp_storage, warp_id, lane_id).ExclusiveScan(input, output, scan_op);
    }


    /**
     * \brief Computes an exclusive prefix scan using the specified binary scan functor in each logical warp.  Because no identity value is supplied, the \p output computed for <em>warp-lane</em><sub>0</sub> is undefined.  Also provides every thread with the warp-wide \p warp_aggregate of all inputs.
     *
     * Supports non-commutative scan operators.
     *
     * \smemreuse
     *
     * The code snippet below illustrates four concurrent warp-wide exclusive prefix max scans within a block of
     * 128 threads (one per each of the 32-thread warps).
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * __global__ void ExampleKernel(...)
     * {
     *     // Specialize WarpScan for 4 warps on type int
     *     typedef cub::WarpScan<int, 4> WarpScan;
     *
     *     // Allocate shared memory for WarpScan
     *     __shared__ typename WarpScan::TempStorage temp_storage;
     *
     *     // Obtain one input item per thread
     *     int thread_data = ...
     *
     *     // Compute exclusive warp-wide prefix max scans
     *     WarpScan(temp_storage).ExclusiveScan(thread_data, thread_data, cub::Max());
     *
     * \endcode
     * \par
     * Suppose the set of input \p thread_data across the block of threads is <tt>0, -1, 2, -3, ..., 126, -127</tt>.
     * The corresponding output \p thread_data in the first warp would be
     * <tt>?, 0, 0, 2, ..., 28, 30</tt>, the output for the second warp would be <tt>?, 32, 32, 34, ..., 60, 62</tt>, etc.
     * (The output \p thread_data in each warp lane0 is undefined.)  Furthermore, \p warp_aggregate would be assigned \p 30 for threads in the first warp, \p 62 for threads
     * in the second warp, etc.
     *
     * \tparam ScanOp     <b>[inferred]</b> Binary scan operator type having member <tt>T operator()(const T &a, const T &b)</tt>
     */
    template <typename ScanOp>
    __device__ __forceinline__ void ExclusiveScan(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        ScanOp          scan_op,            ///< [in] Binary scan operator
        T               &warp_aggregate)    ///< [out] Warp-wide aggregate reduction of input items.
    {
        InternalWarpScan(temp_storage, warp_id, lane_id).ExclusiveScan(input, output, scan_op, warp_aggregate);
    }


    /**
     * \brief Computes an exclusive prefix scan using the specified binary scan functor in each logical warp.  The \p warp_prefix_op value from thread-thread-lane<sub>0</sub> is applied to all scan outputs.  Also computes the warp-wide \p warp_aggregate of all inputs for thread-thread-lane<sub>0</sub>.
     *
     * The \p warp_prefix_op functor must implement a member function <tt>T operator()(T warp_aggregate)</tt>.
     * The functor's input parameter \p warp_aggregate is the same value also returned by the scan operation.
     * The functor will be invoked by the entire warp of threads, however only the return value from
     * <em>lane</em><sub>0</sub> is applied as the threadblock-wide prefix.  Can be stateful.
     *
     * Supports non-commutative scan operators.
     *
     * \smemreuse
     *
     * The code snippet below illustrates a single thread block of 32 threads (one warp) that progressively
     * computes an exclusive prefix max scan over multiple "tiles" of input using a
     * prefix functor to maintain a running total between block-wide scans.  Each tile consists
     * of 32 integer items that are partitioned across the warp.
     * \par
     * \code
     * #include <cub/cub.cuh>
     *
     * // A stateful callback functor that maintains a running prefix to be applied
     * // during consecutive scan operations.
     * struct WarpPrefixOp
     * {
     *     // Running prefix
     *     int running_total;
     *
     *     // Constructor
     *     __device__ WarpPrefixOp(int running_total) : running_total(running_total) {}
     *
     *     // Callback operator to be entered by the entire warp. Lane-0 is responsible
     *     // for returning a value for seeding the warp-wide scan.
     *     __device__ int operator()(int warp_aggregate)
     *     {
     *         int old_prefix = running_total;
     *         running_total = (warp_aggregate > old_prefix) ? warp_aggregate : old_prefix;
     *         return old_prefix;
     *     }
     * };
     *
     * __global__ void ExampleKernel(int *d_data, int num_items, ...)
     * {
     *     // Specialize WarpScan for one warp
     *     typedef cub::WarpScan<int, 1> WarpScan;
     *
     *     // Allocate shared memory for WarpScan
     *     __shared__ typename WarpScan::TempStorage temp_storage;
     *
     *     // Initialize running total
     *     WarpPrefixOp prefix_op(INT_MIN);
     *
     *     // Have the warp iterate over segments of items
     *     for (int block_offset = 0; block_offset < num_items; block_offset += 32)
     *     {
     *         // Load a segment of consecutive items
     *         int thread_data = d_data[block_offset];
     *
     *         // Collectively compute the warp-wide exclusive prefix max scan
     *         int warp_aggregate;
     *         WarpScan(temp_storage).ExclusiveScan(
     *             thread_data, thread_data, INT_MIN, cub::Max(), warp_aggregate, prefix_op);
     *
     *         // Store scanned items to output segment
     *         d_data[block_offset] = thread_data;
     *     }
     * \endcode
     * \par
     * Suppose the input \p d_data is <tt>0, -1, 2, -3, 4, -5, ...</tt>.
     * The corresponding output for the first segment will be <tt>INT_MIN, 0, 0, 2, ..., 28, 30</tt>.
     * The output for the second segment will be <tt>30, 32, 32, 34, ..., 60, 62</tt>.  Furthermore,
     * \p block_aggregate will be assigned \p 30 in all threads after the first scan, assigned \p 62 after the second
     * scan, etc.
     *
     * \tparam ScanOp                       <b>[inferred]</b> Binary scan operator type having member <tt>T operator()(const T &a, const T &b)</tt>
     * \tparam WarpPrefixOp                 <b>[inferred]</b> Call-back functor type having member <tt>T operator()(T warp_aggregate)</tt>
     */
    template <
        typename ScanOp,
        typename WarpPrefixOp>
    __device__ __forceinline__ void ExclusiveScan(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        ScanOp          scan_op,            ///< [in] Binary scan operator
        T               &warp_aggregate,    ///< [out] <b>[<em>warp-lane</em><sub>0</sub> only]</b> Warp-wide aggregate reduction of input items (exclusive of the \p warp_prefix_op value).
        WarpPrefixOp    &warp_prefix_op)    ///< [in-out] <b>[<em>warp-lane</em><sub>0</sub> only]</b> Call-back functor for specifying a warp-wide prefix to be applied to all inputs.
    {
        // Exclusive warp scan
        ExclusiveScan(input, output, scan_op, warp_aggregate);

        // Compute warp-wide prefix from aggregate, then broadcast to other lanes
        T prefix = warp_prefix_op(warp_aggregate);
        prefix = InternalWarpScan(temp_storage, warp_id, lane_id).Broadcast(prefix, 0);

        // Update output with prefix
        output = (lane_id == 0) ?
            prefix :
            scan_op(prefix, output);
    }

    //@}  end member group
};

/** @} */       // end group WarpModule

}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)
