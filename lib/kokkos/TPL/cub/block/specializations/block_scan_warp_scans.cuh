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
 * cub::BlockScanWarpscans provides warpscan-based variants of parallel prefix scan across a CUDA threadblock.
 */

#pragma once

#include "../../util_arch.cuh"
#include "../../warp/warp_scan.cuh"
#include "../../util_namespace.cuh"

/// Optional outer namespace(s)
CUB_NS_PREFIX

/// CUB namespace
namespace cub {

/**
 * \brief BlockScanWarpScans provides warpscan-based variants of parallel prefix scan across a CUDA threadblock.
 */
template <
    typename            T,
    int                 BLOCK_THREADS>
struct BlockScanWarpScans
{
    /// Constants
    enum
    {
        /// Number of active warps
        WARPS = (BLOCK_THREADS + PtxArchProps::WARP_THREADS - 1) / PtxArchProps::WARP_THREADS,
    };

    ///  WarpScan utility type
    typedef WarpScan<T, WARPS, PtxArchProps::WARP_THREADS> WarpScan;

    /// Shared memory storage layout type
    struct _TempStorage
    {
        typename WarpScan::TempStorage      warp_scan;                  ///< Buffer for warp-synchronous scan
        T                                   warp_aggregates[WARPS];     ///< Shared totals from each warp-synchronous scan
        T                                   block_prefix;               ///< Shared prefix for the entire threadblock
    };


    /// Alias wrapper allowing storage to be unioned
    struct TempStorage : Uninitialized<_TempStorage> {};


    // Thread fields
    _TempStorage &temp_storage;
    int linear_tid;
    int warp_id;
    int lane_id;


    /// Constructor
    __device__ __forceinline__ BlockScanWarpScans(
        TempStorage &temp_storage,
        int linear_tid)
    :
        temp_storage(temp_storage.Alias()),
        linear_tid(linear_tid),
        warp_id((BLOCK_THREADS <= PtxArchProps::WARP_THREADS) ?
            0 :
            linear_tid / PtxArchProps::WARP_THREADS),
        lane_id((BLOCK_THREADS <= PtxArchProps::WARP_THREADS) ?
            linear_tid :
            linear_tid % PtxArchProps::WARP_THREADS)
    {}


    /// Update the calling thread's partial reduction with the warp-wide aggregates from preceding warps.  Also returns block-wide aggregate in <em>thread</em><sub>0</sub>.
    template <typename ScanOp>
    __device__ __forceinline__ void ApplyWarpAggregates(
        T               &partial,           ///< [out] The calling thread's partial reduction
        ScanOp          scan_op,            ///< [in] Binary scan operator
        T               warp_aggregate,     ///< [in] <b>[<em>lane</em><sub>0</sub>s only]</b> Warp-wide aggregate reduction of input items
        T               &block_aggregate,   ///< [out] Threadblock-wide aggregate reduction of input items
        bool            lane_valid = true)  ///< [in] Whether or not the partial belonging to the current thread is valid
    {
        // Share lane aggregates
        temp_storage.warp_aggregates[warp_id] = warp_aggregate;

        __syncthreads();

        block_aggregate = temp_storage.warp_aggregates[0];

        #pragma unroll
        for (int WARP = 1; WARP < WARPS; WARP++)
        {
            if (warp_id == WARP)
            {
                partial = (lane_valid) ?
                    scan_op(block_aggregate, partial) :     // fold it in our valid partial
                    block_aggregate;                        // replace our invalid partial with the aggregate
            }

            block_aggregate = scan_op(block_aggregate, temp_storage.warp_aggregates[WARP]);
        }
    }


    /// Computes an exclusive threadblock-wide prefix scan using the specified binary \p scan_op functor.  Each thread contributes one input element.  Also provides every thread with the block-wide \p block_aggregate of all inputs.
    template <typename ScanOp>
    __device__ __forceinline__ void ExclusiveScan(
        T               input,              ///< [in] Calling thread's input items
        T               &output,            ///< [out] Calling thread's output items (may be aliased to \p input)
        const T         &identity,          ///< [in] Identity value
        ScanOp          scan_op,            ///< [in] Binary scan operator
        T               &block_aggregate)   ///< [out] Threadblock-wide aggregate reduction of input items
    {
        T warp_aggregate;
        WarpScan(temp_storage.warp_scan, warp_id, lane_id).ExclusiveScan(input, output, identity, scan_op, warp_aggregate);

        // Update outputs and block_aggregate with warp-wide aggregates
        ApplyWarpAggregates(output, scan_op, warp_aggregate, block_aggregate);
    }


    /// Computes an exclusive threadblock-wide prefix scan using the specified binary \p scan_op functor.  Each thread contributes one input element.  the call-back functor \p block_prefix_op is invoked by the first warp in the block, and the value returned by <em>lane</em><sub>0</sub> in that warp is used as the "seed" value that logically prefixes the threadblock's scan inputs.  Also provides every thread with the block-wide \p block_aggregate of all inputs.
    template <
        typename ScanOp,
        typename BlockPrefixOp>
    __device__ __forceinline__ void ExclusiveScan(
        T               input,                          ///< [in] Calling thread's input item
        T               &output,                        ///< [out] Calling thread's output item (may be aliased to \p input)
        T               identity,                       ///< [in] Identity value
        ScanOp          scan_op,                        ///< [in] Binary scan operator
        T               &block_aggregate,               ///< [out] Threadblock-wide aggregate reduction of input items (exclusive of the \p block_prefix_op value)
        BlockPrefixOp   &block_prefix_op)               ///< [in-out] <b>[<em>warp</em><sub>0</sub> only]</b> Call-back functor for specifying a threadblock-wide prefix to be applied to all inputs.
    {
        ExclusiveScan(input, output, identity, scan_op, block_aggregate);

        // Compute and share threadblock prefix
        if (warp_id == 0)
        {
            temp_storage.block_prefix = block_prefix_op(block_aggregate);
        }

        __syncthreads();

        // Incorporate threadblock prefix into outputs
        output = scan_op(temp_storage.block_prefix, output);
    }


    /// Computes an exclusive threadblock-wide prefix scan using the specified binary \p scan_op functor.  Each thread contributes one input element.  Also provides every thread with the block-wide \p block_aggregate of all inputs.  With no identity value, the output computed for <em>thread</em><sub>0</sub> is undefined.
    template <typename ScanOp>
    __device__ __forceinline__ void ExclusiveScan(
        T               input,                          ///< [in] Calling thread's input item
        T               &output,                        ///< [out] Calling thread's output item (may be aliased to \p input)
        ScanOp          scan_op,                        ///< [in] Binary scan operator
        T               &block_aggregate)               ///< [out] Threadblock-wide aggregate reduction of input items
    {
        T warp_aggregate;
        WarpScan(temp_storage.warp_scan, warp_id, lane_id).ExclusiveScan(input, output, scan_op, warp_aggregate);

        // Update outputs and block_aggregate with warp-wide aggregates
        ApplyWarpAggregates(output, scan_op, warp_aggregate, block_aggregate, (lane_id > 0));
    }


    /// Computes an exclusive threadblock-wide prefix scan using the specified binary \p scan_op functor.  Each thread contributes one input element.  the call-back functor \p block_prefix_op is invoked by the first warp in the block, and the value returned by <em>lane</em><sub>0</sub> in that warp is used as the "seed" value that logically prefixes the threadblock's scan inputs.  Also provides every thread with the block-wide \p block_aggregate of all inputs.
    template <
        typename ScanOp,
        typename BlockPrefixOp>
    __device__ __forceinline__ void ExclusiveScan(
        T               input,                          ///< [in] Calling thread's input item
        T               &output,                        ///< [out] Calling thread's output item (may be aliased to \p input)
        ScanOp          scan_op,                        ///< [in] Binary scan operator
        T               &block_aggregate,               ///< [out] Threadblock-wide aggregate reduction of input items (exclusive of the \p block_prefix_op value)
        BlockPrefixOp   &block_prefix_op)               ///< [in-out] <b>[<em>warp</em><sub>0</sub> only]</b> Call-back functor for specifying a threadblock-wide prefix to be applied to all inputs.
    {
        ExclusiveScan(input, output, scan_op, block_aggregate);

        // Compute and share threadblock prefix
        if (warp_id == 0)
        {
            temp_storage.block_prefix = block_prefix_op(block_aggregate);
        }

        __syncthreads();

        // Incorporate threadblock prefix into outputs
        output = (linear_tid == 0) ?
            temp_storage.block_prefix :
            scan_op(temp_storage.block_prefix, output);
    }


    /// Computes an exclusive threadblock-wide prefix scan using addition (+) as the scan operator.  Each thread contributes one input element.  Also provides every thread with the block-wide \p block_aggregate of all inputs.
    __device__ __forceinline__ void ExclusiveSum(
        T               input,                          ///< [in] Calling thread's input item
        T               &output,                        ///< [out] Calling thread's output item (may be aliased to \p input)
        T               &block_aggregate)               ///< [out] Threadblock-wide aggregate reduction of input items
    {
        T warp_aggregate;
        WarpScan(temp_storage.warp_scan, warp_id, lane_id).ExclusiveSum(input, output, warp_aggregate);

        // Update outputs and block_aggregate with warp-wide aggregates from lane-0s
        ApplyWarpAggregates(output, Sum(), warp_aggregate, block_aggregate);
    }


    /// Computes an exclusive threadblock-wide prefix scan using addition (+) as the scan operator.  Each thread contributes one input element.  Instead of using 0 as the threadblock-wide prefix, the call-back functor \p block_prefix_op is invoked by the first warp in the block, and the value returned by <em>lane</em><sub>0</sub> in that warp is used as the "seed" value that logically prefixes the threadblock's scan inputs.  Also provides every thread with the block-wide \p block_aggregate of all inputs.
    template <typename BlockPrefixOp>
    __device__ __forceinline__ void ExclusiveSum(
        T               input,                          ///< [in] Calling thread's input item
        T               &output,                        ///< [out] Calling thread's output item (may be aliased to \p input)
        T               &block_aggregate,               ///< [out] Threadblock-wide aggregate reduction of input items (exclusive of the \p block_prefix_op value)
        BlockPrefixOp   &block_prefix_op)               ///< [in-out] <b>[<em>warp</em><sub>0</sub> only]</b> Call-back functor for specifying a threadblock-wide prefix to be applied to all inputs.
    {
        ExclusiveSum(input, output, block_aggregate);

        // Compute and share threadblock prefix
        if (warp_id == 0)
        {
            temp_storage.block_prefix = block_prefix_op(block_aggregate);
        }

        __syncthreads();

        // Incorporate threadblock prefix into outputs
        Sum scan_op;
        output = scan_op(temp_storage.block_prefix, output);
    }


    /// Computes an inclusive threadblock-wide prefix scan using the specified binary \p scan_op functor.  Each thread contributes one input element.  Also provides every thread with the block-wide \p block_aggregate of all inputs.
    template <typename ScanOp>
    __device__ __forceinline__ void InclusiveScan(
        T               input,                          ///< [in] Calling thread's input item
        T               &output,                        ///< [out] Calling thread's output item (may be aliased to \p input)
        ScanOp          scan_op,                        ///< [in] Binary scan operator
        T               &block_aggregate)               ///< [out] Threadblock-wide aggregate reduction of input items
    {
        T warp_aggregate;
        WarpScan(temp_storage.warp_scan, warp_id, lane_id).InclusiveScan(input, output, scan_op, warp_aggregate);

        // Update outputs and block_aggregate with warp-wide aggregates from lane-0s
        ApplyWarpAggregates(output, scan_op, warp_aggregate, block_aggregate);

    }


    /// Computes an inclusive threadblock-wide prefix scan using the specified binary \p scan_op functor.  Each thread contributes one input element.  the call-back functor \p block_prefix_op is invoked by the first warp in the block, and the value returned by <em>lane</em><sub>0</sub> in that warp is used as the "seed" value that logically prefixes the threadblock's scan inputs.  Also provides every thread with the block-wide \p block_aggregate of all inputs.
    template <
        typename ScanOp,
        typename BlockPrefixOp>
    __device__ __forceinline__ void InclusiveScan(
        T               input,                          ///< [in] Calling thread's input item
        T               &output,                        ///< [out] Calling thread's output item (may be aliased to \p input)
        ScanOp          scan_op,                        ///< [in] Binary scan operator
        T               &block_aggregate,               ///< [out] Threadblock-wide aggregate reduction of input items (exclusive of the \p block_prefix_op value)
        BlockPrefixOp   &block_prefix_op)               ///< [in-out] <b>[<em>warp</em><sub>0</sub> only]</b> Call-back functor for specifying a threadblock-wide prefix to be applied to all inputs.
    {
        InclusiveScan(input, output, scan_op, block_aggregate);

        // Compute and share threadblock prefix
        if (warp_id == 0)
        {
            temp_storage.block_prefix = block_prefix_op(block_aggregate);
        }

        __syncthreads();

        // Incorporate threadblock prefix into outputs
        output = scan_op(temp_storage.block_prefix, output);
    }


    /// Computes an inclusive threadblock-wide prefix scan using the specified binary \p scan_op functor.  Each thread contributes one input element.  Also provides every thread with the block-wide \p block_aggregate of all inputs.
    __device__ __forceinline__ void InclusiveSum(
        T               input,                          ///< [in] Calling thread's input item
        T               &output,                        ///< [out] Calling thread's output item (may be aliased to \p input)
        T               &block_aggregate)               ///< [out] Threadblock-wide aggregate reduction of input items
    {
        T warp_aggregate;
        WarpScan(temp_storage.warp_scan, warp_id, lane_id).InclusiveSum(input, output, warp_aggregate);

        // Update outputs and block_aggregate with warp-wide aggregates from lane-0s
        ApplyWarpAggregates(output, Sum(), warp_aggregate, block_aggregate);
    }


    /// Computes an inclusive threadblock-wide prefix scan using the specified binary \p scan_op functor.  Each thread contributes one input element.  Instead of using 0 as the threadblock-wide prefix, the call-back functor \p block_prefix_op is invoked by the first warp in the block, and the value returned by <em>lane</em><sub>0</sub> in that warp is used as the "seed" value that logically prefixes the threadblock's scan inputs.  Also provides every thread with the block-wide \p block_aggregate of all inputs.
    template <typename BlockPrefixOp>
    __device__ __forceinline__ void InclusiveSum(
        T               input,                          ///< [in] Calling thread's input item
        T               &output,                        ///< [out] Calling thread's output item (may be aliased to \p input)
        T               &block_aggregate,               ///< [out] Threadblock-wide aggregate reduction of input items (exclusive of the \p block_prefix_op value)
        BlockPrefixOp   &block_prefix_op)               ///< [in-out] <b>[<em>warp</em><sub>0</sub> only]</b> Call-back functor for specifying a threadblock-wide prefix to be applied to all inputs.
    {
        InclusiveSum(input, output, block_aggregate);

        // Compute and share threadblock prefix
        if (warp_id == 0)
        {
            temp_storage.block_prefix = block_prefix_op(block_aggregate);
        }

        __syncthreads();

        // Incorporate threadblock prefix into outputs
        Sum scan_op;
        output = scan_op(temp_storage.block_prefix, output);
    }

};


}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)

