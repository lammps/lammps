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
 * Utility types for device-wide scan
 */

#pragma once

#include <iterator>

#include "../../thread/thread_load.cuh"
#include "../../thread/thread_store.cuh"
#include "../../warp/warp_reduce.cuh"
#include "../../util_namespace.cuh"

/// Optional outer namespace(s)
CUB_NS_PREFIX

/// CUB namespace
namespace cub {


/**
 * Enumerations of tile status
 */
enum ScanTileStatus
{
    SCAN_TILE_OOB,          // Out-of-bounds (e.g., padding)
    SCAN_TILE_INVALID,      // Not yet processed
    SCAN_TILE_PARTIAL,      // Tile aggregate is available
    SCAN_TILE_PREFIX,       // Inclusive tile prefix is available
};


/**
 * Data type of tile status descriptor.
 *
 * Specialized for scan status and value types that can be combined into the same
 * machine word that can be read/written coherently in a single access.
 */
template <
    typename    T,
    bool        SINGLE_WORD = (PowerOfTwo<sizeof(T)>::VALUE && (sizeof(T) <= 8))>
struct ScanTileDescriptor
{
    // Status word type
    typedef typename If<(sizeof(T) == 8),
        long long,
        typename If<(sizeof(T) == 4),
            int,
            typename If<(sizeof(T) == 2),
                short,
                char>::Type>::Type>::Type StatusWord;

    // Vector word type
    typedef typename If<(sizeof(T) == 8),
        longlong2,
        typename If<(sizeof(T) == 4),
            int2,
            typename If<(sizeof(T) == 2),
                int,
                short>::Type>::Type>::Type VectorWord;

    T           value;
    StatusWord  status;

    static __device__ __forceinline__ void SetPrefix(ScanTileDescriptor *ptr, T prefix)
    {
        ScanTileDescriptor tile_descriptor;
        tile_descriptor.status = SCAN_TILE_PREFIX;
        tile_descriptor.value = prefix;

        VectorWord alias;
        *reinterpret_cast<ScanTileDescriptor*>(&alias) = tile_descriptor;
        ThreadStore<STORE_CG>(reinterpret_cast<VectorWord*>(ptr), alias);
    }

    static __device__ __forceinline__ void SetPartial(ScanTileDescriptor *ptr, T partial)
    {
        ScanTileDescriptor tile_descriptor;
        tile_descriptor.status = SCAN_TILE_PARTIAL;
        tile_descriptor.value = partial;

        VectorWord alias;
        *reinterpret_cast<ScanTileDescriptor*>(&alias) = tile_descriptor;
        ThreadStore<STORE_CG>(reinterpret_cast<VectorWord*>(ptr), alias);
    }

    static __device__ __forceinline__ void WaitForValid(
        ScanTileDescriptor    *ptr,
        int                     &status,
        T                       &value)
    {
        ScanTileDescriptor tile_descriptor;
        while (true)
        {
            VectorWord alias = ThreadLoad<LOAD_CG>(reinterpret_cast<VectorWord*>(ptr));

            tile_descriptor = *reinterpret_cast<ScanTileDescriptor*>(&alias);
            if (tile_descriptor.status != SCAN_TILE_INVALID) break;

            __threadfence_block();
        }

        status = tile_descriptor.status;
        value = tile_descriptor.value;
    }

};


/**
 * Data type of tile status descriptor.
 *
 * Specialized for scan status and value types that cannot fused into
 * the same machine word.
 */
template <typename T>
struct ScanTileDescriptor<T, false>
{
    T       prefix_value;
    T       partial_value;

    /// Workaround for the fact that win32 doesn't guarantee 16B alignment 16B values of T
    union
    {
        int                     status;
        Uninitialized<T>        padding;
    };

    static __device__ __forceinline__ void SetPrefix(ScanTileDescriptor *ptr, T prefix)
    {
        ThreadStore<STORE_CG>(&ptr->prefix_value, prefix);
        __threadfence_block();
//        __threadfence();        // __threadfence_block seems sufficient on current architectures to prevent reordeing
        ThreadStore<STORE_CG>(&ptr->status, (int) SCAN_TILE_PREFIX);

    }

    static __device__ __forceinline__ void SetPartial(ScanTileDescriptor *ptr, T partial)
    {
        ThreadStore<STORE_CG>(&ptr->partial_value, partial);
        __threadfence_block();
//        __threadfence();        // __threadfence_block seems sufficient on current architectures to prevent reordeing
        ThreadStore<STORE_CG>(&ptr->status, (int) SCAN_TILE_PARTIAL);
    }

    static __device__ __forceinline__ void WaitForValid(
        ScanTileDescriptor    *ptr,
        int                         &status,
        T                           &value)
    {
        while (true)
        {
            status = ThreadLoad<LOAD_CG>(&ptr->status);
            if (status != SCAN_TILE_INVALID) break;

            __threadfence_block();
        }

        value = (status == SCAN_TILE_PARTIAL) ?
            ThreadLoad<LOAD_CG>(&ptr->partial_value) :
            ThreadLoad<LOAD_CG>(&ptr->prefix_value);
    }
};


/**
 * Stateful prefix functor that provides the the running prefix for
 * the current tile by using the callback warp to wait on on
 * aggregates/prefixes from predecessor tiles to become available
 */
template <
    typename T,
    typename ScanOp>
struct DeviceScanBlockPrefixOp
{
    // Parameterized warp reduce
    typedef WarpReduce<T>                       WarpReduceT;

    // Storage type
    typedef typename WarpReduceT::TempStorage   _TempStorage;

    // Alias wrapper allowing storage to be unioned
    typedef Uninitialized<_TempStorage>         TempStorage;

    // Tile status descriptor type
    typedef ScanTileDescriptor<T>               ScanTileDescriptorT;

    // Fields
    ScanTileDescriptorT         *d_tile_status;     ///< Pointer to array of tile status
    _TempStorage                &temp_storage;      ///< Reference to a warp-reduction instance
    ScanOp                      scan_op;            ///< Binary scan operator
    int                         tile_idx;           ///< The current tile index
    T                           inclusive_prefix;   ///< Inclusive prefix for the tile

    // Constructor
    __device__ __forceinline__
    DeviceScanBlockPrefixOp(
        ScanTileDescriptorT     *d_tile_status,
        TempStorage             &temp_storage,
        ScanOp                  scan_op,
        int                     tile_idx) :
            d_tile_status(d_tile_status),
            temp_storage(temp_storage.Alias()),
            scan_op(scan_op),
            tile_idx(tile_idx) {}


    // Block until all predecessors within the specified window have non-invalid status
    __device__ __forceinline__
    void ProcessWindow(
        int                         predecessor_idx,
        int                         &predecessor_status,
        T                           &window_aggregate)
    {
        T value;
        ScanTileDescriptorT::WaitForValid(d_tile_status + predecessor_idx, predecessor_status, value);

        // Perform a segmented reduction to get the prefix for the current window
        int flag = (predecessor_status != SCAN_TILE_PARTIAL);
        window_aggregate = WarpReduceT(temp_storage).TailSegmentedReduce(value, flag, scan_op);
    }


    // Prefix functor (called by the first warp)
    __device__ __forceinline__
    T operator()(T block_aggregate)
    {
        // Update our status with our tile-aggregate
        if (threadIdx.x == 0)
        {
            ScanTileDescriptorT::SetPartial(d_tile_status + tile_idx, block_aggregate);
        }

        // Wait for the window of predecessor tiles to become valid
        int predecessor_idx = tile_idx - threadIdx.x - 1;
        int predecessor_status;
        T window_aggregate;
        ProcessWindow(predecessor_idx, predecessor_status, window_aggregate);

        // The exclusive tile prefix starts out as the current window aggregate
        T exclusive_prefix = window_aggregate;

        // Keep sliding the window back until we come across a tile whose inclusive prefix is known
        while (WarpAll(predecessor_status != SCAN_TILE_PREFIX))
        {
            predecessor_idx -= PtxArchProps::WARP_THREADS;

            // Update exclusive tile prefix with the window prefix
            ProcessWindow(predecessor_idx, predecessor_status, window_aggregate);
            exclusive_prefix = scan_op(window_aggregate, exclusive_prefix);
        }

        // Compute the inclusive tile prefix and update the status for this tile
        if (threadIdx.x == 0)
        {
            inclusive_prefix = scan_op(exclusive_prefix, block_aggregate);
            ScanTileDescriptorT::SetPrefix(
                d_tile_status + tile_idx,
                inclusive_prefix);
        }

        // Return exclusive_prefix
        return exclusive_prefix;
    }
};


// Running scan prefix callback type for single-block scans.
// Maintains a running prefix that can be applied to consecutive
// scan operations.
template <typename T>
struct RunningBlockPrefixOp
{
    // Running prefix
    T running_total;

    // Callback operator.
    __device__ T operator()(T block_aggregate)
    {
        T old_prefix = running_total;
        running_total += block_aggregate;
        return old_prefix;
    }
};


}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)

