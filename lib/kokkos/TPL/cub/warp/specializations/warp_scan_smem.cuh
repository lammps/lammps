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
 * cub::WarpScanSmem provides smem-based variants of parallel prefix scan across CUDA warps.
 */

#pragma once

#include "../../thread/thread_operators.cuh"
#include "../../thread/thread_load.cuh"
#include "../../thread/thread_store.cuh"
#include "../../util_type.cuh"
#include "../../util_namespace.cuh"

/// Optional outer namespace(s)
CUB_NS_PREFIX

/// CUB namespace
namespace cub {

/**
 * \brief WarpScanSmem provides smem-based variants of parallel prefix scan across CUDA warps.
 */
template <
    typename    T,                      ///< Data type being scanned
    int         LOGICAL_WARPS,          ///< Number of logical warps entrant
    int         LOGICAL_WARP_THREADS>   ///< Number of threads per logical warp
struct WarpScanSmem
{
    /******************************************************************************
     * Constants and typedefs
     ******************************************************************************/

    enum
    {
        /// The number of warp scan steps
        STEPS = Log2<LOGICAL_WARP_THREADS>::VALUE,

        /// The number of threads in half a warp
        HALF_WARP_THREADS = 1 << (STEPS - 1),

        /// The number of shared memory elements per warp
        WARP_SMEM_ELEMENTS =  LOGICAL_WARP_THREADS + HALF_WARP_THREADS,
    };


    /// Shared memory storage layout type (1.5 warps-worth of elements for each warp)
    typedef T _TempStorage[LOGICAL_WARPS][WARP_SMEM_ELEMENTS];

    // Alias wrapper allowing storage to be unioned
    struct TempStorage : Uninitialized<_TempStorage> {};


    /******************************************************************************
     * Thread fields
     ******************************************************************************/

    _TempStorage     &temp_storage;
    unsigned int    warp_id;
    unsigned int    lane_id;


    /******************************************************************************
     * Construction
     ******************************************************************************/

    /// Constructor
    __device__ __forceinline__ WarpScanSmem(
        TempStorage     &temp_storage,
        int             warp_id,
        int             lane_id)
    :
        temp_storage(temp_storage.Alias()),
        warp_id(warp_id),
        lane_id(lane_id)
    {}


    /******************************************************************************
     * Operation
     ******************************************************************************/

    /// Initialize identity padding (specialized for operations that have identity)
    __device__ __forceinline__ void InitIdentity(Int2Type<true> has_identity)
    {
        T identity = T();
        ThreadStore<STORE_VOLATILE>(&temp_storage[warp_id][lane_id], identity);
    }


    /// Initialize identity padding (specialized for operations without identity)
    __device__ __forceinline__ void InitIdentity(Int2Type<false> has_identity)
    {}


    /// Basic inclusive scan iteration(template unrolled, base-case specialization)
    template <
        bool HAS_IDENTITY,
        typename ScanOp>
    __device__ __forceinline__ void ScanStep(
        T               &partial,
        ScanOp          scan_op,
        Int2Type<STEPS>  step)
    {}


    /// Basic inclusive scan iteration (template unrolled, inductive-case specialization)
    template <
        bool        HAS_IDENTITY,
        int         STEP,
        typename    ScanOp>
    __device__ __forceinline__ void ScanStep(
        T               &partial,
        ScanOp          scan_op,
        Int2Type<STEP>  step)
    {
        const int OFFSET = 1 << STEP;

        // Share partial into buffer
        ThreadStore<STORE_VOLATILE>(&temp_storage[warp_id][HALF_WARP_THREADS + lane_id], partial);

        // Update partial if addend is in range
        if (HAS_IDENTITY || (lane_id >= OFFSET))
        {
            T addend = ThreadLoad<LOAD_VOLATILE>(&temp_storage[warp_id][HALF_WARP_THREADS + lane_id - OFFSET]);
            partial = scan_op(addend, partial);
        }

        ScanStep<HAS_IDENTITY>(partial, scan_op, Int2Type<STEP + 1>());
    }


    /// Broadcast
    __device__ __forceinline__ T Broadcast(
        T               input,              ///< [in] The value to broadcast
        unsigned int    src_lane)           ///< [in] Which warp lane is to do the broadcasting
    {
        if (lane_id == src_lane)
        {
            ThreadStore<STORE_VOLATILE>(temp_storage[warp_id], input);
        }

        return ThreadLoad<LOAD_VOLATILE>(temp_storage[warp_id]);
    }


    /// Basic inclusive scan
    template <
        bool        HAS_IDENTITY,
        bool        SHARE_FINAL,
        typename    ScanOp>
    __device__ __forceinline__ T BasicScan(
        T               partial,            ///< Calling thread's input partial reduction
        ScanOp          scan_op)            ///< Binary associative scan functor
    {
        // Iterate scan steps
        ScanStep<HAS_IDENTITY>(partial, scan_op, Int2Type<0>());

        if (SHARE_FINAL)
        {
            // Share partial into buffer
            ThreadStore<STORE_VOLATILE>(&temp_storage[warp_id][HALF_WARP_THREADS + lane_id], partial);
        }

        return partial;
    }


    /// Inclusive prefix sum
    __device__ __forceinline__ void InclusiveSum(
        T               input,              ///< [in] Calling thread's input item.
        T               &output)            ///< [out] Calling thread's output item.  May be aliased with \p input.
    {
        const bool HAS_IDENTITY = Traits<T>::PRIMITIVE;

        // Initialize identity region
        InitIdentity(Int2Type<HAS_IDENTITY>());

        // Compute inclusive warp scan (has identity, don't share final)
        output = BasicScan<HAS_IDENTITY, false>(input, Sum());
    }


    /// Inclusive prefix sum with aggregate
    __device__ __forceinline__ void InclusiveSum(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        T               &warp_aggregate)    ///< [out] Warp-wide aggregate reduction of input items.
    {
        const bool HAS_IDENTITY = Traits<T>::PRIMITIVE;

        // Initialize identity region
        InitIdentity(Int2Type<HAS_IDENTITY>());

        // Compute inclusive warp scan (has identity, share final)
        output = BasicScan<HAS_IDENTITY, true>(input, Sum());

        // Retrieve aggregate in <em>warp-lane</em><sub>0</sub>
        warp_aggregate = ThreadLoad<LOAD_VOLATILE>(&temp_storage[warp_id][WARP_SMEM_ELEMENTS - 1]);
    }


    /// Inclusive scan
    template <typename ScanOp>
    __device__ __forceinline__ void InclusiveScan(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        ScanOp          scan_op)            ///< [in] Binary scan operator
    {
        // Compute inclusive warp scan (no identity, don't share final)
        output = BasicScan<false, false>(input, scan_op);
    }


    /// Inclusive scan with aggregate
    template <typename ScanOp>
    __device__ __forceinline__ void InclusiveScan(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        ScanOp          scan_op,            ///< [in] Binary scan operator
        T               &warp_aggregate)    ///< [out] Warp-wide aggregate reduction of input items.
    {
        // Compute inclusive warp scan (no identity, share final)
        output = BasicScan<false, true>(input, scan_op);

        // Retrieve aggregate
        warp_aggregate = ThreadLoad<LOAD_VOLATILE>(&temp_storage[warp_id][WARP_SMEM_ELEMENTS - 1]);
    }

    /// Exclusive scan
    template <typename ScanOp>
    __device__ __forceinline__ void ExclusiveScan(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        T               identity,           ///< [in] Identity value
        ScanOp          scan_op)            ///< [in] Binary scan operator
    {
        // Initialize identity region
        ThreadStore<STORE_VOLATILE>(&temp_storage[warp_id][lane_id], identity);

        // Compute inclusive warp scan (identity, share final)
        T inclusive = BasicScan<true, true>(input, scan_op);

        // Retrieve exclusive scan
        output = ThreadLoad<LOAD_VOLATILE>(&temp_storage[warp_id][HALF_WARP_THREADS + lane_id - 1]);
    }


    /// Exclusive scan with aggregate
    template <typename ScanOp>
    __device__ __forceinline__ void ExclusiveScan(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        T               identity,           ///< [in] Identity value
        ScanOp          scan_op,            ///< [in] Binary scan operator
        T               &warp_aggregate)    ///< [out] Warp-wide aggregate reduction of input items.
    {
        // Exclusive warp scan (which does share final)
        ExclusiveScan(input, output, identity, scan_op);

        // Retrieve aggregate
        warp_aggregate = ThreadLoad<LOAD_VOLATILE>(&temp_storage[warp_id][WARP_SMEM_ELEMENTS - 1]);
    }


    /// Exclusive scan without identity
    template <typename ScanOp>
    __device__ __forceinline__ void ExclusiveScan(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        ScanOp          scan_op)            ///< [in] Binary scan operator
    {
        // Compute inclusive warp scan (no identity, share final)
        T inclusive = BasicScan<false, true>(input, scan_op);

        // Retrieve exclusive scan
        output = ThreadLoad<LOAD_VOLATILE>(&temp_storage[warp_id][HALF_WARP_THREADS + lane_id - 1]);
    }


    /// Exclusive scan with aggregate, without identity
    template <typename ScanOp>
    __device__ __forceinline__ void ExclusiveScan(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        ScanOp          scan_op,            ///< [in] Binary scan operator
        T               &warp_aggregate)    ///< [out] Warp-wide aggregate reduction of input items.
    {
        // Exclusive warp scan (which does share final)
        ExclusiveScan(input, output, scan_op);

        // Retrieve aggregate
        warp_aggregate = ThreadLoad<LOAD_VOLATILE>(&temp_storage[warp_id][WARP_SMEM_ELEMENTS - 1]);
    }

};


}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)
