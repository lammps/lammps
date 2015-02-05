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
 * cub::WarpScanShfl provides SHFL-based variants of parallel prefix scan across CUDA warps.
 */

#pragma once

#include "../../thread/thread_operators.cuh"
#include "../../util_type.cuh"
#include "../../util_ptx.cuh"
#include "../../util_namespace.cuh"

/// Optional outer namespace(s)
CUB_NS_PREFIX

/// CUB namespace
namespace cub {

/**
 * \brief WarpScanShfl provides SHFL-based variants of parallel prefix scan across CUDA warps.
 */
template <
    typename    T,                      ///< Data type being scanned
    int         LOGICAL_WARPS,          ///< Number of logical warps entrant
    int         LOGICAL_WARP_THREADS>   ///< Number of threads per logical warp
struct WarpScanShfl
{

    /******************************************************************************
     * Constants and typedefs
     ******************************************************************************/

    enum
    {
        /// The number of warp scan steps
        STEPS = Log2<LOGICAL_WARP_THREADS>::VALUE,

        // The 5-bit SHFL mask for logically splitting warps into sub-segments starts 8-bits up
        SHFL_C = ((-1 << STEPS) & 31) << 8,
    };

    /// Shared memory storage layout type
    typedef NullType TempStorage;


    /******************************************************************************
     * Thread fields
     ******************************************************************************/

    int             warp_id;
    int             lane_id;

    /******************************************************************************
     * Construction
     ******************************************************************************/

    /// Constructor
    __device__ __forceinline__ WarpScanShfl(
        TempStorage &temp_storage,
        int warp_id,
        int lane_id)
    :
        warp_id(warp_id),
        lane_id(lane_id)
    {}


    /******************************************************************************
     * Operation
     ******************************************************************************/

    /// Broadcast
    __device__ __forceinline__ T Broadcast(
        T               input,              ///< [in] The value to broadcast
        int             src_lane)           ///< [in] Which warp lane is to do the broadcasting
    {
        typedef typename WordAlignment<T>::ShuffleWord ShuffleWord;

        const int       WORDS           = (sizeof(T) + sizeof(ShuffleWord) - 1) / sizeof(ShuffleWord);
        T               output;
        ShuffleWord     *output_alias   = reinterpret_cast<ShuffleWord *>(&output);
        ShuffleWord     *input_alias    = reinterpret_cast<ShuffleWord *>(&input);

        #pragma unroll
        for (int WORD = 0; WORD < WORDS; ++WORD)
        {
            unsigned int shuffle_word = input_alias[WORD];
            asm("shfl.idx.b32 %0, %1, %2, %3;"
                : "=r"(shuffle_word) : "r"(shuffle_word), "r"(src_lane), "r"(LOGICAL_WARP_THREADS - 1));
            output_alias[WORD] = (ShuffleWord) shuffle_word;
        }

        return output;
    }


    //---------------------------------------------------------------------
    // Inclusive operations
    //---------------------------------------------------------------------

    /// Inclusive prefix sum with aggregate (single-SHFL)
    __device__ __forceinline__ void InclusiveSum(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        T               &warp_aggregate,    ///< [out] Warp-wide aggregate reduction of input items.
        Int2Type<true>  single_shfl)
    {
        unsigned int temp = reinterpret_cast<unsigned int &>(input);

        // Iterate scan steps
        #pragma unroll
        for (int STEP = 0; STEP < STEPS; STEP++)
        {
            // Use predicate set from SHFL to guard against invalid peers
            asm(
                "{"
                "  .reg .u32 r0;"
                "  .reg .pred p;"
                "  shfl.up.b32 r0|p, %1, %2, %3;"
                "  @p add.u32 r0, r0, %4;"
                "  mov.u32 %0, r0;"
                "}"
                : "=r"(temp) : "r"(temp), "r"(1 << STEP), "r"(SHFL_C), "r"(temp));
        }

        output = temp;

        // Grab aggregate from last warp lane
        warp_aggregate = Broadcast(output, LOGICAL_WARP_THREADS - 1);
    }


    /// Inclusive prefix sum with aggregate (multi-SHFL)
    __device__ __forceinline__ void InclusiveSum(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        T               &warp_aggregate,    ///< [out] Warp-wide aggregate reduction of input items.
        Int2Type<false> single_shfl)        ///< [in] Marker type indicating whether only one SHFL instruction is required
    {
        // Delegate to generic scan
        InclusiveScan(input, output, Sum(), warp_aggregate);
    }


    /// Inclusive prefix sum with aggregate (specialized for float)
    __device__ __forceinline__ void InclusiveSum(
        float           input,              ///< [in] Calling thread's input item.
        float           &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        float           &warp_aggregate)    ///< [out] Warp-wide aggregate reduction of input items.
    {
        output = input;

        // Iterate scan steps
        #pragma unroll
        for (int STEP = 0; STEP < STEPS; STEP++)
        {
            // Use predicate set from SHFL to guard against invalid peers
            asm(
                "{"
                "  .reg .f32 r0;"
                "  .reg .pred p;"
                "  shfl.up.b32 r0|p, %1, %2, %3;"
                "  @p add.f32 r0, r0, %4;"
                "  mov.f32 %0, r0;"
                "}"
                : "=f"(output) : "f"(output), "r"(1 << STEP), "r"(SHFL_C), "f"(output));
        }

        // Grab aggregate from last warp lane
        warp_aggregate = Broadcast(output, LOGICAL_WARP_THREADS - 1);
    }


    /// Inclusive prefix sum with aggregate (specialized for unsigned long long)
    __device__ __forceinline__ void InclusiveSum(
        unsigned long long  input,              ///< [in] Calling thread's input item.
        unsigned long long  &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        unsigned long long  &warp_aggregate)    ///< [out] Warp-wide aggregate reduction of input items.
    {
        output = input;

        // Iterate scan steps
        #pragma unroll
        for (int STEP = 0; STEP < STEPS; STEP++)
        {
            // Use predicate set from SHFL to guard against invalid peers
            asm(
                "{"
                "  .reg .u32 r0;"
                "  .reg .u32 r1;"
                "  .reg .u32 lo;"
                "  .reg .u32 hi;"
                "  .reg .pred p;"
                "  mov.b64 {lo, hi}, %1;"
                "  shfl.up.b32 r0|p, lo, %2, %3;"
                "  shfl.up.b32 r1|p, hi, %2, %3;"
                "  @p add.cc.u32 r0, r0, lo;"
                "  @p addc.u32 r1, r1, hi;"
                "  mov.b64 %0, {r0, r1};"
                "}"
                : "=l"(output) : "l"(output), "r"(1 << STEP), "r"(SHFL_C));
        }

        // Grab aggregate from last warp lane
        warp_aggregate = Broadcast(output, LOGICAL_WARP_THREADS - 1);
    }


    /// Inclusive prefix sum with aggregate (generic)
    template <typename _T>
    __device__ __forceinline__ void InclusiveSum(
        _T               input,             ///< [in] Calling thread's input item.
        _T               &output,           ///< [out] Calling thread's output item.  May be aliased with \p input.
        _T               &warp_aggregate)   ///< [out] Warp-wide aggregate reduction of input items.
    {
        // Whether sharing can be done with a single SHFL instruction (vs multiple SFHL instructions)
        Int2Type<(Traits<_T>::PRIMITIVE) && (sizeof(_T) <= sizeof(unsigned int))> single_shfl;

        InclusiveSum(input, output, warp_aggregate, single_shfl);
    }


    /// Inclusive prefix sum
    __device__ __forceinline__ void InclusiveSum(
        T               input,              ///< [in] Calling thread's input item.
        T               &output)            ///< [out] Calling thread's output item.  May be aliased with \p input.
    {
        T warp_aggregate;
        InclusiveSum(input, output, warp_aggregate);
    }


    /// Inclusive scan with aggregate
    template <typename ScanOp>
    __device__ __forceinline__ void InclusiveScan(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        ScanOp          scan_op,            ///< [in] Binary scan operator
        T               &warp_aggregate)    ///< [out] Warp-wide aggregate reduction of input items.
    {
        output = input;

        // Iterate scan steps
        #pragma unroll
        for (int STEP = 0; STEP < STEPS; STEP++)
        {
            // Grab addend from peer
            const int OFFSET = 1 << STEP;
            T temp = ShuffleUp(output, OFFSET);

            // Perform scan op if from a valid peer
            if (lane_id >= OFFSET)
                output = scan_op(temp, output);
        }

        // Grab aggregate from last warp lane
        warp_aggregate = Broadcast(output, LOGICAL_WARP_THREADS - 1);
    }


    /// Inclusive scan
    template <typename ScanOp>
    __device__ __forceinline__ void InclusiveScan(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        ScanOp          scan_op)            ///< [in] Binary scan operator
    {
        T warp_aggregate;
        InclusiveScan(input, output, scan_op, warp_aggregate);
    }


    //---------------------------------------------------------------------
    // Exclusive operations
    //---------------------------------------------------------------------

    /// Exclusive scan with aggregate
    template <typename ScanOp>
    __device__ __forceinline__ void ExclusiveScan(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        T               identity,           ///< [in] Identity value
        ScanOp          scan_op,            ///< [in] Binary scan operator
        T               &warp_aggregate)    ///< [out] Warp-wide aggregate reduction of input items.
    {
        // Compute inclusive scan
        T inclusive;
        InclusiveScan(input, inclusive, scan_op, warp_aggregate);

        // Grab result from predecessor
        T exclusive = ShuffleUp(inclusive, 1);

        output = (lane_id == 0) ?
            identity :
            exclusive;
    }


    /// Exclusive scan
    template <typename ScanOp>
    __device__ __forceinline__ void ExclusiveScan(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        T               identity,           ///< [in] Identity value
        ScanOp          scan_op)            ///< [in] Binary scan operator
    {
        T warp_aggregate;
        ExclusiveScan(input, output, identity, scan_op, warp_aggregate);
    }


    /// Exclusive scan with aggregate, without identity
    template <typename ScanOp>
    __device__ __forceinline__ void ExclusiveScan(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        ScanOp          scan_op,            ///< [in] Binary scan operator
        T               &warp_aggregate)    ///< [out] Warp-wide aggregate reduction of input items.
    {
        // Compute inclusive scan
        T inclusive;
        InclusiveScan(input, inclusive, scan_op, warp_aggregate);

        // Grab result from predecessor
        output = ShuffleUp(inclusive, 1);
    }


    /// Exclusive scan without identity
    template <typename ScanOp>
    __device__ __forceinline__ void ExclusiveScan(
        T               input,              ///< [in] Calling thread's input item.
        T               &output,            ///< [out] Calling thread's output item.  May be aliased with \p input.
        ScanOp          scan_op)            ///< [in] Binary scan operator
    {
        T warp_aggregate;
        ExclusiveScan(input, output, scan_op, warp_aggregate);
    }
};


}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)
