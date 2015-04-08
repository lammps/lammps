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
 * cub::WarpReduceShfl provides SHFL-based variants of parallel reduction across CUDA warps.
 */

#pragma once

#include "../../thread/thread_operators.cuh"
#include "../../util_ptx.cuh"
#include "../../util_type.cuh"
#include "../../util_macro.cuh"
#include "../../util_namespace.cuh"

/// Optional outer namespace(s)
CUB_NS_PREFIX

/// CUB namespace
namespace cub {


/**
 * \brief WarpReduceShfl provides SHFL-based variants of parallel reduction across CUDA warps.
 */
template <
    typename    T,                      ///< Data type being reduced
    int         LOGICAL_WARPS,          ///< Number of logical warps entrant
    int         LOGICAL_WARP_THREADS>   ///< Number of threads per logical warp
struct WarpReduceShfl
{
    /******************************************************************************
     * Constants and typedefs
     ******************************************************************************/

    enum
    {
        /// The number of warp reduction steps
        STEPS = Log2<LOGICAL_WARP_THREADS>::VALUE,

        // The 5-bit SHFL mask for logically splitting warps into sub-segments
        SHFL_MASK = (-1 << STEPS) & 31,

        // The 5-bit SFHL clamp
        SHFL_CLAMP = LOGICAL_WARP_THREADS - 1,

        // The packed C argument (mask starts 8 bits up)
        SHFL_C = (SHFL_MASK << 8) | SHFL_CLAMP,
    };


    /// Shared memory storage layout type
    typedef NullType TempStorage;


    /******************************************************************************
     * Thread fields
     ******************************************************************************/

    int     warp_id;
    int     lane_id;


    /******************************************************************************
     * Construction
     ******************************************************************************/

    /// Constructor
    __device__ __forceinline__ WarpReduceShfl(
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

    /// Summation (single-SHFL)
    template <
        bool                FULL_WARPS,             ///< Whether all lanes in each warp are contributing a valid fold of items
        int                 FOLDED_ITEMS_PER_LANE>  ///< Number of items folded into each lane
    __device__ __forceinline__ T Sum(
        T                   input,                  ///< [in] Calling thread's input
        int                 folded_items_per_warp,  ///< [in] Total number of valid items folded into each logical warp
        Int2Type<true>      single_shfl)            ///< [in] Marker type indicating whether only one SHFL instruction is required
    {
        unsigned int output = reinterpret_cast<unsigned int &>(input);

        // Iterate reduction steps
        #pragma unroll
        for (int STEP = 0; STEP < STEPS; STEP++)
        {
            const int OFFSET = 1 << STEP;

            if (FULL_WARPS)
            {
                // Use predicate set from SHFL to guard against invalid peers
                asm(
                    "{"
                    "  .reg .u32 r0;"
                    "  .reg .pred p;"
                    "  shfl.down.b32 r0|p, %1, %2, %3;"
                    "  @p add.u32 r0, r0, %4;"
                    "  mov.u32 %0, r0;"
                    "}"
                    : "=r"(output) : "r"(output), "r"(OFFSET), "r"(SHFL_C), "r"(output));
            }
            else
            {
                // Set range predicate to guard against invalid peers
                asm(
                    "{"
                    "  .reg .u32 r0;"
                    "  .reg .pred p;"
                    "  shfl.down.b32 r0, %1, %2, %3;"
                    "  setp.lt.u32 p, %5, %6;"
                    "  mov.u32 %0, %1;"
                    "  @p add.u32 %0, %1, r0;"
                    "}"
                    : "=r"(output) : "r"(output), "r"(OFFSET), "r"(SHFL_C), "r"(output), "r"((lane_id + OFFSET) * FOLDED_ITEMS_PER_LANE), "r"(folded_items_per_warp));
            }
        }

        return output;
    }


    /// Summation (multi-SHFL)
    template <
        bool                FULL_WARPS,             ///< Whether all lanes in each warp are contributing a valid fold of items
        int                 FOLDED_ITEMS_PER_LANE>  ///< Number of items folded into each lane
    __device__ __forceinline__ T Sum(
        T                   input,              ///< [in] Calling thread's input
        int                 folded_items_per_warp,        ///< [in] Total number of valid items folded into each logical warp
        Int2Type<false>     single_shfl)        ///< [in] Marker type indicating whether only one SHFL instruction is required
    {
        // Delegate to generic reduce
        return Reduce<FULL_WARPS, FOLDED_ITEMS_PER_LANE>(input, folded_items_per_warp, cub::Sum());
    }


    /// Summation (float)
    template <
        bool                FULL_WARPS,             ///< Whether all lanes in each warp are contributing a valid fold of items
        int                 FOLDED_ITEMS_PER_LANE>  ///< Number of items folded into each lane
    __device__ __forceinline__ float Sum(
        float               input,              ///< [in] Calling thread's input
        int                 folded_items_per_warp)        ///< [in] Total number of valid items folded into each logical warp
    {
        T output = input;

        // Iterate reduction steps
        #pragma unroll
        for (int STEP = 0; STEP < STEPS; STEP++)
        {
            const int OFFSET = 1 << STEP;

            if (FULL_WARPS)
            {
                // Use predicate set from SHFL to guard against invalid peers
                asm(
                    "{"
                    "  .reg .f32 r0;"
                    "  .reg .pred p;"
                    "  shfl.down.b32 r0|p, %1, %2, %3;"
                    "  @p add.f32 r0, r0, %4;"
                    "  mov.f32 %0, r0;"
                    "}"
                    : "=f"(output) : "f"(output), "r"(OFFSET), "r"(SHFL_C), "f"(output));
            }
            else
            {
                // Set range predicate to guard against invalid peers
                asm(
                    "{"
                    "  .reg .f32 r0;"
                    "  .reg .pred p;"
                    "  shfl.down.b32 r0, %1, %2, %3;"
                    "  setp.lt.u32 p, %5, %6;"
                    "  mov.f32 %0, %1;"
                    "  @p add.f32 %0, %0, r0;"
                    "}"
                    : "=f"(output) : "f"(output), "r"(OFFSET), "r"(SHFL_C), "f"(output), "r"((lane_id + OFFSET) * FOLDED_ITEMS_PER_LANE), "r"(folded_items_per_warp));
            }
        }

        return output;
    }

    /// Summation (generic)
    template <
        bool                FULL_WARPS,             ///< Whether all lanes in each warp are contributing a valid fold of items
        int                 FOLDED_ITEMS_PER_LANE,  ///< Number of items folded into each lane
        typename            _T>
    __device__ __forceinline__ _T Sum(
        _T                  input,                  ///< [in] Calling thread's input
        int                 folded_items_per_warp)  ///< [in] Total number of valid items folded into each logical warp
    {
        // Whether sharing can be done with a single SHFL instruction (vs multiple SFHL instructions)
        Int2Type<(Traits<_T>::PRIMITIVE) && (sizeof(_T) <= sizeof(unsigned int))> single_shfl;

        return Sum<FULL_WARPS, FOLDED_ITEMS_PER_LANE>(input, folded_items_per_warp, single_shfl);
    }


    /// Reduction
    template <
        bool            FULL_WARPS,             ///< Whether all lanes in each warp are contributing a valid fold of items
        int             FOLDED_ITEMS_PER_LANE,  ///< Number of items folded into each lane
        typename        ReductionOp>
    __device__ __forceinline__ T Reduce(
        T               input,                  ///< [in] Calling thread's input
        int             folded_items_per_warp,  ///< [in] Total number of valid items folded into each logical warp
        ReductionOp     reduction_op)           ///< [in] Binary reduction operator
    {
        typedef typename WordAlignment<T>::ShuffleWord ShuffleWord;

        const int       WORDS           = (sizeof(T) + sizeof(ShuffleWord) - 1) / sizeof(ShuffleWord);
        T               output          = input;
        T               temp;
        ShuffleWord     *temp_alias     = reinterpret_cast<ShuffleWord *>(&temp);
        ShuffleWord     *output_alias   = reinterpret_cast<ShuffleWord *>(&output);

        // Iterate scan steps
        #pragma unroll
        for (int STEP = 0; STEP < STEPS; STEP++)
        {
            // Grab addend from peer
            const int OFFSET = 1 << STEP;

            #pragma unroll
            for (int WORD = 0; WORD < WORDS; ++WORD)
            {
                unsigned int shuffle_word = output_alias[WORD];
                asm(
                    "  shfl.down.b32 %0, %1, %2, %3;"
                    : "=r"(shuffle_word) : "r"(shuffle_word), "r"(OFFSET), "r"(SHFL_C));
                temp_alias[WORD] = (ShuffleWord) shuffle_word;
            }

            // Perform reduction op if from a valid peer
            if (FULL_WARPS)
            {
                if (lane_id < LOGICAL_WARP_THREADS - OFFSET)
                    output = reduction_op(output, temp);
            }
            else
            {
                if (((lane_id + OFFSET) * FOLDED_ITEMS_PER_LANE) < folded_items_per_warp)
                    output = reduction_op(output, temp);
            }
        }

        return output;
    }


    /// Segmented reduction
    template <
        bool            HEAD_SEGMENTED,     ///< Whether flags indicate a segment-head or a segment-tail
        typename        Flag,
        typename        ReductionOp>
    __device__ __forceinline__ T SegmentedReduce(
        T               input,              ///< [in] Calling thread's input
        Flag            flag,               ///< [in] Whether or not the current lane is a segment head/tail
        ReductionOp     reduction_op)       ///< [in] Binary reduction operator
    {
        typedef typename WordAlignment<T>::ShuffleWord ShuffleWord;

        T output = input;

        const int       WORDS           = (sizeof(T) + sizeof(ShuffleWord) - 1) / sizeof(ShuffleWord);
        T               temp;
        ShuffleWord     *temp_alias     = reinterpret_cast<ShuffleWord *>(&temp);
        ShuffleWord     *output_alias   = reinterpret_cast<ShuffleWord *>(&output);

        // Get the start flags for each thread in the warp.
        int warp_flags = __ballot(flag);

        if (!HEAD_SEGMENTED)
            warp_flags <<= 1;

        // Keep bits above the current thread.
        warp_flags &= LaneMaskGt();

        // Accommodate packing of multiple logical warps in a single physical warp
        if ((LOGICAL_WARPS > 1) && (LOGICAL_WARP_THREADS < 32))
            warp_flags >>= (warp_id * LOGICAL_WARP_THREADS);

        // Find next flag
        int next_flag = __clz(__brev(warp_flags));

        // Clip the next segment at the warp boundary if necessary
        if (LOGICAL_WARP_THREADS != 32)
            next_flag = CUB_MIN(next_flag, LOGICAL_WARP_THREADS);

        // Iterate scan steps
        #pragma unroll
        for (int STEP = 0; STEP < STEPS; STEP++)
        {
            const int OFFSET = 1 << STEP;

            // Grab addend from peer
            #pragma unroll
            for (int WORD = 0; WORD < WORDS; ++WORD)
            {
                unsigned int shuffle_word = output_alias[WORD];

                asm(
                    "  shfl.down.b32 %0, %1, %2, %3;"
                    : "=r"(shuffle_word) : "r"(shuffle_word), "r"(OFFSET), "r"(SHFL_C));
                temp_alias[WORD] = (ShuffleWord) shuffle_word;

            }

            // Perform reduction op if valid
            if (OFFSET < next_flag - lane_id)
                output = reduction_op(output, temp);
        }

        return output;
    }
};


}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)
