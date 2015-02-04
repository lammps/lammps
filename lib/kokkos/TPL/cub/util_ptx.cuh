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
 * PTX intrinsics
 */


#pragma once

#include "util_type.cuh"
#include "util_arch.cuh"
#include "util_namespace.cuh"

/// Optional outer namespace(s)
CUB_NS_PREFIX

/// CUB namespace
namespace cub {


/**
 * \addtogroup UtilModule
 * @{
 */


/******************************************************************************
 * PTX helper macros
 ******************************************************************************/

/**
 * Register modifier for pointer-types (for inlining PTX assembly)
 */
#if defined(_WIN64) || defined(__LP64__)
    #define __CUB_LP64__ 1
    // 64-bit register modifier for inlined asm
    #define _CUB_ASM_PTR_ "l"
    #define _CUB_ASM_PTR_SIZE_ "u64"
#else
    #define __CUB_LP64__ 0
    // 32-bit register modifier for inlined asm
    #define _CUB_ASM_PTR_ "r"
    #define _CUB_ASM_PTR_SIZE_ "u32"
#endif


/******************************************************************************
 * Inlined PTX intrinsics
 ******************************************************************************/

/**
 * Shift-right then add.  Returns (x >> shift) + addend.
 */
__device__ __forceinline__ unsigned int SHR_ADD(
    unsigned int x,
    unsigned int shift,
    unsigned int addend)
{
    unsigned int ret;
#if __CUDA_ARCH__ >= 200
    asm("vshr.u32.u32.u32.clamp.add %0, %1, %2, %3;" :
        "=r"(ret) : "r"(x), "r"(shift), "r"(addend));
#else
    ret = (x >> shift) + addend;
#endif
    return ret;
}


/**
 * Shift-left then add.  Returns (x << shift) + addend.
 */
__device__ __forceinline__ unsigned int SHL_ADD(
    unsigned int x,
    unsigned int shift,
    unsigned int addend)
{
    unsigned int ret;
#if __CUDA_ARCH__ >= 200
    asm("vshl.u32.u32.u32.clamp.add %0, %1, %2, %3;" :
        "=r"(ret) : "r"(x), "r"(shift), "r"(addend));
#else
    ret = (x << shift) + addend;
#endif
    return ret;
}


/**
 * Bitfield-extract.
 */
template <typename UnsignedBits>
__device__ __forceinline__ unsigned int BFE(
    UnsignedBits source,
    unsigned int bit_start,
    unsigned int num_bits)
{
    unsigned int bits;
#if __CUDA_ARCH__ >= 200
    asm("bfe.u32 %0, %1, %2, %3;" : "=r"(bits) : "r"((unsigned int) source), "r"(bit_start), "r"(num_bits));
#else
    const unsigned int MASK = (1 << num_bits) - 1;
    bits = (source >> bit_start) & MASK;
#endif
    return bits;
}


/**
 * Bitfield-extract for 64-bit types.
 */
__device__ __forceinline__ unsigned int BFE(
    unsigned long long source,
    unsigned int bit_start,
    unsigned int num_bits)
{
    const unsigned long long MASK = (1ull << num_bits) - 1;
    return (source >> bit_start) & MASK;
}


/**
 * Bitfield insert.  Inserts the first num_bits of y into x starting at bit_start
 */
__device__ __forceinline__ void BFI(
    unsigned int &ret,
    unsigned int x,
    unsigned int y,
    unsigned int bit_start,
    unsigned int num_bits)
{
#if __CUDA_ARCH__ >= 200
    asm("bfi.b32 %0, %1, %2, %3, %4;" :
        "=r"(ret) : "r"(y), "r"(x), "r"(bit_start), "r"(num_bits));
#else
    // TODO
#endif
}


/**
 * Three-operand add
 */
__device__ __forceinline__ unsigned int IADD3(unsigned int x, unsigned int y, unsigned int z)
{
#if __CUDA_ARCH__ >= 200
    asm("vadd.u32.u32.u32.add %0, %1, %2, %3;" : "=r"(x) : "r"(x), "r"(y), "r"(z));
#else
    x = x + y + z;
#endif
    return x;
}


/**
 * Byte-permute. Pick four arbitrary bytes from two 32-bit registers, and
 * reassemble them into a 32-bit destination register
 */
__device__ __forceinline__ int PRMT(unsigned int a, unsigned int b, unsigned int index)
{
    int ret;
    asm("prmt.b32 %0, %1, %2, %3;" : "=r"(ret) : "r"(a), "r"(b), "r"(index));
    return ret;
}


/**
 * Sync-threads barrier.
 */
__device__ __forceinline__ void BAR(int count)
{
    asm volatile("bar.sync 1, %0;" : : "r"(count));
}


/**
 * Floating point multiply. (Mantissa LSB rounds towards zero.)
 */
__device__ __forceinline__ float FMUL_RZ(float a, float b)
{
    float d;
    asm("mul.rz.f32 %0, %1, %2;" : "=f"(d) : "f"(a), "f"(b));
    return d;
}


/**
 * Floating point multiply-add. (Mantissa LSB rounds towards zero.)
 */
__device__ __forceinline__ float FFMA_RZ(float a, float b, float c)
{
    float d;
    asm("fma.rz.f32 %0, %1, %2, %3;" : "=f"(d) : "f"(a), "f"(b), "f"(c));
    return d;
}


/**
 * Terminates the calling thread
 */
__device__ __forceinline__ void ThreadExit() {
    asm("exit;");
}    


/**
 * Returns the warp lane ID of the calling thread
 */
__device__ __forceinline__ unsigned int LaneId()
{
    unsigned int ret;
    asm("mov.u32 %0, %laneid;" : "=r"(ret) );
    return ret;
}


/**
 * Returns the warp ID of the calling thread
 */
__device__ __forceinline__ unsigned int WarpId()
{
    unsigned int ret;
    asm("mov.u32 %0, %warpid;" : "=r"(ret) );
    return ret;
}

/**
 * Returns the warp lane mask of all lanes less than the calling thread
 */
__device__ __forceinline__ unsigned int LaneMaskLt()
{
    unsigned int ret;
    asm("mov.u32 %0, %lanemask_lt;" : "=r"(ret) );
    return ret;
}

/**
 * Returns the warp lane mask of all lanes less than or equal to the calling thread
 */
__device__ __forceinline__ unsigned int LaneMaskLe()
{
    unsigned int ret;
    asm("mov.u32 %0, %lanemask_le;" : "=r"(ret) );
    return ret;
}

/**
 * Returns the warp lane mask of all lanes greater than the calling thread
 */
__device__ __forceinline__ unsigned int LaneMaskGt()
{
    unsigned int ret;
    asm("mov.u32 %0, %lanemask_gt;" : "=r"(ret) );
    return ret;
}

/**
 * Returns the warp lane mask of all lanes greater than or equal to the calling thread
 */
__device__ __forceinline__ unsigned int LaneMaskGe()
{
    unsigned int ret;
    asm("mov.u32 %0, %lanemask_ge;" : "=r"(ret) );
    return ret;
}

/**
 * Portable implementation of __all
 */
__device__ __forceinline__ int WarpAll(int cond)
{
#if CUB_PTX_ARCH < 120

    __shared__ volatile int warp_signals[PtxArchProps::MAX_SM_THREADS / PtxArchProps::WARP_THREADS];

    if (LaneId() == 0)
        warp_signals[WarpId()] = 1;

    if (cond == 0)
        warp_signals[WarpId()] = 0;

    return warp_signals[WarpId()];

#else

    return __all(cond);

#endif
}


/**
 * Portable implementation of __any
 */
__device__ __forceinline__ int WarpAny(int cond)
{
#if CUB_PTX_ARCH < 120

    __shared__ volatile int warp_signals[PtxArchProps::MAX_SM_THREADS / PtxArchProps::WARP_THREADS];

    if (LaneId() == 0)
        warp_signals[WarpId()] = 0;

    if (cond)
        warp_signals[WarpId()] = 1;

    return warp_signals[WarpId()];

#else

    return __any(cond);

#endif
}


/// Generic shuffle-up
template <typename T>
__device__ __forceinline__ T ShuffleUp(
    T               input,              ///< [in] The value to broadcast
    int             src_offset)         ///< [in] The up-offset of the peer to read from
{
    enum
    {
        SHFL_C = 0,
    };

    typedef typename WordAlignment<T>::ShuffleWord ShuffleWord;

    const int       WORDS           = (sizeof(T) + sizeof(ShuffleWord) - 1) / sizeof(ShuffleWord);
    T               output;
    ShuffleWord     *output_alias   = reinterpret_cast<ShuffleWord *>(&output);
    ShuffleWord     *input_alias    = reinterpret_cast<ShuffleWord *>(&input);

    #pragma unroll
    for (int WORD = 0; WORD < WORDS; ++WORD)
    {
        unsigned int shuffle_word = input_alias[WORD];
        asm(
            "  shfl.up.b32 %0, %1, %2, %3;"
            : "=r"(shuffle_word) : "r"(shuffle_word), "r"(src_offset), "r"(SHFL_C));
        output_alias[WORD] = (ShuffleWord) shuffle_word;
    }

    return output;
}



/** @} */       // end group UtilModule

}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)
