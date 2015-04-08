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
 * Thread utilities for writing memory using PTX cache modifiers.
 */

#pragma once

#include <cuda.h>

#include "../util_ptx.cuh"
#include "../util_type.cuh"
#include "../util_namespace.cuh"

/// Optional outer namespace(s)
CUB_NS_PREFIX

/// CUB namespace
namespace cub {

/**
 * \addtogroup IoModule
 * @{
 */


//-----------------------------------------------------------------------------
// Tags and constants
//-----------------------------------------------------------------------------

/**
 * \brief Enumeration of PTX cache-modifiers for memory store operations.
 */
enum PtxStoreModifier
{
    STORE_DEFAULT,              ///< Default (no modifier)
    STORE_WB,                   ///< Cache write-back all coherent levels
    STORE_CG,                   ///< Cache at global level
    STORE_CS,                   ///< Cache streaming (likely to be accessed once)
    STORE_WT,                   ///< Cache write-through (to system memory)
    STORE_VOLATILE,             ///< Volatile shared (any memory space)
};


/**
 * \name Simple I/O
 * @{
 */

/**
 * \brief Thread utility for writing memory using cub::PtxStoreModifier cache modifiers.
 *
 * Cache modifiers will only be effected for built-in types (i.e., C++
 * primitives and CUDA vector-types).
 *
 * For example:
 * \par
 * \code
 * #include <cub/cub.cuh>
 *
 * // 32-bit store using cache-global modifier:
 * int *d_out;
 * int val;
 * cub::ThreadStore<cub::STORE_CG>(d_out + threadIdx.x, val);
 *
 * // 16-bit store using default modifier
 * short *d_out;
 * short val;
 * cub::ThreadStore<cub::STORE_DEFAULT>(d_out + threadIdx.x, val);
 *
 * // 256-bit store using write-through modifier
 * double4 *d_out;
 * double4 val;
 * cub::ThreadStore<cub::STORE_WT>(d_out + threadIdx.x, val);
 *
 * // 96-bit store using default cache modifier (ignoring STORE_CS)
 * struct TestFoo { bool a; short b; };
 * TestFoo *d_struct;
 * TestFoo val;
 * cub::ThreadStore<cub::STORE_CS>(d_out + threadIdx.x, val);
 * \endcode
 *
 */
template <
    PtxStoreModifier MODIFIER,
    typename OutputIteratorRA,
    typename T>
__device__ __forceinline__ void ThreadStore(OutputIteratorRA itr, T val);


//@}  end member group


#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document


/**
 * Define a int4 (16B) ThreadStore specialization for the given PTX load modifier
 */
#define CUB_STORE_16(cub_modifier, ptx_modifier)                                            \
    template<>                                                                              \
    __device__ __forceinline__ void ThreadStore<cub_modifier, int4*, int4>(int4* ptr, int4 val)              \
    {                                                                                       \
        asm volatile ("st."#ptx_modifier".v4.s32 [%0], {%1, %2, %3, %4};" : :               \
            _CUB_ASM_PTR_(ptr),                                                             \
            "r"(val.x),                                                                     \
            "r"(val.y),                                                                     \
            "r"(val.z),                                                                     \
            "r"(val.w));                                                                    \
    }                                                                                       \
    template<>                                                                              \
    __device__ __forceinline__ void ThreadStore<cub_modifier, longlong2*, longlong2>(longlong2* ptr, longlong2 val)              \
    {                                                                                       \
        asm volatile ("st."#ptx_modifier".v2.s64 [%0], {%1, %2};" : :                       \
            _CUB_ASM_PTR_(ptr),                                                             \
            "l"(val.x),                                                                     \
            "l"(val.y));                                                                    \
    }


/**
 * Define a int2 (8B) ThreadStore specialization for the given PTX load modifier
 */
#define CUB_STORE_8(cub_modifier, ptx_modifier)                                             \
    template<>                                                                              \
    __device__ __forceinline__ void ThreadStore<cub_modifier, short4*, short4>(short4* ptr, short4 val)              \
    {                                                                                       \
        asm volatile ("st."#ptx_modifier".v4.s16 [%0], {%1, %2, %3, %4};" : :               \
            _CUB_ASM_PTR_(ptr),                                                             \
            "h"(val.x),                                                                     \
            "h"(val.y),                                                                     \
            "h"(val.z),                                                                     \
            "h"(val.w));                                                                    \
    }                                                                                       \
    template<>                                                                              \
    __device__ __forceinline__ void ThreadStore<cub_modifier, int2*, int2>(int2* ptr, int2 val)              \
    {                                                                                       \
        asm volatile ("st."#ptx_modifier".v2.s32 [%0], {%1, %2};" : :                       \
            _CUB_ASM_PTR_(ptr),                                                             \
            "r"(val.x),                                                                     \
            "r"(val.y));                                                                    \
    }                                                                                       \
    template<>                                                                              \
    __device__ __forceinline__ void ThreadStore<cub_modifier, long long*, long long>(long long* ptr, long long val)                 \
    {                                                                                       \
        asm volatile ("st."#ptx_modifier".s64 [%0], %1;" : :                                \
            _CUB_ASM_PTR_(ptr),                                                             \
            "l"(val));                                                                      \
    }

/**
 * Define a int (4B) ThreadStore specialization for the given PTX load modifier
 */
#define CUB_STORE_4(cub_modifier, ptx_modifier)                                             \
    template<>                                                                              \
    __device__ __forceinline__ void ThreadStore<cub_modifier, int*, int>(int* ptr, int val)                 \
    {                                                                                       \
        asm volatile ("st."#ptx_modifier".s32 [%0], %1;" : :                                \
            _CUB_ASM_PTR_(ptr),                                                             \
            "r"(val));                                                                      \
    }


/**
 * Define a short (2B) ThreadStore specialization for the given PTX load modifier
 */
#define CUB_STORE_2(cub_modifier, ptx_modifier)                                             \
    template<>                                                                              \
    __device__ __forceinline__ void ThreadStore<cub_modifier, short*, short>(short* ptr, short val)           \
    {                                                                                       \
        asm volatile ("st."#ptx_modifier".s16 [%0], %1;" : :                                \
            _CUB_ASM_PTR_(ptr),                                                             \
            "h"(val));                                                                      \
    }


/**
 * Define a char (1B) ThreadStore specialization for the given PTX load modifier
 */
#define CUB_STORE_1(cub_modifier, ptx_modifier)                                             \
    template<>                                                                              \
    __device__ __forceinline__ void ThreadStore<cub_modifier, char*, char>(char* ptr, char val)              \
    {                                                                                       \
        asm volatile (                                                                      \
        "{"                                                                                 \
        "   .reg .s8 datum;"                                                                \
        "   cvt.s8.s16 datum, %1;"                                                          \
        "   st."#ptx_modifier".s8 [%0], datum;"                                             \
        "}" : :                                                                             \
            _CUB_ASM_PTR_(ptr),                                                             \
            "h"(short(val)));                                                               \
    }

/**
 * Define powers-of-two ThreadStore specializations for the given PTX load modifier
 */
#define CUB_STORE_ALL(cub_modifier, ptx_modifier)                                           \
    CUB_STORE_16(cub_modifier, ptx_modifier)                                                \
    CUB_STORE_8(cub_modifier, ptx_modifier)                                                 \
    CUB_STORE_4(cub_modifier, ptx_modifier)                                                 \
    CUB_STORE_2(cub_modifier, ptx_modifier)                                                 \
    CUB_STORE_1(cub_modifier, ptx_modifier)                                                 \


/**
 * Define ThreadStore specializations for the various PTX load modifiers
 */
#if CUB_PTX_ARCH >= 200
    CUB_STORE_ALL(STORE_WB, ca)
    CUB_STORE_ALL(STORE_CG, cg)
    CUB_STORE_ALL(STORE_CS, cs)
    CUB_STORE_ALL(STORE_WT, cv)
#else
    // STORE_WT on SM10-13 uses "volatile.global" to ensure writes to last level
    CUB_STORE_ALL(STORE_WT, volatile.global)
#endif



/// Helper structure for templated store iteration (inductive case)
template <PtxStoreModifier MODIFIER, int COUNT, int MAX>
struct IterateThreadStore
{
    template <typename T>
    static __device__ __forceinline__ void Store(T *ptr, T *vals)
    {
        ThreadStore<MODIFIER>(ptr + COUNT, vals[COUNT]);
        IterateThreadStore<MODIFIER, COUNT + 1, MAX>::Store(ptr, vals);
    }
};

/// Helper structure for templated store iteration (termination case)
template <PtxStoreModifier MODIFIER, int MAX>
struct IterateThreadStore<MODIFIER, MAX, MAX>
{
    template <typename T>
    static __device__ __forceinline__ void Store(T *ptr, T *vals) {}
};




/**
 * Store with STORE_DEFAULT on iterator types
 */
template <typename OutputIteratorRA, typename T>
__device__ __forceinline__ void ThreadStore(
    OutputIteratorRA            itr,
    T                           val,
    Int2Type<STORE_DEFAULT>     modifier,
    Int2Type<false>             is_pointer)
{
    *itr = val;
}


/**
 * Store with STORE_DEFAULT on pointer types
 */
template <typename T>
__device__ __forceinline__ void ThreadStore(
    T                           *ptr,
    T                           val,
    Int2Type<STORE_DEFAULT>     modifier,
    Int2Type<true>              is_pointer)
{
    *ptr = val;
}


/**
 * Store with STORE_VOLATILE on primitive pointer types
 */
template <typename T>
__device__ __forceinline__ void ThreadStoreVolatile(
    T                           *ptr,
    T                           val,
    Int2Type<true>              is_primitive)
{
    *reinterpret_cast<volatile T*>(ptr) = val;
}


/**
 * Store with STORE_VOLATILE on non-primitive pointer types
 */
template <typename T>
__device__ __forceinline__ void ThreadStoreVolatile(
    T                           *ptr,
    T                           val,
    Int2Type<false>             is_primitive)
{
    typedef typename WordAlignment<T>::VolatileWord VolatileWord;   // Word type for memcopying
    enum { NUM_WORDS = sizeof(T) / sizeof(VolatileWord) };

    // Store into array of uninitialized words
    typename WordAlignment<T>::UninitializedVolatileWords words;
    *reinterpret_cast<T*>(words.buf) = val;

    // Memcopy words to aliased destination
    #pragma unroll
    for (int i = 0; i < NUM_WORDS; ++i)
        reinterpret_cast<volatile VolatileWord*>(ptr)[i] = words.buf[i];
}


/**
 * Store with STORE_VOLATILE on pointer types
 */
template <typename T>
__device__ __forceinline__ void ThreadStore(
    T                           *ptr,
    T                           val,
    Int2Type<STORE_VOLATILE>    modifier,
    Int2Type<true>              is_pointer)
{
    ThreadStoreVolatile(ptr, val, Int2Type<Traits<T>::PRIMITIVE>());
}


#if (CUB_PTX_ARCH <= 350)

/**
 * Store with STORE_CG on pointer types (uses STORE_DEFAULT on current architectures)
 */
template <typename T>
__device__ __forceinline__ void ThreadStore(
    T                           *ptr,
    T                           val,
    Int2Type<STORE_CG>          modifier,
    Int2Type<true>              is_pointer)
{
    ThreadStore<STORE_DEFAULT>(ptr, val);
}

#endif  // (CUB_PTX_ARCH <= 350)


/**
 * Store with arbitrary MODIFIER on pointer types
 */
template <typename T, int MODIFIER>
__device__ __forceinline__ void ThreadStore(
    T                           *ptr,
    T                           val,
    Int2Type<MODIFIER>          modifier,
    Int2Type<true>              is_pointer)
{
    typedef typename WordAlignment<T>::DeviceWord DeviceWord;   // Word type for memcopying
    enum { NUM_WORDS = sizeof(T) / sizeof(DeviceWord) };

    // Store into array of uninitialized words
    typename WordAlignment<T>::UninitializedDeviceWords words;
    *reinterpret_cast<T*>(words.buf) = val;

    // Memcopy words to aliased destination
    IterateThreadStore<PtxStoreModifier(MODIFIER), 0, NUM_WORDS>::Store(
        reinterpret_cast<DeviceWord*>(ptr),
        words.buf);
}


/**
 * Generic ThreadStore definition
 */
template <PtxStoreModifier MODIFIER, typename OutputIteratorRA, typename T>
__device__ __forceinline__ void ThreadStore(OutputIteratorRA itr, T val)
{
    ThreadStore(
        itr,
        val,
        Int2Type<MODIFIER>(),
        Int2Type<IsPointer<OutputIteratorRA>::VALUE>());
}



#endif // DOXYGEN_SHOULD_SKIP_THIS


/** @} */       // end group IoModule


}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)
