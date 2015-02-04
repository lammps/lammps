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
 * Thread utilities for reading memory using PTX cache modifiers.
 */

#pragma once

#include <cuda.h>

#include <iterator>

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
 * \brief Enumeration of PTX cache-modifiers for memory load operations.
 */
enum PtxLoadModifier
{
    LOAD_DEFAULT,       ///< Default (no modifier)
    LOAD_CA,            ///< Cache at all levels
    LOAD_CG,            ///< Cache at global level
    LOAD_CS,            ///< Cache streaming (likely to be accessed once)
    LOAD_CV,            ///< Cache as volatile (including cached system lines)
    LOAD_LDG,           ///< Cache as texture
    LOAD_VOLATILE,      ///< Volatile (any memory space)
};


/**
 * \name Simple I/O
 * @{
 */

/**
 * \brief Thread utility for reading memory using cub::PtxLoadModifier cache modifiers.
 *
 * Cache modifiers will only be effected for built-in types (i.e., C++
 * primitives and CUDA vector-types).
 *
 * For example:
 * \par
 * \code
 * #include <cub/cub.cuh>
 *
 * // 32-bit load using cache-global modifier:
 * int *d_in;
 * int val = cub::ThreadLoad<cub::LOAD_CA>(d_in + threadIdx.x);
 *
 * // 16-bit load using default modifier
 * short *d_in;
 * short val = cub::ThreadLoad<cub::LOAD_DEFAULT>(d_in + threadIdx.x);
 *
 * // 256-bit load using cache-volatile modifier
 * double4 *d_in;
 * double4 val = cub::ThreadLoad<cub::LOAD_CV>(d_in + threadIdx.x);
 *
 * // 96-bit load using default cache modifier (ignoring LOAD_CS)
 * struct TestFoo { bool a; short b; };
 * TestFoo *d_struct;
 * TestFoo val = cub::ThreadLoad<cub::LOAD_CS>(d_in + threadIdx.x);
 * \endcode
 *
 */
template <
    PtxLoadModifier MODIFIER,
    typename InputIteratorRA>
__device__ __forceinline__ typename std::iterator_traits<InputIteratorRA>::value_type ThreadLoad(InputIteratorRA itr);


//@}  end member group



#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document


/**
 * Define a int4 (16B) ThreadLoad specialization for the given PTX load modifier
 */
#define CUB_LOAD_16(cub_modifier, ptx_modifier)                                             \
    template<>                                                                              \
    __device__ __forceinline__ int4 ThreadLoad<cub_modifier, int4*>(int4* ptr)              \
    {                                                                                       \
        int4 retval;                                                                        \
        asm volatile ("ld."#ptx_modifier".v4.s32 {%0, %1, %2, %3}, [%4];" :                 \
            "=r"(retval.x),                                                                 \
            "=r"(retval.y),                                                                 \
            "=r"(retval.z),                                                                 \
            "=r"(retval.w) :                                                                \
            _CUB_ASM_PTR_(ptr));                                                            \
        return retval;                                                                      \
    }                                                                                       \
    template<>                                                                              \
    __device__ __forceinline__ longlong2 ThreadLoad<cub_modifier, longlong2*>(longlong2* ptr)              \
    {                                                                                       \
        longlong2 retval;                                                                   \
        asm volatile ("ld."#ptx_modifier".v2.s64 {%0, %1}, [%2];" :                         \
            "=l"(retval.x),                                                                 \
            "=l"(retval.y) :                                                                \
            _CUB_ASM_PTR_(ptr));                                                            \
        return retval;                                                                      \
    }

/**
 * Define a int2 (8B) ThreadLoad specialization for the given PTX load modifier
 */
#define CUB_LOAD_8(cub_modifier, ptx_modifier)                                              \
    template<>                                                                              \
    __device__ __forceinline__ short4 ThreadLoad<cub_modifier, short4*>(short4* ptr)        \
    {                                                                                       \
        short4 retval;                                                                      \
        asm volatile ("ld."#ptx_modifier".v4.s16 {%0, %1, %2, %3}, [%4];" :                 \
            "=h"(retval.x),                                                                 \
            "=h"(retval.y),                                                                 \
            "=h"(retval.z),                                                                 \
            "=h"(retval.w) :                                                                \
            _CUB_ASM_PTR_(ptr));                                                            \
        return retval;                                                                      \
    }                                                                                       \
    template<>                                                                              \
    __device__ __forceinline__ int2 ThreadLoad<cub_modifier, int2*>(int2* ptr)              \
    {                                                                                       \
        int2 retval;                                                                        \
        asm volatile ("ld."#ptx_modifier".v2.s32 {%0, %1}, [%2];" :                         \
            "=r"(retval.x),                                                                 \
            "=r"(retval.y) :                                                                \
            _CUB_ASM_PTR_(ptr));                                                            \
        return retval;                                                                      \
    }                                                                                       \
    template<>                                                                              \
    __device__ __forceinline__ long long ThreadLoad<cub_modifier, long long*>(long long* ptr)                 \
    {                                                                                       \
        long long retval;                                                                   \
        asm volatile ("ld."#ptx_modifier".s64 %0, [%1];" :                                  \
            "=l"(retval) :                                                                  \
            _CUB_ASM_PTR_(ptr));                                                            \
        return retval;                                                                      \
    }

/**
 * Define a int (4B) ThreadLoad specialization for the given PTX load modifier
 */
#define CUB_LOAD_4(cub_modifier, ptx_modifier)                                              \
    template<>                                                                              \
    __device__ __forceinline__ int ThreadLoad<cub_modifier, int*>(int* ptr)                 \
    {                                                                                       \
        int retval;                                                                         \
        asm volatile ("ld."#ptx_modifier".s32 %0, [%1];" :                                  \
            "=r"(retval) :                                                                  \
            _CUB_ASM_PTR_(ptr));                                                            \
        return retval;                                                                      \
    }


/**
 * Define a short (2B) ThreadLoad specialization for the given PTX load modifier
 */
#define CUB_LOAD_2(cub_modifier, ptx_modifier)                                              \
    template<>                                                                              \
    __device__ __forceinline__ short ThreadLoad<cub_modifier, short*>(short* ptr)           \
    {                                                                                       \
        short retval;                                                                       \
        asm volatile ("ld."#ptx_modifier".s16 %0, [%1];" :                                  \
            "=h"(retval) :                                                                  \
            _CUB_ASM_PTR_(ptr));                                                            \
        return retval;                                                                      \
    }


/**
 * Define a char (1B) ThreadLoad specialization for the given PTX load modifier
 */
#define CUB_LOAD_1(cub_modifier, ptx_modifier)                                              \
    template<>                                                                              \
    __device__ __forceinline__ char ThreadLoad<cub_modifier, char*>(char* ptr)              \
    {                                                                                       \
        short retval;                                                                       \
        asm volatile (                                                                      \
        "{"                                                                                 \
        "   .reg .s8 datum;"                                                                \
        "    ld."#ptx_modifier".s8 datum, [%1];"                                            \
        "    cvt.s16.s8 %0, datum;"                                                         \
        "}" :                                                                               \
            "=h"(retval) :                                                                  \
            _CUB_ASM_PTR_(ptr));                                                            \
        return (char) retval;                                                               \
    }


/**
 * Define powers-of-two ThreadLoad specializations for the given PTX load modifier
 */
#define CUB_LOAD_ALL(cub_modifier, ptx_modifier)                                            \
    CUB_LOAD_16(cub_modifier, ptx_modifier)                                                 \
    CUB_LOAD_8(cub_modifier, ptx_modifier)                                                  \
    CUB_LOAD_4(cub_modifier, ptx_modifier)                                                  \
    CUB_LOAD_2(cub_modifier, ptx_modifier)                                                  \
    CUB_LOAD_1(cub_modifier, ptx_modifier)                                                  \


/**
 * Define ThreadLoad specializations for the various PTX load modifiers
 */
#if CUB_PTX_ARCH >= 200
    CUB_LOAD_ALL(LOAD_CA, ca)
    CUB_LOAD_ALL(LOAD_CG, cg)
    CUB_LOAD_ALL(LOAD_CS, cs)
    CUB_LOAD_ALL(LOAD_CV, cv)
#else
    // LOAD_CV on SM10-13 uses "volatile.global" to ensure reads from last level
    CUB_LOAD_ALL(LOAD_CV, volatile.global)
#endif
#if CUB_PTX_ARCH >= 350
    CUB_LOAD_ALL(LOAD_LDG, global.nc)
#endif


/// Helper structure for templated load iteration (inductive case)
template <PtxLoadModifier MODIFIER, int COUNT, int MAX>
struct IterateThreadLoad
{
    template <typename T>
    static __device__ __forceinline__ void Load(T *ptr, T *vals)
    {
        vals[COUNT] = ThreadLoad<MODIFIER>(ptr + COUNT);
        IterateThreadLoad<MODIFIER, COUNT + 1, MAX>::Load(ptr, vals);
    }
};

/// Helper structure for templated load iteration (termination case)
template <PtxLoadModifier MODIFIER, int MAX>
struct IterateThreadLoad<MODIFIER, MAX, MAX>
{
    template <typename T>
    static __device__ __forceinline__ void Load(T *ptr, T *vals) {}
};



/**
 * Load with LOAD_DEFAULT on iterator types
 */
template <typename InputIteratorRA>
__device__ __forceinline__ typename std::iterator_traits<InputIteratorRA>::value_type ThreadLoad(
    InputIteratorRA         itr,
    Int2Type<LOAD_DEFAULT>  modifier,
    Int2Type<false>         is_pointer)
{
    return *itr;
}


/**
 * Load with LOAD_DEFAULT on pointer types
 */
template <typename T>
__device__ __forceinline__ T ThreadLoad(
    T                       *ptr,
    Int2Type<LOAD_DEFAULT>  modifier,
    Int2Type<true>          is_pointer)
{
    return *ptr;
}


/**
 * Load with LOAD_VOLATILE on primitive pointer types
 */
template <typename T>
__device__ __forceinline__ T ThreadLoadVolatile(
    T                       *ptr,
    Int2Type<true>          is_primitive)
{
    T retval = *reinterpret_cast<volatile T*>(ptr);

#if (CUB_PTX_ARCH <= 130)
    if (sizeof(T) == 1) __threadfence_block();
#endif

    return retval;
}


/**
 * Load with LOAD_VOLATILE on non-primitive pointer types
 */
template <typename T>
__device__ __forceinline__ T ThreadLoadVolatile(
    T                       *ptr,
    Int2Type<false>          is_primitive)
{
    typedef typename WordAlignment<T>::VolatileWord VolatileWord;   // Word type for memcopying
    enum { NUM_WORDS = sizeof(T) / sizeof(VolatileWord) };

    // Memcopy from aliased source into array of uninitialized words
    typename WordAlignment<T>::UninitializedVolatileWords words;

    #pragma unroll
    for (int i = 0; i < NUM_WORDS; ++i)
        words.buf[i] = reinterpret_cast<volatile VolatileWord*>(ptr)[i];

    // Load from words
    return *reinterpret_cast<T*>(words.buf);
}


/**
 * Load with LOAD_VOLATILE on pointer types
 */
template <typename T>
__device__ __forceinline__ T ThreadLoad(
    T                       *ptr,
    Int2Type<LOAD_VOLATILE> modifier,
    Int2Type<true>          is_pointer)
{
    return ThreadLoadVolatile(ptr, Int2Type<Traits<T>::PRIMITIVE>());
}


#if (CUB_PTX_ARCH <= 130)

/**
 * Load with LOAD_CG uses LOAD_CV in pre-SM20 PTX to ensure coherent reads when run on newer architectures with L1
 */
template <typename T>
__device__ __forceinline__ T ThreadLoad(
    T                       *ptr,
    Int2Type<LOAD_CG>       modifier,
    Int2Type<true>          is_pointer)
{
    return ThreadLoad<LOAD_CV>(ptr);
}

#endif  // (CUB_PTX_ARCH <= 130)


/**
 * Load with arbitrary MODIFIER on pointer types
 */
template <typename T, int MODIFIER>
__device__ __forceinline__ T ThreadLoad(
    T                       *ptr,
    Int2Type<MODIFIER>      modifier,
    Int2Type<true>          is_pointer)
{
    typedef typename WordAlignment<T>::DeviceWord DeviceWord;
    enum { NUM_WORDS = sizeof(T) / sizeof(DeviceWord) };

    // Memcopy from aliased source into array of uninitialized words
    typename WordAlignment<T>::UninitializedDeviceWords words;

    IterateThreadLoad<PtxLoadModifier(MODIFIER), 0, NUM_WORDS>::Load(
        reinterpret_cast<DeviceWord*>(ptr),
        words.buf);

    // Load from words
    return *reinterpret_cast<T*>(words.buf);
}


/**
 * Generic ThreadLoad definition
 */
template <
    PtxLoadModifier MODIFIER,
    typename InputIteratorRA>
__device__ __forceinline__ typename std::iterator_traits<InputIteratorRA>::value_type ThreadLoad(InputIteratorRA itr)
{
    return ThreadLoad(
        itr,
        Int2Type<MODIFIER>(),
        Int2Type<IsPointer<InputIteratorRA>::VALUE>());
}



#endif // DOXYGEN_SHOULD_SKIP_THIS


/** @} */       // end group IoModule


}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)
