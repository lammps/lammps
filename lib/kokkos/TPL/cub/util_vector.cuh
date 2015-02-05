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
 * Vector type inference utilities
 */

#pragma once

#include <iostream>

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
 * Vector type inference utilities.  For example:
 *
 * typename VectorHelper<unsigned int, 2>::Type    // Aliases uint2
 *
 ******************************************************************************/

/**
 * \brief Exposes a member typedef \p Type that names the corresponding CUDA vector type if one exists.  Otherwise \p Type refers to the VectorHelper structure itself, which will wrap the corresponding \p x, \p y, etc. vector fields.
 */
template <typename T, int vec_elements> struct VectorHelper;

#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document

enum
{
    /// The maximum number of elements in CUDA vector types
    MAX_VEC_ELEMENTS = 4,
};


/**
 * Generic vector-1 type
 */
template <typename T>
struct VectorHelper<T, 1>
{
    enum { BUILT_IN = false };

    T x;

    typedef VectorHelper<T, 1> Type;
};

/**
 * Generic vector-2 type
 */
template <typename T>
struct VectorHelper<T, 2>
{
    enum { BUILT_IN = false };

    T x;
    T y;

    typedef VectorHelper<T, 2> Type;
};

/**
 * Generic vector-3 type
 */
template <typename T>
struct VectorHelper<T, 3>
{
    enum { BUILT_IN = false };

    T x;
    T y;
    T z;

    typedef VectorHelper<T, 3> Type;
};

/**
 * Generic vector-4 type
 */
template <typename T>
struct VectorHelper<T, 4>
{
    enum { BUILT_IN = false };

    T x;
    T y;
    T z;
    T w;

    typedef VectorHelper<T, 4> Type;
};

/**
 * Macro for expanding partially-specialized built-in vector types
 */
#define CUB_DEFINE_VECTOR_TYPE(base_type,short_type)                                                            \
  template<> struct VectorHelper<base_type, 1> { typedef short_type##1 Type; enum { BUILT_IN = true }; };         \
  template<> struct VectorHelper<base_type, 2> { typedef short_type##2 Type; enum { BUILT_IN = true }; };         \
  template<> struct VectorHelper<base_type, 3> { typedef short_type##3 Type; enum { BUILT_IN = true }; };         \
  template<> struct VectorHelper<base_type, 4> { typedef short_type##4 Type; enum { BUILT_IN = true }; };

// Expand CUDA vector types for built-in primitives
CUB_DEFINE_VECTOR_TYPE(char,               char)
CUB_DEFINE_VECTOR_TYPE(signed char,        char)
CUB_DEFINE_VECTOR_TYPE(short,              short)
CUB_DEFINE_VECTOR_TYPE(int,                int)
CUB_DEFINE_VECTOR_TYPE(long,               long)
CUB_DEFINE_VECTOR_TYPE(long long,          longlong)
CUB_DEFINE_VECTOR_TYPE(unsigned char,      uchar)
CUB_DEFINE_VECTOR_TYPE(unsigned short,     ushort)
CUB_DEFINE_VECTOR_TYPE(unsigned int,       uint)
CUB_DEFINE_VECTOR_TYPE(unsigned long,      ulong)
CUB_DEFINE_VECTOR_TYPE(unsigned long long, ulonglong)
CUB_DEFINE_VECTOR_TYPE(float,              float)
CUB_DEFINE_VECTOR_TYPE(double,             double)
CUB_DEFINE_VECTOR_TYPE(bool,               uchar)

// Undefine macros
#undef CUB_DEFINE_VECTOR_TYPE

#endif // DOXYGEN_SHOULD_SKIP_THIS


/** @} */       // end group UtilModule

}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)
