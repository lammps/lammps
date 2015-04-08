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
 * Random-access iterator types
 */

#pragma once

#include "thread/thread_load.cuh"
#include "util_device.cuh"
#include "util_debug.cuh"
#include "util_namespace.cuh"

/// Optional outer namespace(s)
CUB_NS_PREFIX

/// CUB namespace
namespace cub {


/******************************************************************************
 * Texture references
 *****************************************************************************/

#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document

// Anonymous namespace
namespace {

/// Templated texture reference type
template <typename T>
struct TexIteratorRef
{
    // Texture reference type
    typedef texture<T, cudaTextureType1D, cudaReadModeElementType> TexRef;

    static TexRef ref;

    /**
     * Bind texture
     */
    static cudaError_t BindTexture(void *d_in)
    {
        cudaChannelFormatDesc tex_desc = cudaCreateChannelDesc<T>();
        if (d_in)
            return (CubDebug(cudaBindTexture(NULL, ref, d_in, tex_desc)));

        return cudaSuccess;
    }

    /**
     * Unbind textures
     */
    static cudaError_t UnbindTexture()
    {
        return CubDebug(cudaUnbindTexture(ref));
    }
};

// Texture reference definitions
template <typename Value>
typename TexIteratorRef<Value>::TexRef TexIteratorRef<Value>::ref = 0;

} // Anonymous namespace


#endif // DOXYGEN_SHOULD_SKIP_THIS







/**
 * \addtogroup UtilModule
 * @{
 */


/******************************************************************************
 * Iterators
 *****************************************************************************/

/**
 * \brief A simple random-access iterator pointing to a range of constant values
 *
 * \par Overview
 * ConstantIteratorRA is a random-access iterator that when dereferenced, always
 * returns the supplied constant of type \p OutputType.
 *
 * \tparam OutputType           The value type of this iterator
 */
template <typename OutputType>
class ConstantIteratorRA
{
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document

    typedef ConstantIteratorRA                  self_type;
    typedef OutputType                          value_type;
    typedef OutputType                          reference;
    typedef OutputType*                         pointer;
    typedef std::random_access_iterator_tag     iterator_category;
    typedef int                                 difference_type;

#endif  // DOXYGEN_SHOULD_SKIP_THIS

private:

    OutputType    val;

public:

    /// Constructor
    __host__ __device__ __forceinline__ ConstantIteratorRA(
        const OutputType &val)          ///< Constant value for the iterator instance to report
    :
        val(val)
    {}

#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document

    __host__ __device__ __forceinline__ self_type operator++()
    {
        self_type i = *this;
        return i;
    }

    __host__ __device__ __forceinline__ self_type operator++(int junk)
    {
        return *this;
    }

    __host__ __device__ __forceinline__ reference operator*()
    {
        return val;
    }

    template <typename SizeT>
    __host__ __device__ __forceinline__ self_type operator+(SizeT n)
    {
        return ConstantIteratorRA(val);
    }

    template <typename SizeT>
    __host__ __device__ __forceinline__ self_type operator-(SizeT n)
    {
        return ConstantIteratorRA(val);
    }

    template <typename SizeT>
    __host__ __device__ __forceinline__ reference operator[](SizeT n)
    {
        return ConstantIteratorRA(val);
    }

    __host__ __device__ __forceinline__ pointer operator->()
    {
        return &val;
    }

    __host__ __device__ __forceinline__ bool operator==(const self_type& rhs)
    {
        return (val == rhs.val);
    }

    __host__ __device__ __forceinline__ bool operator!=(const self_type& rhs)
    {
        return (val != rhs.val);
    }

#endif // DOXYGEN_SHOULD_SKIP_THIS

};



/**
 * \brief A simple random-access transform iterator for applying a transformation operator.
 *
 * \par Overview
 * TransformIteratorRA is a random-access iterator that wraps both a native
 * device pointer of type <tt>InputType*</tt> and a unary conversion functor of
 * type \p ConversionOp. \p OutputType references are made by pulling \p InputType
 * values through the \p ConversionOp instance.
 *
 * \tparam InputType            The value type of the pointer being wrapped
 * \tparam ConversionOp         Unary functor type for mapping objects of type \p InputType to type \p OutputType.  Must have member <tt>OutputType operator()(const InputType &datum)</tt>.
 * \tparam OutputType           The value type of this iterator
 */
template <typename OutputType, typename ConversionOp, typename InputType>
class TransformIteratorRA
{
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document

    typedef TransformIteratorRA                 self_type;
    typedef OutputType                          value_type;
    typedef OutputType                          reference;
    typedef OutputType*                         pointer;
    typedef std::random_access_iterator_tag     iterator_category;
    typedef int                                 difference_type;

#endif  // DOXYGEN_SHOULD_SKIP_THIS

private:

    ConversionOp    conversion_op;
    InputType*      ptr;

public:

    /**
     * \brief Constructor
     * @param ptr Native pointer to wrap
     * @param conversion_op Binary transformation functor
     */
    __host__ __device__ __forceinline__ TransformIteratorRA(InputType* ptr, ConversionOp conversion_op) :
        conversion_op(conversion_op),
        ptr(ptr) {}

#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document

    __host__ __device__ __forceinline__ self_type operator++()
    {
        self_type i = *this;
        ptr++;
        return i;
    }

    __host__ __device__ __forceinline__ self_type operator++(int junk)
    {
        ptr++;
        return *this;
    }

    __host__ __device__ __forceinline__ reference operator*()
    {
        return conversion_op(*ptr);
    }

    template <typename SizeT>
    __host__ __device__ __forceinline__ self_type operator+(SizeT n)
    {
        TransformIteratorRA retval(ptr + n, conversion_op);
        return retval;
    }

    template <typename SizeT>
    __host__ __device__ __forceinline__ self_type operator-(SizeT n)
    {
        TransformIteratorRA retval(ptr - n, conversion_op);
        return retval;
    }

    template <typename SizeT>
    __host__ __device__ __forceinline__ reference operator[](SizeT n)
    {
        return conversion_op(ptr[n]);
    }

    __host__ __device__ __forceinline__ pointer operator->()
    {
        return &conversion_op(*ptr);
    }

    __host__ __device__ __forceinline__ bool operator==(const self_type& rhs)
    {
        return (ptr == rhs.ptr);
    }

    __host__ __device__ __forceinline__ bool operator!=(const self_type& rhs)
    {
        return (ptr != rhs.ptr);
    }

#endif // DOXYGEN_SHOULD_SKIP_THIS

};



/**
 * \brief A simple random-access iterator for loading primitive values through texture cache.
 *
 * \par Overview
 * TexIteratorRA is a random-access iterator that wraps a native
 * device pointer of type <tt>T*</tt>. References made through TexIteratorRA
 * causes values to be pulled through texture cache.
 *
 * \par Usage Considerations
 * - Can only be used with primitive types (e.g., \p char, \p int, \p float), with the exception of \p double
 * - Only one TexIteratorRA or TexIteratorRA of a certain \p InputType can be bound at any given time (per host thread)
 *
 * \tparam InputType            The value type of the pointer being wrapped
 * \tparam ConversionOp         Unary functor type for mapping objects of type \p InputType to type \p OutputType.  Must have member <tt>OutputType operator()(const InputType &datum)</tt>.
 * \tparam OutputType           The value type of this iterator
 */
template <typename T>
class TexIteratorRA
{
public:
#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document

    typedef TexIteratorRA                       self_type;
    typedef T                                   value_type;
    typedef T                                   reference;
    typedef T*                                  pointer;
    typedef std::random_access_iterator_tag     iterator_category;
    typedef int                                 difference_type;

#endif // DOXYGEN_SHOULD_SKIP_THIS

    /// Tag identifying iterator type as being texture-bindable
    typedef void TexBindingTag;

private:

    T*                  ptr;
    size_t              tex_align_offset;
    cudaTextureObject_t tex_obj;

public:

    /**
     * \brief Constructor
     */
    __host__ __device__ __forceinline__ TexIteratorRA()
    :
        ptr(NULL),
        tex_align_offset(0),
        tex_obj(0)
    {}

    /// \brief Bind iterator to texture reference
    cudaError_t BindTexture(
        T               *ptr,                   ///< Native pointer to wrap that is aligned to cudaDeviceProp::textureAlignment
        size_t          bytes,                  ///< Number of items
        size_t          tex_align_offset = 0)   ///< Offset (in items) from ptr denoting the position of the iterator
    {
        this->ptr = ptr;
        this->tex_align_offset = tex_align_offset;

        int ptx_version;
        cudaError_t error = cudaSuccess;
        if (CubDebug(error = PtxVersion(ptx_version))) return error;
        if (ptx_version >= 300)
        {
            // Use texture object
            cudaChannelFormatDesc   channel_desc = cudaCreateChannelDesc<T>();
            cudaResourceDesc        res_desc;
            cudaTextureDesc         tex_desc;
            memset(&res_desc, 0, sizeof(cudaResourceDesc));
            memset(&tex_desc, 0, sizeof(cudaTextureDesc));
            res_desc.resType                = cudaResourceTypeLinear;
            res_desc.res.linear.devPtr      = ptr;
            res_desc.res.linear.desc        = channel_desc;
            res_desc.res.linear.sizeInBytes = bytes;
            tex_desc.readMode               = cudaReadModeElementType;
            return cudaCreateTextureObject(&tex_obj, &res_desc, &tex_desc, NULL);
        }
        else
        {
            // Use texture reference
            return TexIteratorRef<T>::BindTexture(ptr);
        }
    }

    /// \brief Unbind iterator to texture reference
    cudaError_t UnbindTexture()
    {
        int ptx_version;
        cudaError_t error = cudaSuccess;
        if (CubDebug(error = PtxVersion(ptx_version))) return error;
        if (ptx_version < 300)
        {
            // Use texture reference
            return TexIteratorRef<T>::UnbindTexture();
        }
        else
        {
            // Use texture object
            return cudaDestroyTextureObject(tex_obj);
        }
    }

#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document

    __host__ __device__ __forceinline__ self_type operator++()
    {
        self_type i = *this;
        ptr++;
        tex_align_offset++;
        return i;
    }

    __host__ __device__ __forceinline__ self_type operator++(int junk)
    {
        ptr++;
        tex_align_offset++;
        return *this;
    }

    __host__ __device__ __forceinline__ reference operator*()
    {
#if (CUB_PTX_ARCH == 0)
        // Simply dereference the pointer on the host
        return *ptr;
#elif (CUB_PTX_ARCH < 300)
        // Use the texture reference
        return tex1Dfetch(TexIteratorRef<T>::ref, tex_align_offset);
#else
        // Use the texture object
        return conversion_op(tex1Dfetch<InputType>(tex_obj, tex_align_offset));
#endif
    }

    template <typename SizeT>
    __host__ __device__ __forceinline__ self_type operator+(SizeT n)
    {
        TexIteratorRA retval;
        retval.ptr = ptr + n;
        retval.tex_align_offset = tex_align_offset + n;
        return retval;
    }

    template <typename SizeT>
    __host__ __device__ __forceinline__ self_type operator-(SizeT n)
    {
        TexIteratorRA retval;
        retval.ptr = ptr - n;
        retval.tex_align_offset = tex_align_offset - n;
        return retval;
    }

    template <typename SizeT>
    __host__ __device__ __forceinline__ reference operator[](SizeT n)
    {
#if (CUB_PTX_ARCH == 0)
        // Simply dereference the pointer on the host
        return ptr[n];
#elif (CUB_PTX_ARCH < 300)
        // Use the texture reference
        return tex1Dfetch(TexIteratorRef<T>::ref, tex_align_offset + n);
#else
        // Use the texture object
        return conversion_op(tex1Dfetch<InputType>(tex_obj, tex_align_offset + n));
#endif
    }

    __host__ __device__ __forceinline__ pointer operator->()
    {
#if (CUB_PTX_ARCH == 0)
        // Simply dereference the pointer on the host
        return &(*ptr);
#elif (CUB_PTX_ARCH < 300)
        // Use the texture reference
        return &(tex1Dfetch(TexIteratorRef<T>::ref, tex_align_offset));
#else
        // Use the texture object
        return conversion_op(tex1Dfetch<InputType>(tex_obj, tex_align_offset));
#endif
    }

    __host__ __device__ __forceinline__ bool operator==(const self_type& rhs)
    {
        return (ptr == rhs.ptr);
    }

    __host__ __device__ __forceinline__ bool operator!=(const self_type& rhs)
    {
        return (ptr != rhs.ptr);
    }

#endif // DOXYGEN_SHOULD_SKIP_THIS

};


/**
 * \brief A simple random-access transform iterator for loading primitive values through texture cache and and subsequently applying a transformation operator.
 *
 * \par Overview
 * TexTransformIteratorRA is a random-access iterator that wraps both a native
 * device pointer of type <tt>InputType*</tt> and a unary conversion functor of
 * type \p ConversionOp. \p OutputType references are made by pulling \p InputType
 * values through the texture cache and then transformed them using the
 * \p ConversionOp instance.
 *
 * \par Usage Considerations
 * - Can only be used with primitive types (e.g., \p char, \p int, \p float), with the exception of \p double
 * - Only one TexIteratorRA or TexTransformIteratorRA of a certain \p InputType can be bound at any given time (per host thread)
 *
 * \tparam InputType            The value type of the pointer being wrapped
 * \tparam ConversionOp         Unary functor type for mapping objects of type \p InputType to type \p OutputType.  Must have member <tt>OutputType operator()(const InputType &datum)</tt>.
 * \tparam OutputType           The value type of this iterator
 */
template <typename OutputType, typename ConversionOp, typename InputType>
class TexTransformIteratorRA
{
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document

    typedef TexTransformIteratorRA              self_type;
    typedef OutputType                          value_type;
    typedef OutputType                          reference;
    typedef OutputType*                         pointer;
    typedef std::random_access_iterator_tag     iterator_category;
    typedef int                                 difference_type;

#endif  // DOXYGEN_SHOULD_SKIP_THIS

    /// Tag identifying iterator type as being texture-bindable
    typedef void TexBindingTag;

private:

    ConversionOp        conversion_op;
    InputType*          ptr;
    size_t              tex_align_offset;
    cudaTextureObject_t tex_obj;

public:

    /**
     * \brief Constructor
     */
    TexTransformIteratorRA(
        ConversionOp    conversion_op)          ///< Binary transformation functor
    :
        conversion_op(conversion_op),
        ptr(NULL),
        tex_align_offset(0),
        tex_obj(0)
    {}

    /// \brief Bind iterator to texture reference
    cudaError_t BindTexture(
        InputType*      ptr,                    ///< Native pointer to wrap that is aligned to cudaDeviceProp::textureAlignment
        size_t          bytes,                  ///< Number of items
        size_t          tex_align_offset = 0)   ///< Offset (in items) from ptr denoting the position of the iterator
    {
        this->ptr = ptr;
        this->tex_align_offset = tex_align_offset;

        int ptx_version;
        cudaError_t error = cudaSuccess;
        if (CubDebug(error = PtxVersion(ptx_version))) return error;
        if (ptx_version >= 300)
        {
            // Use texture object
            cudaChannelFormatDesc   channel_desc = cudaCreateChannelDesc<InputType>();
            cudaResourceDesc        res_desc;
            cudaTextureDesc         tex_desc;
            memset(&res_desc, 0, sizeof(cudaResourceDesc));
            memset(&tex_desc, 0, sizeof(cudaTextureDesc));
            res_desc.resType                = cudaResourceTypeLinear;
            res_desc.res.linear.devPtr      = ptr;
            res_desc.res.linear.desc        = channel_desc;
            res_desc.res.linear.sizeInBytes = bytes;
            tex_desc.readMode               = cudaReadModeElementType;
            return cudaCreateTextureObject(&tex_obj, &res_desc, &tex_desc, NULL);
        }
        else
        {
            // Use texture reference
            return TexIteratorRef<InputType>::BindTexture(ptr);
        }
    }

    /// \brief Unbind iterator to texture reference
    cudaError_t UnbindTexture()
    {
        int ptx_version;
        cudaError_t error = cudaSuccess;
        if (CubDebug(error = PtxVersion(ptx_version))) return error;
        if (ptx_version >= 300)
        {
            // Use texture object
            return cudaDestroyTextureObject(tex_obj);
        }
        else
        {
            // Use texture reference
            return TexIteratorRef<InputType>::UnbindTexture();
        }
    }

#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document

    __host__ __device__ __forceinline__ self_type operator++()
    {
        self_type i = *this;
        ptr++;
        tex_align_offset++;
        return i;
    }

    __host__ __device__ __forceinline__ self_type operator++(int junk)
    {
        ptr++;
        tex_align_offset++;
        return *this;
    }

    __host__ __device__ __forceinline__ reference operator*()
    {
#if (CUB_PTX_ARCH == 0)
        // Simply dereference the pointer on the host
        return conversion_op(*ptr);
#elif (CUB_PTX_ARCH < 300)
        // Use the texture reference
        return conversion_op(tex1Dfetch(TexIteratorRef<InputType>::ref, tex_align_offset));
#else
        // Use the texture object
        return conversion_op(tex1Dfetch<InputType>(tex_obj, tex_align_offset));
#endif
    }

    template <typename SizeT>
    __host__ __device__ __forceinline__ self_type operator+(SizeT n)
    {
        TexTransformIteratorRA retval(conversion_op);
        retval.ptr = ptr + n;
        retval.tex_align_offset = tex_align_offset + n;
        return retval;
    }

    template <typename SizeT>
    __host__ __device__ __forceinline__ self_type operator-(SizeT n)
    {
        TexTransformIteratorRA retval(conversion_op);
        retval.ptr = ptr - n;
        retval.tex_align_offset = tex_align_offset - n;
        return retval;
    }

    template <typename SizeT>
    __host__ __device__ __forceinline__ reference operator[](SizeT n)
    {
#if (CUB_PTX_ARCH == 0)
        // Simply dereference the pointer on the host
        return conversion_op(ptr[n]);
#elif (CUB_PTX_ARCH < 300)
        // Use the texture reference
        return conversion_op(tex1Dfetch(TexIteratorRef<InputType>::ref, tex_align_offset + n));
#else
        // Use the texture object
        return conversion_op(tex1Dfetch<InputType>(tex_obj, tex_align_offset + n));
#endif
    }

    __host__ __device__ __forceinline__ pointer operator->()
    {
#if (CUB_PTX_ARCH == 0)
        // Simply dereference the pointer on the host
        return &conversion_op(*ptr);
#elif (CUB_PTX_ARCH < 300)
        // Use the texture reference
        return &conversion_op(tex1Dfetch(TexIteratorRef<InputType>::ref, tex_align_offset));
#else
        // Use the texture object
        return &conversion_op(tex1Dfetch<InputType>(tex_obj, tex_align_offset));
#endif
    }

    __host__ __device__ __forceinline__ bool operator==(const self_type& rhs)
    {
        return (ptr == rhs.ptr);
    }

    __host__ __device__ __forceinline__ bool operator!=(const self_type& rhs)
    {
        return (ptr != rhs.ptr);
    }

#endif // DOXYGEN_SHOULD_SKIP_THIS

};




/** @} */       // end group UtilModule

}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)
