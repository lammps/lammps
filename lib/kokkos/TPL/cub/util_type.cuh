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
 * Common type manipulation (metaprogramming) utilities
 */

#pragma once

#include <iostream>
#include <limits>

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
 * Type equality
 ******************************************************************************/

/**
 * \brief Type selection (<tt>IF ? ThenType : ElseType</tt>)
 */
template <bool IF, typename ThenType, typename ElseType>
struct If
{
    /// Conditional type result
    typedef ThenType Type;      // true
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document

template <typename ThenType, typename ElseType>
struct If<false, ThenType, ElseType>
{
    typedef ElseType Type;      // false
};

#endif // DOXYGEN_SHOULD_SKIP_THIS


/******************************************************************************
 * Conditional types
 ******************************************************************************/


/**
 * \brief Type equality test
 */
template <typename A, typename B>
struct Equals
{
    enum {
        VALUE = 0,
        NEGATE = 1
    };
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document

template <typename A>
struct Equals <A, A>
{
    enum {
        VALUE = 1,
        NEGATE = 0
    };
};

#endif // DOXYGEN_SHOULD_SKIP_THIS


/******************************************************************************
 * Marker types
 ******************************************************************************/

/**
 * \brief A simple "NULL" marker type
 */
struct NullType
{
#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document
    template <typename T>
    __host__ __device__ __forceinline__ NullType& operator =(const T& b) { return *this; }
#endif // DOXYGEN_SHOULD_SKIP_THIS
};


/**
 * \brief Allows for the treatment of an integral constant as a type at compile-time (e.g., to achieve static call dispatch based on constant integral values)
 */
template <int A>
struct Int2Type
{
   enum {VALUE = A};
};


/******************************************************************************
 * Size and alignment
 ******************************************************************************/

#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document

template <typename T>
struct WordAlignment
{
    struct Pad
    {
        T       val;
        char    byte;
    };

    enum
    {
        /// The alignment of T in bytes
        ALIGN_BYTES = sizeof(Pad) - sizeof(T)
    };

    /// Biggest shuffle word that T is a whole multiple of and is not larger than the alignment of T
    typedef typename If<(ALIGN_BYTES % 4 == 0),
        int,
        typename If<(ALIGN_BYTES % 2 == 0),
            short,
            char>::Type>::Type                  ShuffleWord;

    /// Biggest volatile word that T is a whole multiple of and is not larger than the alignment of T
    typedef typename If<(ALIGN_BYTES % 8 == 0),
        long long,
        ShuffleWord>::Type                      VolatileWord;

    /// Biggest memory-access word that T is a whole multiple of and is not larger than the alignment of T
    typedef typename If<(ALIGN_BYTES % 16 == 0),
        longlong2,
        typename If<(ALIGN_BYTES % 8 == 0),
            long long,                                 // needed to get heterogenous PODs to work on all platforms
            ShuffleWord>::Type>::Type           DeviceWord;

    enum
    {
        DEVICE_MULTIPLE = sizeof(DeviceWord) / sizeof(T)
    };

    struct UninitializedBytes
    {
        char buf[sizeof(T)];
    };

    struct UninitializedShuffleWords
    {
        ShuffleWord buf[sizeof(T) / sizeof(ShuffleWord)];
    };

    struct UninitializedVolatileWords
    {
        VolatileWord buf[sizeof(T) / sizeof(VolatileWord)];
    };

    struct UninitializedDeviceWords
    {
        DeviceWord buf[sizeof(T) / sizeof(DeviceWord)];
    };


};


#endif // DOXYGEN_SHOULD_SKIP_THIS


/******************************************************************************
 * Wrapper types
 ******************************************************************************/

/**
 * \brief A storage-backing wrapper that allows types with non-trivial constructors to be aliased in unions
 */
template <typename T>
struct Uninitialized
{
    /// Biggest memory-access word that T is a whole multiple of and is not larger than the alignment of T
    typedef typename WordAlignment<T>::DeviceWord DeviceWord;

    enum
    {
        WORDS = sizeof(T) / sizeof(DeviceWord)
    };

    /// Backing storage
    DeviceWord storage[WORDS];

    /// Alias
    __host__ __device__ __forceinline__ T& Alias()
    {
        return reinterpret_cast<T&>(*this);
    }
};


/**
 * \brief A wrapper for passing simple static arrays as kernel parameters
 */
template <typename T, int COUNT>
struct ArrayWrapper
{
    /// Static array of type \p T
    T array[COUNT];
};


/**
 * \brief Double-buffer storage wrapper for multi-pass stream transformations that require more than one storage array for streaming intermediate results back and forth.
 *
 * Many multi-pass computations require a pair of "ping-pong" storage
 * buffers (e.g., one for reading from and the other for writing to, and then
 * vice-versa for the subsequent pass).  This structure wraps a set of device
 * buffers and a "selector" member to track which is "current".
 */
template <typename T>
struct DoubleBuffer
{
    /// Pair of device buffer pointers
    T *d_buffers[2];

    ///  Selector into \p d_buffers (i.e., the active/valid buffer)
    int selector;

    /// \brief Constructor
    __host__ __device__ __forceinline__ DoubleBuffer()
    {
        selector = 0;
        d_buffers[0] = NULL;
        d_buffers[1] = NULL;
    }

    /// \brief Constructor
    __host__ __device__ __forceinline__ DoubleBuffer(
        T *d_current,         ///< The currently valid buffer
        T *d_alternate)       ///< Alternate storage buffer of the same size as \p d_current
    {
        selector = 0;
        d_buffers[0] = d_current;
        d_buffers[1] = d_alternate;
    }

    /// \brief Return pointer to the currently valid buffer
    __host__ __device__ __forceinline__ T* Current() { return d_buffers[selector]; }
};



/******************************************************************************
 * Static math
 ******************************************************************************/

/**
 * \brief Statically determine log2(N), rounded up.
 *
 * For example:
 *     Log2<8>::VALUE   // 3
 *     Log2<3>::VALUE   // 2
 */
template <int N, int CURRENT_VAL = N, int COUNT = 0>
struct Log2
{
    /// Static logarithm value
    enum { VALUE = Log2<N, (CURRENT_VAL >> 1), COUNT + 1>::VALUE };         // Inductive case
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document
template <int N, int COUNT>
struct Log2<N, 0, COUNT>
{
    enum {VALUE = (1 << (COUNT - 1) < N) ?                                  // Base case
        COUNT :
        COUNT - 1 };
};
#endif // DOXYGEN_SHOULD_SKIP_THIS


/**
 * \brief Statically determine if N is a power-of-two
 */
template <int N>
struct PowerOfTwo
{
    enum { VALUE = ((N & (N - 1)) == 0) };
};



/******************************************************************************
 * Pointer vs. iterator detection
 ******************************************************************************/


/**
 * \brief Pointer vs. iterator
 */
template <typename Tp>
struct IsPointer
{
    enum { VALUE = 0 };
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document

template <typename Tp>
struct IsPointer<Tp*>
{
    enum { VALUE = 1 };
};

#endif // DOXYGEN_SHOULD_SKIP_THIS



/******************************************************************************
 * Qualifier detection
 ******************************************************************************/

/**
 * \brief Volatile modifier test
 */
template <typename Tp>
struct IsVolatile
{
    enum { VALUE = 0 };
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document

template <typename Tp>
struct IsVolatile<Tp volatile>
{
    enum { VALUE = 1 };
};

#endif // DOXYGEN_SHOULD_SKIP_THIS


/******************************************************************************
 * Qualifier removal
 ******************************************************************************/

/**
 * \brief Removes \p const and \p volatile qualifiers from type \p Tp.
 *
 * For example:
 *     <tt>typename RemoveQualifiers<volatile int>::Type         // int;</tt>
 */
template <typename Tp, typename Up = Tp>
struct RemoveQualifiers
{
    /// Type without \p const and \p volatile qualifiers
    typedef Up Type;
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document

template <typename Tp, typename Up>
struct RemoveQualifiers<Tp, volatile Up>
{
    typedef Up Type;
};

template <typename Tp, typename Up>
struct RemoveQualifiers<Tp, const Up>
{
    typedef Up Type;
};

template <typename Tp, typename Up>
struct RemoveQualifiers<Tp, const volatile Up>
{
    typedef Up Type;
};

#endif // DOXYGEN_SHOULD_SKIP_THIS



/******************************************************************************
 * Typedef-detection
 ******************************************************************************/


/**
 * \brief Defines a structure \p detector_name that is templated on type \p T.  The \p detector_name struct exposes a constant member \p VALUE indicating whether or not parameter \p T exposes a nested type \p nested_type_name
 */
#define CUB_DEFINE_DETECT_NESTED_TYPE(detector_name, nested_type_name)  \
    template <typename T>                                               \
    struct detector_name                                                \
    {                                                                   \
        template <typename C>                                           \
        static char& test(typename C::nested_type_name*);               \
        template <typename>                                             \
        static int& test(...);                                          \
        enum                                                            \
        {                                                               \
            VALUE = sizeof(test<T>(0)) < sizeof(int)                    \
        };                                                              \
    };



/******************************************************************************
 * Simple enable-if (similar to Boost)
 ******************************************************************************/

/**
 * \brief Simple enable-if (similar to Boost)
 */
template <bool Condition, class T = void>
struct EnableIf
{
    /// Enable-if type for SFINAE dummy variables
    typedef T Type;
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document

template <class T>
struct EnableIf<false, T> {};

#endif // DOXYGEN_SHOULD_SKIP_THIS


/******************************************************************************
 * Typedef-detection
 ******************************************************************************/

/**
 * \brief Determine whether or not BinaryOp's functor is of the form <tt>bool operator()(const T& a, const T&b)</tt> or <tt>bool operator()(const T& a, const T&b, unsigned int idx)</tt>
 */
template <typename T, typename BinaryOp>
struct BinaryOpHasIdxParam
{
private:
    template <typename BinaryOpT, bool (BinaryOpT::*)(const T &a, const T &b, unsigned int idx) const>  struct SFINAE1 {};
    template <typename BinaryOpT, bool (BinaryOpT::*)(const T &a, const T &b, unsigned int idx)>        struct SFINAE2 {};
    template <typename BinaryOpT, bool (BinaryOpT::*)(T a, T b, unsigned int idx) const>                struct SFINAE3 {};
    template <typename BinaryOpT, bool (BinaryOpT::*)(T a, T b, unsigned int idx)>                      struct SFINAE4 {};

    template <typename BinaryOpT, bool (BinaryOpT::*)(const T &a, const T &b, int idx) const>           struct SFINAE5 {};
    template <typename BinaryOpT, bool (BinaryOpT::*)(const T &a, const T &b, int idx)>                 struct SFINAE6 {};
    template <typename BinaryOpT, bool (BinaryOpT::*)(T a, T b, int idx) const>                         struct SFINAE7 {};
    template <typename BinaryOpT, bool (BinaryOpT::*)(T a, T b, int idx)>                               struct SFINAE8 {};

    template <typename BinaryOpT> static char Test(SFINAE1<BinaryOpT, &BinaryOpT::operator()> *);
    template <typename BinaryOpT> static char Test(SFINAE2<BinaryOpT, &BinaryOpT::operator()> *);
    template <typename BinaryOpT> static char Test(SFINAE3<BinaryOpT, &BinaryOpT::operator()> *);
    template <typename BinaryOpT> static char Test(SFINAE4<BinaryOpT, &BinaryOpT::operator()> *);

    template <typename BinaryOpT> static char Test(SFINAE5<BinaryOpT, &BinaryOpT::operator()> *);
    template <typename BinaryOpT> static char Test(SFINAE6<BinaryOpT, &BinaryOpT::operator()> *);
    template <typename BinaryOpT> static char Test(SFINAE7<BinaryOpT, &BinaryOpT::operator()> *);
    template <typename BinaryOpT> static char Test(SFINAE8<BinaryOpT, &BinaryOpT::operator()> *);

    template <typename BinaryOpT> static int Test(...);

public:

    /// Whether the functor BinaryOp has a third <tt>unsigned int</tt> index param
    static const bool HAS_PARAM = sizeof(Test<BinaryOp>(NULL)) == sizeof(char);
};



/******************************************************************************
 * Simple type traits utilities.
 *
 * For example:
 *     Traits<int>::CATEGORY             // SIGNED_INTEGER
 *     Traits<NullType>::NULL_TYPE       // true
 *     Traits<uint4>::CATEGORY           // NOT_A_NUMBER
 *     Traits<uint4>::PRIMITIVE;         // false
 *
 ******************************************************************************/

/**
 * \brief Basic type traits categories
 */
enum Category
{
    NOT_A_NUMBER,
    SIGNED_INTEGER,
    UNSIGNED_INTEGER,
    FLOATING_POINT
};


/**
 * \brief Basic type traits
 */
template <Category _CATEGORY, bool _PRIMITIVE, bool _NULL_TYPE, typename _UnsignedBits>
struct BaseTraits
{
    /// Category
    static const Category CATEGORY      = _CATEGORY;
    enum
    {
        PRIMITIVE       = _PRIMITIVE,
        NULL_TYPE       = _NULL_TYPE,
    };
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document

/**
 * Basic type traits (unsigned primitive specialization)
 */
template <typename _UnsignedBits>
struct BaseTraits<UNSIGNED_INTEGER, true, false, _UnsignedBits>
{
    typedef _UnsignedBits       UnsignedBits;

    static const Category       CATEGORY    = UNSIGNED_INTEGER;
    static const UnsignedBits   MIN_KEY     = UnsignedBits(0);
    static const UnsignedBits   MAX_KEY     = UnsignedBits(-1);

    enum
    {
        PRIMITIVE       = true,
        NULL_TYPE       = false,
    };


    static __device__ __forceinline__ UnsignedBits TwiddleIn(UnsignedBits key)
    {
        return key;
    }

    static __device__ __forceinline__ UnsignedBits TwiddleOut(UnsignedBits key)
    {
        return key;
    }
};


/**
 * Basic type traits (signed primitive specialization)
 */
template <typename _UnsignedBits>
struct BaseTraits<SIGNED_INTEGER, true, false, _UnsignedBits>
{
    typedef _UnsignedBits       UnsignedBits;

    static const Category       CATEGORY    = SIGNED_INTEGER;
    static const UnsignedBits   HIGH_BIT    = UnsignedBits(1) << ((sizeof(UnsignedBits) * 8) - 1);
    static const UnsignedBits   MIN_KEY     = HIGH_BIT;
    static const UnsignedBits   MAX_KEY     = UnsignedBits(-1) ^ HIGH_BIT;

    enum
    {
        PRIMITIVE       = true,
        NULL_TYPE       = false,
    };

    static __device__ __forceinline__ UnsignedBits TwiddleIn(UnsignedBits key)
    {
        return key ^ HIGH_BIT;
    };

    static __device__ __forceinline__ UnsignedBits TwiddleOut(UnsignedBits key)
    {
        return key ^ HIGH_BIT;
    };

};


/**
 * Basic type traits (fp primitive specialization)
 */
template <typename _UnsignedBits>
struct BaseTraits<FLOATING_POINT, true, false, _UnsignedBits>
{
    typedef _UnsignedBits       UnsignedBits;

    static const Category       CATEGORY    = FLOATING_POINT;
    static const UnsignedBits   HIGH_BIT    = UnsignedBits(1) << ((sizeof(UnsignedBits) * 8) - 1);
    static const UnsignedBits   MIN_KEY     = UnsignedBits(-1);
    static const UnsignedBits   MAX_KEY     = UnsignedBits(-1) ^ HIGH_BIT;

    static __device__ __forceinline__ UnsignedBits TwiddleIn(UnsignedBits key)
    {
        UnsignedBits mask = (key & HIGH_BIT) ? UnsignedBits(-1) : HIGH_BIT;
        return key ^ mask;
    };

    static __device__ __forceinline__ UnsignedBits TwiddleOut(UnsignedBits key)
    {
        UnsignedBits mask = (key & HIGH_BIT) ? HIGH_BIT : UnsignedBits(-1);
        return key ^ mask;
    };

    enum
    {
        PRIMITIVE       = true,
        NULL_TYPE       = false,
    };
};

#endif // DOXYGEN_SHOULD_SKIP_THIS


/**
 * \brief Numeric type traits
 */
template <typename T> struct NumericTraits :            BaseTraits<NOT_A_NUMBER, false, false, T> {};

#ifndef DOXYGEN_SHOULD_SKIP_THIS    // Do not document

template <> struct NumericTraits<NullType> :            BaseTraits<NOT_A_NUMBER, false, true, NullType> {};

template <> struct NumericTraits<char> :                BaseTraits<(std::numeric_limits<char>::is_signed) ? SIGNED_INTEGER : UNSIGNED_INTEGER, true, false, unsigned char> {};
template <> struct NumericTraits<signed char> :         BaseTraits<SIGNED_INTEGER, true, false, unsigned char> {};
template <> struct NumericTraits<short> :               BaseTraits<SIGNED_INTEGER, true, false, unsigned short> {};
template <> struct NumericTraits<int> :                 BaseTraits<SIGNED_INTEGER, true, false, unsigned int> {};
template <> struct NumericTraits<long> :                BaseTraits<SIGNED_INTEGER, true, false, unsigned long> {};
template <> struct NumericTraits<long long> :           BaseTraits<SIGNED_INTEGER, true, false, unsigned long long> {};

template <> struct NumericTraits<unsigned char> :       BaseTraits<UNSIGNED_INTEGER, true, false, unsigned char> {};
template <> struct NumericTraits<unsigned short> :      BaseTraits<UNSIGNED_INTEGER, true, false, unsigned short> {};
template <> struct NumericTraits<unsigned int> :        BaseTraits<UNSIGNED_INTEGER, true, false, unsigned int> {};
template <> struct NumericTraits<unsigned long> :       BaseTraits<UNSIGNED_INTEGER, true, false, unsigned long> {};
template <> struct NumericTraits<unsigned long long> :  BaseTraits<UNSIGNED_INTEGER, true, false, unsigned long long> {};

template <> struct NumericTraits<float> :               BaseTraits<FLOATING_POINT, true, false, unsigned int> {};
template <> struct NumericTraits<double> :              BaseTraits<FLOATING_POINT, true, false, unsigned long long> {};

#endif // DOXYGEN_SHOULD_SKIP_THIS


/**
 * \brief Type traits
 */
template <typename T>
struct Traits : NumericTraits<typename RemoveQualifiers<T>::Type> {};



/** @} */       // end group UtilModule

}               // CUB namespace
CUB_NS_POSTFIX  // Optional outer namespace(s)
