// -------------------------------------------------------------
// cuDPP -- CUDA Data Parallel Primitives library
// -------------------------------------------------------------
// $Revision: 5632 $
// $Date: 2009-07-01 14:36:01 +1000 (Wed, 01 Jul 2009) $
// ------------------------------------------------------------- 
// This source code is distributed under the terms of license.txt in
// the root directory of this source distribution.
// ------------------------------------------------------------- 

/**
 * @file
 * cudpp_util.h
 *
 * @brief C++ utility functions and classes used internally to cuDPP
 */

#ifndef __CUDPP_UTIL_H__
#define __CUDPP_UTIL_H__

#ifdef WIN32
#include <windows.h>
#endif

#include <cuda.h>
#include <cudpp.h>
#include <limits.h>
#include <float.h>

#if (CUDA_VERSION >= 3000)
#define LAUNCH_BOUNDS(x) __launch_bounds__((x))
#define LAUNCH_BOUNDS_MINBLOCKs(x, y) __launch_bounds__((x),(y))
#else
#define LAUNCH_BOUNDS(x)
#define LAUNCH_BOUNDS_MINBLOCKS(x, y)
#endif


/** @brief Determine if \a n is a power of two.
  * @param n Value to be checked to see if it is a power of two
  * @returns True if \a n is a power of two, false otherwise
  */
inline bool 
isPowerOfTwo(int n)
{
    return ((n&(n-1))==0) ;
}

/** @brief Determine if an integer \a n is a multiple of an integer \a f.
  * @param n Multiple
  * @param f Factor
  * @returns True if \a n is a multiple of \a f, false otherwise
  */
inline bool
isMultiple(int n, int f)
{
    if (isPowerOfTwo(f))
        return ((n&(f-1))==0);
    else
        return (n%f==0);
}

/** @brief Compute the smallest power of two larger than \a n.
  * @param n Input value
  * @returns The smallest power f two larger than \a n
  */
inline int 
ceilPow2(int n) 
{
        double log2n = log2((double)n);
        if (isPowerOfTwo(n))
                return n;
        else
                return 1 << (int)ceil(log2n);
}

/** @brief Compute the largest power of two smaller than \a n.
  * @param n Input value
  * @returns The largest power of two smaller than \a n.
  */
inline int 
floorPow2(int n)
{
#ifdef WIN32
    // method 2
    return 1 << (int)_logb((float)n);
#else
    // method 3
    int exp;
    frexp((float)n, &exp);
    return 1 << (exp - 1);
#endif
}

/** @brief Returns the maximum value for type \a T.
  * 
  * Implemented using template specialization on \a T.
  */
template <class T> 
__host__ __device__ inline T getMax() { return 0; }
/** @brief Returns the minimum value for type \a T.
* 
* Implemented using template specialization on \a T.
*/
template <class T> 
__host__ __device__ inline T getMin() { return 0; }
// type specializations for the above
// getMax
template <> __host__ __device__ inline int getMax() { return INT_MAX; }
template <> __host__ __device__ inline unsigned int getMax() { return INT_MAX; }
template <> __host__ __device__ inline float getMax() { return FLT_MAX; }
template <> __host__ __device__ inline char getMax() { return (char)INT_MAX; }
template <> __host__ __device__ inline unsigned char getMax() { return (unsigned char)INT_MAX; }
// getMin
template <> __host__ __device__ inline int getMin() { return INT_MIN; }
template <> __host__ __device__ inline unsigned int getMin() { return 0; }
template <> __host__ __device__ inline float getMin() { return -FLT_MAX; }
template <> __host__ __device__ inline char getMin() { return (char)INT_MIN; }
template <> __host__ __device__ inline unsigned char getMin() { return (unsigned char)0; }

/** @brief Returns the maximum of three values. 
  * @param a First value. 
  * @param b Second value. 
  * @param c Third value. 
  * @returns The maximum of \a a, \a b and \a c.
  */
template<class T>
inline int max3(T a, T b, T c)
{       
    return (a > b) ? ((a > c)? a : c) : ((b > c) ? b : c);
}

/** @brief Utility template struct for generating small vector types from scalar types
  *
  * Given a base scalar type (\c int, \c float, etc.) and a vector length (1 through 4) as 
  * template parameters, this struct defines a vector type (\c float3, \c int4, etc.) of the 
  * specified length and base type.  For example:
  * \code
  * template <class T>
  * __device__ void myKernel(T *data)
  * {
  *     typeToVector<T,4>::Result myVec4;             // create a vec4 of type T
  *     myVec4 = (typeToVector<T,4>::Result*)data[0]; // load first element of data as a vec4
  * }
  * \endcode
  *
  * This functionality is implemented using template specialization.  Currently specializations
  * for int, float, and unsigned int vectors of lengths 2-4 are defined.  Note that this results 
  * in types being generated at compile time -- there is no runtime cost.  typeToVector is used by 
  * the optimized scan \c __device__ functions in scan_cta.cu.
  */
template <typename T, int N>
struct typeToVector
{
    typedef T Result;
};

template<>
struct typeToVector<int, 4>
{
    typedef int4 Result;
};
template<>
struct typeToVector<unsigned int, 4>
{
    typedef uint4 Result;
};
template<>
struct typeToVector<float, 4>
{
    typedef float4 Result;
};
template<>
struct typeToVector<int, 3>
{
    typedef int3 Result;
};
template<>
struct typeToVector<unsigned int, 3>
{
    typedef uint3 Result;
};
template<>
struct typeToVector<float, 3>
{
    typedef float3 Result;
};
template<>
struct typeToVector<int, 2>
{
    typedef int2 Result;
};
template<>
struct typeToVector<unsigned int, 2>
{
    typedef uint2 Result;
};
template<>
struct typeToVector<float, 2>
{
    typedef float2 Result;
};

/** @brief Templatized operator class used by scan and segmented scan
  * 
  * This Operator class is used to allow generic support of binary 
  * associative operators in scan.  It defines two member functions, 
  * op() and identity(), that are used in place of + and 0 (for 
  * example) in the scan and  segmented scan code. Because this is 
  * template code, all decisions in the code are made at compile 
  * time, resulting in optimal operator code. Currently the operators 
  * CUDPP_ADD, CUDPP_MULTIPLY, CUDPP_MIN, and CUDPP_MAX are supported. 
  * Operator is implemented using template specialization for the 
  * types \c int, \c unsigned int, and \c float.
  */
template <typename T, CUDPPOperator oper>
class Operator
{
public:
    /** Applies the operator to operands \a a and \a b.
      * @param a First operand
      * @param b Second operand
      * @returns a OP b, where OP is defined by ::CUDPPOperator \a oper.
      */
    static __device__ T op(const T a, const T b)
    {
        switch (oper)
        {
        case CUDPP_ADD: 
            return a + b;
        case CUDPP_MULTIPLY:
            return a * b;
        case CUDPP_MIN:
            return min(a, b);
        case CUDPP_MAX: 
            return max(a, b);
        }         
    }

    /** Returns the identity element defined for type \a T */
    static __device__ T identity() { return 0; }
};

// specializations for different types
template <CUDPPOperator oper>
class Operator <int, oper>
{
public:
    static __device__ int op(const int a, const int b)
    {
        switch (oper)
        {
        default:
        case CUDPP_ADD: 
            return a + b;
        case CUDPP_MULTIPLY:
            return a * b;
        case CUDPP_MIN:
            return min(a, b);
        case CUDPP_MAX: 
            return max(a, b);
        }         
    }

    static __device__ int identity()
    {
        switch (oper)
        {
        default:
        case CUDPP_ADD:
            return 0;
        case CUDPP_MULTIPLY:
            return 1;
        case CUDPP_MIN:
            return INT_MAX;
        case CUDPP_MAX:
            return INT_MIN;
        }
    }
};

template <CUDPPOperator oper>
class Operator <unsigned int, oper>
{
public:
    static __device__ unsigned int op(const unsigned int a, const unsigned int b)
    {
        switch (oper)
        {
        default:
        case CUDPP_ADD: 
            return a + b;
        case CUDPP_MULTIPLY:
            return a * b;
        case CUDPP_MIN:
            return min(a, b);
        case CUDPP_MAX: 
            return max(a, b);
        }         
    }

    static __device__ unsigned int identity()
    {
        switch (oper)
        {
        default:
        case CUDPP_ADD:
            return 0;
        case CUDPP_MULTIPLY:
            return 1;
        case CUDPP_MIN:
            return UINT_MAX;
        case CUDPP_MAX:
            return 0;
        }
    }
};


template <CUDPPOperator oper>
class Operator <float, oper>
{
public:
    static __device__ float op(const float a, const float b)
    {
        switch (oper)
        {
        default:
        case CUDPP_ADD: 
            return a + b;
        case CUDPP_MULTIPLY:
            return a * b;
        case CUDPP_MIN:
            return min(a, b);
        case CUDPP_MAX: 
            return max(a, b);
        }         
    }

    static __device__ float identity()
    {
        switch (oper)
        {
        default:
        case CUDPP_ADD:
            return 0.0f;
        case CUDPP_MULTIPLY:
            return 1.0f;
        case CUDPP_MIN:
            return FLT_MAX;
        case CUDPP_MAX:
            return -FLT_MAX;
        }
    }
};

#endif // __CUDPP_UTIL_H__

// Leave this at the end of the file
// Local Variables:
// mode:c++
// c-file-style: "NVIDIA"
// End:
