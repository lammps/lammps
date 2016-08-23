// -------------------------------------------------------------
// cuDPP -- CUDA Data Parallel Primitives library
// -------------------------------------------------------------
// $Revision: 5289 $
// $Date: 2010-11-23 13:04:43 -0700 (Tue, 23 Nov 2010) $
// ------------------------------------------------------------- 
// This source code is distributed under the terms of license.txt 
// in the root directory of this source distribution.
// ------------------------------------------------------------- 
 
/**
 * @file
 * sharedmem.h
 *
 * @brief Shared memory declaration struct for templatized types.
 *
 * Because dynamically sized shared memory arrays are declared "extern" in CUDA,
 * we can't templatize their types directly.  To get around this, we declare a 
 * simple wrapper struct that will declare the extern array with a different 
 * name depending on the type.  This avoids linker errors about multiple
 * definitions.
 * 
 * To use dynamically allocated shared memory in a templatized __global__ or 
 * __device__ function, just replace code like this:
 *
 * <pre>
 *  template<class T>
 *  __global__ void
 *  foo( T* d_out, T* d_in) 
 *  {
 *      // Shared mem size is determined by the host app at run time
 *      extern __shared__  T sdata[];
 *      ...
 *      doStuff(sdata);
 *      ...
 *  }
 * </pre>
 *  
 *  With this
 * <pre>
 *  template<class T>
 *  __global__ void
 *  foo( T* d_out, T* d_in) 
 *  {
 *      // Shared mem size is determined by the host app at run time
 *      SharedMemory<T> smem;
 *      T* sdata = smem.getPointer();
 *      ...
 *      doStuff(sdata);
 *      ...
 *  }
 * </pre>
 */

#ifndef _SHAREDMEM_H_
#define _SHAREDMEM_H_


/** @brief Wrapper class for templatized dynamic shared memory arrays.
  * 
  * This struct uses template specialization on the type \a T to declare
  * a differently named dynamic shared memory array for each type
  * (\code extern __shared__ T s_type[] \endcode).
  * 
  * Currently there are specializations for the following types:
  * \c int, \c uint, \c char, \c uchar, \c short, \c ushort, \c long, 
  * \c unsigned long, \c bool, \c float, and \c double. One can also specialize it
  * for user defined types.
  */
template <typename T>
struct SharedMemory
{
    /** Return a pointer to the runtime-sized shared memory array. **/
    __device__ T* getPointer() 
    { 
        extern __device__ void Error_UnsupportedType(); // Ensure that we won't compile any un-specialized types
        Error_UnsupportedType();
        return (T*)0;
    }
    // TODO: Use operator overloading to make this class look like a regular array
};

// Following are the specializations for the following types.
// int, uint, char, uchar, short, ushort, long, ulong, bool, float, and double
// One could also specialize it for user-defined types.

template <>
struct SharedMemory <int>
{
    __device__ int* getPointer() { extern __shared__ int s_int[]; return s_int; }      
};

template <>
struct SharedMemory <unsigned int>
{
    __device__ unsigned int* getPointer() { extern __shared__ unsigned int s_uint[]; return s_uint; }    
};

template <>
struct SharedMemory <char>
{
    __device__ char* getPointer() { extern __shared__ char s_char[]; return s_char; }    
};

template <>
struct SharedMemory <unsigned char>
{
    __device__ unsigned char* getPointer() { extern __shared__ unsigned char s_uchar[]; return s_uchar; }    
};

template <>
struct SharedMemory <short>
{
    __device__ short* getPointer() { extern __shared__ short s_short[]; return s_short; }    
};

template <>
struct SharedMemory <unsigned short>
{
    __device__ unsigned short* getPointer() { extern __shared__ unsigned short s_ushort[]; return s_ushort; }    
};

template <>
struct SharedMemory <long>
{
    __device__ long* getPointer() { extern __shared__ long s_long[]; return s_long; }    
};

template <>
struct SharedMemory <unsigned long>
{
    __device__ unsigned long* getPointer() { extern __shared__ unsigned long s_ulong[]; return s_ulong; }    
};

template <>
struct SharedMemory <bool>
{
    __device__ bool* getPointer() { extern __shared__ bool s_bool[]; return s_bool; }    
};

template <>
struct SharedMemory <float>
{
    __device__ float* getPointer() { extern __shared__ float s_float[]; return s_float; }    
};

template <>
struct SharedMemory <double>
{
    __device__ double* getPointer() { extern __shared__ double s_double[]; return s_double; }    
};

template <>
struct SharedMemory <uchar4>
{
    __device__ uchar4* getPointer() { extern __shared__ uchar4 s_uchar4[]; return s_uchar4; }    
};


#endif //_SHAREDMEM_H_

// Leave this at the end of the file
// Local Variables:
// mode:c++
// c-file-style: "NVIDIA"
// End:
