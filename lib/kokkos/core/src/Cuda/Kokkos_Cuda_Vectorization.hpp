/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <Kokkos_Cuda.hpp>

namespace Kokkos {


// Shuffle only makes sense on >= Kepler GPUs; it doesn't work on CPUs
// or other GPUs.  We provide a generic definition (which is trivial
// and doesn't do what it claims to do) because we don't actually use
// this function unless we are on a suitable GPU, with a suitable
// Scalar type.  (For example, in the mat-vec, the "ThreadsPerRow"
// internal parameter depends both on the Device and the Scalar type,
// and it controls whether shfl_down() gets called.)
template<typename Scalar>
KOKKOS_INLINE_FUNCTION
Scalar shfl_down(const Scalar &val, const int& delta, const int& width){
  return val;
}

template<>
KOKKOS_INLINE_FUNCTION
unsigned int shfl_down<unsigned int>(const unsigned int &val, const int& delta, const int& width){
#ifdef __CUDA_ARCH__
  #if (__CUDA_ARCH__ >= 300)
    unsigned int tmp1 = val;
    int tmp = *reinterpret_cast<int*>(&tmp1);
    tmp = __shfl_down(tmp,delta,width);
    return *reinterpret_cast<unsigned int*>(&tmp);
  #else
    return val;
  #endif
#else
  return val;
#endif
}

template<>
KOKKOS_INLINE_FUNCTION
int shfl_down<int>(const int &val, const int& delta, const int& width){
#ifdef __CUDA_ARCH__
  #if (__CUDA_ARCH__ >= 300)
    return __shfl_down(val,delta,width);
  #else
    return val;
  #endif
#else
  return val;
#endif
}

template<>
KOKKOS_INLINE_FUNCTION
float shfl_down<float>(const float &val, const int& delta, const int& width){
#ifdef __CUDA_ARCH__
  #if (__CUDA_ARCH__ >= 300)
    return __shfl_down(val,delta,width);
  #else
    return val;
  #endif
#else
  return val;
#endif
}

template<>
KOKKOS_INLINE_FUNCTION
double shfl_down<double>(const double &val, const int& delta, const int& width){
#ifdef __CUDA_ARCH__
  #if (__CUDA_ARCH__ >= 300)
    int lo = __double2loint(val);
    int hi = __double2hiint(val);
    lo = __shfl_down(lo,delta,width);
    hi = __shfl_down(hi,delta,width);
    return __hiloint2double(hi,lo);
  #else
    return val;
  #endif
#else
  return val;
#endif
}

template<>
KOKKOS_INLINE_FUNCTION
long int shfl_down<long int>(const long int &val, const int& delta, const int& width){
#ifdef __CUDA_ARCH__
  #if (__CUDA_ARCH__ >= 300)
    int lo = __double2loint(*reinterpret_cast<const double*>(&val));
    int hi = __double2hiint(*reinterpret_cast<const double*>(&val));
    lo = __shfl_down(lo,delta,width);
    hi = __shfl_down(hi,delta,width);
    const double tmp = __hiloint2double(hi,lo);
    return *(reinterpret_cast<const long int*>(&tmp));
  #else
    return val;
  #endif
#else
  return val;
#endif
}

template<>
KOKKOS_INLINE_FUNCTION
unsigned long shfl_down<unsigned long>(const unsigned long &val, const int& delta, const int& width){
#ifdef __CUDA_ARCH__
  #if (__CUDA_ARCH__ >= 300)
    int lo = __double2loint(*reinterpret_cast<const double*>(&val));
    int hi = __double2hiint(*reinterpret_cast<const double*>(&val));
    lo = __shfl_down(lo,delta,width);
    hi = __shfl_down(hi,delta,width);
    const double tmp = __hiloint2double(hi,lo);
    return *(reinterpret_cast<const unsigned long*>(&tmp));
  #else
    return val;
  #endif
#else
  return val;
#endif
}

template<int N>
struct Vectorization<Cuda,N> {
  enum {increment = N};

#ifdef __CUDA_ARCH__
  KOKKOS_FORCEINLINE_FUNCTION
  static int begin() { return threadIdx.x%N;}
#else
  KOKKOS_FORCEINLINE_FUNCTION
  static int begin() { return 0;}
#endif

  KOKKOS_FORCEINLINE_FUNCTION
  static int thread_rank(const Cuda &dev) {
    return dev.team_rank()/increment;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  static int global_thread_rank(const Cuda &dev) {
    return (dev.league_rank()*dev.team_size()+dev.team_rank())/increment;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  static bool is_lane_0(const Cuda &dev) {
    return (dev.team_rank()%increment)==0;
  }

  template<class Scalar>
  KOKKOS_INLINE_FUNCTION
  static Scalar reduce(const Scalar& val) {
    #ifdef __CUDA_ARCH__
    __shared__ Scalar result[256];
    Scalar myresult;
    for(int k=0;k<blockDim.x;k+=256) {
      const int tid = threadIdx.x - k;
      if(tid > 0 && tid<256) {
        result[tid] = val;
        if ( (N > 1) && (tid%2==0) )
          result[tid] += result[tid+1];
        if ( (N > 2) && (tid%4==0) )
          result[tid] += result[tid+2];
        if ( (N > 4) && (tid%8==0) )
          result[tid] += result[tid+4];
        if ( (N > 8) && (tid%16==0) )
          result[tid] += result[tid+8];
        if ( (N > 16) && (tid%32==0) )
          result[tid] += result[tid+16];
        myresult = result[tid];
      }
      if(blockDim.x>256)
        __syncthreads();
    }
    return myresult;
    #else
    return val;
    #endif
  }

#ifdef __CUDA_ARCH__
  #if (__CUDA_ARCH__ >= 300)
  KOKKOS_INLINE_FUNCTION
  static int reduce(const int& val) {
    int result = val;
    if (N > 1)
      result += shfl_down(result, 1,N);
    if (N > 2)
      result += shfl_down(result, 2,N);
    if (N > 4)
      result += shfl_down(result, 4,N);
    if (N > 8)
      result += shfl_down(result, 8,N);
    if (N > 16)
      result += shfl_down(result, 16,N);
    return result;
  }

  KOKKOS_INLINE_FUNCTION
  static unsigned int reduce(const unsigned int& val) {
    unsigned int result = val;
    if (N > 1)
      result += shfl_down(result, 1,N);
    if (N > 2)
      result += shfl_down(result, 2,N);
    if (N > 4)
      result += shfl_down(result, 4,N);
    if (N > 8)
      result += shfl_down(result, 8,N);
    if (N > 16)
      result += shfl_down(result, 16,N);
    return result;
  }

  KOKKOS_INLINE_FUNCTION
  static long int reduce(const long int& val) {
    long int result = val;
    if (N > 1)
      result += shfl_down(result, 1,N);
    if (N > 2)
      result += shfl_down(result, 2,N);
    if (N > 4)
      result += shfl_down(result, 4,N);
    if (N > 8)
      result += shfl_down(result, 8,N);
    if (N > 16)
      result += shfl_down(result, 16,N);
    return result;
  }

  KOKKOS_INLINE_FUNCTION
  static unsigned long int reduce(const unsigned long int& val) {
    unsigned long int result = val;
    if (N > 1)
      result += shfl_down(result, 1,N);
    if (N > 2)
      result += shfl_down(result, 2,N);
    if (N > 4)
      result += shfl_down(result, 4,N);
    if (N > 8)
      result += shfl_down(result, 8,N);
    if (N > 16)
      result += shfl_down(result, 16,N);
    return result;
  }

  KOKKOS_INLINE_FUNCTION
  static float reduce(const float& val) {
    float result = val;
    if (N > 1)
      result += shfl_down(result, 1,N);
    if (N > 2)
      result += shfl_down(result, 2,N);
    if (N > 4)
      result += shfl_down(result, 4,N);
    if (N > 8)
      result += shfl_down(result, 8,N);
    if (N > 16)
      result += shfl_down(result, 16,N);
    return result;
  }

  KOKKOS_INLINE_FUNCTION
  static double reduce(const double& val) {
    double result = val;
    if (N > 1)
      result += shfl_down(result, 1,N);
    if (N > 2)
      result += shfl_down(result, 2,N);
    if (N > 4)
      result += shfl_down(result, 4,N);
    if (N > 8)
      result += shfl_down(result, 8,N);
    if (N > 16)
      result += shfl_down(result, 16,N);
    return result;
  }
  #endif
#endif

};
}

