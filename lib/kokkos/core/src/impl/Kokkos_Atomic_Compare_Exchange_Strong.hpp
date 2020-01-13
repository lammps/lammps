/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
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
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#if defined( KOKKOS_ENABLE_RFO_PREFETCH )
#include <xmmintrin.h>
#endif

#include <Kokkos_Macros.hpp>
#if defined( KOKKOS_ATOMIC_HPP ) && ! defined( KOKKOS_ATOMIC_COMPARE_EXCHANGE_STRONG_HPP )
#define KOKKOS_ATOMIC_COMPARE_EXCHANGE_STRONG_HPP

#if defined(KOKKOS_ENABLE_CUDA)
#include<Cuda/Kokkos_Cuda_Version_9_8_Compatibility.hpp>
#endif

#include <impl/Kokkos_Atomic_Memory_Order.hpp>
#include <impl/Kokkos_Memory_Fence.hpp>

#if defined(KOKKOS_ENABLE_CUDA)
#include <Cuda/Kokkos_Cuda_Atomic_Intrinsics.hpp>
#endif

namespace Kokkos {

//----------------------------------------------------------------------------
// Cuda native CAS supports int, unsigned int, and unsigned long long int (non-standard type).
// Must cast-away 'volatile' for the CAS call.

#if defined( KOKKOS_ENABLE_CUDA )

#if defined(__CUDA_ARCH__) || defined(KOKKOS_IMPL_CUDA_CLANG_WORKAROUND)
__inline__ __device__
int atomic_compare_exchange( volatile int * const dest, const int compare, const int val)
{ return atomicCAS((int*)dest,compare,val); }

__inline__ __device__
unsigned int atomic_compare_exchange( volatile unsigned int * const dest, const unsigned int compare, const unsigned int val)
{ return atomicCAS((unsigned int*)dest,compare,val); }

__inline__ __device__
unsigned long long int atomic_compare_exchange( volatile unsigned long long int * const dest ,
                                                const unsigned long long int compare ,
                                                const unsigned long long int val )
{ return atomicCAS((unsigned long long int*)dest,compare,val); }

template < typename T >
__inline__ __device__
T atomic_compare_exchange( volatile T * const dest , const T & compare ,
  typename Kokkos::Impl::enable_if< sizeof(T) == sizeof(int) , const T & >::type val )
{
  const int tmp = atomicCAS( (int*) dest , *((int*)&compare) , *((int*)&val) );
  return *((T*)&tmp);
}

template < typename T >
__inline__ __device__
T atomic_compare_exchange( volatile T * const dest , const T & compare ,
  typename Kokkos::Impl::enable_if< sizeof(T) != sizeof(int) &&
                                    sizeof(T) == sizeof(unsigned long long int) , const T & >::type val )
{
  typedef unsigned long long int type ;
  const type tmp = atomicCAS( (type*) dest , *((type*)&compare) , *((type*)&val) );
  return *((T*)&tmp);
}

template < typename T >
__inline__ __device__
T atomic_compare_exchange( volatile T * const dest , const T & compare ,
    typename Kokkos::Impl::enable_if<
                  ( sizeof(T) != 4 )
               && ( sizeof(T) != 8 )
             , const T >::type& val )
{
  T return_val;
  // This is a way to (hopefully) avoid dead lock in a warp
  int done = 0;
#ifdef KOKKOS_IMPL_CUDA_SYNCWARP_NEEDS_MASK
  unsigned int mask = KOKKOS_IMPL_CUDA_ACTIVEMASK;
  unsigned int active = KOKKOS_IMPL_CUDA_BALLOT_MASK(mask,1);
#else
  unsigned int active = KOKKOS_IMPL_CUDA_BALLOT(1);
#endif
  unsigned int done_active = 0;
  while (active!=done_active) {
    if(!done) {
      if( Impl::lock_address_cuda_space( (void*) dest ) ) {
        return_val = *dest;
        if( return_val == compare )
          *dest = val;
        Impl::unlock_address_cuda_space( (void*) dest );
        done = 1;
      }
    }
#ifdef KOKKOS_IMPL_CUDA_SYNCWARP_NEEDS_MASK
    done_active = KOKKOS_IMPL_CUDA_BALLOT_MASK(mask,done);
#else
    done_active = KOKKOS_IMPL_CUDA_BALLOT(done);
#endif
  }
  return return_val;
}
#endif
#endif

//----------------------------------------------------------------------------
// GCC native CAS supports int, long, unsigned int, unsigned long.
// Intel native CAS support int and long with the same interface as GCC.
#if !defined(KOKKOS_ENABLE_ROCM_ATOMICS)
#if !defined(__CUDA_ARCH__) || defined(KOKKOS_IMPL_CUDA_CLANG_WORKAROUND)
#if defined(KOKKOS_ENABLE_GNU_ATOMICS) || defined(KOKKOS_ENABLE_INTEL_ATOMICS)

inline
int atomic_compare_exchange( volatile int * const dest, const int compare, const int val)
{
#if defined( KOKKOS_ENABLE_RFO_PREFETCH )
  _mm_prefetch( (const char*) dest, _MM_HINT_ET0 );
#endif
  return __sync_val_compare_and_swap(dest,compare,val);
}

inline
long atomic_compare_exchange( volatile long * const dest, const long compare, const long val )
{ 
#if defined( KOKKOS_ENABLE_RFO_PREFETCH )
  _mm_prefetch( (const char*) dest, _MM_HINT_ET0 );
#endif
  return __sync_val_compare_and_swap(dest,compare,val);
}

#if defined( KOKKOS_ENABLE_GNU_ATOMICS )

// GCC supports unsigned

inline
unsigned int atomic_compare_exchange( volatile unsigned int * const dest, const unsigned int compare, const unsigned int val )
{ return __sync_val_compare_and_swap(dest,compare,val); }

inline
unsigned long atomic_compare_exchange( volatile unsigned long * const dest ,
                                       const unsigned long compare ,
                                       const unsigned long val )
{ return __sync_val_compare_and_swap(dest,compare,val); }

#endif

template < typename T >
inline
T atomic_compare_exchange( volatile T * const dest, const T & compare,
  typename Kokkos::Impl::enable_if< sizeof(T) == sizeof(int) , const T & >::type val )
{
  union U {
    int i ;
    T t ;
    KOKKOS_INLINE_FUNCTION U() {};
  } tmp ;

#if defined( KOKKOS_ENABLE_RFO_PREFETCH )
  _mm_prefetch( (const char*) dest, _MM_HINT_ET0 );
#endif

  tmp.i = __sync_val_compare_and_swap( (int*) dest , *((int*)&compare) , *((int*)&val) );
  return tmp.t ;
}

template < typename T >
inline
T atomic_compare_exchange( volatile T * const dest, const T & compare,
  typename Kokkos::Impl::enable_if< sizeof(T) != sizeof(int) &&
                                    sizeof(T) == sizeof(long) , const T & >::type val )
{
  union U {
    long i ;
    T t ;
    KOKKOS_INLINE_FUNCTION U() {};
  } tmp ;

#if defined( KOKKOS_ENABLE_RFO_PREFETCH )
  _mm_prefetch( (const char*) dest, _MM_HINT_ET0 );
#endif

  tmp.i = __sync_val_compare_and_swap( (long*) dest , *((long*)&compare) , *((long*)&val) );
  return tmp.t ;
}

#if defined( KOKKOS_ENABLE_ASM) && defined ( KOKKOS_ENABLE_ISA_X86_64 )
template < typename T >
inline
T atomic_compare_exchange( volatile T * const dest, const T & compare,
  typename Kokkos::Impl::enable_if< sizeof(T) != sizeof(int) &&
                                    sizeof(T) != sizeof(long) &&
                                    sizeof(T) == sizeof(Impl::cas128_t), const T & >::type val )
{
  union U {
    Impl::cas128_t i ;
    T t ;
    KOKKOS_INLINE_FUNCTION U() {};
  } tmp ;

#if defined( KOKKOS_ENABLE_RFO_PREFETCH )
  _mm_prefetch( (const char*) dest, _MM_HINT_ET0 );
#endif

  tmp.i = Impl::cas128( (Impl::cas128_t*) dest , *((Impl::cas128_t*)&compare) , *((Impl::cas128_t*)&val) );
  return tmp.t ;
}
#endif

template < typename T >
inline
T atomic_compare_exchange( volatile T * const dest , const T compare ,
    typename Kokkos::Impl::enable_if<
                  ( sizeof(T) != 4 )
               && ( sizeof(T) != 8 )
            #if defined(KOKKOS_ENABLE_ASM) && defined ( KOKKOS_ENABLE_ISA_X86_64 )
               && ( sizeof(T) != 16 )
            #endif
             , const T >::type& val )
{
#if defined( KOKKOS_ENABLE_RFO_PREFETCH )
  _mm_prefetch( (const char*) dest, _MM_HINT_ET0 );
#endif

  while( !Impl::lock_address_host_space( (void*) dest ) );
  T return_val = *dest;
  if( return_val == compare ) {
    // Don't use the following line of code here:
    //
    //const T tmp = *dest = val;
    //
    // Instead, put each assignment in its own statement.  This is
    // because the overload of T::operator= for volatile *this should
    // return void, not volatile T&.  See Kokkos #177:
    //
    // https://github.com/kokkos/kokkos/issues/177
    *dest = val;
    const T tmp = *dest;
    #ifndef KOKKOS_COMPILER_CLANG
    (void) tmp;
    #endif
  }
  Impl::unlock_address_host_space( (void*) dest );
  return return_val;
}
//----------------------------------------------------------------------------

#elif defined( KOKKOS_ENABLE_OPENMP_ATOMICS )

template< typename T >
KOKKOS_INLINE_FUNCTION
T atomic_compare_exchange( volatile T * const dest, const T compare, const T val )
{
  T retval;
#pragma omp critical
  {
    retval = dest[0];
    if ( retval == compare )
        dest[0] = val;
  }
  return retval;
}

#elif defined( KOKKOS_ENABLE_SERIAL_ATOMICS )

template< typename T >
KOKKOS_INLINE_FUNCTION
T atomic_compare_exchange( volatile T * const dest_v, const T compare, const T val )
{
  T* dest = const_cast<T*>(dest_v);
  T retval = *dest;
  if (retval == compare) *dest = val;
  return retval;
}

#endif
#endif
#endif // !defined ROCM_ATOMICS

// dummy for non-CUDA Kokkos headers being processed by NVCC
#if defined(__CUDA_ARCH__) && !defined(KOKKOS_ENABLE_CUDA)
template <typename T>
__inline__ __device__
T atomic_compare_exchange(volatile T * const, const Kokkos::Impl::identity_t<T>, const Kokkos::Impl::identity_t<T>)
{
  return T();
}
#endif

template <typename T>
KOKKOS_INLINE_FUNCTION
bool atomic_compare_exchange_strong(volatile T* const dest, const T compare, const T val)
{
  return compare == atomic_compare_exchange(dest, compare, val);
}
//----------------------------------------------------------------------------

namespace Impl {
// memory-ordered versions are in the Impl namespace

template <class T, class MemoryOrderFailure>
KOKKOS_INLINE_FUNCTION
bool _atomic_compare_exchange_strong_fallback(
  T* dest, T compare, T val, memory_order_seq_cst_t, MemoryOrderFailure
)
{
  Kokkos::memory_fence();
  auto rv = Kokkos::atomic_compare_exchange_strong(
    dest, compare, val
  );
  Kokkos::memory_fence();
  return rv;
}

template <class T, class MemoryOrderFailure>
KOKKOS_INLINE_FUNCTION
bool _atomic_compare_exchange_strong_fallback(
  T* dest, T compare, T val, memory_order_acquire_t, MemoryOrderFailure
)
{
  auto rv = Kokkos::atomic_compare_exchange_strong(
    dest, compare, val
  );
  Kokkos::memory_fence();
  return rv;
}

template <class T, class MemoryOrderFailure>
KOKKOS_INLINE_FUNCTION
bool _atomic_compare_exchange_strong_fallback(
  T* dest, T compare, T val, memory_order_release_t, MemoryOrderFailure
)
{
  Kokkos::memory_fence();
  return Kokkos::atomic_compare_exchange_strong(
    dest, compare, val
  );
}

template <class T, class MemoryOrderFailure>
KOKKOS_INLINE_FUNCTION
bool _atomic_compare_exchange_strong_fallback(
  T* dest, T compare, T val, memory_order_relaxed_t, MemoryOrderFailure
)
{
  return Kokkos::atomic_compare_exchange_strong(
    dest, compare, val
  );
}

#if (defined(KOKKOS_ENABLE_GNU_ATOMICS) && !defined(__CUDA_ARCH__)) \
    || (defined(KOKKOS_ENABLE_INTEL_ATOMICS) && !defined(__CUDA_ARCH__)) \
    || defined(KOKKOS_ENABLE_CUDA_ASM_ATOMICS)

#if defined(__CUDA_ARCH__)
  #define KOKKOS_INTERNAL_INLINE_DEVICE_IF_CUDA_ARCH __inline__ __device__
#else
  #define KOKKOS_INTERNAL_INLINE_DEVICE_IF_CUDA_ARCH inline
#endif

template <class T, class MemoryOrderSuccess, class MemoryOrderFailure>
KOKKOS_INTERNAL_INLINE_DEVICE_IF_CUDA_ARCH
bool _atomic_compare_exchange_strong(
  T* dest, T compare, T val,
  MemoryOrderSuccess,
  MemoryOrderFailure,
  typename std::enable_if<
    (
      sizeof(T) == 1
      || sizeof(T) == 2
      || sizeof(T) == 4
      || sizeof(T) == 8
    )
    && std::is_same<
      typename MemoryOrderSuccess::memory_order,
      typename std::remove_cv<MemoryOrderSuccess>::type
    >::value
    && std::is_same<
      typename MemoryOrderFailure::memory_order,
      typename std::remove_cv<MemoryOrderFailure>::type
    >::value,
    void const**
  >::type = nullptr
) {
  return __atomic_compare_exchange_n(
    dest, &compare, val, /* weak = */ false,
    MemoryOrderSuccess::gnu_constant,
    MemoryOrderFailure::gnu_constant
  );
}

template <class T, class MemoryOrderSuccess, class MemoryOrderFailure>
KOKKOS_INTERNAL_INLINE_DEVICE_IF_CUDA_ARCH
bool _atomic_compare_exchange_strong(
  T* dest, T compare, T val,
  MemoryOrderSuccess order_success,
  MemoryOrderFailure order_failure,
  typename std::enable_if<
    !(
      sizeof(T) == 1
      || sizeof(T) == 2
      || sizeof(T) == 4
      || sizeof(T) == 8
    )
    && std::is_same<
      typename MemoryOrderSuccess::memory_order,
      typename std::remove_cv<MemoryOrderSuccess>::type
    >::value
    && std::is_same<
      typename MemoryOrderFailure::memory_order,
      typename std::remove_cv<MemoryOrderFailure>::type
    >::value,
    void const**
  >::type = nullptr
) {
  return _atomic_compare_exchange_fallback(
    dest, compare, val,
    order_success, order_failure
  );
}

#else

template <class T, class MemoryOrderSuccess, class MemoryOrderFailure>
KOKKOS_INLINE_FUNCTION
bool _atomic_compare_exchange_strong(
  T* dest, T compare, T val,
  MemoryOrderSuccess order_success,
  MemoryOrderFailure order_failure
) {
  return _atomic_compare_exchange_strong_fallback(
    dest, compare, val, order_success, order_failure
  );
}

#endif

// TODO static asserts in overloads that don't make sense (as listed in https://gcc.gnu.org/onlinedocs/gcc-5.2.0/gcc/_005f_005fatomic-Builtins.html)
template <class T, class MemoryOrderSuccess, class MemoryOrderFailure>
KOKKOS_FORCEINLINE_FUNCTION
bool atomic_compare_exchange_strong(
  T* dest, T compare, T val,
  MemoryOrderSuccess order_success,
  MemoryOrderFailure order_failure
) {
  return _atomic_compare_exchange_strong(dest, compare, val, order_success, order_failure);
}


} // end namespace Impl

} // namespace Kokkos

#if defined(KOKKOS_ENABLE_CUDA)
#include <Cuda/Kokkos_Cuda_Atomic_Intrinsics_Restore_Builtins.hpp>
#endif

#endif

