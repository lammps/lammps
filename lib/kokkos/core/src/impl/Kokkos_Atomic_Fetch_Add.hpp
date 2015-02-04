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

#if defined( KOKKOS_ATOMIC_HPP ) && ! defined( KOKKOS_ATOMIC_FETCH_ADD_HPP )
#define KOKKOS_ATOMIC_FETCH_ADD_HPP

namespace Kokkos {

//----------------------------------------------------------------------------

#if defined( KOKKOS_ATOMICS_USE_CUDA )

// Support for int, unsigned int, unsigned long long int, and float

__inline__ __device__
int atomic_fetch_add( volatile int * const dest , const int val )
{ return atomicAdd((int*)dest,val); }

__inline__ __device__
unsigned int atomic_fetch_add( volatile unsigned int * const dest , const unsigned int val )
{ return atomicAdd((unsigned int*)dest,val); }

__inline__ __device__
unsigned long long int atomic_fetch_add( volatile unsigned long long int * const dest ,
                                         const unsigned long long int val )
{ return atomicAdd((unsigned long long int*)dest,val); }

__inline__ __device__
float atomic_fetch_add( volatile float * const dest , const float val )
{ return atomicAdd((float*)dest,val); }

template < typename T >
__inline__ __device__
T atomic_fetch_add( volatile T * const dest ,
  typename Kokkos::Impl::enable_if< sizeof(T) == sizeof(int) , const T >::type val )
{
#ifdef KOKKOS_HAVE_CXX11
  union U {
    int i ;
    T t ;
    KOKKOS_INLINE_FUNCTION U() {};
  } assume , oldval , newval ;
#else
  union U {
    int i ;
    T t ;
  } assume , oldval , newval ;
#endif

  oldval.t = *dest ;

  do {
    assume.i = oldval.i ;
    newval.t = assume.t + val ;
    oldval.i = atomicCAS( (int*)dest , assume.i , newval.i );
  } while ( assumed.i != oldval.i );

  return oldval.t ;
}

template < typename T >
__inline__ __device__
T atomic_fetch_add( volatile T * const dest ,
  typename Kokkos::Impl::enable_if< sizeof(T) != sizeof(int) &&
                                    sizeof(T) == sizeof(unsigned long long int) , const T >::type val )
{
#ifdef KOKKOS_HAVE_CXX11
  union U {
    unsigned long long int i ;
    T t ;
    KOKKOS_INLINE_FUNCTION U() {};
  } assume , oldval , newval ;
#else
  union U {
    unsigned long long int i ;
    T t ;
  } assume , oldval , newval ;
#endif

  oldval.t = *dest ;

  do {
    assume.i = oldval.i ;
    newval.t = assume.t + val ;
    oldval.i = atomicCAS( (unsigned long long int*)dest , assume.i , newval.i );
  } while ( assume.i != oldval.i );

  return oldval.t ;
}

template < typename T >
__inline__ __device__
T atomic_fetch_add( volatile T * const dest ,
  typename Kokkos::Impl::enable_if< sizeof(T) != sizeof(int) &&
                                    sizeof(T) != sizeof(unsigned long long int) &&
                                    sizeof(T) == sizeof(Impl::cas128_t), const T >::type val )
{
  Kokkos::abort("Error: calling atomic_fetch_add with 128bit type is not supported on CUDA execution space.");
  return T();
}

//----------------------------------------------------------------------------

#elif defined(KOKKOS_ATOMICS_USE_GCC) || defined(KOKKOS_ATOMICS_USE_INTEL)

KOKKOS_INLINE_FUNCTION
int atomic_fetch_add( volatile int * const dest , const int val )
{ return __sync_fetch_and_add(dest,val); }

KOKKOS_INLINE_FUNCTION
long int atomic_fetch_add( volatile long int * const dest , const long int val )
{ return __sync_fetch_and_add(dest,val); }

#if defined( KOKKOS_ATOMICS_USE_GCC )

KOKKOS_INLINE_FUNCTION
unsigned int atomic_fetch_add( volatile unsigned int * const dest , const unsigned int val )
{ return __sync_fetch_and_add(dest,val); }

KOKKOS_INLINE_FUNCTION
unsigned long int atomic_fetch_add( volatile unsigned long int * const dest , const unsigned long int val )
{ return __sync_fetch_and_add(dest,val); }

#endif

template < typename T >
KOKKOS_INLINE_FUNCTION
T atomic_fetch_add( volatile T * const dest ,
  typename Kokkos::Impl::enable_if< sizeof(T) == sizeof(int) , const T >::type val )
{
#ifdef KOKKOS_HAVE_CXX11
  union U {
    int i ;
    T t ;
    KOKKOS_INLINE_FUNCTION U() {};
  } assume , oldval , newval ;
#else
  union U {
    int i ;
    T t ;
  } assume , oldval , newval ;
#endif

  oldval.t = *dest ;

  do {
    assume.i = oldval.i ;
    newval.t = assume.t + val ;
    oldval.i = __sync_val_compare_and_swap( (int*) dest , assume.i , newval.i );
  } while ( assume.i != oldval.i );

  return oldval.t ;
}

template < typename T >
KOKKOS_INLINE_FUNCTION
T atomic_fetch_add( volatile T * const dest ,
  typename Kokkos::Impl::enable_if< sizeof(T) != sizeof(int) &&
                                    sizeof(T) == sizeof(long) , const T >::type val )
{
#ifdef KOKKOS_HAVE_CXX11
  union U {
    long i ;
    T t ;
    KOKKOS_INLINE_FUNCTION U() {};
  } assume , oldval , newval ;
#else
  union U {
    long i ;
    T t ;
  } assume , oldval , newval ;
#endif

  oldval.t = *dest ;

  do {
    assume.i = oldval.i ;
    newval.t = assume.t + val ;
    oldval.i = __sync_val_compare_and_swap( (long*) dest , assume.i , newval.i );
  } while ( assume.i != oldval.i );

  return oldval.t ;
}

template < typename T >
KOKKOS_INLINE_FUNCTION
T atomic_fetch_add( volatile T * const dest ,
  typename Kokkos::Impl::enable_if< sizeof(T) != sizeof(int) &&
                                    sizeof(T) != sizeof(long) &&
                                    sizeof(T) == sizeof(Impl::cas128_t) , const T >::type val )
{
#ifdef KOKKOS_HAVE_CXX11
  union U {
    Impl::cas128_t i ;
    T t ;
    KOKKOS_INLINE_FUNCTION U() {};
  } assume , oldval , newval ;
#else
  union U {
    Impl::cas128_t i ;
    T t ;
  } assume , oldval , newval ;
#endif

  oldval.t = *dest ;

  do {
    assume.i = oldval.i ;
    newval.t = assume.t + val ;
    oldval.i = Impl::cas128( (volatile Impl::cas128_t*) dest , assume.i , newval.i );
  } while ( assume.i != oldval.i );

  return oldval.t ;
}
//----------------------------------------------------------------------------

#elif defined( KOKKOS_ATOMICS_USE_OMP31 )

template< typename T >
T atomic_fetch_add( volatile T * const dest , const T val )
{
  T retval;
#pragma omp atomic capture
  {
    retval = dest[0];
    dest[0] += val;
  }
  return retval;
}

#endif

//----------------------------------------------------------------------------

// Simpler version of atomic_fetch_add without the fetch
template <typename T>
KOKKOS_INLINE_FUNCTION
void atomic_add(volatile T * const dest, const T src) {
  atomic_fetch_add(dest,src);
}

// Atomic increment
template<typename T>
KOKKOS_INLINE_FUNCTION
void atomic_increment(volatile T* a) {
  Kokkos::atomic_fetch_add(a,1);
}

template<typename T>
KOKKOS_INLINE_FUNCTION
void atomic_decrement(volatile T* a) {
  Kokkos::atomic_fetch_add(a,-1);
}

}
#endif

