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

#if defined( KOKKOS_ATOMIC_HPP ) && ! defined( KOKKOS_ATOMIC_EXCHANGE_HPP )
#define KOKKOS_ATOMIC_EXCHANGE_HPP

namespace Kokkos {

//----------------------------------------------------------------------------

#if defined( KOKKOS_ATOMICS_USE_CUDA )

__inline__ __device__
int atomic_exchange( volatile int * const dest , const int val )
{
  // return __iAtomicExch( (int*) dest , val );
  return atomicExch( (int*) dest , val );
}

__inline__ __device__
unsigned int atomic_exchange( volatile unsigned int * const dest , const unsigned int val )
{
  // return __uAtomicExch( (unsigned int*) dest , val );
  return atomicExch( (unsigned int*) dest , val );
}

__inline__ __device__
unsigned long long int atomic_exchange( volatile unsigned long long int * const dest , const unsigned long long int val )
{
  // return __ullAtomicExch( (unsigned long long*) dest , val );
  return atomicExch( (unsigned long long*) dest , val );
}

/** \brief  Atomic exchange for any type with compatible size */
template< typename T >
__inline__ __device__
T atomic_exchange(
  volatile T * const dest ,
  typename Kokkos::Impl::enable_if< sizeof(T) == sizeof(int) , const T & >::type val )
{
  // int tmp = __ullAtomicExch( (int*) dest , *((int*)&val) );
  int tmp = atomicExch( ((int*)dest) , *((int*)&val) );
  return *((T*)&tmp);
}

template< typename T >
__inline__ __device__
T atomic_exchange(
  volatile T * const dest ,
  typename Kokkos::Impl::enable_if< sizeof(T) != sizeof(int) &&
                                    sizeof(T) == sizeof(unsigned long long int) , const T & >::type val )
{
  typedef unsigned long long int type ;
  // type tmp = __ullAtomicExch( (type*) dest , *((type*)&val) );
  type tmp = atomicExch( ((type*)dest) , *((type*)&val) );
  return *((T*)&tmp);
}

/** \brief  Atomic exchange for any type with compatible size */
template< typename T >
__inline__ __device__
void atomic_assign(
  volatile T * const dest ,
  typename Kokkos::Impl::enable_if< sizeof(T) == sizeof(int) , const T & >::type val )
{
  // (void) __ullAtomicExch( (int*) dest , *((int*)&val) );
  (void) atomicExch( ((int*)dest) , *((int*)&val) );
}

template< typename T >
__inline__ __device__
void atomic_assign(
  volatile T * const dest ,
  typename Kokkos::Impl::enable_if< sizeof(T) != sizeof(int) &&
                                    sizeof(T) == sizeof(unsigned long long int) , const T & >::type val )
{
  typedef unsigned long long int type ;
  // (void) __ullAtomicExch( (type*) dest , *((type*)&val) );
  (void) atomicExch( ((type*)dest) , *((type*)&val) );
}

//----------------------------------------------------------------------------

#elif defined(KOKKOS_ATOMICS_USE_GCC) || defined(KOKKOS_ATOMICS_USE_INTEL)

template< typename T >
KOKKOS_INLINE_FUNCTION
T atomic_exchange( volatile T * const dest ,
  typename Kokkos::Impl::enable_if< sizeof(T) == sizeof(int) || sizeof(T) == sizeof(long)
                                  , const T & >::type val )
{
  typedef typename Kokkos::Impl::if_c< sizeof(T) == sizeof(int) , int , long >::type type ;

  const type v = *((type*)&val); // Extract to be sure the value doesn't change
  
  type assumed ;

  union { T val_T ; type val_type ; } old ;

  old.val_T = *dest ;

  do {
    assumed = old.val_type ;
    old.val_type = __sync_val_compare_and_swap( (volatile type *) dest , assumed , v );
  } while ( assumed != old.val_type );

  return old.val_T ;
}

template< typename T >
KOKKOS_INLINE_FUNCTION
void atomic_assign( volatile T * const dest ,
  typename Kokkos::Impl::enable_if< sizeof(T) == sizeof(int) || sizeof(T) == sizeof(long)
                                  , const T & >::type val )
{
  typedef typename Kokkos::Impl::if_c< sizeof(T) == sizeof(int) , int , long >::type type ;

  const type v = *((type*)&val); // Extract to be sure the value doesn't change
  
  type assumed ;

  union { T val_T ; type val_type ; } old ;

  old.val_T = *dest ;

  do {
    assumed = old.val_type ;
    old.val_type = __sync_val_compare_and_swap( (volatile type *) dest , assumed , v );
  } while ( assumed != old.val_type );
}

//----------------------------------------------------------------------------

#elif defined( KOKKOS_ATOMICS_USE_OMP31 )

template < typename T >
KOKKOS_INLINE_FUNCTION
T atomic_exchange( volatile T * const dest , const T val )
{
  T retval;
#pragma omp critical
  {
    retval = dest[0];
    dest[0] = val;
  }
  return retval;
}

template < typename T >
KOKKOS_INLINE_FUNCTION
void atomic_assign( volatile T * const dest , const T val )
{
#pragma omp critical
  {
    dest[0] = val;
  }
}

#endif

//----------------------------------------------------------------------------

} // namespace Kokkos

#endif

//----------------------------------------------------------------------------

