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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <memory.h>
#include <stddef.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <cstring>

#include <Kokkos_HostSpace.hpp>
#include <impl/Kokkos_MemoryTracking.hpp>
#include <impl/Kokkos_Error.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {
namespace {

Impl::MemoryTracking<> & host_space_singleton()
{
  static Impl::MemoryTracking<> self("Kokkos::HostSpace");
  return self ;
}

} // namespace <blank>
} // namespace Impl
} // namespade Kokkos

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

void * host_allocate_not_thread_safe( const std::string & label , const size_t size )
{
  void * ptr = 0 ;

  if ( size ) {
    size_t size_padded = size ;
    void * ptr_alloc = 0 ;

#if defined( __INTEL_COMPILER ) && !defined ( KOKKOS_HAVE_CUDA )

    ptr = ptr_alloc = _mm_malloc( size , MEMORY_ALIGNMENT );

#elif ( defined( _POSIX_C_SOURCE ) && _POSIX_C_SOURCE >= 200112L ) || \
      ( defined( _XOPEN_SOURCE )   && _XOPEN_SOURCE   >= 600 )

    posix_memalign( & ptr_alloc , MEMORY_ALIGNMENT , size );
    ptr = ptr_alloc ;

#else

    {
      // Over-allocate to and round up to guarantee proper alignment.

      size_padded = ( size + MEMORY_ALIGNMENT - 1 );

      ptr_alloc = malloc( size_padded );

      const size_t rem = reinterpret_cast<ptrdiff_t>(ptr_alloc) % MEMORY_ALIGNMENT ;

      ptr = static_cast<unsigned char *>(ptr_alloc) + ( rem ? MEMORY_ALIGNMENT - rem : 0 );
    }

#endif

    if ( ptr_alloc && ptr_alloc <= ptr &&
         0 == ( reinterpret_cast<ptrdiff_t>(ptr) % MEMORY_ALIGNMENT ) ) {
      // Insert allocated pointer and allocation count
      Impl::host_space_singleton().insert( label , ptr_alloc , size_padded );
    }
    else {
      std::ostringstream msg ;
      msg << "Kokkos::Impl::host_allocate_not_thread_safe( "
          << label
          << " , " << size
          << " ) FAILED aligned memory allocation" ;
      Kokkos::Impl::throw_runtime_exception( msg.str() );
    }
  }

  return ptr ;
}

void host_decrement_not_thread_safe( const void * ptr )
{
  void * ptr_alloc = Impl::host_space_singleton().decrement( ptr );

  if ( ptr_alloc ) {
#if defined( __INTEL_COMPILER ) && !defined ( KOKKOS_HAVE_CUDA )
     _mm_free( ptr_alloc );
#else
     free( ptr_alloc );
#endif
  }
}

DeepCopy<HostSpace,HostSpace>::DeepCopy( void * dst , const void * src , size_t n )
{
  memcpy( dst , src , n );
}

}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace {

static const int QUERY_SPACE_IN_PARALLEL_MAX = 16 ;

typedef int (* QuerySpaceInParallelPtr )();

QuerySpaceInParallelPtr s_in_parallel_query[ QUERY_SPACE_IN_PARALLEL_MAX ] ;
int s_in_parallel_query_count = 0 ;

} // namespace <empty>

void HostSpace::register_in_parallel( int (*device_in_parallel)() )
{
  if ( 0 == device_in_parallel ) {
    Kokkos::Impl::throw_runtime_exception( std::string("Kokkos::HostSpace::register_in_parallel ERROR : given NULL" ) );
  }

  int i = -1 ;

  if ( ! (device_in_parallel)() ) {
    for ( i = 0 ; i < s_in_parallel_query_count && ! (*(s_in_parallel_query[i]))() ; ++i );
  }

  if ( i < s_in_parallel_query_count ) {
    Kokkos::Impl::throw_runtime_exception( std::string("Kokkos::HostSpace::register_in_parallel_query ERROR : called in_parallel" ) );

  }

  if ( QUERY_SPACE_IN_PARALLEL_MAX <= i ) {
    Kokkos::Impl::throw_runtime_exception( std::string("Kokkos::HostSpace::register_in_parallel_query ERROR : exceeded maximum" ) );

  }

  for ( i = 0 ; i < s_in_parallel_query_count && s_in_parallel_query[i] != device_in_parallel ; ++i );

  if ( i == s_in_parallel_query_count ) {
    s_in_parallel_query[s_in_parallel_query_count++] = device_in_parallel ;
  }
}

int HostSpace::in_parallel()
{
  const int n = s_in_parallel_query_count ;

  int i = 0 ;

  while ( i < n && ! (*(s_in_parallel_query[i]))() ) { ++i ; }

  return i < n ;
}

} // namespace Kokkos

/*--------------------------------------------------------------------------*/

namespace Kokkos {

void * HostSpace::allocate( const std::string & label , const size_t size )
{
  void * ptr = 0 ;

  if ( ! HostSpace::in_parallel() ) {
    ptr = Impl::host_allocate_not_thread_safe( label , size );
  }
  else {
    Kokkos::Impl::throw_runtime_exception( std::string("Kokkos::HostSpace::allocate called within a parallel functor") );
  }

  return ptr ;
}

void HostSpace::increment( const void * ptr )
{
  if ( ! HostSpace::in_parallel() ) {
    Impl::host_space_singleton().increment( ptr );
  }
  else {
    Kokkos::Impl::throw_runtime_exception( std::string("Kokkos::HostSpace::increment called within a parallel functor") );
  }
}

void HostSpace::decrement( const void * ptr )
{
  if ( ! HostSpace::in_parallel() ) {
    Impl::host_decrement_not_thread_safe( ptr );
  }
  else {
    Kokkos::Impl::throw_runtime_exception( std::string("Kokkos::HostSpace::decrement called within a parallel functor") );
  }
}

int HostSpace::count( const void * ptr ) {
  if ( ! HostSpace::in_parallel() ) {
    Impl::MemoryTracking<>::Entry * const entry =
        Impl::host_space_singleton().query(ptr);
    return entry != NULL?entry->count():0;
  }
  else {
    Kokkos::Impl::throw_runtime_exception( std::string("Kokkos::HostSpace::count called within a parallel functor") );
    return -1;
  }
}

void HostSpace::print_memory_view( std::ostream & o )
{
  Impl::host_space_singleton().print( o , std::string("  ") );
}

std::string HostSpace::query_label( const void * p )
{
  Impl::MemoryTracking<>::Entry * const entry = Impl::host_space_singleton().query(p);
  return std::string( entry ? entry->label() : "<NOT ALLOCATED>" );
}

} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

