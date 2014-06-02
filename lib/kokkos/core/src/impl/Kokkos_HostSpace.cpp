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
#include <impl/Kokkos_Error.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace {

class HostMemoryTrackingEntry : public Impl::MemoryTrackingEntry
{
public:

  void * const ptr_alloc ;

  HostMemoryTrackingEntry( const std::string & arg_label ,
                           const std::type_info & arg_info ,
                           void * const           arg_ptr ,
                           const unsigned         arg_size )
    : Impl::MemoryTrackingEntry( arg_label , arg_info , arg_ptr , arg_size )
    , ptr_alloc( arg_ptr )
    {}

  ~HostMemoryTrackingEntry();
};

HostMemoryTrackingEntry::~HostMemoryTrackingEntry()
{
#if defined( __INTEL_COMPILER ) && !defined ( KOKKOS_HAVE_CUDA )
   _mm_free( ptr_alloc );
#else
   free( ptr_alloc );
#endif
}

Impl::MemoryTracking & host_space_singleton()
{
  static Impl::MemoryTracking self("Kokkos::HostSpace");
  return self ;
}

bool host_space_verify_modifiable( const char * const label )
{
  static const char error_in_parallel[] = "Called with HostSpace::in_parallel()" ;
  static const char error_not_exists[]  = "Called after return from main()" ;

  const char * const error_msg =
    HostSpace::in_parallel() ? error_in_parallel : (
    ! host_space_singleton().exists() ? error_not_exists : (const char *) 0 );

  if ( error_msg ) {
    std::cerr << "Kokkos::HostSpace::" << label << " ERROR : " << error_msg << std::endl ;
  }

  return error_msg == 0  ;
}

} // namespace <blank>
} // namespade Kokkos

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

void * host_allocate_not_thread_safe(
  const std::string    & label ,
  const std::type_info & scalar_type ,
  const size_t           scalar_size ,
  const size_t           scalar_count )
{
  void * ptr = 0 ;

  if ( 0 < scalar_size && 0 < scalar_count ) {
    void * ptr_alloc = 0 ;
    size_t count_alloc = scalar_count ;

#if defined( __INTEL_COMPILER ) && !defined ( KOKKOS_HAVE_CUDA )

    ptr = ptr_alloc = _mm_malloc( scalar_size * count_alloc , MEMORY_ALIGNMENT );
   
#elif ( defined( _POSIX_C_SOURCE ) && _POSIX_C_SOURCE >= 200112L ) || \
      ( defined( _XOPEN_SOURCE )   && _XOPEN_SOURCE   >= 600 )

    posix_memalign( & ptr_alloc , MEMORY_ALIGNMENT , scalar_size * count_alloc );
    ptr = ptr_alloc ;

#else

    // Over-allocate to guarantee enough aligned space.

    count_alloc += ( MEMORY_ALIGNMENT + scalar_size - 1 ) / scalar_size ;

    ptr_alloc = malloc( scalar_size * count_alloc );

    ptr = static_cast<unsigned char *>(ptr_alloc) + 
          ( MEMORY_ALIGNMENT - reinterpret_cast<ptrdiff_t>(ptr_alloc) % MEMORY_ALIGNMENT );

#endif

    if ( ptr_alloc && ptr_alloc <= ptr &&
         0 == ( reinterpret_cast<ptrdiff_t>(ptr) % MEMORY_ALIGNMENT ) ) {
      host_space_singleton().insert(
        new HostMemoryTrackingEntry( label , scalar_type , ptr_alloc , scalar_size * count_alloc ) );
    }
    else {
      std::ostringstream msg ;
      msg << "Kokkos::Impl::host_allocate_not_thread_safe( "
          << label
          << " , " << scalar_type.name()
          << " , " << scalar_size
          << " , " << scalar_count
          << " ) FAILED aligned memory allocation" ;
      Kokkos::Impl::throw_runtime_exception( msg.str() );
    }
  }

  return ptr ;
}

void host_decrement_not_thread_safe( const void * ptr )
{
  host_space_singleton().decrement( ptr );
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

static const int QUERY_DEVICE_IN_PARALLEL_MAX = 16 ;

typedef int (* QueryDeviceInParallelPtr )();

QueryDeviceInParallelPtr s_in_parallel_query[ QUERY_DEVICE_IN_PARALLEL_MAX ] ;
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

  if ( QUERY_DEVICE_IN_PARALLEL_MAX <= i ) {
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

void * HostSpace::allocate(
  const std::string    & label ,
  const std::type_info & scalar_type ,
  const size_t           scalar_size ,
  const size_t           scalar_count )
{
  void * ptr = 0 ;

  if ( host_space_verify_modifiable("allocate") ) {
    ptr = Impl::host_allocate_not_thread_safe( label , scalar_type , scalar_size , scalar_count );
  }

  return ptr ;
}

void HostSpace::increment( const void * ptr )
{
  if ( host_space_verify_modifiable("increment") ) {
    host_space_singleton().increment( ptr );
  }
}

void HostSpace::decrement( const void * ptr )
{
  if ( host_space_verify_modifiable("decrement") ) {
    Impl::host_decrement_not_thread_safe( ptr );
  }
}

void HostSpace::print_memory_view( std::ostream & o )
{
  host_space_singleton().print( o , std::string("  ") );
}

std::string HostSpace::query_label( const void * p )
{
  const Impl::MemoryTrackingEntry * const info = 
    host_space_singleton().query( p );

  return 0 != info ? info->label : std::string("ERROR NOT DEFINED");
}

} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

