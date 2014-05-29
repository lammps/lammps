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

#include <stdio.h>
#include <limits>
#include <iostream>
#include <Kokkos_OpenMP.hpp>
#include <Kokkos_hwloc.hpp>
#include <iostream>

namespace Kokkos {
namespace Impl {
namespace {

KOKKOS_INLINE_FUNCTION
int kokkos_omp_in_parallel();

int kokkos_omp_in_critical_region = ( Kokkos::HostSpace::register_in_parallel( kokkos_omp_in_parallel ) , 0 );

KOKKOS_INLINE_FUNCTION
int kokkos_omp_in_parallel()
{
#ifndef __CUDA_ARCH__
  return omp_in_parallel() && ! kokkos_omp_in_critical_region ;
#else
  return 0;
#endif
}

unsigned s_threads_per_core = 0 ;
unsigned s_threads_per_numa = 0 ;
bool s_using_hwloc = false;

KOKKOS_INLINE_FUNCTION
unsigned fan_size( const unsigned rank , const unsigned size )
{
  const unsigned rank_rev = size - ( rank + 1 );
  unsigned count = 0 ;
  for ( unsigned n = 1 ; ( rank_rev + n < size ) && ! ( rank_rev & n ) ; n <<= 1 ) { ++count ; }
  return count ;
}

} // namespace
} // namespace Impl
} // namespace Kokkos


namespace Kokkos {
namespace Impl {

OpenMPexec * OpenMPexec::m_thread[ OpenMPexec::MAX_THREAD_COUNT ] = { 0 }; // Indexed by omp_get_thread_num()
OpenMPexec * OpenMPexec::m_pool[   OpenMPexec::MAX_THREAD_COUNT ] = { 0 }; // Indexed by OpenMPexec::m_pool_rank

OpenMPexec::OpenMPexec( const unsigned pool_rank )
  : m_team_base(0)
  , m_alloc_reduce(0)
  , m_alloc_shared(0)
  , m_team_shared(0)
  , m_alloc_shared_size(0)
  , m_pool_rank( pool_rank )
  , m_team_shared_end(0)
  , m_team_shared_iter(0)
  , m_team_rank(0)
  , m_team_size(0)
  , m_team_fan_size(0)
  , m_league_rank(0)
  , m_league_end(0)
  , m_league_size(0)
  , m_barrier_state( OpenMPexec::Active )
  , m_scan_state( OpenMPexec::Active )
{}

OpenMPexec::~OpenMPexec() {}

void OpenMPexec::team_work_init( size_t league_size , size_t team_size )
{
  m_team_base        = 0 ;
  m_team_shared      = 0 ;
  m_team_shared_end  = 0 ;
  m_team_size        = 0 ;
  m_team_rank        = 0 ;
  m_team_fan_size    = 0 ;
  m_league_size      = 0 ;
  m_league_rank      = 0 ;
  m_league_end       = 0 ;

  if ( league_size ) {

    if ( s_threads_per_numa < team_size ) { team_size = s_threads_per_numa ; }

    // Execution is using device-team interface:

    const unsigned pool_size     = omp_get_num_threads();
    const unsigned team_alloc    = s_threads_per_core * ( ( team_size + s_threads_per_core - 1 ) / s_threads_per_core );
    const unsigned pool_rank_rev = pool_size - ( m_pool_rank + 1 );
    const unsigned team_rank_rev = pool_rank_rev % team_alloc ;

    // May be using fewer threads per team than a multiple of threads per core,
    // some threads will idle.

    if ( team_rank_rev < team_size ) {
      const size_t pool_league_size     = pool_size     / team_alloc ;
      const size_t pool_league_rank_rev = pool_rank_rev / team_alloc ;
      const size_t pool_league_rank     = pool_league_size - ( pool_league_rank_rev + 1 );

      m_team_base        = m_pool + team_alloc * pool_league_rank_rev ;
      m_team_shared      = (*m_team_base)->m_alloc_shared ;
      m_team_shared_end  = (*m_team_base)->m_alloc_shared_size ;
      m_team_size        = team_size ;
      m_team_rank        = team_size - ( team_rank_rev + 1 );
      m_team_fan_size    = fan_size( m_team_rank , team_size );
      m_league_size      = league_size ;
      m_league_rank      = ( league_size *  pool_league_rank    ) / pool_league_size ;
      m_league_end       = ( league_size * (pool_league_rank+1) ) / pool_league_size ;
    }
  }
}


void OpenMPexec::verify_is_process( const char * const label )
{
  if ( omp_in_parallel() ) {
    std::string msg( label );
    msg.append( " ERROR: in parallel" );
    Kokkos::Impl::throw_runtime_exception( msg );
  }
}

void OpenMPexec::verify_initialized( const char * const label )
{
  if ( 0 == m_thread[0] ) {
    std::string msg( label );
    msg.append( " ERROR: not initialized" );
    Kokkos::Impl::throw_runtime_exception( msg );
  }
}

void OpenMPexec::resize_reduce_scratch( size_t size )
{
  static size_t s_size = 0 ;

  verify_initialized( "OpenMP::resize_reduce_scratch" );
  verify_is_process( "OpenMP::resize_reduce_scratch" );

  if ( size ) { size += REDUCE_TEAM_BASE ; }

  const size_t rem = size % Kokkos::Impl::MEMORY_ALIGNMENT ;

  if ( rem ) size += Kokkos::Impl::MEMORY_ALIGNMENT - rem ;

  if ( ( 0 == size && 0 != s_size ) || s_size < size ) {

#pragma omp parallel
    {
      OpenMPexec & th = * m_thread[ omp_get_thread_num() ];

#pragma omp critical
      {
        kokkos_omp_in_critical_region = 1 ;

        if ( th.m_alloc_reduce ) {
          HostSpace::decrement( th.m_alloc_reduce );
          th.m_alloc_reduce = 0 ;
        }

        if ( size ) {
          th.m_alloc_reduce = HostSpace::allocate( "openmp_reduce_scratch" , typeid(unsigned char) , 1 , size );
        }
        kokkos_omp_in_critical_region = 0 ;
      }
/* END #pragma omp critical */
    }
/* END #pragma omp parallel */
  }

  s_size = size ;
}

void OpenMPexec::resize_shared_scratch( size_t size )
{
  static size_t s_size = 0 ;

  verify_initialized( "OpenMP::resize_shared_scratch" );
  verify_is_process( "OpenMP::resize_shared_scratch" );

  const size_t rem = size % Kokkos::Impl::MEMORY_ALIGNMENT ;

  if ( rem ) size += Kokkos::Impl::MEMORY_ALIGNMENT - rem ;

  if ( ( 0 == size && 0 != s_size ) || s_size < size ) {

#pragma omp parallel
    {
      OpenMPexec & th = * m_thread[ omp_get_thread_num() ];

      const unsigned rank_rev = omp_get_num_threads() - ( th.m_pool_rank + 1 );

      if ( ! ( rank_rev % s_threads_per_core ) ) {
#pragma omp critical
        {
          kokkos_omp_in_critical_region = 1 ;

          if ( th.m_alloc_shared ) {
            HostSpace::decrement( th.m_alloc_shared );
            th.m_alloc_shared = 0 ;
          }

          if ( size ) {
            th.m_alloc_shared = HostSpace::allocate( "openmp_shared_scratch" , typeid(unsigned char) , 1 , size );
            th.m_alloc_shared_size = size ;
          }

          kokkos_omp_in_critical_region = 0 ;
        }
/* END #pragma omp critical */
      }
    }
/* END #pragma omp parallel */
  }

  s_size = size ;
}


KOKKOS_FUNCTION
void * OpenMPexec::get_shmem( const int size )
{
#ifndef __CUDA_ARCH__
  // m_shared_iter is in bytes, convert to integer offsets
  const int offset = m_team_shared_iter >> power_of_two<sizeof(int)>::value ;

  m_team_shared_iter += size ;

  if ( m_team_shared_end < m_team_shared_iter ) {
    Kokkos::Impl::throw_runtime_exception( std::string("OpenMPexec::get_shmem FAILED : exceeded shared memory size" ) );
  }

  return ((int*)m_team_shared) + offset ;
#else
  return NULL;
#endif
}

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

KOKKOS_FUNCTION
unsigned OpenMP::league_max()
{
#ifndef __CUDA_ARCH__
  Impl::OpenMPexec::verify_initialized("Kokkos::OpenMP::league_max" );
  Impl::OpenMPexec::verify_is_process("Kokkos::OpenMP::league_max" );

  return unsigned( std::numeric_limits<int>::max() );
#else
  return 0;
#endif
}

KOKKOS_FUNCTION
unsigned OpenMP::team_max()
{
#ifndef __CUDA_ARCH__
  Impl::OpenMPexec::verify_initialized("Kokkos::OpenMP::team_max" );
  Impl::OpenMPexec::verify_is_process("Kokkos::OpenMP::team_max" );

  return Impl::s_threads_per_numa ;
#else
  return 0;
#endif
}

//----------------------------------------------------------------------------

int OpenMP::is_initialized()
{ return 0 != Impl::OpenMPexec::m_thread[0]; }

void OpenMP::initialize( unsigned thread_count ,
                         unsigned use_numa_count ,
                         unsigned use_cores_per_numa )
{
  if(thread_count==0) thread_count = omp_get_max_threads();
  const bool is_initialized = 0 != Impl::OpenMPexec::m_thread[0] ;

  bool thread_spawn_failed = false ;

  if ( ! is_initialized ) {

    Impl::s_using_hwloc = hwloc::available() && (use_cores_per_numa > 0);

    std::pair<unsigned,unsigned> threads_coord[ Impl::OpenMPexec::MAX_THREAD_COUNT ];

    if(Impl::s_using_hwloc)
      hwloc::thread_mapping( "Kokkos::OpenMP::initialize" ,
                           false /* do not allow asynchronous */ ,
                           thread_count ,
                           use_numa_count ,
                           use_cores_per_numa ,
                           threads_coord );

    // Spawn threads:

    omp_set_num_threads( thread_count );

    // Verify OMP interaction:
    if ( int(thread_count) != omp_get_max_threads() ) {
      thread_spawn_failed = true ;
    }

    // Verify spawning and bind threads:
#pragma omp parallel
    {
#pragma omp critical
      {
        if ( int(thread_count) != omp_get_num_threads() ) {
          thread_spawn_failed = true ;
        }

        // Call to 'bind_this_thread' is not thread safe so place this whole block in a critical region.
        // Call to 'new' may not be thread safe as well.

        // Reverse the rank for threads so that the scan operation reduces to the highest rank thread.

        const unsigned omp_rank    = omp_get_thread_num();
        const unsigned thread_r    = Impl::s_using_hwloc ? Kokkos::hwloc::bind_this_thread( thread_count , threads_coord ) : omp_rank ;
        const unsigned thread_rank = thread_count - ( thread_r + 1 );

        Impl::OpenMPexec::m_thread[ omp_rank ] = new Impl::OpenMPexec( thread_rank );

        Impl::OpenMPexec::m_pool[ thread_r ] = Impl::OpenMPexec::m_thread[ omp_rank ] ;
      }
/* END #pragma omp critical */
    }
/* END #pragma omp parallel */

    if ( ! thread_spawn_failed ) {
      Impl::s_threads_per_numa = Impl::s_using_hwloc ? thread_count / use_numa_count : thread_count;
      Impl::s_threads_per_core = Impl::s_using_hwloc ? thread_count / ( use_numa_count * use_cores_per_numa ) : 1;

      Impl::OpenMPexec::resize_reduce_scratch( 4096 - Impl::OpenMPexec::REDUCE_TEAM_BASE );
      Impl::OpenMPexec::resize_shared_scratch( 4096 );
    }
  }

  if ( is_initialized || thread_spawn_failed ) {
    std::string msg("Kokkos::OpenMP::initialize ERROR");

    if ( is_initialized ) { msg.append(" : already initialized"); }
    if ( thread_spawn_failed ) { msg.append(" : failed spawning threads"); }

    Kokkos::Impl::throw_runtime_exception(msg);
  }
}

//----------------------------------------------------------------------------

void OpenMP::finalize()
{
  Impl::OpenMPexec::verify_initialized( "OpenMP::finalize" );
  Impl::OpenMPexec::verify_is_process( "OpenMP::finalize" );

  Impl::OpenMPexec::resize_reduce_scratch(0);
  Impl::OpenMPexec::resize_shared_scratch(0);

  for ( int i = 0 ; i < Impl::OpenMPexec::MAX_THREAD_COUNT ; ++i ) {
    if ( Impl::OpenMPexec::m_thread[i] ) { delete Impl::OpenMPexec::m_thread[i] ; }
    Impl::OpenMPexec::m_thread[i] = 0 ;
  }
  for ( int i = 0 ; i < Impl::OpenMPexec::MAX_THREAD_COUNT ; ++i ) {
    Impl::OpenMPexec::m_pool[i] = 0 ;
  }

  omp_set_num_threads(0);

  if(Impl::s_using_hwloc)
    hwloc::unbind_this_thread();
}

} // namespace Kokkos

