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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <Kokkos_Macros.hpp>
#if defined( KOKKOS_ENABLE_OPENMP )

#include <cstdio>
#include <limits>
#include <iostream>
#include <vector>
#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Error.hpp>
#include <iostream>
#include <impl/Kokkos_CPUDiscovery.hpp>
#include <impl/Kokkos_Profiling_Interface.hpp>


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

bool s_using_hwloc = false;

} // namespace
} // namespace Impl
} // namespace Kokkos


namespace Kokkos {
namespace Impl {

int OpenMPExec::m_map_rank[ OpenMPExec::MAX_THREAD_COUNT ] = { 0 };

int OpenMPExec::m_pool_topo[ 4 ] = { 0 };

HostThreadTeamData * OpenMPExec::m_pool[ OpenMPExec::MAX_THREAD_COUNT ] = { 0 };

void OpenMPExec::verify_is_process( const char * const label )
{
  if ( omp_in_parallel() ) {
    std::string msg( label );
    msg.append( " ERROR: in parallel" );
    Kokkos::Impl::throw_runtime_exception( msg );
  }
}

void OpenMPExec::verify_initialized( const char * const label )
{
  if ( 0 == m_pool[0] ) {
    std::string msg( label );
    msg.append( " ERROR: not initialized" );
    Kokkos::Impl::throw_runtime_exception( msg );
  }

  if ( omp_get_max_threads() != Kokkos::OpenMP::thread_pool_size(0) ) {
    std::string msg( label );
    msg.append( " ERROR: Initialized but threads modified inappropriately" );
    Kokkos::Impl::throw_runtime_exception( msg );
  }

}

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

void OpenMPExec::clear_thread_data()
{
  const size_t member_bytes =
    sizeof(int64_t) *
    HostThreadTeamData::align_to_int64( sizeof(HostThreadTeamData) );

  const int old_alloc_bytes =
    m_pool[0] ? ( member_bytes + m_pool[0]->scratch_bytes() ) : 0 ;

  Kokkos::HostSpace space ;

#pragma omp parallel
  {
    const int rank = m_map_rank[ omp_get_thread_num() ];

    if ( 0 != m_pool[rank] ) {

      m_pool[rank]->disband_pool();

      space.deallocate( m_pool[rank] , old_alloc_bytes );

      m_pool[rank] = 0 ;
    }
  }
/* END #pragma omp parallel */
}

void OpenMPExec::resize_thread_data( size_t pool_reduce_bytes
                                   , size_t team_reduce_bytes
                                   , size_t team_shared_bytes
                                   , size_t thread_local_bytes )
{
  const size_t member_bytes =
    sizeof(int64_t) *
    HostThreadTeamData::align_to_int64( sizeof(HostThreadTeamData) );

  HostThreadTeamData * root = m_pool[0] ;

  const size_t old_pool_reduce  = root ? root->pool_reduce_bytes() : 0 ;
  const size_t old_team_reduce  = root ? root->team_reduce_bytes() : 0 ;
  const size_t old_team_shared  = root ? root->team_shared_bytes() : 0 ;
  const size_t old_thread_local = root ? root->thread_local_bytes() : 0 ;
  const size_t old_alloc_bytes  = root ? ( member_bytes + root->scratch_bytes() ) : 0 ;

  // Allocate if any of the old allocation is tool small:

  const bool allocate = ( old_pool_reduce  < pool_reduce_bytes ) ||
                        ( old_team_reduce  < team_reduce_bytes ) ||
                        ( old_team_shared  < team_shared_bytes ) ||
                        ( old_thread_local < thread_local_bytes );

  if ( allocate ) {

    if ( pool_reduce_bytes < old_pool_reduce ) { pool_reduce_bytes = old_pool_reduce ; }
    if ( team_reduce_bytes < old_team_reduce ) { team_reduce_bytes = old_team_reduce ; }
    if ( team_shared_bytes < old_team_shared ) { team_shared_bytes = old_team_shared ; }
    if ( thread_local_bytes < old_thread_local ) { thread_local_bytes = old_thread_local ; }

    const size_t alloc_bytes =
      member_bytes +
      HostThreadTeamData::scratch_size( pool_reduce_bytes
                                      , team_reduce_bytes
                                      , team_shared_bytes
                                      , thread_local_bytes );

    const int pool_size = omp_get_max_threads();

    Kokkos::HostSpace space ;

#pragma omp parallel
    {
      const int rank = m_map_rank[ omp_get_thread_num() ];

      if ( 0 != m_pool[rank] ) {

        m_pool[rank]->disband_pool();

        space.deallocate( m_pool[rank] , old_alloc_bytes );
      }

      void * const ptr = space.allocate( alloc_bytes );

      m_pool[ rank ] = new( ptr ) HostThreadTeamData();

      m_pool[ rank ]->
        scratch_assign( ((char *)ptr) + member_bytes
                      , alloc_bytes
                      , pool_reduce_bytes
                      , team_reduce_bytes
                      , team_shared_bytes
                      , thread_local_bytes );
    }
/* END #pragma omp parallel */

    HostThreadTeamData::organize_pool( m_pool , pool_size );
  }
}

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

//----------------------------------------------------------------------------

int OpenMP::is_initialized()
{ return 0 != Impl::OpenMPExec::m_pool[0]; }

void OpenMP::initialize( unsigned thread_count ,
                         unsigned use_numa_count ,
                         unsigned use_cores_per_numa )
{
  // Before any other call to OMP query the maximum number of threads
  // and save the value for re-initialization unit testing.

  // Using omp_get_max_threads(); is problematic in conjunction with
  // Hwloc on Intel (essentially an initial call to the OpenMP runtime
  // without a parallel region before will set a process mask for a single core
  // The runtime will than bind threads for a parallel region to other cores on the
  // entering the first parallel region and make the process mask the aggregate of
  // the thread masks. The intend seems to be to make serial code run fast, if you
  // compile with OpenMP enabled but don't actually use parallel regions or so
  // static int omp_max_threads = omp_get_max_threads();
  int nthreads = 0;
  #pragma omp parallel
  {
    #pragma omp atomic
    nthreads++;
  }

  static int omp_max_threads = nthreads;

  const bool is_initialized = 0 != Impl::OpenMPExec::m_pool[0] ;

  bool thread_spawn_failed = false ;

  if ( ! is_initialized ) {

    // Use hwloc thread pinning if concerned with locality.
    // If spreading threads across multiple NUMA regions.
    // If hyperthreading is enabled.
    Impl::s_using_hwloc = hwloc::available() && (
                            ( 1 < Kokkos::hwloc::get_available_numa_count() ) ||
                            ( 1 < Kokkos::hwloc::get_available_threads_per_core() ) );

    std::pair<unsigned,unsigned> threads_coord[ Impl::OpenMPExec::MAX_THREAD_COUNT ];

    // If hwloc available then use it's maximum value.

    if ( thread_count == 0 ) {
      thread_count = Impl::s_using_hwloc
      ? Kokkos::hwloc::get_available_numa_count() *
        Kokkos::hwloc::get_available_cores_per_numa() *
        Kokkos::hwloc::get_available_threads_per_core()
      : omp_max_threads ;
    }

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

        const unsigned omp_rank    = omp_get_thread_num();
        const unsigned thread_r    = Impl::s_using_hwloc && Kokkos::hwloc::can_bind_threads()
                                   ? Kokkos::hwloc::bind_this_thread( thread_count , threads_coord )
                                   : omp_rank ;

        Impl::OpenMPExec::m_map_rank[ omp_rank ] = thread_r ;
      }
/* END #pragma omp critical */
    }
/* END #pragma omp parallel */

    if ( ! thread_spawn_failed ) {
      Impl::OpenMPExec::m_pool_topo[0] = thread_count ;
      Impl::OpenMPExec::m_pool_topo[1] = Impl::s_using_hwloc ? thread_count / use_numa_count : thread_count;
      Impl::OpenMPExec::m_pool_topo[2] = Impl::s_using_hwloc ? thread_count / ( use_numa_count * use_cores_per_numa ) : 1;

      // New, unified host thread team data:
      {
        size_t pool_reduce_bytes  =   32 * thread_count ;
        size_t team_reduce_bytes  =   32 * thread_count ;
        size_t team_shared_bytes  = 1024 * thread_count ;
        size_t thread_local_bytes = 1024 ;

        Impl::OpenMPExec::resize_thread_data( pool_reduce_bytes
                                            , team_reduce_bytes
                                            , team_shared_bytes
                                            , thread_local_bytes
                                            );
      }
    }
  }

  if ( is_initialized || thread_spawn_failed ) {
    std::string msg("Kokkos::OpenMP::initialize ERROR");

    if ( is_initialized ) { msg.append(" : already initialized"); }
    if ( thread_spawn_failed ) { msg.append(" : failed spawning threads"); }

    Kokkos::Impl::throw_runtime_exception(msg);
  }

  // Check for over-subscription
  //if( Impl::mpi_ranks_per_node() * long(thread_count) > Impl::processors_per_node() ) {
  //  std::cout << "Kokkos::OpenMP::initialize WARNING: You are likely oversubscribing your CPU cores." << std::endl;
  //  std::cout << "                                    Detected: " << Impl::processors_per_node() << " cores per node." << std::endl;
  //  std::cout << "                                    Detected: " << Impl::mpi_ranks_per_node() << " MPI_ranks per node." << std::endl;
  //  std::cout << "                                    Requested: " << thread_count << " threads per process." << std::endl;
  //}
  // Init the array for used for arbitrarily sized atomics
  Impl::init_lock_array_host_space();

  #if defined(KOKKOS_ENABLE_PROFILING)
    Kokkos::Profiling::initialize();
  #endif
}

//----------------------------------------------------------------------------

void OpenMP::finalize()
{
  Impl::OpenMPExec::verify_initialized( "OpenMP::finalize" );
  Impl::OpenMPExec::verify_is_process( "OpenMP::finalize" );

  // New, unified host thread team data:
  Impl::OpenMPExec::clear_thread_data();

  Impl::OpenMPExec::m_pool_topo[0] = 0 ;
  Impl::OpenMPExec::m_pool_topo[1] = 0 ;
  Impl::OpenMPExec::m_pool_topo[2] = 0 ;

  omp_set_num_threads(1);

  if ( Impl::s_using_hwloc && Kokkos::hwloc::can_bind_threads() ) {
    hwloc::unbind_this_thread();
  }

  #if defined(KOKKOS_ENABLE_PROFILING)
    Kokkos::Profiling::finalize();
  #endif
}

//----------------------------------------------------------------------------

void OpenMP::print_configuration( std::ostream & s , const bool detail )
{
  Impl::OpenMPExec::verify_is_process( "OpenMP::print_configuration" );

  s << "Kokkos::OpenMP" ;

#if defined( KOKKOS_ENABLE_OPENMP )
  s << " KOKKOS_ENABLE_OPENMP" ;
#endif
#if defined( KOKKOS_ENABLE_HWLOC )

  const unsigned numa_count_       = Kokkos::hwloc::get_available_numa_count();
  const unsigned cores_per_numa   = Kokkos::hwloc::get_available_cores_per_numa();
  const unsigned threads_per_core = Kokkos::hwloc::get_available_threads_per_core();

  s << " hwloc[" << numa_count_ << "x" << cores_per_numa << "x" << threads_per_core << "]"
    << " hwloc_binding_" << ( Impl::s_using_hwloc ? "enabled" : "disabled" )
    ;
#endif

  const bool is_initialized = 0 != Impl::OpenMPExec::m_pool[0] ;

  if ( is_initialized ) {
    const int numa_count      = Kokkos::Impl::OpenMPExec::m_pool_topo[0] / Kokkos::Impl::OpenMPExec::m_pool_topo[1] ;
    const int core_per_numa   = Kokkos::Impl::OpenMPExec::m_pool_topo[1] / Kokkos::Impl::OpenMPExec::m_pool_topo[2] ;
    const int thread_per_core = Kokkos::Impl::OpenMPExec::m_pool_topo[2] ;

    s << " thread_pool_topology[ " << numa_count
      << " x " << core_per_numa
      << " x " << thread_per_core
      << " ]"
      << std::endl ;

    if ( detail ) {
      std::vector< std::pair<unsigned,unsigned> > coord( Kokkos::Impl::OpenMPExec::m_pool_topo[0] );

#pragma omp parallel
      {
#pragma omp critical
        {
          coord[ omp_get_thread_num() ] = hwloc::get_this_thread_coordinate();
        }
/* END #pragma omp critical */
      }
/* END #pragma omp parallel */

      for ( unsigned i = 0 ; i < coord.size() ; ++i ) {
        s << "  thread omp_rank[" << i << "]"
          << " kokkos_rank[" << Impl::OpenMPExec::m_map_rank[ i ] << "]"
          << " hwloc_coord[" << coord[i].first << "." << coord[i].second << "]"
          << std::endl ;
      }
    }
  }
  else {
    s << " not initialized" << std::endl ;
  }
}

int OpenMP::concurrency() {
  return thread_pool_size(0);
}

const char* OpenMP::name() { return "OpenMP"; }

} // namespace Kokkos

#else
void KOKKOS_CORE_SRC_OPENMP_EXEC_PREVENT_LINK_ERROR() {}
#endif //KOKKOS_ENABLE_OPENMP

