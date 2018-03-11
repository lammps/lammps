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
#include <cstdlib>

#include <limits>
#include <iostream>
#include <vector>

#include <Kokkos_Core.hpp>

#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_CPUDiscovery.hpp>
#include <impl/Kokkos_Profiling_Interface.hpp>


namespace Kokkos {
namespace Impl {

int g_openmp_hardware_max_threads = 1;

__thread int t_openmp_hardware_id = 0;
__thread Impl::OpenMPExec * t_openmp_instance = nullptr;

void OpenMPExec::validate_partition( const int nthreads
                                   , int & num_partitions
                                   , int & partition_size
                                  )
{
  if (nthreads == 1) {
    num_partitions = 1;
    partition_size = 1;
  }
  else if( num_partitions < 1 && partition_size < 1) {
    int idle = nthreads;
    for (int np = 2; np <= nthreads ; ++np) {
      for (int ps = 1; ps <= nthreads/np; ++ps) {
        if (nthreads - np*ps < idle) {
          idle = nthreads - np*ps;
          num_partitions = np;
          partition_size = ps;
        }
        if (idle == 0) {
          break;
        }
      }
    }
  }
  else if( num_partitions < 1 && partition_size > 0 ) {
    if ( partition_size <= nthreads ) {
      num_partitions = nthreads / partition_size;
    }
    else {
      num_partitions = 1;
      partition_size = nthreads;
    }
  }
  else if( num_partitions > 0 && partition_size < 1 ) {
    if ( num_partitions <= nthreads ) {
      partition_size = nthreads / num_partitions;
    }
    else {
      num_partitions = nthreads;
      partition_size = 1;
    }
  }
  else if ( num_partitions * partition_size > nthreads ) {
    int idle = nthreads;
    const int NP = num_partitions;
    const int PS = partition_size;
    for (int np = NP; np > 0; --np) {
      for (int ps = PS; ps > 0; --ps) {
        if (  (np*ps <= nthreads)
           && (nthreads - np*ps < idle) ) {
          idle = nthreads - np*ps;
          num_partitions = np;
          partition_size = ps;
        }
        if (idle == 0) {
          break;
        }
      }
    }
  }

}

void OpenMPExec::verify_is_master( const char * const label )
{
  if ( !t_openmp_instance )
  {
    std::string msg( label );
    msg.append( " ERROR: in parallel or not initialized" );
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

  OpenMP::memory_space space ;

  #pragma omp parallel num_threads( m_pool_size )
  {
    const int rank = omp_get_thread_num();

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

    OpenMP::memory_space space ;

    memory_fence();

    #pragma omp parallel num_threads(m_pool_size)
    {
      const int rank = omp_get_thread_num();

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
                      , thread_local_bytes
                      );

      memory_fence();
    }
/* END #pragma omp parallel */

    HostThreadTeamData::organize_pool( m_pool , m_pool_size );
  }
}

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

//----------------------------------------------------------------------------

int OpenMP::get_current_max_threads() noexcept
{
  // Using omp_get_max_threads(); is problematic in conjunction with
  // Hwloc on Intel (essentially an initial call to the OpenMP runtime
  // without a parallel region before will set a process mask for a single core
  // The runtime will than bind threads for a parallel region to other cores on the
  // entering the first parallel region and make the process mask the aggregate of
  // the thread masks. The intend seems to be to make serial code run fast, if you
  // compile with OpenMP enabled but don't actually use parallel regions or so
  // static int omp_max_threads = omp_get_max_threads();

  int count = 0;
  #pragma omp parallel
  {
    #pragma omp atomic
     ++count;
  }
  return count;
}


void OpenMP::initialize( int thread_count )
{
  if ( omp_in_parallel() ) {
    std::string msg("Kokkos::OpenMP::initialize ERROR : in parallel");
    Kokkos::Impl::throw_runtime_exception(msg);
  }

  if ( Impl::t_openmp_instance )
  {
    finalize();
  }

  {
    if ( Kokkos::show_warnings() && nullptr == std::getenv("OMP_PROC_BIND") ) {
      printf("Kokkos::OpenMP::initialize WARNING: OMP_PROC_BIND environment variable not set\n");
      printf("  In general, for best performance with OpenMP 4.0 or better set OMP_PROC_BIND=spread and OMP_PLACES=threads\n");
      printf("  For best performance with OpenMP 3.1 set OMP_PROC_BIND=true\n");
      printf("  For unit testing set OMP_PROC_BIND=false\n");
    }

    OpenMP::memory_space space ;

    // Before any other call to OMP query the maximum number of threads
    // and save the value for re-initialization unit testing.

    Impl::g_openmp_hardware_max_threads = get_current_max_threads();

    int process_num_threads = Impl::g_openmp_hardware_max_threads;

    if ( Kokkos::hwloc::available() ) {
      process_num_threads = Kokkos::hwloc::get_available_numa_count()
                          * Kokkos::hwloc::get_available_cores_per_numa()
                          * Kokkos::hwloc::get_available_threads_per_core();
    }

    // if thread_count  < 0, use g_openmp_hardware_max_threads;
    // if thread_count == 0, set g_openmp_hardware_max_threads to process_num_threads
    // if thread_count  > 0, set g_openmp_hardware_max_threads to thread_count
    if (thread_count < 0 ) {
      thread_count = Impl::g_openmp_hardware_max_threads;
    }
    else if( thread_count == 0 && Impl::g_openmp_hardware_max_threads != process_num_threads ) {
      Impl::g_openmp_hardware_max_threads = process_num_threads;
      omp_set_num_threads(Impl::g_openmp_hardware_max_threads);
    }
    else {
      if( Kokkos::show_warnings() && thread_count > process_num_threads ) {
        printf( "Kokkos::OpenMP::initialize WARNING: You are likely oversubscribing your CPU cores.\n");
        printf( "  process threads available : %3d,  requested thread : %3d\n", process_num_threads, thread_count );
      }
      Impl::g_openmp_hardware_max_threads = thread_count;
      omp_set_num_threads(Impl::g_openmp_hardware_max_threads);
    }

    // setup thread local
    #pragma omp parallel num_threads(Impl::g_openmp_hardware_max_threads)
    {
      Impl::t_openmp_instance = nullptr;
      Impl::t_openmp_hardware_id = omp_get_thread_num();
      Impl::SharedAllocationRecord< void, void >::tracking_enable();
    }

    void * const ptr = space.allocate( sizeof(Impl::OpenMPExec) );

    Impl::t_openmp_instance = new (ptr) Impl::OpenMPExec( Impl::g_openmp_hardware_max_threads );

    // New, unified host thread team data:
    {
      size_t pool_reduce_bytes  =   32 * thread_count ;
      size_t team_reduce_bytes  =   32 * thread_count ;
      size_t team_shared_bytes  = 1024 * thread_count ;
      size_t thread_local_bytes = 1024 ;

      Impl::t_openmp_instance->resize_thread_data( pool_reduce_bytes
                                                 , team_reduce_bytes
                                                 , team_shared_bytes
                                                 , thread_local_bytes
                                                 );
    }
  }


  // Check for over-subscription
  if( Kokkos::show_warnings() && (Impl::mpi_ranks_per_node() * long(thread_count) > Impl::processors_per_node()) ) {
    std::cout << "Kokkos::OpenMP::initialize WARNING: You are likely oversubscribing your CPU cores." << std::endl;
    std::cout << "                                    Detected: " << Impl::processors_per_node() << " cores per node." << std::endl;
    std::cout << "                                    Detected: " << Impl::mpi_ranks_per_node() << " MPI_ranks per node." << std::endl;
    std::cout << "                                    Requested: " << thread_count << " threads per process." << std::endl;
  }
  // Init the array for used for arbitrarily sized atomics
  Impl::init_lock_array_host_space();

  #if defined(KOKKOS_ENABLE_PROFILING)
    Kokkos::Profiling::initialize();
  #endif
}

//----------------------------------------------------------------------------

void OpenMP::finalize()
{
  if ( omp_in_parallel() )
  {
    std::string msg("Kokkos::OpenMP::finalize ERROR ");
    if( !Impl::t_openmp_instance ) msg.append(": not initialized");
    if( omp_in_parallel() ) msg.append(": in parallel");
    Kokkos::Impl::throw_runtime_exception(msg);
  }

  if ( Impl::t_openmp_instance ) {

    const int nthreads = Impl::t_openmp_instance->m_pool_size <= Impl::g_openmp_hardware_max_threads
                       ? Impl::g_openmp_hardware_max_threads
                       : Impl::t_openmp_instance->m_pool_size;

    using Exec = Impl::OpenMPExec;
    Exec * instance = Impl::t_openmp_instance;
    instance->~Exec();

    OpenMP::memory_space space;
    space.deallocate( instance, sizeof(Exec) );

    #pragma omp parallel num_threads(nthreads)
    {
      Impl::t_openmp_hardware_id = 0;
      Impl::t_openmp_instance    = nullptr;
      Impl::SharedAllocationRecord< void, void >::tracking_disable();
    }

    // allow main thread to track
    Impl::SharedAllocationRecord< void, void >::tracking_enable();

    Impl::g_openmp_hardware_max_threads = 1;
  }

  #if defined(KOKKOS_ENABLE_PROFILING)
    Kokkos::Profiling::finalize();
  #endif
}

//----------------------------------------------------------------------------

void OpenMP::print_configuration( std::ostream & s , const bool verbose )
{
  s << "Kokkos::OpenMP" ;

  const bool is_initialized =  Impl::t_openmp_instance != nullptr;

  if ( is_initialized ) {
    Impl::OpenMPExec::verify_is_master( "OpenMP::print_configuration" );

    const int numa_count      = 1;
    const int core_per_numa   = Impl::g_openmp_hardware_max_threads;
    const int thread_per_core = 1;

    s << " thread_pool_topology[ " << numa_count
      << " x " << core_per_numa
      << " x " << thread_per_core
      << " ]"
      << std::endl ;
  }
  else {
    s << " not initialized" << std::endl ;
  }
}

std::vector<OpenMP> OpenMP::partition(...)
{ return std::vector<OpenMP>(1); }

OpenMP OpenMP::create_instance(...) { return OpenMP(); }


#if !defined( KOKKOS_DISABLE_DEPRECATED )

int OpenMP::concurrency() {
  return Impl::g_openmp_hardware_max_threads;
}

void OpenMP::initialize( int thread_count , int, int )
{
  initialize(thread_count);
}

#endif

} // namespace Kokkos

#else
void KOKKOS_CORE_SRC_OPENMP_EXEC_PREVENT_LINK_ERROR() {}
#endif //KOKKOS_ENABLE_OPENMP

