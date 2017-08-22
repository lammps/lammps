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
#if defined( KOKKOS_ENABLE_SERIAL )

#include <cstdlib>
#include <sstream>
#include <Kokkos_Serial.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Error.hpp>

#include <impl/Kokkos_SharedAlloc.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {
namespace {

HostThreadTeamData g_serial_thread_team_data ;

}

// Resize thread team data scratch memory
void serial_resize_thread_team_data( size_t pool_reduce_bytes
                                   , size_t team_reduce_bytes
                                   , size_t team_shared_bytes
                                   , size_t thread_local_bytes )
{
  if ( pool_reduce_bytes < 512 ) pool_reduce_bytes = 512 ;
  if ( team_reduce_bytes < 512 ) team_reduce_bytes = 512 ;

  const size_t old_pool_reduce  = g_serial_thread_team_data.pool_reduce_bytes();
  const size_t old_team_reduce  = g_serial_thread_team_data.team_reduce_bytes();
  const size_t old_team_shared  = g_serial_thread_team_data.team_shared_bytes();
  const size_t old_thread_local = g_serial_thread_team_data.thread_local_bytes();
  const size_t old_alloc_bytes  = g_serial_thread_team_data.scratch_bytes();

  // Allocate if any of the old allocation is tool small:

  const bool allocate = ( old_pool_reduce  < pool_reduce_bytes ) ||
                        ( old_team_reduce  < team_reduce_bytes ) ||
                        ( old_team_shared  < team_shared_bytes ) ||
                        ( old_thread_local < thread_local_bytes );

  if ( allocate ) {

    Kokkos::HostSpace space ;

    if ( old_alloc_bytes ) {
      g_serial_thread_team_data.disband_team();
      g_serial_thread_team_data.disband_pool();

      space.deallocate( g_serial_thread_team_data.scratch_buffer()
                      , g_serial_thread_team_data.scratch_bytes() );
    }

    if ( pool_reduce_bytes < old_pool_reduce ) { pool_reduce_bytes = old_pool_reduce ; }
    if ( team_reduce_bytes < old_team_reduce ) { team_reduce_bytes = old_team_reduce ; }
    if ( team_shared_bytes < old_team_shared ) { team_shared_bytes = old_team_shared ; }
    if ( thread_local_bytes < old_thread_local ) { thread_local_bytes = old_thread_local ; }

    const size_t alloc_bytes =
      HostThreadTeamData::scratch_size( pool_reduce_bytes
                                      , team_reduce_bytes
                                      , team_shared_bytes
                                      , thread_local_bytes );

    void * const ptr = space.allocate( alloc_bytes );

    g_serial_thread_team_data.
      scratch_assign( ((char *)ptr)
                    , alloc_bytes
                    , pool_reduce_bytes
                    , team_reduce_bytes
                    , team_shared_bytes
                    , thread_local_bytes );

    HostThreadTeamData * pool[1] = { & g_serial_thread_team_data };

    g_serial_thread_team_data.organize_pool( pool , 1 );
    g_serial_thread_team_data.organize_team(1);
  }
}

HostThreadTeamData * serial_get_thread_team_data()
{
  return & g_serial_thread_team_data ;
}

} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/

namespace Kokkos {

int Serial::is_initialized()
{
  return 1 ;
}

void Serial::initialize( unsigned threads_count
                       , unsigned use_numa_count
                       , unsigned use_cores_per_numa
                       , bool allow_asynchronous_threadpool )
{
  (void) threads_count;
  (void) use_numa_count;
  (void) use_cores_per_numa;
  (void) allow_asynchronous_threadpool;

  Impl::SharedAllocationRecord< void, void >::tracking_enable();

  // Init the array of locks used for arbitrarily sized atomics
  Impl::init_lock_array_host_space();
  #if defined(KOKKOS_ENABLE_PROFILING)
    Kokkos::Profiling::initialize();
  #endif
}

void Serial::finalize()
{
  if ( Impl::g_serial_thread_team_data.scratch_buffer() ) {
    Impl::g_serial_thread_team_data.disband_team();
    Impl::g_serial_thread_team_data.disband_pool();

    Kokkos::HostSpace space ;

    space.deallocate( Impl::g_serial_thread_team_data.scratch_buffer()
                    , Impl::g_serial_thread_team_data.scratch_bytes() );

    Impl::g_serial_thread_team_data.scratch_assign( (void*) 0, 0, 0, 0, 0, 0 );
  }

  #if defined(KOKKOS_ENABLE_PROFILING)
    Kokkos::Profiling::finalize();
  #endif
}

const char* Serial::name() { return "Serial"; }

} // namespace Kokkos

#else
void KOKKOS_CORE_SRC_IMPL_SERIAL_PREVENT_LINK_ERROR() {}
#endif // defined( KOKKOS_ENABLE_SERIAL )

