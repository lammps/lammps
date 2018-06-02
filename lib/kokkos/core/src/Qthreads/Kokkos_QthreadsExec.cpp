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

#include <Kokkos_Macros.hpp>
#if defined( KOKKOS_ENABLE_QTHREADS )

#include <Kokkos_Core_fwd.hpp>

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <utility>

#include <Kokkos_Qthreads.hpp>
#include <Kokkos_Atomic.hpp>
#include <impl/Kokkos_Error.hpp>

// Defines to enable experimental Qthreads functionality.
//#define QTHREAD_LOCAL_PRIORITY
//#define CLONED_TASKS

//#include <qthread.h>

//----------------------------------------------------------------------------

namespace Kokkos {

namespace Impl {

namespace {

enum { MAXIMUM_QTHREADS_WORKERS = 1024 };

/** s_exec is indexed by the reverse rank of the workers
 *  for faster fan-in / fan-out lookups
 *  [ n - 1, n - 2, ..., 0 ]
 */
QthreadsExec * s_exec[ MAXIMUM_QTHREADS_WORKERS ];

int  s_number_shepherds            = 0;
int  s_number_workers_per_shepherd = 0;
int  s_number_workers              = 0;

inline
QthreadsExec ** worker_exec()
{
  return s_exec + s_number_workers - ( qthread_shep() * s_number_workers_per_shepherd + qthread_worker_local( NULL ) + 1 );
}

const int s_base_size = QthreadsExec::align_alloc( sizeof(QthreadsExec) );

int s_worker_reduce_end   = 0;  // End of worker reduction memory.
int s_worker_shared_end   = 0;  // Total of worker scratch memory.
int s_worker_shared_begin = 0;  // Beginning of worker shared memory.

QthreadsExecFunctionPointer volatile s_active_function     = 0;
const void                * volatile s_active_function_arg = 0;

} // namespace

} // namespace Impl

} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

int Qthreads::is_initialized()
{
  return Impl::s_number_workers != 0;
}

int Qthreads::concurrency()
{
  return Impl::s_number_workers_per_shepherd;
}

int Qthreads::in_parallel()
{
  return Impl::s_active_function != 0;
}

void Qthreads::initialize( int thread_count )
{
  // Environment variable: QTHREAD_NUM_SHEPHERDS
  // Environment variable: QTHREAD_NUM_WORKERS_PER_SHEP
  // Environment variable: QTHREAD_HWPAR

  {
    char buffer[256];
    snprintf( buffer, sizeof(buffer), "QTHREAD_HWPAR=%d", thread_count );
    putenv( buffer );
  }

  const bool ok_init = ( QTHREAD_SUCCESS == qthread_initialize() ) &&
                       ( thread_count    == qthread_num_shepherds() * qthread_num_workers_local( NO_SHEPHERD ) ) &&
                       ( thread_count    == qthread_num_workers() );

  bool ok_symmetry = true;

  if ( ok_init ) {
    Impl::s_number_shepherds            = qthread_num_shepherds();
    Impl::s_number_workers_per_shepherd = qthread_num_workers_local( NO_SHEPHERD );
    Impl::s_number_workers              = Impl::s_number_shepherds * Impl::s_number_workers_per_shepherd;

    for ( int i = 0; ok_symmetry && i < Impl::s_number_shepherds; ++i ) {
      ok_symmetry = ( Impl::s_number_workers_per_shepherd == qthread_num_workers_local( i ) );
    }
  }

  if ( ! ok_init || ! ok_symmetry ) {
    std::ostringstream msg;

    msg << "Kokkos::Qthreads::initialize(" << thread_count << ") FAILED";
    msg << " : qthread_num_shepherds = " << qthread_num_shepherds();
    msg << " : qthread_num_workers_per_shepherd = " << qthread_num_workers_local( NO_SHEPHERD );
    msg << " : qthread_num_workers = " << qthread_num_workers();

    if ( ! ok_symmetry ) {
      msg << " : qthread_num_workers_local = {";
      for ( int i = 0; i < Impl::s_number_shepherds; ++i ) {
        msg << " " << qthread_num_workers_local( i );
      }
      msg << " }";
    }

    Impl::s_number_workers              = 0;
    Impl::s_number_shepherds            = 0;
    Impl::s_number_workers_per_shepherd = 0;

    if ( ok_init ) { qthread_finalize(); }

    Kokkos::Impl::throw_runtime_exception( msg.str() );
  }

  Impl::QthreadsExec::resize_worker_scratch( 256, 256 );

  // Init the array for used for arbitrarily sized atomics.
  Impl::init_lock_array_host_space();

}

void Qthreads::finalize()
{
  Impl::QthreadsExec::clear_workers();

  if ( Impl::s_number_workers ) {
    qthread_finalize();
  }

  Impl::s_number_workers              = 0;
  Impl::s_number_shepherds            = 0;
  Impl::s_number_workers_per_shepherd = 0;
}

void Qthreads::print_configuration( std::ostream & s, const bool detail )
{
  s << "Kokkos::Qthreads {"
    << " num_shepherds(" << Impl::s_number_shepherds << ")"
    << " num_workers_per_shepherd(" << Impl::s_number_workers_per_shepherd << ")"
    << " }" << std::endl;
}

Qthreads & Qthreads::instance( int )
{
  static Qthreads q;
  return q;
}

void Qthreads::fence()
{
}

int Qthreads::shepherd_size() const { return Impl::s_number_shepherds; }
int Qthreads::shepherd_worker_size() const { return Impl::s_number_workers_per_shepherd; }

const char* Qthreads::name() { return "Qthreads"; }

} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

namespace Impl {

namespace {

aligned_t driver_exec_all( void * arg )
{
  QthreadsExec & exec = **worker_exec();

  (*s_active_function)( exec, s_active_function_arg );

/*
  fprintf( stdout
         , "QthreadsExec driver worker(%d:%d) shepherd(%d:%d) shepherd_worker(%d:%d) done\n"
         , exec.worker_rank()
         , exec.worker_size()
         , exec.shepherd_rank()
         , exec.shepherd_size()
         , exec.shepherd_worker_rank()
         , exec.shepherd_worker_size()
         );
  fflush(stdout);
*/

  return 0;
}

aligned_t driver_resize_worker_scratch( void * arg )
{
  static volatile int lock_begin = 0;
  static volatile int lock_end   = 0;

  QthreadsExec ** const exec = worker_exec();

  //----------------------------------------
  // Serialize allocation for thread safety.

  while ( ! atomic_compare_exchange_strong( & lock_begin, 0, 1 ) ); // Spin wait to claim lock.

  const bool ok = 0 == *exec;

  if ( ok ) { *exec = (QthreadsExec *) malloc( s_base_size + s_worker_shared_end ); }

  lock_begin = 0; // Release lock.

  if ( ok ) { new( *exec ) QthreadsExec(); }

  //----------------------------------------
  // Wait for all calls to complete to insure that each worker has executed.

  if ( s_number_workers == 1 + atomic_fetch_add( & lock_end, 1 ) ) { lock_end = 0; }

  while ( lock_end );

/*
  fprintf( stdout
         , "QthreadsExec resize worker(%d:%d) shepherd(%d:%d) shepherd_worker(%d:%d) done\n"
         , (**exec).worker_rank()
         , (**exec).worker_size()
         , (**exec).shepherd_rank()
         , (**exec).shepherd_size()
         , (**exec).shepherd_worker_rank()
         , (**exec).shepherd_worker_size()
         );
  fflush(stdout);
*/

  //----------------------------------------

  if ( ! ok ) {
    fprintf( stderr, "Kokkos::QthreadsExec resize failed\n" );
    fflush( stderr );
  }

  return 0;
}

void verify_is_process( const char * const label, bool not_active = false )
{
  const bool not_process = 0 != qthread_shep() || 0 != qthread_worker_local( NULL );
  const bool is_active   = not_active && ( s_active_function || s_active_function_arg );

  if ( not_process || is_active ) {
    std::string msg( label );
    msg.append( " : FAILED" );
    if ( not_process ) msg.append(" : not called by main process");
    if ( is_active )   msg.append(" : parallel execution in progress");
    Kokkos::Impl::throw_runtime_exception( msg );
  }
}

} // namespace

int QthreadsExec::worker_per_shepherd()
{
  return s_number_workers_per_shepherd;
}

QthreadsExec::QthreadsExec()
{
  const int shepherd_rank        = qthread_shep();
  const int shepherd_worker_rank = qthread_worker_local( NULL );
  const int worker_rank          = shepherd_rank * s_number_workers_per_shepherd + shepherd_worker_rank;

  m_worker_base          = s_exec;
  m_shepherd_base        = s_exec + s_number_workers_per_shepherd * ( ( s_number_shepherds - ( shepherd_rank + 1 ) ) );
  m_scratch_alloc        = ( (unsigned char *) this ) + s_base_size;
  m_reduce_end           = s_worker_reduce_end;
  m_shepherd_rank        = shepherd_rank;
  m_shepherd_size        = s_number_shepherds;
  m_shepherd_worker_rank = shepherd_worker_rank;
  m_shepherd_worker_size = s_number_workers_per_shepherd;
  m_worker_rank          = worker_rank;
  m_worker_size          = s_number_workers;
  m_worker_state         = QthreadsExec::Active;
}

void QthreadsExec::clear_workers()
{
  for ( int iwork = 0; iwork < s_number_workers; ++iwork ) {
    QthreadsExec * const exec = s_exec[iwork];
    s_exec[iwork] = 0;
    free( exec );
  }
}

void QthreadsExec::shared_reset( Qthreads::scratch_memory_space & space )
{
  new( & space )
    Qthreads::scratch_memory_space(
      ((unsigned char *) (**m_shepherd_base).m_scratch_alloc ) + s_worker_shared_begin,
      s_worker_shared_end - s_worker_shared_begin
    );
}

void QthreadsExec::resize_worker_scratch( const int reduce_size, const int shared_size )
{
  const int exec_all_reduce_alloc = align_alloc( reduce_size );
  const int shepherd_scan_alloc   = align_alloc( 8 );
  const int shepherd_shared_end   = exec_all_reduce_alloc + shepherd_scan_alloc + align_alloc( shared_size );

  if ( s_worker_reduce_end < exec_all_reduce_alloc ||
       s_worker_shared_end < shepherd_shared_end ) {

/*
  fprintf( stdout, "QthreadsExec::resize\n");
  fflush(stdout);
*/

    // Clear current worker memory before allocating new worker memory.
    clear_workers();

    // Increase the buffers to an aligned allocation.
    s_worker_reduce_end   = exec_all_reduce_alloc;
    s_worker_shared_begin = exec_all_reduce_alloc + shepherd_scan_alloc;
    s_worker_shared_end   = shepherd_shared_end;

    // Need to query which shepherd this main 'process' is running.

    const int main_shep = qthread_shep();

    // Have each worker resize its memory for proper first-touch.
#if 0
    for ( int jshep = 0; jshep < s_number_shepherds; ++jshep ) {
      for ( int i = jshep != main_shep ? 0 : 1; i < s_number_workers_per_shepherd; ++i ) {
        qthread_fork_to( driver_resize_worker_scratch, NULL, NULL, jshep );
      }
    }
#else
    // If this function is used before the 'qthreads.task_policy' unit test,
    // the 'qthreads.task_policy' unit test fails with a seg-fault within libqthread.so.
    for ( int jshep = 0; jshep < s_number_shepherds; ++jshep ) {
      const int num_clone = jshep != main_shep ? s_number_workers_per_shepherd : s_number_workers_per_shepherd - 1;

      if ( num_clone ) {
        const int ret = qthread_fork_clones_to_local_priority
          ( driver_resize_worker_scratch   // Function
          , NULL                           // Function data block
          , NULL                           // Pointer to return value feb
          , jshep                          // Shepherd number
          , num_clone - 1                  // Number of instances - 1
          );

        assert( ret == QTHREAD_SUCCESS );
      }
    }
#endif

    driver_resize_worker_scratch( NULL );

    // Verify all workers allocated.

    bool ok = true;
    for ( int iwork = 0; ok && iwork < s_number_workers; ++iwork ) { ok = 0 != s_exec[iwork]; }

    if ( ! ok ) {
      std::ostringstream msg;
      msg << "Kokkos::Impl::QthreadsExec::resize : FAILED for workers {";
      for ( int iwork = 0; iwork < s_number_workers; ++iwork ) {
         if ( 0 == s_exec[iwork] ) { msg << " " << ( s_number_workers - ( iwork + 1 ) ); }
      }
      msg << " }";
      Kokkos::Impl::throw_runtime_exception( msg.str() );
    }
  }
}

void QthreadsExec::exec_all( Qthreads &, QthreadsExecFunctionPointer func, const void * arg )
{
  verify_is_process("QthreadsExec::exec_all(...)",true);

/*
  fprintf( stdout, "QthreadsExec::exec_all\n");
  fflush(stdout);
*/

  s_active_function     = func;
  s_active_function_arg = arg;

  // Need to query which shepherd this main 'process' is running.

  const int main_shep = qthread_shep();

#if 0
  for ( int jshep = 0, iwork = 0; jshep < s_number_shepherds; ++jshep ) {
    for ( int i = jshep != main_shep ? 0 : 1; i < s_number_workers_per_shepherd; ++i, ++iwork ) {
      qthread_fork_to( driver_exec_all, NULL, NULL, jshep );
    }
  }
#else
  // If this function is used before the 'qthreads.task_policy' unit test,
  // the 'qthreads.task_policy' unit test fails with a seg-fault within libqthread.so.
  for ( int jshep = 0; jshep < s_number_shepherds; ++jshep ) {
    const int num_clone = jshep != main_shep ? s_number_workers_per_shepherd : s_number_workers_per_shepherd - 1;

    if ( num_clone ) {
      const int ret = qthread_fork_clones_to_local_priority
        ( driver_exec_all   // Function
        , NULL              // Function data block
        , NULL              // Pointer to return value feb
        , jshep             // Shepherd number
        , num_clone - 1     // Number of instances - 1
        );

      assert(ret == QTHREAD_SUCCESS);
    }
  }
#endif

  driver_exec_all( NULL );

  s_active_function     = 0;
  s_active_function_arg = 0;
}

void * QthreadsExec::exec_all_reduce_result()
{
  return s_exec[0]->m_scratch_alloc;
}

} // namespace Impl

} // namespace Kokkos

namespace Kokkos {

namespace Impl {

QthreadsTeamPolicyMember::QthreadsTeamPolicyMember()
  : m_exec( **worker_exec() )
  , m_team_shared( 0, 0 )
  , m_team_size( 1 )
  , m_team_rank( 0 )
  , m_league_size( 1 )
  , m_league_end( 1 )
  , m_league_rank( 0 )
{
  m_exec.shared_reset( m_team_shared );
}

QthreadsTeamPolicyMember::QthreadsTeamPolicyMember( const QthreadsTeamPolicyMember::TaskTeam & )
  : m_exec( **worker_exec() )
  , m_team_shared( 0, 0 )
  , m_team_size( s_number_workers_per_shepherd )
  , m_team_rank( m_exec.shepherd_worker_rank() )
  , m_league_size( 1 )
  , m_league_end( 1 )
  , m_league_rank( 0 )
{
  m_exec.shared_reset( m_team_shared );
}

} // namespace Impl

} // namespace Kokkos

#else
void KOKKOS_SRC_QTHREADS_EXEC_PREVENT_LINK_ERROR() {}
#endif // #if defined( KOKKOS_ENABLE_QTHREADS )

