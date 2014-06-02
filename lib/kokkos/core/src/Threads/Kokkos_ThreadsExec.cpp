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

#include <limits>
#include <utility>
#include <iostream>
#include <Kokkos_Threads.hpp>
#include <Kokkos_hwloc.hpp>
#include <Kokkos_Atomic.hpp>

#include <stdint.h>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {
namespace {

ThreadsExec                  s_threads_process ;
ThreadsExec                * s_threads_exec[  ThreadsExec::MAX_THREAD_COUNT ];
std::pair<unsigned,unsigned> s_threads_coord[ ThreadsExec::MAX_THREAD_COUNT ];

unsigned s_threads_count       = 0 ;
unsigned s_threads_per_numa    = 0 ;
unsigned s_threads_per_core    = 0 ;

unsigned s_current_reduce_size = 0 ;
unsigned s_current_shared_size = 0 ;
unsigned s_current_team_alloc  = 0 ;
unsigned s_current_team_size   = 0 ;
unsigned s_current_league_size = 0 ;

void (* volatile s_current_function)( ThreadsExec & , const void * );
const void * volatile s_current_function_arg = 0 ;

struct Sentinel {
  Sentinel()
  {
    HostSpace::register_in_parallel( ThreadsExec::in_parallel );
  }

  ~Sentinel()
  {
    if ( s_threads_count ||
         s_threads_per_numa ||
         s_threads_per_core ||
         s_current_reduce_size ||
         s_current_shared_size ||
         s_current_function ||
         s_current_function_arg ||
         s_threads_exec[0] ) {
      std::cerr << "ERROR : Process exiting without calling Kokkos::Threads::terminate()" << std::endl ;
    }
  }
};

inline
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

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

void execute_function_noop( ThreadsExec & , const void * ) {}

void ThreadsExec::driver(void)
{
  ThreadsExec this_thread ;

  while ( ThreadsExec::Active == this_thread.m_pool_state ) {

    this_thread.set_team_relations();

    (*s_current_function)( this_thread , s_current_function_arg );

    // Deactivate thread and wait for reactivation
    this_thread.m_pool_state = ThreadsExec::Inactive ;

    wait_yield( this_thread.m_pool_state , ThreadsExec::Inactive );
  }
}

void ThreadsExec::set_team_relations()
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

  const size_t league_size = s_current_league_size ;

  if ( league_size ) {
    // Execution is using device-team interface:

    const unsigned team_alloc    = s_current_team_alloc ;
    const unsigned team_size     = s_current_team_size ;
    const unsigned pool_rank_rev = m_pool_size - ( m_pool_rank + 1 );
    const unsigned team_rank_rev = pool_rank_rev % team_alloc ;

    // May be using fewer threads per team than a multiple of threads per core,
    // some threads will idle.

    if ( team_rank_rev < team_size ) {
      const size_t pool_league_size     = m_pool_size   / team_alloc ;
      const size_t pool_league_rank_rev = pool_rank_rev / team_alloc ;
      const size_t pool_league_rank     = pool_league_size - ( pool_league_rank_rev + 1 );

      m_team_base        = m_pool_base + team_alloc * pool_league_rank_rev ;
      m_team_shared      = (*m_team_base)->m_alloc_shared ;
      m_team_shared_end  = s_current_shared_size ;
      m_team_size        = team_size ;
      m_team_rank        = team_size - ( team_rank_rev + 1 );
      m_team_fan_size    = fan_size( m_team_rank , team_size );
      m_league_size      = league_size ;
      m_league_rank      = ( league_size *  pool_league_rank    ) / pool_league_size ;
      m_league_end       = ( league_size * (pool_league_rank+1) ) / pool_league_size ;
    }
  }
}

ThreadsExec::ThreadsExec()
  : m_pool_base(0)
  , m_team_base(0)
  , m_alloc_reduce(0)
  , m_alloc_shared(0)
  , m_team_shared(0)
  , m_team_shared_end(0)
  , m_team_shared_iter(0)

  , m_pool_rank(0)
  , m_pool_size(0)
  , m_pool_fan_size(0)

  , m_team_rank(0)
  , m_team_size(0)
  , m_team_fan_size(0)

  , m_league_rank(0)
  , m_league_end(0)
  , m_league_size(0)

  , m_pool_state( ThreadsExec::Terminating )
  , m_team_state( ThreadsExec::Inactive )
{
  if ( & s_threads_process != this ) {

    // A spawned thread

    ThreadsExec * const nil = 0 ;

    // Which entry in 's_threads_exec', possibly determined from hwloc binding
    const unsigned entry = ((size_t)s_current_function_arg) < s_threads_count
                         ? ((size_t)s_current_function_arg)
                         : size_t(Kokkos::hwloc::bind_this_thread( s_threads_count , s_threads_coord ));

    // Given a good entry set this thread in the 's_threads_exec' array
    if ( entry < s_threads_count &&
         nil == atomic_compare_exchange( s_threads_exec + entry , nil , this ) ) {

      m_pool_base     = s_threads_exec ;
      m_pool_rank     = s_threads_count - ( entry + 1 );
      m_pool_size     = s_threads_count ;
      m_pool_fan_size = fan_size( m_pool_rank , m_pool_size );
      m_pool_state    = ThreadsExec::Active ;

      // Inform spawning process that the threads_exec entry has been set.
      s_threads_process.m_pool_state = ThreadsExec::Active ;
    }
    else {
      // Inform spawning process that the threads_exec entry could not be set.
      s_threads_process.m_pool_state = ThreadsExec::Terminating ;
    }
  }
  else {
    // Enables 'parallel_for' to execute on unitialized Threads device
    m_pool_rank  = 0 ;
    m_pool_size  = 1 ;
    m_pool_state = ThreadsExec::Inactive ;
  }
}

ThreadsExec::~ThreadsExec()
{
  const unsigned entry = m_pool_size - ( m_pool_rank + 1 );

  m_pool_base   = 0 ;
  m_team_base   = 0 ;

  m_alloc_reduce     = 0 ;
  m_alloc_shared     = 0 ;
  m_team_shared      = 0 ;
  m_team_shared_end  = 0 ;
  m_team_shared_iter = 0 ;

  m_pool_rank     = 0 ;
  m_pool_size     = 0 ;
  m_pool_fan_size = 0 ;
  m_team_rank     = 0 ;
  m_team_size     = 0 ;
  m_team_fan_size = 0 ;
  m_league_rank   = 0 ;
  m_league_end    = 0 ;
  m_league_size   = 0 ;

  m_pool_state  = ThreadsExec::Terminating ;
  m_team_state  = ThreadsExec::Inactive ;

  if ( & s_threads_process != this && entry < MAX_THREAD_COUNT ) {
    ThreadsExec * const nil = 0 ;

    atomic_compare_exchange( s_threads_exec + entry , this , nil );

    s_threads_process.m_pool_state = ThreadsExec::Terminating ;
  }
}


int ThreadsExec::get_thread_count()
{
  return s_threads_count ;
}

ThreadsExec * ThreadsExec::get_thread( const int init_thread_rank )
{
  ThreadsExec * const th =
    unsigned(init_thread_rank) < s_threads_count
    ? s_threads_exec[ s_threads_count - ( init_thread_rank + 1 ) ] : 0 ;

  if ( 0 == th || th->m_pool_rank != init_thread_rank ) {
    std::ostringstream msg ;
    msg << "Kokkos::Impl::ThreadsExec::get_thread ERROR : "
        << "thread " << init_thread_rank << " of " << s_threads_count ;
    if ( 0 == th ) {
      msg << " does not exist" ;
    }
    else {
      msg << " has wrong thread_rank " << th->m_pool_rank ;
    }
    Kokkos::Impl::throw_runtime_exception( msg.str() );
  }

  return th ;
}

//----------------------------------------------------------------------------

void ThreadsExec::execute_get_binding( ThreadsExec & exec , const void * )
{
  s_threads_coord[ exec.m_pool_rank ] = Kokkos::hwloc::get_this_thread_coordinate();
}

void ThreadsExec::execute_sleep( ThreadsExec & exec , const void * )
{
  ThreadsExec::global_lock();
  ThreadsExec::global_unlock();

  const int n = exec.m_pool_fan_size ;
  const int rank_rev = exec.m_pool_size - ( exec.m_pool_rank + 1 );

  for ( int i = 0 ; i < n ; ++i ) {
    Impl::spinwait( exec.m_pool_base[ rank_rev + (1<<i) ]->m_pool_state , ThreadsExec::Active );
  }

  exec.m_pool_state = ThreadsExec::Inactive ;
}

void ThreadsExec::execute_reduce_resize( ThreadsExec & exec , const void * )
{
  if ( exec.m_alloc_reduce ) {
    HostSpace::decrement( exec.m_alloc_reduce );
    exec.m_alloc_reduce = 0 ;
  }

  if ( s_current_reduce_size ) {

    exec.m_alloc_reduce =
      HostSpace::allocate( "reduce_scratch_space" , typeid(unsigned char) , 1 , s_current_reduce_size );

    // Guaranteed multiple of 'unsigned'

    unsigned * ptr = (unsigned *)( exec.m_alloc_reduce );
    unsigned * const end = ptr + s_current_reduce_size / sizeof(unsigned);

    // touch on this thread
    while ( ptr < end ) *ptr++ = 0 ;
  }
}

void ThreadsExec::execute_shared_resize( ThreadsExec & exec , const void * )
{
  // First thread pinned to a core allocates shared memory
  const int rank_rev = exec.m_pool_size - ( exec.m_pool_rank + 1 );

  if ( ! ( rank_rev % s_threads_per_core ) ) {

    if ( exec.m_alloc_shared ) {
      HostSpace::decrement( exec.m_alloc_shared );
      exec.m_alloc_shared = 0 ;
    }

    if ( s_current_shared_size ) {
      exec.m_alloc_shared =
        HostSpace::allocate( "shared_scratch_space" , typeid(unsigned char) , 1 , s_current_shared_size );

      // Guaranteed multiple of 'unsigned'

      unsigned * ptr = (unsigned *)( exec.m_alloc_shared );
      unsigned * const end = ptr + s_current_shared_size / sizeof(unsigned);

      // touch on this thread
      while ( ptr < end ) *ptr++ = 0 ;
    }
  }
  else {
    exec.m_alloc_shared = 0 ;
  }
}

void * ThreadsExec::get_shmem( const int size )
{
  // m_team_shared_iter is in bytes, convert to integer offsets
  const int offset = m_team_shared_iter >> power_of_two<sizeof(int)>::value ;

  m_team_shared_iter += size ;

  if ( m_team_shared_end < m_team_shared_iter ) {
    Kokkos::Impl::throw_runtime_exception( std::string("ThreadsExec::get_shmem FAILED : exceeded shared memory size" ) );
  }

  return ((int*)m_team_shared) + offset ;
}

}
}

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

void ThreadsExec::verify_is_process( const std::string & name , const bool initialized )
{
  if ( ! is_process() ) {
    std::string msg( name );
    msg.append( " FAILED : Called by a worker thread, can only be called by the master process." );
    Kokkos::Impl::throw_runtime_exception( msg );
  }

  if ( initialized && 0 == s_threads_count ) {
    std::string msg( name );
    msg.append( " FAILED : Threads not initialized." );
    Kokkos::Impl::throw_runtime_exception( msg );
  }
}

int ThreadsExec::in_parallel()
{
  // A thread function is in execution and
  // the function argument is not the special threads process argument and
  // the master process is a worker or is not the master process.
  return s_current_function &&
         ( & s_threads_process != s_current_function_arg ) &&
         ( s_threads_process.m_pool_base || ! is_process() );
}

// Wait for root thread to become inactive
void ThreadsExec::fence()
{
  if ( s_threads_count ) {
    // Wait for the root thread to complete:
    Impl::spinwait( s_threads_exec[0]->m_pool_state , ThreadsExec::Active );
  }

  s_current_function     = 0 ;
  s_current_function_arg = 0 ;
  s_current_team_size    = 0 ;
  s_current_team_alloc   = 0 ;
  s_current_league_size  = 0 ;
}

/** \brief  Begin execution of the asynchronous functor */
void ThreadsExec::start( void (*func)( ThreadsExec & , const void * ) , const void * arg ,
                         int work_league_size ,
                         int work_team_size )
{
  const bool work_spec = work_league_size || work_team_size ;

  verify_is_process("ThreadsExec::start" , work_spec );

  if ( s_current_function || s_current_function_arg ) {
    Kokkos::Impl::throw_runtime_exception( std::string( "ThreadsExec::start() FAILED : already executing" ) );
  }

  if ( work_spec ) {
    s_current_team_size    = work_team_size ? std::min( s_threads_per_numa , unsigned(work_team_size) ) : s_threads_per_numa ;
    s_current_team_alloc   = s_threads_per_core * ( ( s_current_team_size + s_threads_per_core - 1 ) / s_threads_per_core );
    s_current_league_size  = work_league_size ;
  }

  s_current_function     = func ;
  s_current_function_arg = arg ;

  // Activate threads:
  for ( int i = s_threads_count ; 0 < i-- ; ) {
    s_threads_exec[i]->m_pool_state = ThreadsExec::Active ;
  }

  if ( s_threads_process.m_pool_size ) {
    // Master process is the root thread, run it:
    s_threads_process.set_team_relations();
    (*func)( s_threads_process , arg );
    s_threads_process.m_pool_state = ThreadsExec::Inactive ;
  }
}

//----------------------------------------------------------------------------

bool ThreadsExec::sleep()
{
  verify_is_process("ThreadsExec::sleep", true );

  if ( & execute_sleep == s_current_function ) return false ;

  fence();

  ThreadsExec::global_lock();

  s_current_function = & execute_sleep ;

  // Activate threads:
  for ( unsigned i = s_threads_count ; 0 < i ; ) {
    s_threads_exec[--i]->m_pool_state = ThreadsExec::Active ;
  }

  return true ;
}

bool ThreadsExec::wake()
{
  verify_is_process("ThreadsExec::wake", true );

  if ( & execute_sleep != s_current_function ) return false ;

  ThreadsExec::global_unlock();

  if ( s_threads_process.m_pool_base ) {
    execute_sleep( s_threads_process , 0 );
    s_threads_process.m_pool_state = ThreadsExec::Inactive ;
  }

  fence();

  return true ;
}

//----------------------------------------------------------------------------

void ThreadsExec::execute_serial( void (*func)( ThreadsExec & , const void * ) )
{
  s_current_function = func ;
  s_current_function_arg = & s_threads_process ;

  const unsigned begin = s_threads_process.m_pool_base ? 1 : 0 ;

  for ( unsigned i = s_threads_count ; begin < i ; ) {
    ThreadsExec & th = * s_threads_exec[ --i ];

    th.m_pool_state = ThreadsExec::Active ;

    wait_yield( th.m_pool_state , ThreadsExec::Active );
  }

  if ( s_threads_process.m_pool_base ) {
    s_threads_process.m_pool_state = ThreadsExec::Active ;
    (*func)( s_threads_process , 0 );
    s_threads_process.m_pool_state = ThreadsExec::Inactive ;
  }

  s_current_function_arg = 0 ;
  s_current_function = 0 ;
}

//----------------------------------------------------------------------------

void * ThreadsExec::root_reduce_scratch()
{
  return s_threads_process.reduce_base();
}

void ThreadsExec::resize_reduce_scratch( size_t size )
{
  fence();

  if ( size ) { size += REDUCE_TEAM_BASE ; }

  const size_t rem = size % Kokkos::Impl::MEMORY_ALIGNMENT ;

  if ( rem ) size += Kokkos::Impl::MEMORY_ALIGNMENT - rem ;

  if ( ( s_current_reduce_size < size ) ||
       ( 0 == size && s_current_reduce_size ) ) {

    verify_is_process( "ThreadsExec::resize_reduce_scratch" , true );

    s_current_reduce_size = size ;

    execute_serial( & execute_reduce_resize );

    s_threads_process.m_alloc_reduce = s_threads_exec[0]->m_alloc_reduce ;
  }
}

void ThreadsExec::resize_shared_scratch( size_t size )
{
  fence();

  const size_t rem = size % Kokkos::Impl::MEMORY_ALIGNMENT ;

  if ( rem ) size += Kokkos::Impl::MEMORY_ALIGNMENT - rem ;

  if ( s_current_shared_size < size || ( 0 == size && s_current_shared_size ) ) {

    verify_is_process( "ThreadsExec::resize_shared_scratch" , true );

    s_current_shared_size = size ;

    execute_serial( & execute_shared_resize );
  }
}

//----------------------------------------------------------------------------

void ThreadsExec::print_configuration( std::ostream & s , const bool detail )
{
  verify_is_process("ThreadsExec::print_configuration",false);

  fence();

  const unsigned numa_count       = Kokkos::hwloc::get_available_numa_count();
  const unsigned cores_per_numa   = Kokkos::hwloc::get_available_cores_per_numa();
  const unsigned threads_per_core = Kokkos::hwloc::get_available_threads_per_core();

  // Forestall compiler warnings for unused variables.
  (void) numa_count;
  (void) cores_per_numa;
  (void) threads_per_core;

  s << "Kokkos::Threads" ;

#if defined( KOKKOS_HAVE_PTHREAD )
  s << " KOKKOS_HAVE_PTHREAD" ;
#endif
#if defined( KOKKOS_HAVE_HWLOC )
  s << " hwloc[" << numa_count << "x" << cores_per_numa << "x" << threads_per_core << "]" ;
#endif

  if ( s_threads_count ) {
    s << " threads[" << s_threads_count << "]"
      << " threads_per_numa[" << s_threads_per_numa << "]"
      << " threads_per_core[" << s_threads_per_core << "]"
      ;
    if ( 0 == s_threads_process.m_pool_base ) { s << " Asynchronous" ; }
    s << " ReduceScratch[" << s_current_reduce_size << "]"
      << " SharedScratch[" << s_current_shared_size << "]" ;
    s << std::endl ;

    if ( detail ) {

      execute_serial( & execute_get_binding );

      for ( unsigned i = 0 ; i < s_threads_count ; ++i ) {
        ThreadsExec * const th = s_threads_exec[i] ;
        s << "  Thread hwloc("
          << s_threads_coord[i].first << "."
          << s_threads_coord[i].second << ")" ;

        s_threads_coord[i].first  = ~0u ;
        s_threads_coord[i].second = ~0u ;

        if ( th ) {
          const int rank_rev = th->m_pool_size - ( th->m_pool_rank + 1 );

          s << " rank(" << th->m_pool_rank << ")" ;

          if ( th->m_pool_fan_size ) {
            s << " Fan{" ;
            for ( int j = 0 ; j < th->m_pool_fan_size ; ++j ) {
              s << " " << th->m_pool_base[rank_rev+(1<<j)]->m_pool_rank ;
            }
            s << " }" ;
          }

          if ( th->m_team_base && th->m_team_size ) {
            s << " Team[ " << th->m_team_base[0]->m_pool_rank
              << " .. " << th->m_team_base[ th->m_team_size - 1 ]->m_pool_rank
              << " ]" ;
          }

          if ( th == & s_threads_process ) {
            s << " is_process" ;
          }
        }
        s << std::endl ;
      }
    }
  }
  else {
    s << " not initialized" << std::endl ;
  }
}

//----------------------------------------------------------------------------

int ThreadsExec::league_max()
{ return std::numeric_limits<int>::max(); }

int ThreadsExec::team_max()
{ return s_threads_per_numa ; }

//----------------------------------------------------------------------------

int ThreadsExec::is_initialized()
{ return 0 != s_threads_exec[0] ; }

void ThreadsExec::initialize( unsigned thread_count ,
                              unsigned use_numa_count ,
                              unsigned use_cores_per_numa ,
                              bool allow_asynchronous_threadpool )
{
  static const Sentinel sentinel ;

  const bool is_initialized = 0 != s_threads_count ;

  unsigned thread_spawn_failed = 0 ;

  if ( ! is_initialized ) {

    // If thread_count, use_numa_count, or use_cores_per_numa are zero
    // then they will be given default values based upon hwloc detection
    // and allowed asynchronous execution.

    const bool hwloc_avail = hwloc::available();

    const unsigned thread_spawn_begin =
      hwloc::thread_mapping( "Kokkos::Threads::initialize" ,
                             allow_asynchronous_threadpool ,
                             thread_count ,
                             use_numa_count ,
                             use_cores_per_numa ,
                             s_threads_coord );

    const std::pair<unsigned,unsigned> proc_coord = s_threads_coord[0] ;

    if ( thread_spawn_begin ) {
      // Synchronous with s_threads_coord[0] as the process core
      // Claim entry #0 for binding the process core.
      s_threads_coord[0] = std::pair<unsigned,unsigned>(~0u,~0u);
    }

    s_threads_count    = thread_count ;
    s_threads_per_numa = s_threads_count / use_numa_count ;
    s_threads_per_core = s_threads_per_numa / use_cores_per_numa ;
    s_current_function = & execute_function_noop ; // Initialization work function

    for ( unsigned ith = thread_spawn_begin ; ith < thread_count ; ++ith ) {

      s_threads_process.m_pool_state = ThreadsExec::Inactive ;

      // If hwloc available then spawned thread will
      // choose its own entry in 's_threads_coord'
      // otherwise specify the entry.
      s_current_function_arg = (void*)static_cast<uintptr_t>( hwloc_avail ? ~0u : ith );

      // Spawn thread executing the 'driver()' function.
      // Wait until spawned thread has attempted to initialize.
      // If spawning and initialization is successfull then
      // an entry in 's_threads_exec' will be assigned.
      if ( ThreadsExec::spawn() ) {
        wait_yield( s_threads_process.m_pool_state , ThreadsExec::Inactive );
      }
      if ( s_threads_process.m_pool_state == ThreadsExec::Terminating ) break ;
    }

    // Wait for all spawned threads to deactivate before zeroing the function.

    for ( unsigned ith = thread_spawn_begin ; ith < thread_count ; ++ith ) {
      // Try to protect against cache coherency failure by casting to volatile.
      ThreadsExec * const th = ((ThreadsExec * volatile *)s_threads_exec)[ith] ;
      if ( th ) {
        wait_yield( th->m_pool_state , ThreadsExec::Active );
      }
      else {
        ++thread_spawn_failed ;
      }
    }

    s_current_function     = 0 ;
    s_current_function_arg = 0 ;
    s_threads_process.m_pool_state = ThreadsExec::Inactive ;

    if ( ! thread_spawn_failed ) {
      // Bind process to the core on which it was located before spawning occured
      Kokkos::hwloc::bind_this_thread( proc_coord );

      if ( thread_spawn_begin ) { // Include process in pool.
        s_threads_exec[0]                 = & s_threads_process ;
        s_threads_process.m_pool_base     = s_threads_exec ;
        s_threads_process.m_pool_rank     = thread_count - 1 ; // Reversed for scan-compatible reductions
        s_threads_process.m_pool_size     = thread_count ;
        s_threads_process.m_pool_fan_size = fan_size( s_threads_process.m_pool_rank , s_threads_process.m_pool_size );
      }
      else {
        s_threads_process.m_pool_base = 0 ;
        s_threads_process.m_pool_rank = 0 ;
        s_threads_process.m_pool_size = 0 ;
        s_threads_process.m_pool_fan_size = 0 ;
      }

      // Initial allocations:
      ThreadsExec::resize_reduce_scratch( 4096 - REDUCE_TEAM_BASE );
      ThreadsExec::resize_shared_scratch( 4096 );
    }
    else {
      s_threads_count    = 0 ;
      s_threads_per_numa = 0 ;
      s_threads_per_core = 0 ;
    }
  }

  if ( is_initialized || thread_spawn_failed ) {

    std::ostringstream msg ;

    msg << "Kokkos::Threads::initialize ERROR" ;

    if ( is_initialized ) {
      msg << " : already initialized" ;
    }
    if ( thread_spawn_failed ) {
      msg << " : failed to spawn " << thread_spawn_failed << " threads" ;
    }

    Kokkos::Impl::throw_runtime_exception( msg.str() );
  }
}

//----------------------------------------------------------------------------

void ThreadsExec::finalize()
{
  verify_is_process("ThreadsExec::finalize",false);

  fence();

  resize_reduce_scratch(0);
  resize_shared_scratch(0);

  const unsigned begin = s_threads_process.m_pool_base ? 1 : 0 ;

  for ( unsigned i = s_threads_count ; begin < i-- ; ) {

    if ( s_threads_exec[i] ) {

      s_threads_exec[i]->m_pool_state = ThreadsExec::Terminating ;

      wait_yield( s_threads_process.m_pool_state , ThreadsExec::Inactive );

      s_threads_process.m_pool_state = ThreadsExec::Inactive ;
    }
  }

  if ( s_threads_process.m_pool_base ) {
    ( & s_threads_process )->~ThreadsExec();
    s_threads_exec[0] = 0 ;
  }

  Kokkos::hwloc::unbind_this_thread();

  s_threads_count    = 0 ;
  s_threads_per_numa = 0 ;
  s_threads_per_core = 0 ;

  // Reset master thread to run solo.
  s_threads_process.m_pool_base     = 0 ;
  s_threads_process.m_pool_rank     = 0 ;
  s_threads_process.m_pool_size     = 1 ;
  s_threads_process.m_pool_fan_size = 0 ;
  s_threads_process.m_pool_state = ThreadsExec::Inactive ;
}

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace Kokkos */


