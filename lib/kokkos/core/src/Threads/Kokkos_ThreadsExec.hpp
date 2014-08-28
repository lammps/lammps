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

#ifndef KOKKOS_THREADSEXEC_HPP
#define KOKKOS_THREADSEXEC_HPP

#include <stdio.h>

#include <utility>
#include <impl/Kokkos_spinwait.hpp>

#include <Kokkos_Atomic.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template< class > struct ThreadsExecAdapter ;

//----------------------------------------------------------------------------

class ThreadsExecTeamMember ;

class ThreadsExec {
public:

  // Fan array has log_2(NT) reduction threads plus 2 scan threads
  // Currently limited to 16k threads.
  enum { MAX_FAN_COUNT    = 16 };
  enum { MAX_THREAD_COUNT = 1 << ( MAX_FAN_COUNT - 2 ) };
  enum { VECTOR_LENGTH    = 8 };

  /** \brief States of a worker thread */
  enum { Terminating ///<  Termination in progress
       , Inactive    ///<  Exists, waiting for work
       , Active      ///<  Exists, performing work
       , Rendezvous  ///<  Exists, waiting in a barrier or reduce

       , ScanCompleted
       , ScanAvailable
       , ReductionAvailable
       };

private:

  friend class ThreadsExecTeamMember ;
  friend class Kokkos::Threads ;

  // Fan-in operations' root is the highest ranking thread
  // to place the 'scan' reduction intermediate values on
  // the threads that need them.
  // For a simple reduction the thread location is arbitrary.

  /** \brief  Reduction memory reserved for team reductions */
  enum { REDUCE_TEAM_BASE = 512 };

  ThreadsExec * const * m_pool_base ; ///< Base for pool fan-in

  void        * m_scratch ;
  int           m_scratch_reduce_end ;
  int           m_scratch_thread_end ;
  int           m_pool_rank ;
  int           m_pool_size ;
  int           m_pool_fan_size ;
  int volatile  m_pool_state ;  ///< State for global synchronizations


  static void global_lock();
  static void global_unlock();
  static bool spawn();

  static void execute_resize_scratch( ThreadsExec & , const void * );
  static void execute_sleep(          ThreadsExec & , const void * );
  static void execute_get_binding(    ThreadsExec & , const void * );

  ThreadsExec( const ThreadsExec & );
  ThreadsExec & operator = ( const ThreadsExec & );

  static void execute_serial( void (*)( ThreadsExec & , const void * ) );

public:

  KOKKOS_INLINE_FUNCTION int pool_size() const { return m_pool_size ; }
  KOKKOS_INLINE_FUNCTION int pool_rank() const { return m_pool_rank ; }

  static int get_thread_count();
  static ThreadsExec * get_thread( const int init_thread_rank );

  inline void * reduce_memory() const { return ((unsigned char *) m_scratch ); }
  inline void * scratch_memory() const { return ((unsigned char *) m_scratch ) + m_scratch_reduce_end ; }

  static void driver(void);

  ~ThreadsExec();
  ThreadsExec();

  static void * resize_scratch( size_t reduce_size , size_t thread_size );

  static void * root_reduce_scratch();

  static bool is_process();

  static void verify_is_process( const std::string & , const bool initialized );

  static int is_initialized();

  static void initialize( unsigned thread_count ,
                          unsigned use_numa_count ,
                          unsigned use_cores_per_numa ,
                          bool allow_asynchronous_threadpool );

  static void finalize();

  /* Given a requested team size, return valid team size */
  static unsigned team_size_valid( unsigned );

  static void print_configuration( std::ostream & , const bool detail = false );

  //------------------------------------

  static void wait_yield( volatile int & , const int );

  //------------------------------------
  // All-thread functions:

  template< class Functor >
  inline
  void fan_in_reduce( const Functor & f ) const
    {
      typedef ReduceAdapter< Functor > Reduce ;

      const int rev_rank  = m_pool_size - ( m_pool_rank + 1 );

      for ( int i = 0 ; i < m_pool_fan_size ; ++i ) {

        ThreadsExec & fan = *m_pool_base[ rev_rank + ( 1 << i ) ] ;

        Impl::spinwait( fan.m_pool_state , ThreadsExec::Active );

        Reduce::join( f , reduce_memory() , fan.reduce_memory() );
      }

      if ( ! rev_rank ) {
        Reduce::final( f , reduce_memory() );
      }
    }

  inline
  void fan_in() const
    {
      const int rev_rank = m_pool_size - ( m_pool_rank + 1 );

      for ( int i = 0 ; i < m_pool_fan_size ; ++i ) {
        Impl::spinwait( m_pool_base[rev_rank+(1<<i)]->m_pool_state , ThreadsExec::Active );
      }
    }

  template< class FunctorType >
  inline
  void scan_large( const FunctorType & f )
    {
      // Sequence of states:
      //  0) Active             : entry and exit state
      //  1) ReductionAvailable : reduction value available
      //  2) ScanAvailable      : inclusive scan value available
      //  3) Rendezvous         : All threads inclusive scan value are available
      //  4) ScanCompleted      : exclusive scan value copied

      typedef ReduceAdapter< FunctorType > Reduce ;
      typedef typename Reduce::scalar_type scalar_type ;

      const int      rev_rank = m_pool_size - ( m_pool_rank + 1 );
      const unsigned count    = Reduce::value_count( f );

      scalar_type * const work_value = (scalar_type *) reduce_memory();

      //--------------------------------
      // Fan-in reduction with highest ranking thread as the root
      for ( int i = 0 ; i < m_pool_fan_size ; ++i ) {
        ThreadsExec & fan = *m_pool_base[ rev_rank + (1<<i) ];

        // Wait: Active -> ReductionAvailable (or ScanAvailable)
        Impl::spinwait( fan.m_pool_state , ThreadsExec::Active );
        Reduce::join( f , work_value , fan.reduce_memory() );
      }

      // Copy reduction value to scan value before releasing from this phase.
      for ( unsigned i = 0 ; i < count ; ++i ) { work_value[i+count] = work_value[i] ; }

      if ( rev_rank ) {

        // Set: Active -> ReductionAvailable
        m_pool_state = ThreadsExec::ReductionAvailable ;

        // Wait for contributing threads' scan value to be available.
        if ( ( 1 << m_pool_fan_size ) < ( m_pool_rank + 1 ) ) {
          ThreadsExec & th = *m_pool_base[ rev_rank + ( 1 << m_pool_fan_size ) ] ;

          // Wait: Active             -> ReductionAvailable
          // Wait: ReductionAvailable -> ScanAvailable
          Impl::spinwait( th.m_pool_state , ThreadsExec::Active );
          Impl::spinwait( th.m_pool_state , ThreadsExec::ReductionAvailable );

          Reduce::join( f , work_value + count , ((scalar_type *)th.reduce_memory()) + count );
        }

        // This thread has completed inclusive scan
        // Set: ReductionAvailable -> ScanAvailable
        m_pool_state = ThreadsExec::ScanAvailable ;

        // Wait for all threads to complete inclusive scan
        // Wait: ScanAvailable -> Rendezvous
        Impl::spinwait( m_pool_state , ThreadsExec::ScanAvailable );
      }

      //--------------------------------

      for ( int i = 0 ; i < m_pool_fan_size ; ++i ) {
        ThreadsExec & fan = *m_pool_base[ rev_rank + (1<<i) ];
        // Wait: ReductionAvailable -> ScanAvailable
        Impl::spinwait( fan.m_pool_state , ThreadsExec::ReductionAvailable );
        // Set: ScanAvailable -> Rendezvous
        fan.m_pool_state = ThreadsExec::Rendezvous ;
      }

      // All threads have completed the inclusive scan.
      // All non-root threads are in the Rendezvous state.
      // Threads are free to overwrite their reduction value.
      //--------------------------------

      if ( ( rev_rank + 1 ) < m_pool_size ) {
        // Exclusive scan: copy the previous thread's inclusive scan value

        ThreadsExec & th = *m_pool_base[ rev_rank + 1 ] ; // Not the root thread

        const scalar_type * const src_value = ((scalar_type *)th.reduce_memory()) + count ;

        for ( unsigned j = 0 ; j < count ; ++j ) { work_value[j] = src_value[j]; }
      }
      else {
        (void) Reduce::init( f , work_value );
      }

      //--------------------------------
      // Wait for all threads to copy previous thread's inclusive scan value
      // Wait for all threads: Rendezvous -> ScanCompleted
      for ( int i = 0 ; i < m_pool_fan_size ; ++i ) {
        Impl::spinwait( m_pool_base[ rev_rank + (1<<i) ]->m_pool_state , ThreadsExec::Rendezvous );
      }
      if ( rev_rank ) {
        // Set: ScanAvailable -> ScanCompleted
        m_pool_state = ThreadsExec::ScanCompleted ;
        // Wait: ScanCompleted -> Active
        Impl::spinwait( m_pool_state , ThreadsExec::ScanCompleted );
      }
      // Set: ScanCompleted -> Active
      for ( int i = 0 ; i < m_pool_fan_size ; ++i ) {
        m_pool_base[ rev_rank + (1<<i) ]->m_pool_state = ThreadsExec::Active ;
      }
    }

  template< class FunctorType >
  inline
  void scan_small( const FunctorType & f )
    {
      typedef ReduceAdapter< FunctorType > Reduce ;
      typedef typename Reduce::scalar_type scalar_type ;

      const int      rev_rank = m_pool_size - ( m_pool_rank + 1 );
      const unsigned count    = Reduce::value_count( f );

      scalar_type * const work_value = (scalar_type *) reduce_memory();

      //--------------------------------
      // Fan-in reduction with highest ranking thread as the root
      for ( int i = 0 ; i < m_pool_fan_size ; ++i ) {
        // Wait: Active -> Rendezvous
        Impl::spinwait( m_pool_base[ rev_rank + (1<<i) ]->m_pool_state , ThreadsExec::Active );
      }

      for ( unsigned i = 0 ; i < count ; ++i ) { work_value[i+count] = work_value[i]; }

      if ( rev_rank ) {
        m_pool_state = ThreadsExec::Rendezvous ;
        // Wait: Rendezvous -> Active
        Impl::spinwait( m_pool_state , ThreadsExec::Rendezvous );
      }
      else {
        // Root thread does the thread-scan before releasing threads

        scalar_type * ptr_prev = 0 ;

        for ( int rank = 0 ; rank < m_pool_size ; ++rank ) {
          scalar_type * const ptr = (scalar_type *) get_thread( rank )->reduce_memory();
          if ( rank ) {
            for ( unsigned i = 0 ; i < count ; ++i ) { ptr[i] = ptr_prev[ i + count ]; }
            Reduce::join( f , ptr + count , ptr );
          }
          else {
            (void) Reduce::init( f , ptr );
          }
          ptr_prev = ptr ;
        }
      }

      for ( int i = 0 ; i < m_pool_fan_size ; ++i ) {
        m_pool_base[ rev_rank + (1<<i) ]->m_pool_state = ThreadsExec::Active ;
      }
    }

  //------------------------------------
  /** \brief  Wait for previous asynchronous functor to
   *          complete and release the Threads device.
   *          Acquire the Threads device and start this functor.
   */
  static void start( void (*)( ThreadsExec & , const void * ) , const void * );

/*
  static unsigned team_max();
  static unsigned team_recommended();
  static unsigned hardware_thread_id();
  static unsigned max_hardware_threads();
*/

  static int  in_parallel();
  static void fence();
  static bool sleep();
  static bool wake();
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class ThreadsExecTeamMember {
private:

  enum { TEAM_REDUCE_SIZE = 512 };

  typedef Kokkos::Threads execution_space ;
  typedef execution_space::scratch_memory_space space ;

  Impl::ThreadsExec   & m_exec ;
  space                 m_team_shared ;
  ThreadsExec * const * m_team_base ; ///< Base for team fan-in
  int                   m_team_shared_size ;
  int                   m_team_size ;
  int                   m_team_rank ;
  int                   m_team_rank_rev ;
  int                   m_league_size ;
  int                   m_league_end ;
  int                   m_league_rank ;

  inline
  void set_team_shared()
    { new( & m_team_shared ) space( ((char *) (*m_team_base)->scratch_memory()) + TEAM_REDUCE_SIZE , m_team_shared_size ); }
  
  // Fan-in and wait until the matching fan-out is called.
  // The root thread which does not wait will return true.
  // All other threads will return false during the fan-out.
  KOKKOS_INLINE_FUNCTION bool team_fan_in() const
    {
      int n , j ;

      // Wait for fan-in threads
      for ( n = 1 ; ( ! ( m_team_rank_rev & n ) ) && ( ( j = m_team_rank_rev + n ) < m_team_size ) ; n <<= 1 ) {
        Impl::spinwait( m_team_base[j]->m_pool_state , ThreadsExec::Active );
      }

      // If not root then wait for release
      if ( m_team_rank_rev ) {
        m_exec.m_pool_state = ThreadsExec::Rendezvous ;
        Impl::spinwait( m_exec.m_pool_state , ThreadsExec::Rendezvous );
      }

      return ! m_team_rank_rev ;
    }

  KOKKOS_INLINE_FUNCTION void team_fan_out() const
    {
      int n , j ;
      for ( n = 1 ; ( ! ( m_team_rank_rev & n ) ) && ( ( j = m_team_rank_rev + n ) < m_team_size ) ; n <<= 1 ) {
        m_team_base[j]->m_pool_state = ThreadsExec::Active ;
      }
    }

public:

  KOKKOS_INLINE_FUNCTION static int team_reduce_size() { return TEAM_REDUCE_SIZE ; }

  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space & team_shmem() const
    { return m_team_shared ; }

  KOKKOS_INLINE_FUNCTION int league_rank() const { return m_league_rank ; }
  KOKKOS_INLINE_FUNCTION int league_size() const { return m_league_size ; }
  KOKKOS_INLINE_FUNCTION int team_rank() const { return m_team_rank ; }
  KOKKOS_INLINE_FUNCTION int team_size() const { return m_team_size ; }

  KOKKOS_INLINE_FUNCTION void team_barrier() const
    {
      team_fan_in();
      team_fan_out();
    }

  template< typename Type >
  KOKKOS_INLINE_FUNCTION Type team_reduce( const Type & value ) const
    {
      // Make sure there is enough scratch space:
      typedef typename if_c< sizeof(Type) < ThreadsExec::REDUCE_TEAM_BASE , Type , void >::type type ;

      *((volatile type*) m_exec.scratch_memory() ) = value ;

      memory_fence();

      type & accum = *((type *) m_team_base[0]->scratch_memory() );

      if ( team_fan_in() ) {
        for ( int i = 1 ; i < m_team_size ; ++i ) {
          accum += *((type *) m_team_base[i]->scratch_memory() );
        }
        memory_fence();
      }

      team_fan_out();

      return accum ;
    }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering
   *          with intra-team non-deterministic ordering accumulation.
   *
   *  The global inter-team accumulation value will, at the end of the
   *  league's parallel execution, be the scan's total.
   *  Parallel execution ordering of the league's teams is non-deterministic.
   *  As such the base value for each team's scan operation is similarly
   *  non-deterministic.
   */
  template< typename ArgType >
  KOKKOS_INLINE_FUNCTION ArgType team_scan( const ArgType & value , ArgType * const global_accum ) const
    {
      // Make sure there is enough scratch space:
      typedef typename if_c< sizeof(ArgType) < ThreadsExec::REDUCE_TEAM_BASE , ArgType , void >::type type ;

      volatile type * const work_value  = ((type*) m_exec.scratch_memory());

      *work_value = value ;

      memory_fence();

      if ( team_fan_in() ) {
        // The last thread to synchronize returns true, all other threads wait for team_fan_out()
        // m_team_base[0]                 == highest ranking team member
        // m_team_base[ m_team_size - 1 ] == lowest ranking team member
        //
        // 1) copy from lower to higher rank, initialize lowest rank to zero
        // 2) prefix sum from lowest to highest rank, skipping lowest rank

        type accum = 0 ;

        if ( global_accum ) {
          for ( int i = m_team_size ; i-- ; ) {
            type & val = *((type*) m_team_base[i]->scratch_memory());
            accum += val ;
          }
          accum = atomic_fetch_add( global_accum , accum );
        }

        for ( int i = m_team_size ; i-- ; ) {
          type & val = *((type*) m_team_base[i]->scratch_memory());
          const type offset = accum ;
          accum += val ;
          val = offset ;
        }

        memory_fence();
      }

      team_fan_out();

      return *work_value ;
    }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
   *
   *  The highest rank thread can compute the reduction total as
   *    reduction_total = dev.team_scan( value ) + value ;
   */
  template< typename ArgType >
  KOKKOS_INLINE_FUNCTION ArgType team_scan( const ArgType & value ) const
    { return this-> template team_scan<ArgType>( value , 0 ); }

  //----------------------------------------
  // Private for the driver

  template< class WorkArgTag >
  ThreadsExecTeamMember( Impl::ThreadsExec & exec
                       , const TeamPolicy< execution_space , WorkArgTag > & team 
                       , const int shared_size )
    : m_exec( exec )
    , m_team_shared(0,0)
    , m_team_base(0)
    , m_team_shared_size( shared_size )
    , m_team_size(0)
    , m_team_rank(0)
    , m_team_rank_rev(0)
    , m_league_size(0)
    , m_league_end(0)
    , m_league_rank(0)
    {
      if ( team.league_size() ) {
        // Execution is using device-team interface:

        const int pool_rank_rev = exec.pool_size() - ( exec.pool_rank() + 1 );
        const int team_rank_rev = pool_rank_rev % team.team_alloc();

        // May be using fewer threads per team than a multiple of threads per core,
        // some threads will idle.

        if ( team_rank_rev < team.team_size() ) {
          const size_t pool_league_size     = exec.pool_size() / team.team_alloc() ;
          const size_t pool_league_rank_rev = pool_rank_rev / team.team_alloc() ;
          const size_t pool_league_rank     = pool_league_size - ( pool_league_rank_rev + 1 );

          m_team_base        = exec.m_pool_base + team.team_alloc() * pool_league_rank_rev ;
          m_team_size        = team.team_size() ;
          m_team_rank        = team.team_size() - ( team_rank_rev + 1 );
          m_team_rank_rev    = team_rank_rev ;
          m_league_size      = team.league_size();
          m_league_rank      = ( team.league_size() *  pool_league_rank    ) / pool_league_size ;
          m_league_end       = ( team.league_size() * (pool_league_rank+1) ) / pool_league_size ;

          set_team_shared();
        }
      }
    }

  bool valid() const
    { return m_league_rank < m_league_end ; }

  void next()
    {
      if ( ++m_league_rank < m_league_end ) {
        team_barrier();
        set_team_shared();
      }
    }
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

inline int Threads::in_parallel()
{ return Impl::ThreadsExec::in_parallel(); }

inline int Threads::is_initialized()
{ return Impl::ThreadsExec::is_initialized(); }

inline void Threads::initialize(
  unsigned threads_count ,
  unsigned use_numa_count ,
  unsigned use_cores_per_numa ,
  bool allow_asynchronous_threadpool )
{
  Impl::ThreadsExec::initialize( threads_count , use_numa_count , use_cores_per_numa , allow_asynchronous_threadpool );
}

inline void Threads::finalize()
{
  Impl::ThreadsExec::finalize();
}

inline void Threads::print_configuration( std::ostream & s , const bool detail )
{
  Impl::ThreadsExec::print_configuration( s , detail );
}

inline bool Threads::sleep()
{ return Impl::ThreadsExec::sleep() ; }

inline bool Threads::wake()
{ return Impl::ThreadsExec::wake() ; }

inline void Threads::fence()
{ Impl::ThreadsExec::fence() ; }

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template < class WorkArgTag >
class TeamPolicy< Kokkos::Threads , WorkArgTag > {
private:

  int m_league_size ;
  int m_team_size ;
  int m_team_alloc ;

  inline
  void init( const int league_size_request 
           , const int team_size_request )
   {
      const int pool_size  = execution_space::thread_pool_size(0);
      const int team_max   = execution_space::thread_pool_size(1);
      const int team_grain = execution_space::thread_pool_size(2);

      m_league_size = league_size_request ;

      m_team_size = team_size_request < team_max ?
                    team_size_request : team_max ;

      // Round team size up to a multiple of 'team_gain'
      const int team_size_grain = team_grain * ( ( m_team_size + team_grain - 1 ) / team_grain );
      const int team_count      = pool_size / team_size_grain ;

      // Constraint : pool_size = m_team_alloc * team_count
      m_team_alloc = pool_size / team_count ;
   }


public:

  typedef Impl::ExecutionPolicyTag   kokkos_tag ;      ///< Concept tag
  typedef Kokkos::Threads            execution_space ; ///< Execution space

  inline int team_size() const { return m_team_size ; }
  inline int team_alloc() const { return m_team_alloc ; }
  inline int league_size() const { return m_league_size ; }

  /** \brief  Specify league size, request team size */
  TeamPolicy( execution_space & , int league_size_request , int team_size_request )
    : m_league_size(0)
    , m_team_size(0)
    , m_team_alloc(0)
    { init(league_size_request,team_size_request); }

  TeamPolicy( int league_size_request , int team_size_request )
    : m_league_size(0)
    , m_team_size(0)
    , m_team_alloc(0)
    { init(league_size_request,team_size_request); }

  template< class FunctorType >
  inline static
  int team_size_max( const FunctorType & )
    { return execution_space::thread_pool_size(1); }

  typedef Impl::ThreadsExecTeamMember member_type ;

  friend class Impl::ThreadsExecTeamMember ;
};

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #define KOKKOS_THREADSEXEC_HPP */

