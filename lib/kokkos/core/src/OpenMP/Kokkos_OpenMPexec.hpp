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

#ifndef KOKKOS_OPENMPEXEC_HPP
#define KOKKOS_OPENMPEXEC_HPP

#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_spinwait.hpp>

#include <Kokkos_Atomic.hpp>

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
/** \brief  Data for OpenMP thread execution */

class OpenMPexec {
public:

  enum { MAX_THREAD_COUNT = 4096 };

private:

  static int          m_pool_topo[ 4 ];
  static int          m_map_rank[ MAX_THREAD_COUNT ];
  static OpenMPexec * m_pool[ MAX_THREAD_COUNT ]; // Indexed by: m_pool_rank_rev

  friend class Kokkos::OpenMP ;

  int const  m_pool_rank ;
  int const  m_pool_rank_rev ;
  int const  m_scratch_exec_end ;
  int const  m_scratch_reduce_end ;
  int const  m_scratch_thread_end ;

  int volatile  m_barrier_state ;

  OpenMPexec();
  OpenMPexec( const OpenMPexec & );
  OpenMPexec & operator = ( const OpenMPexec & );

  static void clear_scratch();

public:

  // Topology of a cache coherent thread pool:
  //   TOTAL = NUMA x GRAIN
  //   pool_size( depth = 0 )
  //   pool_size(0) = total number of threads
  //   pool_size(1) = number of threads per NUMA
  //   pool_size(2) = number of threads sharing finest grain memory hierarchy

  inline static
  int pool_size( int depth = 0 ) { return m_pool_topo[ depth ]; }

  inline static
  OpenMPexec * pool_rev( int pool_rank_rev ) { return m_pool[ pool_rank_rev ]; }

  inline int pool_rank() const { return m_pool_rank ; }
  inline int pool_rank_rev() const { return m_pool_rank_rev ; }

  inline void * scratch_reduce() const { return ((char *) this) + m_scratch_exec_end ; }
  inline void * scratch_thread() const { return ((char *) this) + m_scratch_reduce_end ; }

  inline
  void state_wait( int state )
    { Impl::spinwait( m_barrier_state , state ); }

  inline
  void state_set( int state ) { m_barrier_state = state ; }

  ~OpenMPexec() {}

  OpenMPexec( const int poolRank 
            , const int scratch_exec_size
            , const int scratch_reduce_size
            , const int scratch_thread_size )
    : m_pool_rank( poolRank )
    , m_pool_rank_rev( pool_size() - ( poolRank + 1 ) )
    , m_scratch_exec_end( scratch_exec_size )
    , m_scratch_reduce_end( m_scratch_exec_end   + scratch_reduce_size )
    , m_scratch_thread_end( m_scratch_reduce_end + scratch_thread_size )
    , m_barrier_state(0)
    {}

  static void finalize();

  static void initialize( const unsigned  team_count ,
                          const unsigned threads_per_team ,
                          const unsigned numa_count ,
                          const unsigned cores_per_numa );

  static void verify_is_process( const char * const );
  static void verify_initialized( const char * const );

  static void resize_scratch( size_t reduce_size , size_t thread_size );

  inline static
  OpenMPexec * get_thread_omp() { return m_pool[ m_map_rank[ omp_get_thread_num() ] ]; }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

class OpenMPexecTeamMember {
private:

  enum { TEAM_REDUCE_SIZE = 512 };

  /** \brief  Thread states for team synchronization */
  enum { Active = 0 , Rendezvous = 1 };

  typedef Kokkos::OpenMP                         execution_space ;
  typedef execution_space::scratch_memory_space  scratch_memory_space ;

  Impl::OpenMPexec    & m_exec ;
  scratch_memory_space  m_team_shared ;
  int                   m_team_shmem ;
  int                   m_team_base_rev ;
  int                   m_team_rank_rev ;
  int                   m_team_rank ;
  int                   m_team_size ;
  int                   m_league_rank ;
  int                   m_league_end ;
  int                   m_league_size ;

  // Fan-in team threads, root of the fan-in which does not block returns true
  inline
  bool team_fan_in() const
    {
      for ( int n = 1 , j ; ( ( j = m_team_rank_rev + n ) < m_team_size ) && ! ( m_team_rank_rev & n ) ; n <<= 1 ) {
        m_exec.pool_rev( m_team_base_rev + j )->state_wait( Active );
      }

      if ( m_team_rank_rev ) {
        m_exec.state_set( Rendezvous );
        m_exec.state_wait( Rendezvous );
      }

      return 0 == m_team_rank_rev ;
    }

  inline
  void team_fan_out() const
    {
      for ( int n = 1 , j ; ( ( j = m_team_rank_rev + n ) < m_team_size ) && ! ( m_team_rank_rev & n ) ; n <<= 1 ) {
        m_exec.pool_rev( m_team_base_rev + j )->state_set( Active );
      }
    }

public:

  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space & team_shmem() const
    { return m_team_shared ; }

  KOKKOS_INLINE_FUNCTION int league_rank() const { return m_league_rank ; }
  KOKKOS_INLINE_FUNCTION int league_size() const { return m_league_size ; }
  KOKKOS_INLINE_FUNCTION int team_rank() const { return m_team_rank ; }
  KOKKOS_INLINE_FUNCTION int team_size() const { return m_team_size ; }

  KOKKOS_INLINE_FUNCTION void team_barrier() const
#if ! defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
    {}
#else
    {
      if ( 1 < m_team_size ) {
        team_fan_in();
        team_fan_out();
      }
    }
#endif

  template<class ValueType>
  KOKKOS_INLINE_FUNCTION
  void team_broadcast(ValueType& value, const int& thread_id) const
  {
#if ! defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
    { }
#else
    // Make sure there is enough scratch space:
    typedef typename if_c< sizeof(ValueType) < TEAM_REDUCE_SIZE
                         , ValueType , void >::type type ;

    type * const local_value = ((type*) m_exec.scratch_thread());
    if(team_rank() == thread_id)
      *local_value = value;
    memory_fence();
    team_barrier();
    value = *local_value;
#endif
  }

#ifdef KOKKOS_HAVE_CXX11
  template< class ValueType, class JoinOp >
  KOKKOS_INLINE_FUNCTION ValueType
    team_reduce( const ValueType & value
               , const JoinOp & op_in ) const
  #if ! defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
    { return ValueType(); }
  #else
    {
      typedef ValueType value_type;
      const JoinLambdaAdapter<value_type,JoinOp> op(op_in);
  #endif
#else // KOKKOS_HAVE_CXX11
  template< class JoinOp >
  KOKKOS_INLINE_FUNCTION typename JoinOp::value_type
    team_reduce( const typename JoinOp::value_type & value
               , const JoinOp & op ) const
  #if ! defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
    { return typename JoinOp::value_type(); }
  #else
    {
      typedef typename JoinOp::value_type value_type;
  #endif
#endif // KOKKOS_HAVE_CXX11
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      // Make sure there is enough scratch space:
      typedef typename if_c< sizeof(value_type) < TEAM_REDUCE_SIZE
                           , value_type , void >::type type ;

      type * const local_value = ((type*) m_exec.scratch_thread());

      // Set this thread's contribution
      *local_value = value ;

      // Fence to make sure the base team member has access:
      memory_fence();

      if ( team_fan_in() ) {
        // The last thread to synchronize returns true, all other threads wait for team_fan_out()
        type * const team_value  = ((type*) m_exec.pool_rev( m_team_base_rev )->scratch_thread());

        // Join to the team value:
        for ( int i = 1 ; i < m_team_size ; ++i ) {
          op.join( *team_value , *((type*) m_exec.pool_rev( m_team_base_rev + i )->scratch_thread()) );
        }

        // The base team member may "lap" the other team members,
        // copy to their local value before proceeding.
        for ( int i = 1 ; i < m_team_size ; ++i ) {
          *((type*) m_exec.pool_rev( m_team_base_rev + i )->scratch_thread()) = *team_value ;
        }

        // Fence to make sure all team members have access
        memory_fence();
      }

      team_fan_out();

      return *((type volatile const *)local_value);
    }
#endif
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
#if ! defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
    { return ArgType(); }
#else
    {
      // Make sure there is enough scratch space:
      typedef typename if_c< sizeof(ArgType) < TEAM_REDUCE_SIZE , ArgType , void >::type type ;

      volatile type * const work_value  = ((type*) m_exec.scratch_thread());

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
            type & val = *((type*) m_exec.pool_rev( m_team_base_rev + i )->scratch_thread());
            accum += val ;
          }
          accum = atomic_fetch_add( global_accum , accum );
        }

        for ( int i = m_team_size ; i-- ; ) {
          type & val = *((type*) m_exec.pool_rev( m_team_base_rev + i )->scratch_thread());
          const type offset = accum ;  
          accum += val ;
          val = offset ;
        }

        memory_fence();
      }

      team_fan_out();

      return *work_value ;
    }
#endif

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
   *
   *  The highest rank thread can compute the reduction total as
   *    reduction_total = dev.team_scan( value ) + value ;
   */
  template< typename Type >
  KOKKOS_INLINE_FUNCTION Type team_scan( const Type & value ) const
    { return this-> template team_scan<Type>( value , 0 ); }

#ifdef KOKKOS_HAVE_CXX11

  /** \brief  Inter-thread parallel for. Executes op(iType i) for each i=0..N-1.
   *
   * The range i=0..N-1 is mapped to all threads of the the calling thread team.
   * This functionality requires C++11 support.*/
  template< typename iType, class Operation>
  KOKKOS_INLINE_FUNCTION void team_par_for(const iType n, const Operation & op) const {
    const int chunk = ((n+m_team_size-1)/m_team_size);
    const int start = chunk*m_team_rank;
    const int end = start+chunk<n?start+chunk:n;
    for(int i=start; i<end ; i++) {
      op(i);
    }
  }
#endif

  //----------------------------------------
  // Private for the driver

private:

  typedef execution_space::scratch_memory_space space ;

public:

  template< class Arg0 , class Arg1 >
  inline
  OpenMPexecTeamMember( Impl::OpenMPexec & exec
                      , const TeamPolicy< Arg0 , Arg1 , Kokkos::OpenMP > & team
                      , const int shmem_size
                      )
    : m_exec( exec )
    , m_team_shared(0,0)
    , m_team_shmem( shmem_size )
    , m_team_base_rev(0)
    , m_team_rank_rev(0)
    , m_team_rank(0)
    , m_team_size( team.team_size() )
    , m_league_rank(0)
    , m_league_end(0)
    , m_league_size( team.league_size() )
    {
      const int pool_rank_rev        = m_exec.pool_rank_rev();
      const int pool_team_rank_rev   = pool_rank_rev % team.team_alloc();
      const int pool_league_rank_rev = pool_rank_rev / team.team_alloc();
      const int league_iter_end      = team.league_size() - pool_league_rank_rev * team.team_iter();

      if ( pool_team_rank_rev < m_team_size && 0 < league_iter_end ) {
        m_team_base_rev  = team.team_alloc() * pool_league_rank_rev ;
        m_team_rank_rev  = pool_team_rank_rev ;
        m_team_rank      = m_team_size - ( m_team_rank_rev + 1 );
        m_league_end     = league_iter_end ;
        m_league_rank    = league_iter_end > team.team_iter() ? league_iter_end - team.team_iter() : 0 ;
        new( (void*) &m_team_shared ) space( ( (char*) m_exec.pool_rev(m_team_base_rev)->scratch_thread() ) + TEAM_REDUCE_SIZE , m_team_shmem );
      }
    }

  bool valid() const
    { return m_league_rank < m_league_end ; }

  void next()
    {
      if ( ++m_league_rank < m_league_end ) {
        team_barrier();
        new( (void*) &m_team_shared ) space( ( (char*) m_exec.pool_rev(m_team_base_rev)->scratch_thread() ) + TEAM_REDUCE_SIZE , m_team_shmem );
      }
    }

  static inline int team_reduce_size() { return TEAM_REDUCE_SIZE ; }
};



} // namespace Impl

template< class Arg0 , class Arg1 >
class TeamPolicy< Arg0 , Arg1 , Kokkos::OpenMP >
{
public:

  //! Tag this class as a kokkos execution policy
  typedef TeamPolicy      execution_policy ;

  //! Execution space of this execution policy.
  typedef Kokkos::OpenMP  execution_space ;

  typedef typename
    Impl::if_c< ! Impl::is_same< Kokkos::OpenMP , Arg0 >::value , Arg0 , Arg1 >::type
      work_tag ;

  //----------------------------------------

  template< class FunctorType >
  inline static
  int team_size_max( const FunctorType & )
    { return execution_space::thread_pool_size(1); }

  template< class FunctorType >
  inline static
  int team_size_recommended( const FunctorType & )
    { return execution_space::thread_pool_size(2); }

  //----------------------------------------

private:

  int m_league_size ;
  int m_team_size ;
  int m_team_alloc ;
  int m_team_iter ;

  inline void init( const int league_size_request
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

      // Maxumum number of iterations each team will take:
      m_team_iter  = ( m_league_size + team_count - 1 ) / team_count ;
    }

public:

  inline int team_size()   const { return m_team_size ; }
  inline int league_size() const { return m_league_size ; }

  /** \brief  Specify league size, request team size */
  TeamPolicy( execution_space & , int league_size_request , int team_size_request , int vector_length_request = 1)
    { init( league_size_request , team_size_request ); (void) vector_length_request; }

  TeamPolicy( int league_size_request , int team_size_request , int vector_length_request = 1 )
    { init( league_size_request , team_size_request ); (void) vector_length_request; }

  inline int team_alloc() const { return m_team_alloc ; }
  inline int team_iter()  const { return m_team_iter ; }

  typedef Impl::OpenMPexecTeamMember member_type ;
};

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

inline
int OpenMP::thread_pool_size( int depth )
{
  return Impl::OpenMPexec::pool_size(depth);
}

KOKKOS_INLINE_FUNCTION
int OpenMP::thread_pool_rank()
{
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
  return Impl::OpenMPexec::m_map_rank[ omp_get_thread_num() ];
#else
  return -1 ;
#endif
}

} // namespace Kokkos


#ifdef KOKKOS_HAVE_CXX11

namespace Kokkos {

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::TeamThreadLoopBoundariesStruct<iType,Impl::OpenMPexecTeamMember>
  TeamThreadLoop(const Impl::OpenMPexecTeamMember& thread, const iType& count) {
  return Impl::TeamThreadLoopBoundariesStruct<iType,Impl::OpenMPexecTeamMember>(thread,count);
}

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::ThreadVectorLoopBoundariesStruct<iType,Impl::OpenMPexecTeamMember >
  ThreadVectorLoop(const Impl::OpenMPexecTeamMember& thread, const iType& count) {
  return Impl::ThreadVectorLoopBoundariesStruct<iType,Impl::OpenMPexecTeamMember >(thread,count);
}

KOKKOS_INLINE_FUNCTION
Impl::ThreadSingleStruct<Impl::OpenMPexecTeamMember> PerTeam(const Impl::OpenMPexecTeamMember& thread) {
  return Impl::ThreadSingleStruct<Impl::OpenMPexecTeamMember>(thread);
}

KOKKOS_INLINE_FUNCTION
Impl::VectorSingleStruct<Impl::OpenMPexecTeamMember> PerThread(const Impl::OpenMPexecTeamMember& thread) {
  return Impl::VectorSingleStruct<Impl::OpenMPexecTeamMember>(thread);
}
} // namespace Kokkos

namespace Kokkos {

  /** \brief  Inter-thread parallel_for. Executes lambda(iType i) for each i=0..N-1.
   *
   * The range i=0..N-1 is mapped to all threads of the the calling thread team.
   * This functionality requires C++11 support.*/
template<typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION
void parallel_for(const Impl::TeamThreadLoopBoundariesStruct<iType,Impl::OpenMPexecTeamMember>& loop_boundaries, const Lambda& lambda) {
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment)
    lambda(i);
}

/** \brief  Inter-thread vector parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team and a summation of
 * val is performed and put into result. This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ValueType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::TeamThreadLoopBoundariesStruct<iType,Impl::OpenMPexecTeamMember>& loop_boundaries,
                     const Lambda & lambda, ValueType& result) {

  result = ValueType();

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    ValueType tmp = ValueType();
    lambda(i,tmp);
    result+=tmp;
  }

  result = loop_boundaries.thread.team_reduce(result,Impl::JoinAdd<ValueType>());
}

/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread and a reduction of
 * val is performed using JoinType(ValueType& val, const ValueType& update) and put into init_result.
 * The input value of init_result is used as initializer for temporary variables of ValueType. Therefore
 * the input value should be the neutral element with respect to the join operation (e.g. '0 for +-' or
 * '1 for *'). This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ValueType, class JoinType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::TeamThreadLoopBoundariesStruct<iType,Impl::OpenMPexecTeamMember>& loop_boundaries,
                     const Lambda & lambda, const JoinType& join, ValueType& init_result) {

  ValueType result = init_result;

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    ValueType tmp = ValueType();
    lambda(i,tmp);
    join(result,tmp);
  }

  init_result = loop_boundaries.thread.team_reduce(result,join);
}

} //namespace Kokkos


namespace Kokkos {
/** \brief  Intra-thread vector parallel_for. Executes lambda(iType i) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread.
 * This functionality requires C++11 support.*/
template<typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION
void parallel_for(const Impl::ThreadVectorLoopBoundariesStruct<iType,Impl::OpenMPexecTeamMember >&
    loop_boundaries, const Lambda& lambda) {
  #ifdef KOKKOS_HAVE_PRAGMA_IVDEP
  #pragma ivdep
  #endif
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment)
    lambda(i);
}

/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread and a summation of
 * val is performed and put into result. This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ValueType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::ThreadVectorLoopBoundariesStruct<iType,Impl::OpenMPexecTeamMember >&
      loop_boundaries, const Lambda & lambda, ValueType& result) {
  result = ValueType();
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    ValueType tmp = ValueType();
    lambda(i,tmp);
    result+=tmp;
  }
}

/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread and a reduction of
 * val is performed using JoinType(ValueType& val, const ValueType& update) and put into init_result.
 * The input value of init_result is used as initializer for temporary variables of ValueType. Therefore
 * the input value should be the neutral element with respect to the join operation (e.g. '0 for +-' or
 * '1 for *'). This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ValueType, class JoinType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::ThreadVectorLoopBoundariesStruct<iType,Impl::OpenMPexecTeamMember >&
      loop_boundaries, const Lambda & lambda, const JoinType& join, ValueType& init_result) {

  ValueType result = init_result;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    ValueType tmp = ValueType();
    lambda(i,tmp);
    join(result,tmp);
  }
  init_result = result;
}

/** \brief  Intra-thread vector parallel exclusive prefix sum. Executes lambda(iType i, ValueType & val, bool final)
 *          for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes in the thread and a scan operation is performed.
 * Depending on the target execution space the operator might be called twice: once with final=false
 * and once with final=true. When final==true val contains the prefix sum value. The contribution of this
 * "i" needs to be added to val no matter whether final==true or not. In a serial execution
 * (i.e. team_size==1) the operator is only called once with final==true. Scan_val will be set
 * to the final sum value over all vector lanes.
 * This functionality requires C++11 support.*/
template< typename iType, class FunctorType >
KOKKOS_INLINE_FUNCTION
void parallel_scan(const Impl::ThreadVectorLoopBoundariesStruct<iType,Impl::OpenMPexecTeamMember >&
      loop_boundaries, const FunctorType & lambda) {

  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , void > ValueTraits ;
  typedef typename ValueTraits::value_type value_type ;

  value_type scan_val = value_type();

#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    lambda(i,scan_val,true);
  }
}

} // namespace Kokkos

namespace Kokkos {

template<class FunctorType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::VectorSingleStruct<Impl::OpenMPexecTeamMember>& single_struct, const FunctorType& lambda) {
  lambda();
}

template<class FunctorType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::ThreadSingleStruct<Impl::OpenMPexecTeamMember>& single_struct, const FunctorType& lambda) {
  if(single_struct.team_member.team_rank()==0) lambda();
}

template<class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::VectorSingleStruct<Impl::OpenMPexecTeamMember>& single_struct, const FunctorType& lambda, ValueType& val) {
  lambda(val);
}

template<class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::ThreadSingleStruct<Impl::OpenMPexecTeamMember>& single_struct, const FunctorType& lambda, ValueType& val) {
  if(single_struct.team_member.team_rank()==0) {
    lambda(val);
  }
  single_struct.team_member.team_broadcast(val,0);
}
}

#endif // KOKKOS_HAVE_CXX11

#endif /* #ifndef KOKKOS_OPENMPEXEC_HPP */

