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

#ifndef KOKKOS_THREADSTEAM_HPP
#define KOKKOS_THREADSTEAM_HPP

#include <stdio.h>

#include <utility>
#include <impl/Kokkos_spinwait.hpp>
#include <impl/Kokkos_FunctorAdapter.hpp>

#include <Kokkos_Atomic.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template< class > struct ThreadsExecAdapter ;

//----------------------------------------------------------------------------

class ThreadsExecTeamMember {
private:

  enum { TEAM_REDUCE_SIZE = 512 };

  typedef Kokkos::Threads execution_space ;
  typedef execution_space::scratch_memory_space space ;

  ThreadsExec * const   m_exec ;
  ThreadsExec * const * m_team_base ; ///< Base for team fan-in
  space                 m_team_shared ;
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
  
public:

  // Fan-in and wait until the matching fan-out is called.
  // The root thread which does not wait will return true.
  // All other threads will return false during the fan-out.
  KOKKOS_INLINE_FUNCTION bool team_fan_in() const
    {
      int n , j ;

      // Wait for fan-in threads
      for ( n = 1 ; ( ! ( m_team_rank_rev & n ) ) && ( ( j = m_team_rank_rev + n ) < m_team_size ) ; n <<= 1 ) {
        Impl::spinwait( m_team_base[j]->state() , ThreadsExec::Active );
      }

      // If not root then wait for release
      if ( m_team_rank_rev ) {
        m_exec->state() = ThreadsExec::Rendezvous ;
        Impl::spinwait( m_exec->state() , ThreadsExec::Rendezvous );
      }

      return ! m_team_rank_rev ;
    }

  KOKKOS_INLINE_FUNCTION void team_fan_out() const
    {
      int n , j ;
      for ( n = 1 ; ( ! ( m_team_rank_rev & n ) ) && ( ( j = m_team_rank_rev + n ) < m_team_size ) ; n <<= 1 ) {
        m_team_base[j]->state() = ThreadsExec::Active ;
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

    if ( m_team_base ) {
      type * const local_value = ((type*) m_team_base[0]->scratch_memory());
      if(team_rank() == thread_id) *local_value = value;
      memory_fence();
      team_barrier();
      value = *local_value;
    }
#endif
  }

  template< typename Type >
  KOKKOS_INLINE_FUNCTION Type team_reduce( const Type & value ) const
#if ! defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
    { return Type(); }
#else
    {
      // Make sure there is enough scratch space:
      typedef typename if_c< sizeof(Type) < TEAM_REDUCE_SIZE , Type , void >::type type ;

      if ( 0 == m_exec ) return value ;

      *((volatile type*) m_exec->scratch_memory() ) = value ;

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
#endif

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

      if ( 0 == m_exec ) return value ;

      type * const local_value = ((type*) m_exec->scratch_memory());

      // Set this thread's contribution
      *local_value = value ;

      // Fence to make sure the base team member has access:
      memory_fence();

      if ( team_fan_in() ) {
        // The last thread to synchronize returns true, all other threads wait for team_fan_out()
        type * const team_value = ((type*) m_team_base[0]->scratch_memory());

        // Join to the team value:
        for ( int i = 1 ; i < m_team_size ; ++i ) {
          op.join( *team_value , *((type*) m_team_base[i]->scratch_memory()) );
        }

        // Team base thread may "lap" member threads so copy out to their local value.
        for ( int i = 1 ; i < m_team_size ; ++i ) {
          *((type*) m_team_base[i]->scratch_memory()) = *team_value ;
        }

        // Fence to make sure all team members have access
        memory_fence();
      }

      team_fan_out();

      // Value was changed by the team base
      return *((type volatile const *) local_value);
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

      if ( 0 == m_exec ) return type(0);

      volatile type * const work_value  = ((type*) m_exec->scratch_memory());

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
#endif

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

  template< class Arg0 , class Arg1 >
  ThreadsExecTeamMember( Impl::ThreadsExec * exec
                       , const TeamPolicy< Arg0 , Arg1 , Kokkos::Threads > & team 
                       , const int shared_size )
    : m_exec( exec )
    , m_team_base(0)
    , m_team_shared(0,0)
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

        const int pool_rank_rev = m_exec->pool_size() - ( m_exec->pool_rank() + 1 );
        const int team_rank_rev = pool_rank_rev % team.team_alloc();

        // May be using fewer threads per team than a multiple of threads per core,
        // some threads will idle.

        if ( team_rank_rev < team.team_size() ) {
          const size_t pool_league_size     = m_exec->pool_size() / team.team_alloc() ;
          const size_t pool_league_rank_rev = pool_rank_rev / team.team_alloc() ;
          const size_t pool_league_rank     = pool_league_size - ( pool_league_rank_rev + 1 );

          m_team_base        = m_exec->pool_base() + team.team_alloc() * pool_league_rank_rev ;
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

  ThreadsExecTeamMember()
    : m_exec(0)
    , m_team_base(0)
    , m_team_shared(0,0)
    , m_team_shared_size(0)
    , m_team_size(1)
    , m_team_rank(0)
    , m_team_rank_rev(0)
    , m_league_size(1)
    , m_league_end(0)
    , m_league_rank(0)
    {}

  inline
  ThreadsExec & threads_exec_team_base() const { return m_team_base ? **m_team_base : *m_exec ; }

  bool valid() const
    { return m_league_rank < m_league_end ; }

  void next()
    {
      if ( ++m_league_rank < m_league_end ) {
        team_barrier();
        set_team_shared();
      }
    }

  void set_league_shmem( const int arg_league_rank
                       , const int arg_league_size
                       , const int arg_shmem_size
                       )
    {
      m_league_rank = arg_league_rank ;
      m_league_size = arg_league_size ;
      m_team_shared_size = arg_shmem_size ;
      set_team_shared();
    }
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template< class Arg0 , class Arg1 >
class TeamPolicy< Arg0 , Arg1 , Kokkos::Threads >
{
private:

  int m_league_size ;
  int m_team_size ;
  int m_team_alloc ;

  size_t m_scratch_size;

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

  //! Tag this class as a kokkos execution policy
  typedef TeamPolicy       execution_policy ; 
  typedef Kokkos::Threads  execution_space ;

  typedef typename
    Impl::if_c< ! Impl::is_same< Kokkos::Threads , Arg0 >::value , Arg0 , Arg1 >::type
      work_tag ;

  //----------------------------------------

  template< class FunctorType >
  inline static
  int team_size_max( const FunctorType & )
    { return execution_space::thread_pool_size(1); }

  template< class FunctorType >
  static int team_size_recommended( const FunctorType & )
    { return execution_space::thread_pool_size(2); }


  template< class FunctorType >
  inline static
  int team_size_recommended( const FunctorType &, const int& )
    { return execution_space::thread_pool_size(2); }

  //----------------------------------------

  inline int team_size() const { return m_team_size ; }
  inline int team_alloc() const { return m_team_alloc ; }
  inline int league_size() const { return m_league_size ; }
  inline size_t scratch_size() const { return m_scratch_size ; }

  /** \brief  Specify league size, request team size */
  TeamPolicy( execution_space &
            , int league_size_request
            , int team_size_request
            , int vector_length_request = 1 )
    : m_league_size(0)
    , m_team_size(0)
    , m_team_alloc(0)
    , m_scratch_size ( 0 )
    { init(league_size_request,team_size_request); (void) vector_length_request; }

  /** \brief  Specify league size, request team size */
  TeamPolicy( execution_space &
            , int league_size_request
            , const Kokkos::AUTO_t & /* team_size_request */
            , int /* vector_length_request */ = 1 )
    : m_league_size(0)
    , m_team_size(0)
    , m_team_alloc(0)
    , m_scratch_size ( 0 )
    { init(league_size_request,execution_space::thread_pool_size(2)); }

  TeamPolicy( int league_size_request
            , int team_size_request
            , int /* vector_length_request */ = 1 )
    : m_league_size(0)
    , m_team_size(0)
    , m_team_alloc(0)
    , m_scratch_size ( 0 )
    { init(league_size_request,team_size_request); }

  TeamPolicy( int league_size_request
            , const Kokkos::AUTO_t & /* team_size_request */
            , int /* vector_length_request */ = 1 )
    : m_league_size(0)
    , m_team_size(0)
    , m_team_alloc(0)
    , m_scratch_size ( 0 )
    { init(league_size_request,execution_space::thread_pool_size(2)); }

  template<class MemorySpace>
  TeamPolicy( int league_size_request
            , int team_size_request
            , const Experimental::TeamScratchRequest<MemorySpace> & scratch_request )
    : m_league_size(0)
    , m_team_size(0)
    , m_team_alloc(0)
    , m_scratch_size(scratch_request.total(team_size_request))
    { init(league_size_request,team_size_request); }


  template<class MemorySpace>
  TeamPolicy( int league_size_request
            , const Kokkos::AUTO_t & /* team_size_request */
            , const Experimental::TeamScratchRequest<MemorySpace> & scratch_request )
    : m_league_size(0)
    , m_team_size(0)
    , m_team_alloc(0)
    , m_scratch_size(scratch_request.total(execution_space::thread_pool_size(2)))
    { init(league_size_request,execution_space::thread_pool_size(2)); }

  typedef Impl::ThreadsExecTeamMember member_type ;

  friend class Impl::ThreadsExecTeamMember ;
};


} /* namespace Kokkos */


namespace Kokkos {

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::TeamThreadRangeBoundariesStruct<iType,Impl::ThreadsExecTeamMember>
TeamThreadRange(const Impl::ThreadsExecTeamMember& thread, const iType& count)
{
  return Impl::TeamThreadRangeBoundariesStruct<iType,Impl::ThreadsExecTeamMember>(thread,count);
}

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::TeamThreadRangeBoundariesStruct<iType,Impl::ThreadsExecTeamMember>
TeamThreadRange( const Impl::ThreadsExecTeamMember& thread
               , const iType & begin
               , const iType & end
               )
{
  return Impl::TeamThreadRangeBoundariesStruct<iType,Impl::ThreadsExecTeamMember>(thread,begin,end);
}


template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::ThreadsExecTeamMember >
  ThreadVectorRange(const Impl::ThreadsExecTeamMember& thread, const iType& count) {
  return Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::ThreadsExecTeamMember >(thread,count);
}


KOKKOS_INLINE_FUNCTION
Impl::ThreadSingleStruct<Impl::ThreadsExecTeamMember> PerTeam(const Impl::ThreadsExecTeamMember& thread) {
  return Impl::ThreadSingleStruct<Impl::ThreadsExecTeamMember>(thread);
}

KOKKOS_INLINE_FUNCTION
Impl::VectorSingleStruct<Impl::ThreadsExecTeamMember> PerThread(const Impl::ThreadsExecTeamMember& thread) {
  return Impl::VectorSingleStruct<Impl::ThreadsExecTeamMember>(thread);
}
} // namespace Kokkos

namespace Kokkos {

  /** \brief  Inter-thread parallel_for. Executes lambda(iType i) for each i=0..N-1.
   *
   * The range i=0..N-1 is mapped to all threads of the the calling thread team.
   * This functionality requires C++11 support.*/
template<typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION
void parallel_for(const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::ThreadsExecTeamMember>& loop_boundaries, const Lambda& lambda) {
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment)
    lambda(i);
}

/** \brief  Inter-thread vector parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team and a summation of
 * val is performed and put into result. This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ValueType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::ThreadsExecTeamMember>& loop_boundaries,
                     const Lambda & lambda, ValueType& result) {

  result = ValueType();

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    ValueType tmp = ValueType();
    lambda(i,tmp);
    result+=tmp;
  }

  result = loop_boundaries.thread.team_reduce(result,Impl::JoinAdd<ValueType>());
}

#if defined( KOKKOS_HAVE_CXX11 )

/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread and a reduction of
 * val is performed using JoinType(ValueType& val, const ValueType& update) and put into init_result.
 * The input value of init_result is used as initializer for temporary variables of ValueType. Therefore
 * the input value should be the neutral element with respect to the join operation (e.g. '0 for +-' or
 * '1 for *'). This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ValueType, class JoinType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::ThreadsExecTeamMember>& loop_boundaries,
                     const Lambda & lambda, const JoinType& join, ValueType& init_result) {

  ValueType result = init_result;

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    ValueType tmp = ValueType();
    lambda(i,tmp);
    join(result,tmp);
  }

  init_result = loop_boundaries.thread.team_reduce(result,Impl::JoinLambdaAdapter<ValueType,JoinType>(join));
}

#endif /* #if defined( KOKKOS_HAVE_CXX11 ) */

} //namespace Kokkos


namespace Kokkos {
/** \brief  Intra-thread vector parallel_for. Executes lambda(iType i) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread.
 * This functionality requires C++11 support.*/
template<typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION
void parallel_for(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::ThreadsExecTeamMember >&
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
void parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::ThreadsExecTeamMember >&
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
void parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::ThreadsExecTeamMember >&
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
void parallel_scan(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::ThreadsExecTeamMember >&
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
void single(const Impl::VectorSingleStruct<Impl::ThreadsExecTeamMember>& single_struct, const FunctorType& lambda) {
  lambda();
}

template<class FunctorType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::ThreadSingleStruct<Impl::ThreadsExecTeamMember>& single_struct, const FunctorType& lambda) {
  if(single_struct.team_member.team_rank()==0) lambda();
}

template<class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::VectorSingleStruct<Impl::ThreadsExecTeamMember>& single_struct, const FunctorType& lambda, ValueType& val) {
  lambda(val);
}

template<class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::ThreadSingleStruct<Impl::ThreadsExecTeamMember>& single_struct, const FunctorType& lambda, ValueType& val) {
  if(single_struct.team_member.team_rank()==0) {
    lambda(val);
  }
  single_struct.team_member.team_broadcast(val,0);
}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #define KOKKOS_THREADSTEAM_HPP */

