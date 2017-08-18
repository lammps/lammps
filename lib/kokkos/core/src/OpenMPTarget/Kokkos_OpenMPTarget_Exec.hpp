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

#ifndef KOKKOS_OPENMPTARGETEXEC_HPP
#define KOKKOS_OPENMPTARGETEXEC_HPP

#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Spinwait.hpp>

#include <Kokkos_Atomic.hpp>
#include <iostream>
#include <sstream>
#include <fstream>
namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
/** \brief  Data for OpenMPTarget thread execution */


class OpenMPTargetExec {
public:
  enum { MAX_ACTIVE_THREADS = 256*8*56*4 };
  enum { MAX_ACTIVE_TEAMS = MAX_ACTIVE_THREADS/32 };

private:
  static void* scratch_ptr;

public:
  static void verify_is_process( const char * const );
  static void verify_initialized( const char * const );

  static void* get_scratch_ptr();
  static void clear_scratch();
  static void resize_scratch( int64_t reduce_bytes , int64_t team_reduce_bytes, int64_t team_shared_bytes, int64_t thread_local_bytes );

  static void* m_scratch_ptr;
  static int64_t m_scratch_size;
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

class OpenMPTargetExecTeamMember {
public:

  enum { TEAM_REDUCE_SIZE = 512 };

  /** \brief  Thread states for team synchronization */
  enum { Active = 0 , Rendezvous = 1 };

  typedef Kokkos::Experimental::OpenMPTarget                         execution_space ;
  typedef execution_space::scratch_memory_space  scratch_memory_space ;

  scratch_memory_space  m_team_shared ;
  int                   m_team_scratch_size[2] ;
  int                   m_team_rank ;
  int                   m_team_size ;
  int                   m_league_rank ;
  int                   m_league_size ;
  int                   m_vector_length ;
  int                   m_vector_lane ;
  void* 		m_glb_scratch ;

  /*
  // Fan-in team threads, root of the fan-in which does not block returns true
  inline
  bool team_fan_in() const
    {
      memory_fence();
      for ( int n = 1 , j ; ( ( j = m_team_rank_rev + n ) < m_team_size ) && ! ( m_team_rank_rev & n ) ; n <<= 1 ) {

        m_exec.pool_rev( m_team_base_rev + j )->state_wait( Active );
      }

      if ( m_team_rank_rev ) {
        m_exec.state_set( Rendezvous );
        memory_fence();
        m_exec.state_wait( Rendezvous );
      }

      return 0 == m_team_rank_rev ;
    }

  inline
  void team_fan_out() const
    {
      memory_fence();
      for ( int n = 1 , j ; ( ( j = m_team_rank_rev + n ) < m_team_size ) && ! ( m_team_rank_rev & n ) ; n <<= 1 ) {
        m_exec.pool_rev( m_team_base_rev + j )->state_set( Active );
        memory_fence();
      }
    }
  */
public:

  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space& team_shmem() const
    { return m_team_shared.set_team_thread_mode(0,1,0) ; }

  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space& team_scratch(int) const
    { return m_team_shared.set_team_thread_mode(0,1,0) ; }

  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space& thread_scratch(int) const
    { return m_team_shared.set_team_thread_mode(0,team_size(),team_rank()) ; }

  KOKKOS_INLINE_FUNCTION int league_rank() const { return m_league_rank ; }
  KOKKOS_INLINE_FUNCTION int league_size() const { return m_league_size ; }
  KOKKOS_INLINE_FUNCTION int team_rank() const { return m_team_rank ; }
  KOKKOS_INLINE_FUNCTION int team_size() const { return m_team_size ; }

  KOKKOS_INLINE_FUNCTION void team_barrier() const
    {
      #pragma omp barrier
    }

  template<class ValueType>
  KOKKOS_INLINE_FUNCTION
  void team_broadcast(ValueType& value, const int& thread_id) const
  {
/*#if ! defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
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
#endif*/
  }

  template< class ValueType, class JoinOp >
  KOKKOS_INLINE_FUNCTION ValueType
    team_reduce( const ValueType & value
               , const JoinOp & op_in ) const {

      #pragma omp barrier

      typedef ValueType value_type;
      const JoinLambdaAdapter<value_type,JoinOp> op(op_in);

      // Make sure there is enough scratch space:
      typedef typename if_c< sizeof(value_type) < TEAM_REDUCE_SIZE
                           , value_type , void >::type type ;

      const int n_values = TEAM_REDUCE_SIZE/sizeof(value_type);
      type * team_scratch = (type*) ((char*)m_glb_scratch + TEAM_REDUCE_SIZE*omp_get_team_num());
      for(int i = m_team_rank; i < n_values; i+= m_team_size) {
        team_scratch[i] = value_type();
      }

      #pragma omp barrier

      for(int k=0; k<m_team_size; k+=n_values) {
        if((k <= m_team_rank) && (k+n_values > m_team_rank))
          team_scratch[m_team_rank%n_values]+=value;
        #pragma omp barrier
      }

      for(int d = 1; d<n_values;d*=2) {
        if((m_team_rank+d<n_values) && (m_team_rank%(2*d)==0)) {
          team_scratch[m_team_rank] += team_scratch[m_team_rank+d];
        }
        #pragma omp barrier
      }
      return team_scratch[0];
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
    /*  // Make sure there is enough scratch space:
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

      return *work_value ;*/
      return ArgType();
    }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
   *
   *  The highest rank thread can compute the reduction total as
   *    reduction_total = dev.team_scan( value ) + value ;
   */
  template< typename Type >
  KOKKOS_INLINE_FUNCTION Type team_scan( const Type & value ) const
    { return this-> template team_scan<Type>( value , 0 ); }

  //----------------------------------------
  // Private for the driver

private:

  typedef execution_space::scratch_memory_space space ;

public:

  inline
  OpenMPTargetExecTeamMember( const int league_rank, const int league_size, const int team_size, const int vector_length //const TeamPolicyInternal< OpenMPTarget, Properties ...> & team
                      , void* const glb_scratch
                      , const int shmem_size_L1
                      , const int shmem_size_L2
                      )
    : m_team_shared(0,0)
    , m_team_scratch_size{ shmem_size_L1 , shmem_size_L2 }
    , m_team_rank(0)
    , m_vector_length( vector_length )
    , m_team_size( team_size )
    , m_league_rank( league_rank )
    , m_league_size( league_size )
    , m_glb_scratch( glb_scratch )
    {
      const int omp_tid = omp_get_thread_num();
      m_league_rank = league_rank;
      m_team_rank = omp_tid/m_vector_length;
      m_vector_lane = omp_tid%m_vector_length;
    }

  static inline int team_reduce_size() { return TEAM_REDUCE_SIZE ; }
};



template< class ... Properties >
class TeamPolicyInternal< Kokkos::Experimental::OpenMPTarget, Properties ... >: public PolicyTraits<Properties ...>
{
public:

  //! Tag this class as a kokkos execution policy
  typedef TeamPolicyInternal      execution_policy ;

  typedef PolicyTraits<Properties ... > traits;

  TeamPolicyInternal& operator = (const TeamPolicyInternal& p) {
    m_league_size = p.m_league_size;
    m_team_size = p.m_team_size;
    m_vector_length = p.m_vector_length;
    m_team_alloc = p.m_team_alloc;
    m_team_iter = p.m_team_iter;
    m_team_scratch_size[0] = p.m_team_scratch_size[0];
    m_thread_scratch_size[0] = p.m_thread_scratch_size[0];
    m_team_scratch_size[1] = p.m_team_scratch_size[1];
    m_thread_scratch_size[1] = p.m_thread_scratch_size[1];
    m_chunk_size = p.m_chunk_size;
    return *this;
  }

  //----------------------------------------

  template< class FunctorType >
  inline static
  int team_size_max( const FunctorType & )
    { return 1024; }

  template< class FunctorType >
  inline static
  int team_size_recommended( const FunctorType & )
    { return 256; }

  template< class FunctorType >
  inline static
  int team_size_recommended( const FunctorType &, const int& vector_length)
    { return 256/vector_length; }

  //----------------------------------------

private:

  int m_league_size ;
  int m_team_size ;
  int m_vector_length;
  int m_team_alloc ;
  int m_team_iter ;

  size_t m_team_scratch_size[2];
  size_t m_thread_scratch_size[2];

  int m_chunk_size;

  inline void init( const int league_size_request
                  , const int team_size_request
                  , const int vector_length_request )
    {
      m_league_size = league_size_request ;

      m_team_size = team_size_request;

      m_vector_length = vector_length_request;

      set_auto_chunk_size();
    }

public:

  inline int vector_length() const { return m_vector_length ; }
  inline int team_size()   const { return m_team_size ; }
  inline int league_size() const { return m_league_size ; }
  inline size_t scratch_size(const int& level, int team_size_ = -1) const {
    if(team_size_ < 0)
      team_size_ = m_team_size;
    return m_team_scratch_size[level] + team_size_*m_thread_scratch_size[level] ;
  }

  /** \brief  Specify league size, request team size */
  TeamPolicyInternal( typename traits::execution_space &
            , int league_size_request
            , int team_size_request
            , int vector_length_request = 1 )
            : m_team_scratch_size { 0 , 0 }
            , m_thread_scratch_size { 0 , 0 }
            , m_chunk_size(0)
    { init( league_size_request , team_size_request , vector_length_request); }

  TeamPolicyInternal( typename traits::execution_space &
            , int league_size_request
            , const Kokkos::AUTO_t & /* team_size_request */
            , int vector_length_request = 1)
            : m_team_scratch_size { 0 , 0 }
            , m_thread_scratch_size { 0 , 0 }
            , m_chunk_size(0)
    { init( league_size_request , 256/vector_length_request , vector_length_request ); }

  TeamPolicyInternal( int league_size_request
            , int team_size_request
            , int vector_length_request = 1 )
            : m_team_scratch_size { 0 , 0 }
            , m_thread_scratch_size { 0 , 0 }
            , m_chunk_size(0)
    { init( league_size_request , team_size_request , vector_length_request); }

  TeamPolicyInternal( int league_size_request
            , const Kokkos::AUTO_t & /* team_size_request */
            , int vector_length_request = 1 )
            : m_team_scratch_size { 0 , 0 }
            , m_thread_scratch_size { 0 , 0 }
            , m_chunk_size(0)
    { init( league_size_request , 256/vector_length_request , vector_length_request ); }

  inline int team_alloc() const { return m_team_alloc ; }
  inline int team_iter()  const { return m_team_iter ; }

  inline int chunk_size() const { return m_chunk_size ; }

  /** \brief set chunk_size to a discrete value*/
  inline TeamPolicyInternal set_chunk_size(typename traits::index_type chunk_size_) const {
    TeamPolicyInternal p = *this;
    p.m_chunk_size = chunk_size_;
    return p;
  }

  inline TeamPolicyInternal set_scratch_size(const int& level, const PerTeamValue& per_team) const {
    TeamPolicyInternal p = *this;
    p.m_team_scratch_size[level] = per_team.value;
    return p;
  };

  inline TeamPolicyInternal set_scratch_size(const int& level, const PerThreadValue& per_thread) const {
    TeamPolicyInternal p = *this;
    p.m_thread_scratch_size[level] = per_thread.value;
    return p;
  };

  inline TeamPolicyInternal set_scratch_size(const int& level, const PerTeamValue& per_team, const PerThreadValue& per_thread) const {
    TeamPolicyInternal p = *this;
    p.m_team_scratch_size[level] = per_team.value;
    p.m_thread_scratch_size[level] = per_thread.value;
    return p;
  };

private:
  /** \brief finalize chunk_size if it was set to AUTO*/
  inline void set_auto_chunk_size() {

    int concurrency = traits::execution_space::thread_pool_size(0)/m_team_alloc;
    if( concurrency==0 ) concurrency=1;

    if(m_chunk_size > 0) {
      if(!Impl::is_integral_power_of_two( m_chunk_size ))
        Kokkos::abort("TeamPolicy blocking granularity must be power of two" );
    }

    int new_chunk_size = 1;
    while(new_chunk_size*100*concurrency < m_league_size)
      new_chunk_size *= 2;
    if(new_chunk_size < 128) {
      new_chunk_size = 1;
      while( (new_chunk_size*40*concurrency < m_league_size ) && (new_chunk_size<128) )
        new_chunk_size*=2;
    }
    m_chunk_size = new_chunk_size;
  }

public:
  typedef Impl::OpenMPTargetExecTeamMember member_type ;
};
} // namespace Impl


} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {

inline
int OpenMPTarget::thread_pool_size( int depth )
{
  //return Impl::OpenMPTargetExec::pool_size(depth);
  return omp_get_max_threads();
}

KOKKOS_INLINE_FUNCTION
int OpenMPTarget::thread_pool_rank()
{
  return omp_get_thread_num();
}

} // namespace Experimental
} // namespace Kokkos


namespace Kokkos {

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::TeamThreadRangeBoundariesStruct<iType,Impl::OpenMPTargetExecTeamMember>
  TeamThreadRange(const Impl::OpenMPTargetExecTeamMember& thread, const iType& count) {
  return Impl::TeamThreadRangeBoundariesStruct<iType,Impl::OpenMPTargetExecTeamMember>(thread,count);
}

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::TeamThreadRangeBoundariesStruct<iType,Impl::OpenMPTargetExecTeamMember>
  TeamThreadRange(const Impl::OpenMPTargetExecTeamMember& thread, const iType& begin, const iType& end) {
  return Impl::TeamThreadRangeBoundariesStruct<iType,Impl::OpenMPTargetExecTeamMember>(thread,begin,end);
}

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::OpenMPTargetExecTeamMember >
  ThreadVectorRange(const Impl::OpenMPTargetExecTeamMember& thread, const iType& count) {
  return Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::OpenMPTargetExecTeamMember >(thread,count);
}

KOKKOS_INLINE_FUNCTION
Impl::ThreadSingleStruct<Impl::OpenMPTargetExecTeamMember> PerTeam(const Impl::OpenMPTargetExecTeamMember& thread) {
  return Impl::ThreadSingleStruct<Impl::OpenMPTargetExecTeamMember>(thread);
}

KOKKOS_INLINE_FUNCTION
Impl::VectorSingleStruct<Impl::OpenMPTargetExecTeamMember> PerThread(const Impl::OpenMPTargetExecTeamMember& thread) {
  return Impl::VectorSingleStruct<Impl::OpenMPTargetExecTeamMember>(thread);
}
} // namespace Kokkos

namespace Kokkos {

  /** \brief  Inter-thread parallel_for. Executes lambda(iType i) for each i=0..N-1.
   *
   * The range i=0..N-1 is mapped to all threads of the the calling thread team.
   * This functionality requires C++11 support.*/
template<typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION
void parallel_for(const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::OpenMPTargetExecTeamMember>& loop_boundaries, const Lambda& lambda) {
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment)
    lambda(i);
}

/** \brief  Inter-thread vector parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team and a summation of
 * val is performed and put into result. This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ValueType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
                     const Lambda & lambda, ValueType& result) {

  result = ValueType();

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    ValueType tmp = ValueType();
    lambda(i,tmp);
    result+=tmp;
  }

  //result = loop_boundaries.thread.team_reduce(result,Impl::JoinAdd<ValueType>());
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
void parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
                     const Lambda & lambda, const JoinType& join, ValueType& init_result) {

  ValueType result = init_result;

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    ValueType tmp = ValueType();
    lambda(i,tmp);
    join(result,tmp);
  }

  //init_result = loop_boundaries.thread.team_reduce(result,join);
}

} //namespace Kokkos


namespace Kokkos {
/** \brief  Intra-thread vector parallel_for. Executes lambda(iType i) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread.
 * This functionality requires C++11 support.*/
template<typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION
void parallel_for(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::OpenMPTargetExecTeamMember >&
    loop_boundaries, const Lambda& lambda) {
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment)
    lambda(i);
}

/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread and a summation of
 * val is performed and put into result. This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ValueType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::OpenMPTargetExecTeamMember >&
      loop_boundaries, const Lambda & lambda, ValueType& result) {
  result = ValueType();
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
void parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::OpenMPTargetExecTeamMember >&
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
void parallel_scan(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::OpenMPTargetExecTeamMember >&
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
void single(const Impl::VectorSingleStruct<Impl::OpenMPTargetExecTeamMember>& single_struct, const FunctorType& lambda) {
  lambda();
}

template<class FunctorType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::ThreadSingleStruct<Impl::OpenMPTargetExecTeamMember>& single_struct, const FunctorType& lambda) {
  if(single_struct.team_member.team_rank()==0) lambda();
}

template<class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::VectorSingleStruct<Impl::OpenMPTargetExecTeamMember>& single_struct, const FunctorType& lambda, ValueType& val) {
  lambda(val);
}

template<class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::ThreadSingleStruct<Impl::OpenMPTargetExecTeamMember>& single_struct, const FunctorType& lambda, ValueType& val) {
  if(single_struct.team_member.team_rank()==0) {
    lambda(val);
  }
  single_struct.team_member.team_broadcast(val,0);
}
}

#endif /* #ifndef KOKKOS_OPENMPTARGETEXEC_HPP */

