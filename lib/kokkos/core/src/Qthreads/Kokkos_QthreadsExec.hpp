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

#ifndef KOKKOS_QTHREADSEXEC_HPP
#define KOKKOS_QTHREADSEXEC_HPP

#include <Kokkos_Macros.hpp>
#if defined( KOKKOS_ENABLE_QTHREADS )

#include <impl/Kokkos_Spinwait.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {

namespace Impl {

class QthreadsExec;

typedef void (*QthreadsExecFunctionPointer)( QthreadsExec &, const void * );

class QthreadsExec {
private:
  enum { Inactive = 0, Active = 1 };

  const QthreadsExec * const * m_worker_base;
  const QthreadsExec * const * m_shepherd_base;

  void  * m_scratch_alloc;  ///< Scratch memory [ reduce, team, shared ]
  int     m_reduce_end;     ///< End of scratch reduction memory

  int     m_shepherd_rank;
  int     m_shepherd_size;

  int     m_shepherd_worker_rank;
  int     m_shepherd_worker_size;

  /*
   *  m_worker_rank = m_shepherd_rank * m_shepherd_worker_size + m_shepherd_worker_rank
   *  m_worker_size = m_shepherd_size * m_shepherd_worker_size
   */
  int     m_worker_rank;
  int     m_worker_size;

  int mutable volatile m_worker_state;

  friend class Kokkos::Qthreads;

  ~QthreadsExec();
  QthreadsExec( const QthreadsExec & );
  QthreadsExec & operator = ( const QthreadsExec & );

public:
  QthreadsExec();

  /** Execute the input function on all available Qthreads workers. */
  static void exec_all( Qthreads &, QthreadsExecFunctionPointer, const void * );

  /** Barrier across all workers participating in the 'exec_all'. */
  void exec_all_barrier() const
  {
    const int rev_rank = m_worker_size - ( m_worker_rank + 1 );

    int n, j;

    for ( n = 1; ( ! ( rev_rank & n ) ) && ( ( j = rev_rank + n ) < m_worker_size ); n <<= 1 ) {
      Impl::spinwait_while_equal( m_worker_base[j]->m_worker_state, QthreadsExec::Active );
    }

    if ( rev_rank ) {
      m_worker_state = QthreadsExec::Inactive;
      Impl::spinwait_while_equal( m_worker_state, QthreadsExec::Inactive );
    }

    for ( n = 1; ( ! ( rev_rank & n ) ) && ( ( j = rev_rank + n ) < m_worker_size ); n <<= 1 ) {
      m_worker_base[j]->m_worker_state = QthreadsExec::Active;
    }
  }

  /** Barrier across workers within the shepherd with rank < team_rank. */
  void shepherd_barrier( const int team_size ) const
  {
    if ( m_shepherd_worker_rank < team_size ) {

      const int rev_rank = team_size - ( m_shepherd_worker_rank + 1 );

      int n, j;

      for ( n = 1; ( ! ( rev_rank & n ) ) && ( ( j = rev_rank + n ) < team_size ); n <<= 1 ) {
        Impl::spinwait_while_equal( m_shepherd_base[j]->m_worker_state, QthreadsExec::Active );
      }

      if ( rev_rank ) {
        m_worker_state = QthreadsExec::Inactive;
        Impl::spinwait_while_equal( m_worker_state, QthreadsExec::Inactive );
      }

      for ( n = 1; ( ! ( rev_rank & n ) ) && ( ( j = rev_rank + n ) < team_size ); n <<= 1 ) {
        m_shepherd_base[j]->m_worker_state = QthreadsExec::Active;
      }
    }
  }

  /** Reduce across all workers participating in the 'exec_all'. */
  template< class FunctorType, class ReducerType, class ArgTag >
  inline
  void exec_all_reduce( const FunctorType & func, const ReducerType & reduce ) const
  {
    typedef Kokkos::Impl::if_c< std::is_same<InvalidType, ReducerType>::value, FunctorType, ReducerType > ReducerConditional;
    typedef typename ReducerConditional::type ReducerTypeFwd;
    typedef Kokkos::Impl::FunctorValueJoin< ReducerTypeFwd, ArgTag > ValueJoin;

    const int rev_rank = m_worker_size - ( m_worker_rank + 1 );

    int n, j;

    for ( n = 1; ( ! ( rev_rank & n ) ) && ( ( j = rev_rank + n ) < m_worker_size ); n <<= 1 ) {
      const QthreadsExec & fan = *m_worker_base[j];

      Impl::spinwait_while_equal( fan.m_worker_state, QthreadsExec::Active );

      ValueJoin::join( ReducerConditional::select( func, reduce ), m_scratch_alloc, fan.m_scratch_alloc );
    }

    if ( rev_rank ) {
      m_worker_state = QthreadsExec::Inactive;
      Impl::spinwait_while_equal( m_worker_state, QthreadsExec::Inactive );
    }

    for ( n = 1; ( ! ( rev_rank & n ) ) && ( ( j = rev_rank + n ) < m_worker_size ); n <<= 1 ) {
      m_worker_base[j]->m_worker_state = QthreadsExec::Active;
    }
  }

  /** Scan across all workers participating in the 'exec_all'. */
  template< class FunctorType, class ArgTag >
  inline
  void exec_all_scan( const FunctorType & func ) const
  {
    typedef Kokkos::Impl::FunctorValueInit< FunctorType, ArgTag > ValueInit;
    typedef Kokkos::Impl::FunctorValueJoin< FunctorType, ArgTag > ValueJoin;
    typedef Kokkos::Impl::FunctorValueOps<  FunctorType, ArgTag > ValueOps;

    const int rev_rank = m_worker_size - ( m_worker_rank + 1 );

    int n, j;

    for ( n = 1; ( ! ( rev_rank & n ) ) && ( ( j = rev_rank + n ) < m_worker_size ); n <<= 1 ) {
      Impl::spinwait_while_equal( m_worker_base[j]->m_worker_state, QthreadsExec::Active );
    }

    if ( rev_rank ) {
      m_worker_state = QthreadsExec::Inactive;
      Impl::spinwait_while_equal( m_worker_state, QthreadsExec::Inactive );
    }
    else {
      // Root thread scans across values before releasing threads.
      // Worker data is in reverse order, so m_worker_base[0] is the
      // highest ranking thread.

      // Copy from lower ranking to higher ranking worker.
      for ( int i = 1; i < m_worker_size; ++i ) {
        ValueOps::copy( func
                      , m_worker_base[i-1]->m_scratch_alloc
                      , m_worker_base[i]->m_scratch_alloc
                      );
      }

      ValueInit::init( func, m_worker_base[m_worker_size-1]->m_scratch_alloc );

      // Join from lower ranking to higher ranking worker.
      // Value at m_worker_base[n-1] is zero so skip adding it to m_worker_base[n-2].
      for ( int i = m_worker_size - 1; --i > 0; ) {
        ValueJoin::join( func, m_worker_base[i-1]->m_scratch_alloc, m_worker_base[i]->m_scratch_alloc );
      }
    }

    for ( n = 1; ( ! ( rev_rank & n ) ) && ( ( j = rev_rank + n ) < m_worker_size ); n <<= 1 ) {
      m_worker_base[j]->m_worker_state = QthreadsExec::Active;
    }
  }

  //----------------------------------------

  template< class Type >
  inline
  volatile Type * shepherd_team_scratch_value() const
  { return (volatile Type*)( ( (unsigned char *) m_scratch_alloc ) + m_reduce_end ); }

  template< class Type >
  inline
  void shepherd_broadcast( Type & value, const int team_size, const int team_rank ) const
  {
    if ( m_shepherd_base ) {
      Type * const shared_value = m_shepherd_base[0]->shepherd_team_scratch_value<Type>();
      if ( m_shepherd_worker_rank == team_rank ) { *shared_value = value; }
      memory_fence();
      shepherd_barrier( team_size );
      value = *shared_value;
    }
  }

  template< class Type >
  inline
  Type shepherd_reduce( const int team_size, const Type & value ) const
  {
    volatile Type * const shared_value = shepherd_team_scratch_value<Type>();
    *shared_value = value;
//    *shepherd_team_scratch_value<Type>() = value;

    memory_fence();

    const int rev_rank = team_size - ( m_shepherd_worker_rank + 1 );

    int n, j;

    for ( n = 1; ( ! ( rev_rank & n ) ) && ( ( j = rev_rank + n ) < team_size ); n <<= 1 ) {
      Impl::spinwait_while_equal( m_shepherd_base[j]->m_worker_state, QthreadsExec::Active );
    }

    if ( rev_rank ) {
      m_worker_state = QthreadsExec::Inactive;
      Impl::spinwait_while_equal( m_worker_state, QthreadsExec::Inactive );
    }
    else {
      Type & accum = *m_shepherd_base[0]->shepherd_team_scratch_value<Type>();
      for ( int i = 1; i < n; ++i ) {
        accum += *m_shepherd_base[i]->shepherd_team_scratch_value<Type>();
      }
      for ( int i = 1; i < n; ++i ) {
        *m_shepherd_base[i]->shepherd_team_scratch_value<Type>() = accum;
      }

      memory_fence();
    }

    for ( n = 1; ( ! ( rev_rank & n ) ) && ( ( j = rev_rank + n ) < team_size ); n <<= 1 ) {
      m_shepherd_base[j]->m_worker_state = QthreadsExec::Active;
    }

    return *shepherd_team_scratch_value<Type>();
  }

  template< class JoinOp >
  inline
  typename JoinOp::value_type
  shepherd_reduce( const int team_size
                 , const typename JoinOp::value_type & value
                 , const JoinOp & op ) const
  {
    typedef typename JoinOp::value_type Type;

    volatile Type * const shared_value = shepherd_team_scratch_value<Type>();
    *shared_value = value;
//    *shepherd_team_scratch_value<Type>() = value;

    memory_fence();

    const int rev_rank = team_size - ( m_shepherd_worker_rank + 1 );

    int n, j;

    for ( n = 1; ( ! ( rev_rank & n ) ) && ( ( j = rev_rank + n ) < team_size ); n <<= 1 ) {
      Impl::spinwait_while_equal( m_shepherd_base[j]->m_worker_state, QthreadsExec::Active );
    }

    if ( rev_rank ) {
      m_worker_state = QthreadsExec::Inactive;
      Impl::spinwait_while_equal( m_worker_state, QthreadsExec::Inactive );
    }
    else {
      volatile Type & accum = *m_shepherd_base[0]->shepherd_team_scratch_value<Type>();
      for ( int i = 1; i < team_size; ++i ) {
        op.join( accum, *m_shepherd_base[i]->shepherd_team_scratch_value<Type>() );
      }
      for ( int i = 1; i < team_size; ++i ) {
        *m_shepherd_base[i]->shepherd_team_scratch_value<Type>() = accum;
      }

      memory_fence();
    }

    for ( n = 1; ( ! ( rev_rank & n ) ) && ( ( j = rev_rank + n ) < team_size ); n <<= 1 ) {
      m_shepherd_base[j]->m_worker_state = QthreadsExec::Active;
    }

    return *shepherd_team_scratch_value<Type>();
  }

  template< class Type >
  inline
  Type shepherd_scan( const int team_size
                    , const Type & value
                    ,       Type * const global_value = 0 ) const
  {
    *shepherd_team_scratch_value<Type>() = value;

    memory_fence();

    const int rev_rank = team_size - ( m_shepherd_worker_rank + 1 );

    int n, j;

    for ( n = 1; ( ! ( rev_rank & n ) ) && ( ( j = rev_rank + n ) < team_size ); n <<= 1 ) {
      Impl::spinwait_while_equal( m_shepherd_base[j]->m_worker_state, QthreadsExec::Active );
    }

    if ( rev_rank ) {
      m_worker_state = QthreadsExec::Inactive;
      Impl::spinwait_while_equal( m_worker_state, QthreadsExec::Inactive );
    }
    else {
      // Root thread scans across values before releasing threads.
      // Worker data is in reverse order, so m_shepherd_base[0] is the
      // highest ranking thread.

      // Copy from lower ranking to higher ranking worker.

      Type accum = *m_shepherd_base[0]->shepherd_team_scratch_value<Type>();
      for ( int i = 1; i < team_size; ++i ) {
        const Type tmp = *m_shepherd_base[i]->shepherd_team_scratch_value<Type>();
        accum += tmp;
        *m_shepherd_base[i-1]->shepherd_team_scratch_value<Type>() = tmp;
      }

      *m_shepherd_base[team_size-1]->shepherd_team_scratch_value<Type>() =
        global_value ? atomic_fetch_add( global_value, accum ) : 0;

      // Join from lower ranking to higher ranking worker.
      for ( int i = team_size; --i; ) {
        *m_shepherd_base[i-1]->shepherd_team_scratch_value<Type>() += *m_shepherd_base[i]->shepherd_team_scratch_value<Type>();
      }

      memory_fence();
    }

    for ( n = 1; ( ! ( rev_rank & n ) ) && ( ( j = rev_rank + n ) < team_size ); n <<= 1 ) {
      m_shepherd_base[j]->m_worker_state = QthreadsExec::Active;
    }

    return *shepherd_team_scratch_value<Type>();
  }

  //----------------------------------------

  static inline
  int align_alloc( int size )
  {
    enum { ALLOC_GRAIN = 1 << 6 /* power of two, 64bytes */ };
    enum { ALLOC_GRAIN_MASK = ALLOC_GRAIN - 1 };
    return ( size + ALLOC_GRAIN_MASK ) & ~ALLOC_GRAIN_MASK;
  }

  void shared_reset( Qthreads::scratch_memory_space & );

  void * exec_all_reduce_value() const { return m_scratch_alloc; }

  static void * exec_all_reduce_result();

  static void resize_worker_scratch( const int reduce_size, const int shared_size );
  static void clear_workers();

  //----------------------------------------

  inline int worker_rank() const { return m_worker_rank; }
  inline int worker_size() const { return m_worker_size; }
  inline int shepherd_worker_rank() const { return m_shepherd_worker_rank; }
  inline int shepherd_worker_size() const { return m_shepherd_worker_size; }
  inline int shepherd_rank() const { return m_shepherd_rank; }
  inline int shepherd_size() const { return m_shepherd_size; }

  static int worker_per_shepherd();
};

} // namespace Impl

} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

namespace Impl {

class QthreadsTeamPolicyMember {
private:
  typedef Kokkos::Qthreads                       execution_space;
  typedef execution_space::scratch_memory_space  scratch_memory_space;

  Impl::QthreadsExec   & m_exec;
  scratch_memory_space   m_team_shared;
  const int              m_team_size;
  const int              m_team_rank;
  const int              m_league_size;
  const int              m_league_end;
        int              m_league_rank;

public:
  KOKKOS_INLINE_FUNCTION
  const scratch_memory_space & team_shmem() const { return m_team_shared; }

  KOKKOS_INLINE_FUNCTION int league_rank() const { return m_league_rank; }
  KOKKOS_INLINE_FUNCTION int league_size() const { return m_league_size; }
  KOKKOS_INLINE_FUNCTION int team_rank() const { return m_team_rank; }
  KOKKOS_INLINE_FUNCTION int team_size() const { return m_team_size; }

  KOKKOS_INLINE_FUNCTION void team_barrier() const
#if ! defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
  {}
#else
  { m_exec.shepherd_barrier( m_team_size ); }
#endif

  template< typename Type >
  KOKKOS_INLINE_FUNCTION Type team_broadcast( const Type & value, int rank ) const
#if ! defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
  { return Type(); }
#else
  { return m_exec.template shepherd_broadcast<Type>( value, m_team_size, rank ); }
#endif

  template< typename Type >
  KOKKOS_INLINE_FUNCTION Type team_reduce( const Type & value ) const
#if ! defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
  { return Type(); }
#else
  { return m_exec.template shepherd_reduce<Type>( m_team_size, value ); }
#endif

  template< typename JoinOp >
  KOKKOS_INLINE_FUNCTION typename JoinOp::value_type
  team_reduce( const typename JoinOp::value_type & value
             , const JoinOp & op ) const
#if ! defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
  { return typename JoinOp::value_type(); }
#else
  { return m_exec.template shepherd_reduce<JoinOp>( m_team_size, value, op ); }
#endif

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
   *
   *  The highest rank thread can compute the reduction total as
   *    reduction_total = dev.team_scan( value ) + value;
   */
  template< typename Type >
  KOKKOS_INLINE_FUNCTION Type team_scan( const Type & value ) const
#if ! defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
  { return Type(); }
#else
  { return m_exec.template shepherd_scan<Type>( m_team_size, value ); }
#endif

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering
   *          with intra-team non-deterministic ordering accumulation.
   *
   *  The global inter-team accumulation value will, at the end of the league's
   *  parallel execution, be the scan's total.  Parallel execution ordering of
   *  the league's teams is non-deterministic.  As such the base value for each
   *  team's scan operation is similarly non-deterministic.
   */
  template< typename Type >
  KOKKOS_INLINE_FUNCTION Type team_scan( const Type & value, Type * const global_accum ) const
#if ! defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
  { return Type(); }
#else
  { return m_exec.template shepherd_scan<Type>( m_team_size, value, global_accum ); }
#endif

  //----------------------------------------
  // Private driver for task-team parallel.

  struct TaskTeam {};

  QthreadsTeamPolicyMember();
  explicit QthreadsTeamPolicyMember( const TaskTeam & );

  //----------------------------------------
  // Private for the driver ( for ( member_type i( exec, team ); i; i.next_team() ) { ... }

  // Initialize.
  template< class ... Properties >
  QthreadsTeamPolicyMember( Impl::QthreadsExec & exec
                          , const Kokkos::Impl::TeamPolicyInternal< Qthreads, Properties... > & team )
    : m_exec( exec )
    , m_team_shared( 0, 0 )
    , m_team_size( team.m_team_size )
    , m_team_rank( exec.shepherd_worker_rank() )
    , m_league_size( team.m_league_size )
    , m_league_end( team.m_league_size - team.m_shepherd_iter * ( exec.shepherd_size() - ( exec.shepherd_rank() + 1 ) ) )
    , m_league_rank( m_league_end > team.m_shepherd_iter ? m_league_end - team.m_shepherd_iter : 0 )
  {
    m_exec.shared_reset( m_team_shared );
  }

  // Continue.
  operator bool () const { return m_league_rank < m_league_end; }

  // Iterate.
  void next_team() { ++m_league_rank; m_exec.shared_reset( m_team_shared ); }
};

template< class ... Properties >
class TeamPolicyInternal< Kokkos::Qthreads, Properties ... >
  : public PolicyTraits< Properties... >
{
private:
  const int m_league_size;
  const int m_team_size;
  const int m_shepherd_iter;

public:
  //! Tag this class as a kokkos execution policy.
  typedef TeamPolicyInternal              execution_policy;
  typedef Qthreads                        execution_space;
  typedef PolicyTraits< Properties ... >  traits;

  //----------------------------------------

  template< class FunctorType >
  inline static
  int team_size_max( const FunctorType & )
  { return Qthreads::instance().shepherd_worker_size(); }

  template< class FunctorType >
  static int team_size_recommended( const FunctorType & f )
  { return team_size_max( f ); }

  template< class FunctorType >
  inline static
  int team_size_recommended( const FunctorType & f, const int& )
  { return team_size_max( f ); }

  //----------------------------------------

  inline int team_size()   const { return m_team_size; }
  inline int league_size() const { return m_league_size; }

  // One active team per shepherd.
  TeamPolicyInternal( Kokkos::Qthreads & q
                    , const int league_size
                    , const int team_size
                    , const int /* vector_length */ = 0
                    )
    : m_league_size( league_size )
    , m_team_size( team_size < q.shepherd_worker_size()
                 ? team_size : q.shepherd_worker_size() )
    , m_shepherd_iter( ( league_size + q.shepherd_size() - 1 ) / q.shepherd_size() )
  {}

  // TODO: Make sure this is correct.
  // One active team per shepherd.
  TeamPolicyInternal( Kokkos::Qthreads & q
                    , const int league_size
                    , const Kokkos::AUTO_t & /* team_size_request */
                    , const int /* vector_length */ = 0
                    )
    : m_league_size( league_size )
    , m_team_size( q.shepherd_worker_size() )
    , m_shepherd_iter( ( league_size + q.shepherd_size() - 1 ) / q.shepherd_size() )
  {}

  // One active team per shepherd.
  TeamPolicyInternal( const int league_size
                    , const int team_size
                    , const int /* vector_length */ = 0
                    )
    : m_league_size( league_size )
    , m_team_size( team_size < Qthreads::instance().shepherd_worker_size()
                 ? team_size : Qthreads::instance().shepherd_worker_size() )
    , m_shepherd_iter( ( league_size + Qthreads::instance().shepherd_size() - 1 ) / Qthreads::instance().shepherd_size() )
  {}

  // TODO: Make sure this is correct.
  // One active team per shepherd.
  TeamPolicyInternal( const int league_size
                    , const Kokkos::AUTO_t & /* team_size_request */
                    , const int /* vector_length */ = 0
                    )
    : m_league_size( league_size )
    , m_team_size( Qthreads::instance().shepherd_worker_size() )
    , m_shepherd_iter( ( league_size + Qthreads::instance().shepherd_size() - 1 ) / Qthreads::instance().shepherd_size() )
  {}

  // TODO: Doesn't do anything yet.  Fix this.
  /** \brief set chunk_size to a discrete value*/
  inline TeamPolicyInternal set_chunk_size(typename traits::index_type chunk_size_) const {
    TeamPolicyInternal p = *this;
//    p.m_chunk_size = chunk_size_;
    return p;
  }

  typedef Impl::QthreadsTeamPolicyMember member_type;

  friend class Impl::QthreadsTeamPolicyMember;
};

} // namespace Impl

} // namespace Kokkos

//----------------------------------------------------------------------------

#endif
#endif // #define KOKKOS_QTHREADSEXEC_HPP

