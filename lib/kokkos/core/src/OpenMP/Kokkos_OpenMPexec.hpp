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

  // Fan array has log_2(NT) reduction threads plus 2 scan threads
  // Currently limited to 16k threads.
  enum { MAX_FAN_COUNT    = 16 };
  enum { MAX_THREAD_COUNT = 1 << ( MAX_FAN_COUNT - 2 ) };
  enum { VECTOR_LENGTH    = 8 };
  enum { REDUCE_TEAM_BASE = 512 };

  /** \brief  Thread states for team synchronization */
  enum { Active , Rendezvous , ReductionAvailable , ScanAvailable };

private:

  friend class Kokkos::OpenMP ;

  OpenMPexec * const * m_team_base ;

  void        * m_alloc_reduce ;    ///< Reduction memory
  void        * m_alloc_shared ;    ///< Shared memory
  void        * m_team_shared ;     ///< Shared memory
  int           m_alloc_shared_size ;
  int const     m_pool_rank ;
  int           m_team_shared_end ;
  int           m_team_shared_iter ;
  int           m_team_rank ;
  int           m_team_size ;
  int           m_team_fan_size ;
  int           m_league_rank ;
  int           m_league_end ;
  int           m_league_size ;

  int volatile  m_barrier_state ;
  int volatile  m_scan_state ;

  static OpenMPexec * m_thread[ MAX_THREAD_COUNT ]; // Indexed by 'omp_get_thread_num()'
  static OpenMPexec * m_pool[   MAX_THREAD_COUNT ]; // Indexed by 'm_pool_rank'

  OpenMPexec();
  OpenMPexec( const OpenMPexec & );
  OpenMPexec & operator = ( const OpenMPexec & );

public:

  void * reduce_team() const { return m_alloc_reduce ; }
  void * reduce_base() const { return ((unsigned char *)m_alloc_reduce) + REDUCE_TEAM_BASE ; }

  ~OpenMPexec();

  explicit
  OpenMPexec( const unsigned pool_rank );

  static void finalize();

  static void initialize( const unsigned  team_count ,
                          const unsigned threads_per_team ,
                          const unsigned numa_count ,
                          const unsigned cores_per_numa );

  static void verify_is_process( const char * const );
  static void verify_initialized( const char * const );

  static void resize_reduce_scratch( size_t );
  static void resize_shared_scratch( size_t );

  inline static
  OpenMPexec * get_thread_omp() { return m_thread[ omp_get_thread_num() ]; }

  inline static
  OpenMPexec * get_thread_rank_rev( const int rank_rev ) { return m_pool[ rank_rev ]; }

  //----------------------------------------------------------------------
  /** \brief  Compute a range of work for this thread's rank */

  inline
  std::pair< size_t , size_t >
  work_range( const size_t work_count ) const
  {
    typedef integral_constant< size_t , VECTOR_LENGTH - 1 > work_mask ;

    const size_t thread_size = omp_get_num_threads();

    // work per thread rounded up and aligned to vector length:

    const size_t work_per_thread =
      ( ( ( work_count + thread_size - 1 ) / thread_size ) + work_mask::value ) & ~(work_mask::value);

    const size_t work_begin = std::min( work_count , work_per_thread * m_pool_rank );
    const size_t work_end   = std::min( work_count , work_per_thread + work_begin );

    return std::pair< size_t , size_t >( work_begin , work_end );
  }

  //----------------------------------------------------------------------

  KOKKOS_FUNCTION
  void * get_shmem( const int );

  KOKKOS_INLINE_FUNCTION
  void team_barrier()
    {
      if(m_team_size==1) return;
      const int rank_rev = m_team_size - ( m_team_rank + 1 );

      for ( int i = 0 ; i < m_team_fan_size ; ++i ) {
        Impl::spinwait( m_team_base[ rank_rev + (1<<i) ]->m_barrier_state , OpenMPexec::Active );
      }
      if ( rank_rev ) {
        m_barrier_state = Rendezvous ;
        Impl::spinwait( m_barrier_state , OpenMPexec::Rendezvous );
      }
      for ( int i = 0 ; i < m_team_fan_size ; ++i ) {
        m_team_base[ rank_rev + (1<<i) ]->m_barrier_state = OpenMPexec::Active ;
      }
    }

  template< class ArgType >
  KOKKOS_INLINE_FUNCTION
  ArgType team_scan( const ArgType & value , ArgType * const global_accum = 0 )
    {
      // Sequence of m_scan_state states:
      //  0) Active              : entry and exit state
      //  1) ReductionAvailable  : reduction value available, waiting for scan value
      //  2) ScanAvailable       : reduction value available, scan value available
      //  3) Rendezvous          : broadcasting global inter-team accumulation value

      // Make sure there is enough scratch space:
      typedef typename if_c< 2 * sizeof(ArgType) < REDUCE_TEAM_BASE , ArgType , void >::type type ;

      const int rank_rev = m_team_size - ( m_team_rank + 1 );

      type * const work_value = (type*) reduce_team();

      // OpenMPexec::Active == m_scan_state

      work_value[0] = value ;
      memory_fence();

      // Fan-in reduction, wait for source thread to complete it's fan-in reduction.
      for ( int i = 0 ; i < m_team_fan_size ; ++i ) {
        OpenMPexec & th = *m_team_base[ rank_rev + (1<<i) ];

        // Wait for source thread to exit Active state.
        Impl::spinwait( th.m_scan_state , OpenMPexec::Active );
        // Source thread is 'ReductionAvailable' or 'ScanAvailable'
        work_value[0] += ((volatile type*)th.reduce_team())[0];
        memory_fence();
      }

      work_value[1] = work_value[0] ;
      memory_fence();

      if ( rank_rev ) {

        m_scan_state = OpenMPexec::ReductionAvailable ; // Reduction value is available.

        // Wait for contributing threads' scan value to be available.
        if ( ( 1 << m_team_fan_size ) < ( m_team_rank + 1 ) ) {
          OpenMPexec & th = *m_team_base[ rank_rev + ( 1 << m_team_fan_size ) ];

          // Wait: Active -> ReductionAvailable
          Impl::spinwait( th.m_scan_state , OpenMPexec::Active );
          // Wait: ReductionAvailable -> ScanAvailable:
          Impl::spinwait( th.m_scan_state , OpenMPexec::ReductionAvailable );

          work_value[1] += ((volatile type*)th.reduce_team())[1] ;
          memory_fence();
        }

        m_scan_state = OpenMPexec::ScanAvailable ; // Scan value is available.
      }
      else {
         // Root thread add team's total to global inter-team accumulation
        work_value[0] = global_accum ? atomic_fetch_add( global_accum , work_value[0] ) : 0 ;
      }

      for ( int i = 0 ; i < m_team_fan_size ; ++i ) {
        OpenMPexec & th = *m_team_base[ rank_rev + (1<<i) ];
        // Wait: ReductionAvailable -> ScanAvailable
        Impl::spinwait( th.m_scan_state , OpenMPexec::ReductionAvailable );
        // Wait: ScanAvailable -> Rendezvous
        Impl::spinwait( th.m_scan_state , OpenMPexec::ScanAvailable );
      }

      // All fan-in threads are in the ScanAvailable state
      if ( rank_rev ) {
        m_scan_state = OpenMPexec::Rendezvous ;
        Impl::spinwait( m_scan_state , OpenMPexec::Rendezvous );
      }

      // Broadcast global inter-team accumulation value
      volatile type & global_val = work_value[0] ;
      for ( int i = 0 ; i < m_team_fan_size ; ++i ) {
        OpenMPexec & th = *m_team_base[ rank_rev + (1<<i) ];
        ((volatile type*)th.reduce_team())[0] = global_val ;
        memory_fence();
        th.m_scan_state = OpenMPexec::Active ;
      }
      // Exclusive scan, subtract contributed value
      return global_val + work_value[1] - value ;
    }

  void team_work_init( size_t league_size , size_t team_size );

  inline
  bool team_work_avail()
    { m_team_shared_iter = 0 ; return m_league_rank < m_league_end ; }

  inline
  void team_work_next()
    { if ( ++m_league_rank < m_league_end ) team_barrier(); }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

inline OpenMP::OpenMP( Impl::OpenMPexec & e ) : m_exec(e) {}

KOKKOS_INLINE_FUNCTION
int OpenMP::league_rank() const { return m_exec.m_league_rank ; }
KOKKOS_INLINE_FUNCTION
int OpenMP::league_size() const { return m_exec.m_league_size ; }
KOKKOS_INLINE_FUNCTION
int OpenMP::team_rank() const { return m_exec.m_team_rank ; }
KOKKOS_INLINE_FUNCTION
int OpenMP::team_size() const { return m_exec.m_team_size ; }

KOKKOS_INLINE_FUNCTION
void OpenMP::team_barrier() { m_exec.team_barrier() ; }

KOKKOS_INLINE_FUNCTION
void * OpenMP::get_shmem( const int size ) { return m_exec.get_shmem(size) ; }

template< typename Type >
KOKKOS_INLINE_FUNCTION
Type OpenMP::team_scan( const Type & value )
{ return m_exec.team_scan( value ); }

template< typename TypeLocal , typename TypeGlobal >
KOKKOS_INLINE_FUNCTION
TypeGlobal OpenMP::team_scan( const TypeLocal & value , TypeGlobal * const global_accum )
{ return m_exec.template team_scan< TypeGlobal >( value , global_accum ); }

} // namespace Kokkos

#endif /* #ifndef KOKKOS_OPENMPEXEC_HPP */

