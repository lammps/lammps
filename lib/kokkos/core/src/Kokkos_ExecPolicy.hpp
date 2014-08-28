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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_EXECPOLICY_HPP
#define KOKKOS_EXECPOLICY_HPP

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_Tags.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief  Execution policy for work over a range of an integral type.
 */
template< class ExecSpace  = Kokkos::DefaultExecutionSpace
        , class WorkArgTag = void
        , typename IntType = int
        , unsigned GranularityPowerOfTwo = 3 /* Chunk size 8 */
        >
class RangePolicy {
private:

  enum { Granularity     = IntType(1) << GranularityPowerOfTwo };
  enum { GranularityMask = IntType(Granularity) - 1 };

  IntType m_begin ;
  IntType m_end ;

public:

  typedef Impl::ExecutionPolicyTag   kokkos_tag ;      ///< Concept tag
  typedef ExecSpace                  execution_space ; ///< Execution type
  typedef IntType                    member_type ;

  KOKKOS_INLINE_FUNCTION member_type begin() const { return m_begin ; }
  KOKKOS_INLINE_FUNCTION member_type end()   const { return m_end ; }

  KOKKOS_INLINE_FUNCTION RangePolicy() : m_begin(0), m_end(0) {}

  /** \brief  Total range */
  KOKKOS_INLINE_FUNCTION
  RangePolicy( const member_type work_begin
             , const member_type work_end
             )
    : m_begin( work_begin < work_end ? work_begin : 0 )
    , m_end(   work_begin < work_end ? work_end : 0 )
    {}

  /** \brief  Subrange for a partition's rank and size.
   *
   *  Typically used to partition a range over a group of threads.
   */
  KOKKOS_INLINE_FUNCTION
  RangePolicy( const RangePolicy & range
             , const int part_rank
             , const int part_size
             )
    : m_begin(0), m_end(0)
    {
      if ( part_size ) {

        // Split evenly among partitions, then round up to the granularity.
        const member_type work_part =
          ( ( ( ( range.m_end - range.m_begin ) + ( part_size - 1 ) ) / part_size ) + GranularityMask ) & ~member_type(GranularityMask);

        m_begin = range.m_begin + work_part * part_rank ;
        m_end   = m_begin       + work_part ;

        if ( range.m_end < m_begin ) m_begin = range.m_end ;
        if ( range.m_end < m_end )   m_end   = range.m_end ;
      }
    }
};

} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief  Execution policy for parallel work over a league of teams of threads.
 *
 *  The work functor is called for each thread of each team such that
 *  the team's member threads are guaranteed to be concurrent.
 *
 *  The team's threads have access to team shared scratch memory and
 *  team collective operations.
 *
 *  If the WorkArgTag is non-void then the first calling argument of the
 *  work functor's parentheses operator is 'const WorkArgTag &'.
 *  This allows a functor to have multiple work member functions.
 */
template< class ExecSpace  = DefaultExecutionSpace
        , class WorkArgTag = void >
class TeamPolicy {
public:

  typedef Impl::ExecutionPolicyTag   kokkos_tag ;      ///< Concept tag
  typedef ExecSpace                  execution_space ; ///< Execution space

  /** \brief  Query maximum team size for a given functor.
   *
   *  This size takes into account execution space concurrency limitations and
   *  scratch memory space limitations for reductions, team reduce/scan, and
   *  team shared memory.
   */
  template< class FunctorType >
  static int team_size_max( const FunctorType & );

  /** \brief  Construct policy with the given instance of the execution space */
  TeamPolicy( execution_space & , int league_size_request , int team_size_request );

  /** \brief  Construct policy with the default instance of the execution space */
  TeamPolicy( int league_size_request , int team_size_request );

  /** \brief  The actual league size (number of teams) of the policy.
   *
   *  This may be smaller than the requested league size due to limitations
   *  of the execution space.
   */
  KOKKOS_INLINE_FUNCTION int league_size() const ;

  /** \brief  The actual team size (number of threads per team) of the policy.
   *
   *  This may be smaller than the requested team size due to limitations
   *  of the execution space.
   */
  KOKKOS_INLINE_FUNCTION int team_size() const ;

  /** \brief  Parallel execution of a functor calls the functor once with
   *          each member of the execution policy.
   */
  struct member_type {

    /** \brief  Handle to the currently executing team shared scratch memory */
    KOKKOS_INLINE_FUNCTION
    typename execution_space::scratch_memory_space team_shmem() const ;

    /** \brief  Rank of this team within the league of teams */
    KOKKOS_INLINE_FUNCTION int league_rank() const ;

    /** \brief  Number of teams in the league */
    KOKKOS_INLINE_FUNCTION int league_size() const ;

    /** \brief  Rank of this thread within this team */
    KOKKOS_INLINE_FUNCTION int team_rank() const ;

    /** \brief  Number of threads in this team */
    KOKKOS_INLINE_FUNCTION int team_size() const ;

    /** \brief  Barrier among the threads of this team */
    KOKKOS_INLINE_FUNCTION void team_barrier();

    /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
     *
     *  The highest rank thread can compute the reduction total as
     *    reduction_total = dev.team_scan( value ) + value ;
     */
    template< typename Type >
    KOKKOS_INLINE_FUNCTION Type team_scan( const Type & value );

    /** \brief  Intra-team exclusive prefix sum with team_rank() ordering
     *          with intra-team non-deterministic ordering accumulation.
     *
     *  The global inter-team accumulation value will, at the end of the
     *  league's parallel execution, be the scan's total.
     *  Parallel execution ordering of the league's teams is non-deterministic.
     *  As such the base value for each team's scan operation is similarly
     *  non-deterministic.
     */
    template< typename Type >
    KOKKOS_INLINE_FUNCTION Type team_scan( const Type & value , Type * const global_accum );
  };
};

} // namespace Kokkos

#endif /* #define KOKKOS_EXECPOLICY_HPP */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

