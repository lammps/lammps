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

#ifndef KOKKOS_QTHREAD_PARALLEL_HPP
#define KOKKOS_QTHREAD_PARALLEL_HPP

#include <vector>

#include <Kokkos_Parallel.hpp>

#include <impl/Kokkos_StaticAssert.hpp>

#include <Qthread/Kokkos_QthreadExec.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template< class FunctorType , typename IntType , unsigned P >
class ParallelFor< FunctorType
                 , Kokkos::RangePolicy< Kokkos::Qthread , void , IntType , P >
                 , Kokkos::Qthread
                 >
{
public:
  typedef Kokkos::RangePolicy< Kokkos::Qthread , void , IntType , P >  Policy ;

  const FunctorType  m_func ;
  const Policy       m_policy ;

  // Function is called once by every concurrent thread.
  static void execute( QthreadExec & exec , const void * arg )
  {

    const ParallelFor & self = * ((const ParallelFor *) arg );
    const Policy range( self.m_policy , exec.worker_rank() , exec.worker_size() );

    const typename Policy::member_type work_end = range.end();
    for ( typename Policy::member_type iwork = range.begin() ; iwork < work_end ; ++iwork ) {
      self.m_func( iwork );
    }

    // All threads wait for completion.
    exec.exec_all_barrier();
  }

  ParallelFor( const FunctorType & functor
             , const Policy      & policy
             )
    : m_func( functor )
    , m_policy( policy )
    {
      Impl::QthreadExec::exec_all( Qthread::instance() , & ParallelFor::execute , this );
    }
};

//----------------------------------------------------------------------------

template< class FunctorType , typename IntType , unsigned P >
class ParallelReduce< FunctorType
                    , Kokkos::RangePolicy< Kokkos::Qthread , void , IntType , P >
                    , Kokkos::Qthread
                    >
{
public:

  typedef ReduceAdapter< FunctorType >   Reduce ;
  typedef typename Reduce::pointer_type  pointer_type ;
  typedef Kokkos::RangePolicy< Kokkos::Qthread , void , IntType , P >  Policy ;

  const FunctorType  m_func ;
  const Policy       m_policy ;

  static void execute( QthreadExec & exec , const void * arg )
  {
    const ParallelReduce & self = * ((const ParallelReduce *) arg );
    const Policy range( self.m_policy , exec.worker_rank() , exec.worker_size() );

    // Initialize thread-local value
    typename Reduce::reference_type update = Reduce::init( self.m_func , exec.exec_all_reduce_value() );

    const typename Policy::member_type work_end = range.end();
    for ( typename Policy::member_type iwork = range.begin() ; iwork < work_end ; ++iwork ) {
      self.m_func( iwork , update );
    }

    exec.exec_all_reduce( self.m_func );
  }

  template< class HostViewType >
  ParallelReduce( const FunctorType  & functor
                , const Policy       & policy
                , const HostViewType & result_view )
    : m_func( functor )
    , m_policy( policy )
    {
      QthreadExec::resize_worker_scratch( Reduce::value_size( m_func ) , 0 );

      Impl::QthreadExec::exec_all( Qthread::instance() , & ParallelReduce::execute , this );

      const pointer_type data = (pointer_type) QthreadExec::exec_all_reduce_result();

      Reduce::final( m_func , data );

      if ( result_view.ptr_on_device() ) {
        const unsigned n = Reduce::value_count( m_func );
        for ( unsigned i = 0 ; i < n ; ++i ) { result_view.ptr_on_device()[i] = data[i]; }
      }
    }
};

//----------------------------------------------------------------------------

template< class FunctorType >
class ParallelReduce< FunctorType , TeamPolicy< Kokkos::Qthread > , Kokkos::Qthread >
{
public:

  typedef ReduceAdapter< FunctorType >   Reduce ;
  typedef typename Reduce::pointer_type  pointer_type ;
  typedef TeamPolicy< Kokkos::Qthread >  policy_team ;

  const FunctorType  m_func ;
  const policy_team  m_team ;

  static void execute( QthreadExec & exec , const void * arg )
  {
    const ParallelReduce & self = * ((const ParallelReduce *) arg );

    // Initialize thread-local value
    typename Reduce::reference_type update = Reduce::init( self.m_func , exec.exec_all_reduce_value() );

    typename policy_team::member_type team_index( exec , self.m_team );

    while ( team_index ) {
      // Reset shared memory offset to beginning of reduction range.
      exec.shared_reset();
      self.m_func( team_index , update );
      team_index.team_barrier();
      team_index.next_team();
    }

    exec.exec_all_reduce( self.m_func );
  }

  template< class ViewType >
  ParallelReduce( const FunctorType & functor ,
                  const policy_team & policy ,
                  const ViewType    & result )
    : m_func( functor )
    , m_team( policy )
    {
      QthreadExec::resize_worker_scratch
        ( /* reduction   memory */ Reduce::value_size( functor )
        , /* team shared memory */ FunctorTeamShmemSize< FunctorType >::value( functor ) , policy.team_size() );

      Impl::QthreadExec::exec_all( Qthread::instance() , & ParallelReduce::execute , this );

      const pointer_type data = (pointer_type) QthreadExec::exec_all_reduce_result();

      Reduce::final( m_func , data );

      const unsigned n = Reduce::value_count( m_func );
      for ( unsigned i = 0 ; i < n ; ++i ) { result.ptr_on_device()[i] = data[i]; }
    }

  inline void wait() {}

  inline ~ParallelReduce() { wait(); }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class FunctorType , typename IntType , unsigned P >
class ParallelScan< FunctorType
                  , Kokkos::RangePolicy< Kokkos::Qthread , void , IntType , P >
                  , Kokkos::Qthread
                  >
{
public:

  typedef ReduceAdapter< FunctorType >   Reduce ;
  typedef typename Reduce::pointer_type  pointer_type ;
  typedef Kokkos::RangePolicy< Kokkos::Qthread , void , IntType , P > Policy ;

  const FunctorType  m_func ;
  const Policy       m_policy ;

  static void execute( QthreadExec & exec , const void * arg )
  {
    const ParallelScan & self = * ((const ParallelScan *) arg );
    const Policy range( self.m_policy , exec.worker_rank() , exec.worker_size() );

    // Initialize thread-local value
    typename Reduce::reference_type update = Reduce::init( self.m_func , exec.exec_all_reduce_value() );

    const typename Policy::member_type work_end = range.end();
    for ( typename Policy::member_type iwork = range.begin() ; iwork < work_end ; ++iwork ) {
      self.m_func( iwork , update , false );
    }

    exec.exec_all_scan( self.m_func );

    for ( typename Policy::member_type iwork = range.begin() ; iwork < work_end ; ++iwork ) {
      self.m_func( iwork , update , true );
    }

    exec.exec_all_barrier();
  }

  ParallelScan( const FunctorType & functor
              , const Policy      & policy
              )
    : m_func( functor )
    , m_policy( policy )
    {
      QthreadExec::resize_worker_scratch( Reduce::value_size( m_func ) , 0 );

      Impl::QthreadExec::exec_all( Qthread::instance() , & ParallelScan::execute , this );
    }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #define KOKKOS_QTHREAD_PARALLEL_HPP */

