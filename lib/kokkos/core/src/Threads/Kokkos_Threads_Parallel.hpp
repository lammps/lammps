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

#ifndef KOKKOS_THREADS_PARALLEL_HPP
#define KOKKOS_THREADS_PARALLEL_HPP

#include <vector>
#include <iostream> 

#include <Kokkos_Parallel.hpp>

#include <impl/Kokkos_StaticAssert.hpp>
#include <impl/Kokkos_FunctorAdapter.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/* ParallelFor Kokkos::Threads with RangePolicy */

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelFor< FunctorType
                 , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Threads >
                 >
{
private:

  typedef Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Threads > Policy ;
  typedef typename Policy::work_tag    WorkTag ;
  typedef typename Policy::WorkRange   WorkRange ;
  typedef typename Policy::member_type Member ;

  const FunctorType  m_functor ;
  const Policy       m_policy ;

  template< class TagType >
  inline static
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_range( const FunctorType & functor
            , const Member ibeg , const Member iend )
    {
      #if defined( KOKKOS_OPT_RANGE_AGGRESSIVE_VECTORIZATION ) && \
          defined( KOKKOS_HAVE_PRAGMA_IVDEP )
      #pragma ivdep
      #endif
      for ( Member i = ibeg ; i < iend ; ++i ) {
        functor( i );
      }
    }

  template< class TagType >
  inline static
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec_range( const FunctorType & functor
            , const Member ibeg , const Member iend )
    {
      const TagType t{} ;
      #if defined( KOKKOS_OPT_RANGE_AGGRESSIVE_VECTORIZATION ) && \
          defined( KOKKOS_HAVE_PRAGMA_IVDEP )
      #pragma ivdep
      #endif
      for ( Member i = ibeg ; i < iend ; ++i ) {
        functor( t , i );
      }
    }

  static void exec( ThreadsExec & exec , const void * arg )
  {
    const ParallelFor & self = * ((const ParallelFor *) arg );

    WorkRange range( self.m_policy , exec.pool_rank() , exec.pool_size() );

    ParallelFor::template exec_range< WorkTag >
      ( self.m_functor , range.begin() , range.end() );

    exec.fan_in();
  }

public:

  inline
  void execute() const
    {
      ThreadsExec::start( & ParallelFor::exec , this );
      ThreadsExec::fence();
    }

  ParallelFor( const FunctorType & arg_functor
             , const Policy      & arg_policy )
    : m_functor( arg_functor )
    , m_policy( arg_policy )
    {}
};

//----------------------------------------------------------------------------
/* ParallelFor Kokkos::Threads with TeamPolicy */

template< class FunctorType , class Arg0 , class Arg1 >
class ParallelFor< FunctorType
                 , Kokkos::TeamPolicy< Arg0 , Arg1 , Kokkos::Threads >
                 >
{
private:

  typedef TeamPolicy< Arg0 , Arg1 , Kokkos::Threads >  Policy ;
  typedef typename Policy::work_tag                    WorkTag ;
  typedef typename Policy::member_type                 Member ;

  const FunctorType  m_functor ;
  const Policy       m_policy ;
  const int          m_shared ;

  template< class TagType >
  inline static
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_team( const FunctorType & functor , Member member )
    {
      for ( ; member.valid() ; member.next() ) {
        functor( member );
      }
    }

  template< class TagType >
  inline static
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec_team( const FunctorType & functor , Member member )
    {
      const TagType t{} ;
      for ( ; member.valid() ; member.next() ) {
        functor( t , member );
      }
    }

  static void exec( ThreadsExec & exec , const void * arg )
  {
    const ParallelFor & self = * ((const ParallelFor *) arg );

    ParallelFor::exec_team< WorkTag >
      ( self.m_functor , Member( & exec , self.m_policy , self.m_shared ) );

    exec.fan_in();
  }

public:

  inline
  void execute() const
    {
      ThreadsExec::resize_scratch( 0 , Policy::member_type::team_reduce_size() + m_shared );

      ThreadsExec::start( & ParallelFor::exec , this );

      ThreadsExec::fence();
    }

  ParallelFor( const FunctorType & arg_functor
             , const Policy      & arg_policy )
    : m_functor( arg_functor )
    , m_policy(  arg_policy )
    , m_shared( arg_policy.scratch_size() + FunctorTeamShmemSize< FunctorType >::value( arg_functor , arg_policy.team_size() ) )
    { }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/* ParallelReduce with Kokkos::Threads and RangePolicy */

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelReduce< FunctorType
                    , Kokkos::RangePolicy< Arg0, Arg1, Arg2, Kokkos::Threads >
                    >
{
private:

  typedef Kokkos::RangePolicy< Arg0 , Arg1, Arg2, Kokkos::Threads > Policy ;

  typedef typename Policy::work_tag    WorkTag ;
  typedef typename Policy::WorkRange   WorkRange ;
  typedef typename Policy::member_type Member ;

  typedef Kokkos::Impl::FunctorValueTraits< FunctorType, WorkTag > ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   FunctorType, WorkTag > ValueInit ;

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;

  const FunctorType  m_functor ;
  const Policy       m_policy ;
  const pointer_type m_result_ptr ;

  template< class TagType >
  inline static
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_range( const FunctorType & functor
            , const Member & ibeg , const Member & iend
            , reference_type update )
    {
      #if defined( KOKKOS_OPT_RANGE_AGGRESSIVE_VECTORIZATION ) && \
          defined( KOKKOS_HAVE_PRAGMA_IVDEP )
      #pragma ivdep
      #endif
      for ( Member i = ibeg ; i < iend ; ++i ) {
        functor( i , update );
      }
    }

  template< class TagType >
  inline static
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec_range( const FunctorType & functor
            , const Member & ibeg , const Member & iend
            , reference_type update )
    {
      const TagType t{} ;
      #if defined( KOKKOS_OPT_RANGE_AGGRESSIVE_VECTORIZATION ) && \
          defined( KOKKOS_HAVE_PRAGMA_IVDEP )
      #pragma ivdep
      #endif
      for ( Member i = ibeg ; i < iend ; ++i ) {
        functor( t , i , update );
      }
    }

  static void exec( ThreadsExec & exec , const void * arg )
  {
    const ParallelReduce & self = * ((const ParallelReduce *) arg );
    const WorkRange range( self.m_policy, exec.pool_rank(), exec.pool_size() );

    ParallelReduce::template exec_range< WorkTag >
      ( self.m_functor , range.begin() , range.end() 
      , ValueInit::init( self.m_functor , exec.reduce_memory() ) );

    exec.template fan_in_reduce< FunctorType , WorkTag >( self.m_functor );
  }

public:

  inline
  void execute() const
    {
      ThreadsExec::resize_scratch( ValueTraits::value_size( m_functor ) , 0 );

      ThreadsExec::start( & ParallelReduce::exec , this );

      ThreadsExec::fence();

      if ( m_result_ptr ) {

        const pointer_type data =
          (pointer_type) ThreadsExec::root_reduce_scratch();

        const unsigned n = ValueTraits::value_count( m_functor );
        for ( unsigned i = 0 ; i < n ; ++i ) { m_result_ptr[i] = data[i]; }
      }
    }

  template< class HostViewType >
  ParallelReduce( const FunctorType  & arg_functor ,
                  const Policy       & arg_policy ,
                  const HostViewType & arg_result_view )
    : m_functor( arg_functor )
    , m_policy( arg_policy )
    , m_result_ptr( arg_result_view.ptr_on_device() )
    {
      static_assert( Kokkos::is_view< HostViewType >::value
        , "Kokkos::Threads reduce result must be a View" );

      static_assert( std::is_same< typename HostViewType::memory_space , HostSpace >::value
        , "Kokkos::Threads reduce result must be a View in HostSpace" );
    }
};

//----------------------------------------------------------------------------
/* ParallelReduce with Kokkos::Threads and TeamPolicy */

template< class FunctorType , class Arg0 , class Arg1 >
class ParallelReduce< FunctorType
                    , Kokkos::TeamPolicy< Arg0 , Arg1 , Kokkos::Threads >
                    >
{
private:

  typedef TeamPolicy< Arg0 , Arg1 , Kokkos::Threads >              Policy ;
  typedef typename Policy::work_tag                                WorkTag ;
  typedef typename Policy::member_type                             Member ;
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType, WorkTag > ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   FunctorType, WorkTag > ValueInit ;

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;

  const FunctorType  m_functor ;
  const Policy       m_policy ;
  const pointer_type m_result_ptr ;
  const int          m_shared ;

  template< class TagType >
  inline static
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_team( const FunctorType & functor , Member member , reference_type update )
    {
      for ( ; member.valid() ; member.next() ) {
        functor( member , update );
      }
    }

  template< class TagType >
  inline static
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec_team( const FunctorType & functor , Member member , reference_type update )
    {
      const TagType t{} ;
      for ( ; member.valid() ; member.next() ) {
        functor( t , member , update );
      }
    }

  static void exec( ThreadsExec & exec , const void * arg )
  {
    const ParallelReduce & self = * ((const ParallelReduce *) arg );

    ParallelReduce::template exec_team< WorkTag >
      ( self.m_functor , Member( & exec , self.m_policy , self.m_shared )
      , ValueInit::init( self.m_functor , exec.reduce_memory() ) );

    exec.template fan_in_reduce< FunctorType , WorkTag >( self.m_functor );
  }

public:

  inline
  void execute() const
    {
      ThreadsExec::resize_scratch( ValueTraits::value_size( m_functor ) , Policy::member_type::team_reduce_size() + m_shared );

      ThreadsExec::start( & ParallelReduce::exec , this );

      ThreadsExec::fence();

      if ( m_result_ptr ) {

        const pointer_type data = (pointer_type) ThreadsExec::root_reduce_scratch();

        const unsigned n = ValueTraits::value_count( m_functor );
        for ( unsigned i = 0 ; i < n ; ++i ) { m_result_ptr[i] = data[i]; }
      }
    }

  template< class ViewType >
  ParallelReduce( const FunctorType & arg_functor
                , const Policy      & arg_policy
                , const ViewType    & arg_result )
    : m_functor( arg_functor )
    , m_policy( arg_policy )
    , m_result_ptr( arg_result.ptr_on_device() )
    , m_shared( arg_policy.scratch_size() + FunctorTeamShmemSize< FunctorType >::value( arg_functor , arg_policy.team_size() ) )
    { }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/* ParallelScan with Kokkos::Threads and RangePolicy */

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelScan< FunctorType
                  , Kokkos::RangePolicy< Arg0, Arg1, Arg2, Kokkos::Threads >
                  >
{
private:

  typedef Kokkos::RangePolicy< Arg0, Arg1, Arg2, Kokkos::Threads > Policy ;
  typedef typename Policy::WorkRange                               WorkRange ;
  typedef typename Policy::work_tag                                WorkTag ;
  typedef typename Policy::member_type                             Member ;
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType, WorkTag > ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   FunctorType, WorkTag > ValueInit ;

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;

  const FunctorType  m_functor ;
  const Policy       m_policy ;

  template< class TagType >
  inline static
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_range( const FunctorType & functor
            , const Member & ibeg , const Member & iend
            , reference_type update , const bool final )
    {
      #if defined( KOKKOS_OPT_RANGE_AGGRESSIVE_VECTORIZATION ) && \
          defined( KOKKOS_HAVE_PRAGMA_IVDEP )
      #pragma ivdep
      #endif
      for ( Member i = ibeg ; i < iend ; ++i ) {
        functor( i , update , final );
      }
    }

  template< class TagType >
  inline static
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec_range( const FunctorType & functor
            , const Member & ibeg , const Member & iend
            , reference_type update , const bool final )
    {
      const TagType t{} ;
      #if defined( KOKKOS_OPT_RANGE_AGGRESSIVE_VECTORIZATION ) && \
          defined( KOKKOS_HAVE_PRAGMA_IVDEP )
      #pragma ivdep
      #endif
      for ( Member i = ibeg ; i < iend ; ++i ) {
        functor( t , i , update , final );
      }
    }

  static void exec( ThreadsExec & exec , const void * arg )
  {
    const ParallelScan & self = * ((const ParallelScan *) arg );

    const WorkRange range( self.m_policy, exec.pool_rank(), exec.pool_size() );

    reference_type update =
      ValueInit::init( self.m_functor , exec.reduce_memory() );

    ParallelScan::template exec_range< WorkTag >
      ( self.m_functor , range.begin(), range.end(), update, false );

    //  exec.template scan_large<FunctorType,WorkTag>( self.m_functor );
    exec.template scan_small<FunctorType,WorkTag>( self.m_functor );

    ParallelScan::template exec_range< WorkTag >
      ( self.m_functor , range.begin(), range.end(), update, true );

    exec.fan_in();
  }

public:

  inline
  void execute() const
    {
      ThreadsExec::resize_scratch( 2 * ValueTraits::value_size( m_functor ) , 0 );
      ThreadsExec::start( & ParallelScan::exec , this );
      ThreadsExec::fence();
    }

  ParallelScan( const FunctorType & arg_functor
              , const Policy      & arg_policy )
    : m_functor( arg_functor )
    , m_policy( arg_policy )
    { }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #define KOKKOS_THREADS_PARALLEL_HPP */

