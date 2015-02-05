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

#ifndef KOKKOS_THREADS_PARALLEL_HPP
#define KOKKOS_THREADS_PARALLEL_HPP

#include <vector>

#include <Kokkos_Parallel.hpp>

#include <impl/Kokkos_StaticAssert.hpp>
#include <impl/Kokkos_FunctorAdapter.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelFor< FunctorType , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Threads > >
{
private:

  typedef Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Threads > Policy ;

  const FunctorType  m_func ;
  const Policy       m_policy ;

  template< class PType >
  KOKKOS_FORCEINLINE_FUNCTION static
  void driver( typename Impl::enable_if<
                 ( Impl::is_same< typename PType::work_tag , void >::value )
                 , const FunctorType & >::type functor
             , const PType & range )
    {
      const typename PType::member_type e = range.end();
      for ( typename PType::member_type i = range.begin() ; i < e ; ++i ) {
        functor( i );
      }
    }

  template< class PType >
  KOKKOS_FORCEINLINE_FUNCTION static
  void driver( typename Impl::enable_if<
                 ( ! Impl::is_same< typename PType::work_tag , void >::value )
                 , const FunctorType & >::type functor
             , const PType & range )
    {
      const typename PType::member_type e = range.end();
      for ( typename PType::member_type i = range.begin() ; i < e ; ++i ) {
        functor( typename PType::work_tag() , i );
      }
    }

  static void execute( ThreadsExec & exec , const void * arg )
  {
    const ParallelFor & self = * ((const ParallelFor *) arg );

    driver( self.m_func , typename Policy::WorkRange( self.m_policy , exec.pool_rank() , exec.pool_size() ) );

    exec.fan_in();
  }

public:

  ParallelFor( const FunctorType & functor
             , const Policy      & policy )
    : m_func( functor )
    , m_policy( policy )
    {
      ThreadsExec::start( & ParallelFor::execute , this );

      ThreadsExec::fence();
    }
};

template< class FunctorType , class Arg0 , class Arg1 >
class ParallelFor< FunctorType , Kokkos::TeamPolicy< Arg0 , Arg1 , Kokkos::Threads > >
{
private:

  typedef TeamPolicy< Arg0 , Arg1 , Kokkos::Threads >  Policy ;

  const FunctorType  m_func ;
  const Policy       m_policy ;
  const int          m_shared ;

  template< class TagType >
  KOKKOS_FORCEINLINE_FUNCTION
  void driver( typename Impl::enable_if< Impl::is_same< TagType , void >::value ,
                 const typename Policy::member_type & >::type member ) const
    { m_func( member ); }

  template< class TagType >
  KOKKOS_FORCEINLINE_FUNCTION
  void driver( typename Impl::enable_if< ! Impl::is_same< TagType , void >::value ,
                 const typename Policy::member_type & >::type member ) const
    { m_func( TagType() , member ); }

  static void execute( ThreadsExec & exec , const void * arg )
  {
    const ParallelFor & self = * ((const ParallelFor *) arg );

    typename Policy::member_type member( exec , self.m_policy , self.m_shared );

    for ( ; member.valid() ; member.next() ) {
      self.ParallelFor::template driver< typename Policy::work_tag >( member );
    }

    exec.fan_in();
  }

public:

  ParallelFor( const FunctorType & functor
              , const Policy      & policy )
    : m_func( functor )
    , m_policy( policy )
    , m_shared( FunctorTeamShmemSize< FunctorType >::value( functor , policy.team_size() ) )
    {
      ThreadsExec::resize_scratch( 0 , Policy::member_type::team_reduce_size() + m_shared );

      ThreadsExec::start( & ParallelFor::execute , this );

      ThreadsExec::fence();
    }
};



//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelReduce< FunctorType , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Threads > >
{
private:

  typedef Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Threads > Policy ;
  typedef typename Policy::work_tag                                   work_tag ;
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , work_tag >  ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   FunctorType , work_tag >  ValueInit ;

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;

  const FunctorType  m_func ;
  const Policy       m_policy ;

  template< class PType >
  KOKKOS_FORCEINLINE_FUNCTION static
  void driver( typename Impl::enable_if<
                 ( Impl::is_same< typename PType::work_tag , void >::value )
                 , const FunctorType & >::type functor
             , reference_type update
             , const PType & range )
    {
      const typename PType::member_type e = range.end();
      for ( typename PType::member_type i = range.begin() ; i < e ; ++i ) {
        functor( i , update );
      }
    }

  template< class PType >
  KOKKOS_FORCEINLINE_FUNCTION static
  void driver( typename Impl::enable_if<
                 ( ! Impl::is_same< typename PType::work_tag , void >::value )
                 , const FunctorType & >::type functor
             , reference_type update
             , const PType & range )
    {
      const typename PType::member_type e = range.end();
      for ( typename PType::member_type i = range.begin() ; i < e ; ++i ) {
        functor( typename PType::work_tag() , i , update );
      }
    }

  static void execute( ThreadsExec & exec , const void * arg )
  {
    const ParallelReduce & self = * ((const ParallelReduce *) arg );

    driver( self.m_func
          , ValueInit::init( self.m_func , exec.reduce_memory() )
          , typename Policy::WorkRange( self.m_policy , exec.pool_rank() , exec.pool_size() )
          );

    exec.template fan_in_reduce< FunctorType , work_tag >( self.m_func );
  }

public:

  template< class HostViewType >
  ParallelReduce( const FunctorType  & functor ,
                  const Policy       & policy ,
                  const HostViewType & result_view )
    : m_func( functor )
    , m_policy( policy )
    {
      ThreadsExec::resize_scratch( ValueTraits::value_size( m_func ) , 0 );

      ThreadsExec::start( & ParallelReduce::execute , this );

      const pointer_type data = (pointer_type) ThreadsExec::root_reduce_scratch();

      ThreadsExec::fence();

      if ( result_view.ptr_on_device() ) {
        const unsigned n = ValueTraits::value_count( m_func );
        for ( unsigned i = 0 ; i < n ; ++i ) { result_view.ptr_on_device()[i] = data[i]; }
      }
    }
};

//----------------------------------------------------------------------------

template< class FunctorType , class Arg0 , class Arg1 >
class ParallelReduce< FunctorType , Kokkos::TeamPolicy< Arg0 , Arg1 , Kokkos::Threads > >
{
private:

  typedef TeamPolicy< Arg0 , Arg1 , Kokkos::Threads >                 Policy ;
  typedef typename Policy::work_tag                                   work_tag ;
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , work_tag >  ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   FunctorType , work_tag >  ValueInit ;

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;

  const FunctorType  m_func ;
  const Policy       m_policy ;
  const int          m_shared ;

  template< class TagType >
  KOKKOS_FORCEINLINE_FUNCTION
  void driver( typename Impl::enable_if< Impl::is_same< TagType , void >::value ,
                 const typename Policy::member_type & >::type member
             , reference_type update ) const
    { m_func( member , update ); }

  template< class TagType >
  KOKKOS_FORCEINLINE_FUNCTION
  void driver( typename Impl::enable_if< ! Impl::is_same< TagType , void >::value ,
                 const typename Policy::member_type & >::type member
             , reference_type update ) const
    { m_func( TagType() , member , update ); }

  static void execute( ThreadsExec & exec , const void * arg )
  {
    const ParallelReduce & self = * ((const ParallelReduce *) arg );

    // Initialize thread-local value
    reference_type update = ValueInit::init( self.m_func , exec.reduce_memory() );

    typename Policy::member_type member( exec , self.m_policy , self.m_shared );
    for ( ; member.valid() ; member.next() ) {
      self.ParallelReduce::template driver< work_tag >( member , update );
    }

    exec.template fan_in_reduce< FunctorType , work_tag >( self.m_func );
  }

public:

  ParallelReduce( const FunctorType & functor
                , const Policy      & policy )
    : m_func( functor )
    , m_policy( policy )
    , m_shared( FunctorTeamShmemSize< FunctorType >::value( functor , policy.team_size() ) )
    {
      ThreadsExec::resize_scratch( ValueTraits::value_size( m_func ) , Policy::member_type::team_reduce_size() + m_shared );

      ThreadsExec::start( & ParallelReduce::execute , this );

      ThreadsExec::fence();
    }

  template< class ViewType >
  ParallelReduce( const FunctorType & functor
                , const Policy      & policy
                , const ViewType    & result )
    : m_func( functor )
    , m_policy( policy )
    , m_shared( FunctorTeamShmemSize< FunctorType >::value( functor , policy.team_size() ) )
    {
      ThreadsExec::resize_scratch( ValueTraits::value_size( m_func ) , Policy::member_type::team_reduce_size() + m_shared );

      ThreadsExec::start( & ParallelReduce::execute , this );

      const pointer_type data = (pointer_type) ThreadsExec::root_reduce_scratch();

      ThreadsExec::fence();

      const unsigned n = ValueTraits::value_count( m_func );
      for ( unsigned i = 0 ; i < n ; ++i ) { result.ptr_on_device()[i] = data[i]; }
    }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelScan< FunctorType , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Threads > >
{
private:

  typedef Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Threads > Policy ;
  typedef typename Policy::work_tag                                   work_tag ;
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , work_tag >  ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   FunctorType , work_tag >  ValueInit ;

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;

  const FunctorType  m_func ;
  const Policy       m_policy ;

  template< class PType >
  KOKKOS_FORCEINLINE_FUNCTION static
  void driver( typename Impl::enable_if<
                 ( Impl::is_same< typename PType::work_tag , void >::value )
                 , const FunctorType & >::type functor
             , reference_type update
             , const bool    final
             , const PType & range )
    {
      const typename PType::member_type e = range.end();
      for ( typename PType::member_type i = range.begin() ; i < e ; ++i ) {
        functor( i , update , final );
      }
    }

  template< class PType >
  KOKKOS_FORCEINLINE_FUNCTION static
  void driver( typename Impl::enable_if<
                 ( ! Impl::is_same< typename PType::work_tag , void >::value )
                 , const FunctorType & >::type functor
             , reference_type update
             , const bool    final
             , const PType & range )
    {
      const typename PType::member_type e = range.end();
      for ( typename PType::member_type i = range.begin() ; i < e ; ++i ) {
        functor( typename PType::work_tag() , i , update , final );
      }
    }

  static void execute( ThreadsExec & exec , const void * arg )
  {
    const ParallelScan & self = * ((const ParallelScan *) arg );

    const typename Policy::WorkRange range( self.m_policy , exec.pool_rank() , exec.pool_size() );

    reference_type update = ValueInit::init( self.m_func , exec.reduce_memory() );

    driver( self.m_func , update , false , range );

    //  exec.<FunctorType,work_tag>scan_large( self.m_func );
    exec.template scan_small<FunctorType,work_tag>( self.m_func );

    driver( self.m_func , update , true , range );

    exec.fan_in();
  }

public:

  ParallelScan( const FunctorType & functor , const Policy & policy )
    : m_func( functor )
    , m_policy( policy )
    {
      ThreadsExec::resize_scratch( 2 * ValueTraits::value_size( m_func ) , 0 );
      ThreadsExec::start( & ParallelScan::execute , this );
      ThreadsExec::fence();
    }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #define KOKKOS_THREADS_PARALLEL_HPP */

