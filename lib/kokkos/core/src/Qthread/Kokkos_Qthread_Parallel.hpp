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

#ifndef KOKKOS_QTHREAD_PARALLEL_HPP
#define KOKKOS_QTHREAD_PARALLEL_HPP

#include <vector>

#include <Kokkos_Parallel.hpp>

#include <impl/Kokkos_StaticAssert.hpp>
#include <impl/Kokkos_FunctorAdapter.hpp>

#include <Qthread/Kokkos_QthreadExec.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelFor< FunctorType , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Qthread > >
{
private:

  typedef Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Qthread >  Policy ;

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

  // Function is called once by every concurrent thread.
  static void execute( QthreadExec & exec , const void * arg )
  {
    const ParallelFor & self = * ((const ParallelFor *) arg );

    driver( self.m_func , typename Policy::WorkRange( self.m_policy , exec.worker_rank() , exec.worker_size() ) );

    // All threads wait for completion.
    exec.exec_all_barrier();
  }

public:

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

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelReduce< FunctorType , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Qthread > >
{
private:

  typedef Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Qthread >  Policy ;
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , typename Policy::work_tag > ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   FunctorType , typename Policy::work_tag > ValueInit ;

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

  static void execute( QthreadExec & exec , const void * arg )
  {
    const ParallelReduce & self = * ((const ParallelReduce *) arg );

    driver( self.m_func
          , ValueInit::init( self.m_func , exec.exec_all_reduce_value() )
          , typename Policy::WorkRange( self.m_policy , exec.worker_rank() , exec.worker_size() )
          );

    exec.template exec_all_reduce<FunctorType, typename Policy::work_tag >( self.m_func );
  }

public:

  template< class HostViewType >
  ParallelReduce( const FunctorType  & functor
                , const Policy       & policy
                , const HostViewType & result_view )
    : m_func( functor )
    , m_policy( policy )
    {
      QthreadExec::resize_worker_scratch( ValueTraits::value_size( m_func ) , 0 );

      Impl::QthreadExec::exec_all( Qthread::instance() , & ParallelReduce::execute , this );

      const pointer_type data = (pointer_type) QthreadExec::exec_all_reduce_result();

      Kokkos::Impl::FunctorFinal< FunctorType , typename Policy::work_tag >::final( m_func , data );

      if ( result_view.ptr_on_device() ) {
        const unsigned n = ValueTraits::value_count( m_func );
        for ( unsigned i = 0 ; i < n ; ++i ) { result_view.ptr_on_device()[i] = data[i]; }
      }
    }
};

//----------------------------------------------------------------------------

template< class FunctorType , class Arg0 , class Arg1 >
class ParallelFor< FunctorType , TeamPolicy< Arg0 , Arg1 , Kokkos::Qthread > >
{
private:

  typedef TeamPolicy< Arg0 , Arg1 , Kokkos::Qthread >  Policy ;

  const FunctorType  m_func ;
  const Policy       m_team ;

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

  static void execute( QthreadExec & exec , const void * arg )
  {
    const ParallelFor & self = * ((const ParallelFor *) arg );

    typename Policy::member_type member( exec , self.m_team );

    while ( member ) {
      self.ParallelFor::template driver< typename Policy::work_tag >( member );
      member.team_barrier();
      member.next_team();
    }

    exec.exec_all_barrier();
  }

public:

  ParallelFor( const FunctorType & functor ,
               const Policy      & policy )
    : m_func( functor )
    , m_team( policy )
    {
      QthreadExec::resize_worker_scratch
        ( /* reduction   memory */ 0
        , /* team shared memory */ FunctorTeamShmemSize< FunctorType >::value( functor , policy.team_size() ) );

      Impl::QthreadExec::exec_all( Qthread::instance() , & ParallelFor::execute , this );
    }
};

//----------------------------------------------------------------------------

template< class FunctorType , class Arg0 , class Arg1 >
class ParallelReduce< FunctorType , TeamPolicy< Arg0 , Arg1 , Kokkos::Qthread > >
{
private:

  typedef TeamPolicy< Arg0 , Arg1 , Kokkos::Qthread >  Policy ;

  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , typename Policy::work_tag > ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   FunctorType , typename Policy::work_tag > ValueInit ;

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;

  const FunctorType  m_func ;
  const Policy       m_team ;

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

  static void execute( QthreadExec & exec , const void * arg )
  {
    const ParallelReduce & self = * ((const ParallelReduce *) arg );

    // Initialize thread-local value
    reference_type update = ValueInit::init( self.m_func , exec.exec_all_reduce_value() );

    typename Policy::member_type member( exec , self.m_team );

    while ( member ) {
      self.ParallelReduce::template driver< typename Policy::work_tag >( member , update );
      member.team_barrier();
      member.next_team();
    }

    exec.template exec_all_reduce< FunctorType , typename Policy::work_tag >( self.m_func );
  }

public:

  template< class ViewType >
  ParallelReduce( const FunctorType & functor ,
                  const Policy      & policy ,
                  const ViewType    & result )
    : m_func( functor )
    , m_team( policy )
    {
      QthreadExec::resize_worker_scratch
        ( /* reduction   memory */ ValueTraits::value_size( functor )
        , /* team shared memory */ FunctorTeamShmemSize< FunctorType >::value( functor , policy.team_size() ) );

      Impl::QthreadExec::exec_all( Qthread::instance() , & ParallelReduce::execute , this );

      const pointer_type data = (pointer_type) QthreadExec::exec_all_reduce_result();

      Kokkos::Impl::FunctorFinal< FunctorType , typename Policy::work_tag >::final( m_func , data );

      const unsigned n = ValueTraits::value_count( m_func );
      for ( unsigned i = 0 ; i < n ; ++i ) { result.ptr_on_device()[i] = data[i]; }
    }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelScan< FunctorType , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Qthread > >
{
private:

  typedef Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Qthread >  Policy ;
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , typename Policy::work_tag > ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   FunctorType , typename Policy::work_tag > ValueInit ;

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

  static void execute( QthreadExec & exec , const void * arg )
  {
    const ParallelScan & self = * ((const ParallelScan *) arg );

    const typename Policy::WorkRange range( self.m_policy , exec.worker_rank() , exec.worker_size() );

    // Initialize thread-local value
    reference_type update = ValueInit::init( self.m_func , exec.exec_all_reduce_value() );

    driver( self.m_func , update , false , range );

    exec.template exec_all_scan< FunctorType , typename Policy::work_tag >( self.m_func );

    driver( self.m_func , update , true , range );

    exec.exec_all_barrier();
  }

public:

  ParallelScan( const FunctorType & functor
              , const Policy      & policy
              )
    : m_func( functor )
    , m_policy( policy )
    {
      QthreadExec::resize_worker_scratch( ValueTraits::value_size( m_func ) , 0 );

      Impl::QthreadExec::exec_all( Qthread::instance() , & ParallelScan::execute , this );
    }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::TeamThreadRangeBoundariesStruct<iType,Impl::QthreadTeamPolicyMember>
TeamThreadRange(const Impl::QthreadTeamPolicyMember& thread, const iType& count)
{
  return Impl::TeamThreadRangeBoundariesStruct<iType,Impl::QthreadTeamPolicyMember>(thread,count);
}

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::TeamThreadRangeBoundariesStruct<iType,Impl::QthreadTeamPolicyMember>
TeamThreadRange( const Impl::QthreadTeamPolicyMember& thread
               , const iType & begin
               , const iType & end
               )
{
  return Impl::TeamThreadRangeBoundariesStruct<iType,Impl::QthreadTeamPolicyMember>(thread,begin,end);
}


template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::QthreadTeamPolicyMember >
  ThreadVectorRange(const Impl::QthreadTeamPolicyMember& thread, const iType& count) {
  return Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::QthreadTeamPolicyMember >(thread,count);
}


KOKKOS_INLINE_FUNCTION
Impl::ThreadSingleStruct<Impl::QthreadTeamPolicyMember> PerTeam(const Impl::QthreadTeamPolicyMember& thread) {
  return Impl::ThreadSingleStruct<Impl::QthreadTeamPolicyMember>(thread);
}

KOKKOS_INLINE_FUNCTION
Impl::VectorSingleStruct<Impl::QthreadTeamPolicyMember> PerThread(const Impl::QthreadTeamPolicyMember& thread) {
  return Impl::VectorSingleStruct<Impl::QthreadTeamPolicyMember>(thread);
}

} // namespace Kokkos

namespace Kokkos {

  /** \brief  Inter-thread parallel_for. Executes lambda(iType i) for each i=0..N-1.
   *
   * The range i=0..N-1 is mapped to all threads of the the calling thread team.
   * This functionality requires C++11 support.*/
template<typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION
void parallel_for(const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::QthreadTeamPolicyMember>& loop_boundaries, const Lambda& lambda) {
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment)
    lambda(i);
}

/** \brief  Inter-thread vector parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team and a summation of
 * val is performed and put into result. This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ValueType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::QthreadTeamPolicyMember>& loop_boundaries,
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
void parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::QthreadTeamPolicyMember>& loop_boundaries,
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

} // namespace Kokkos

namespace Kokkos {
/** \brief  Intra-thread vector parallel_for. Executes lambda(iType i) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread.
 * This functionality requires C++11 support.*/
template<typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION
void parallel_for(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::QthreadTeamPolicyMember >&
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
void parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::QthreadTeamPolicyMember >&
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
void parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::QthreadTeamPolicyMember >&
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
void parallel_scan(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::QthreadTeamPolicyMember >&
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
void single(const Impl::VectorSingleStruct<Impl::QthreadTeamPolicyMember>& single_struct, const FunctorType& lambda) {
  lambda();
}

template<class FunctorType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::ThreadSingleStruct<Impl::QthreadTeamPolicyMember>& single_struct, const FunctorType& lambda) {
  if(single_struct.team_member.team_rank()==0) lambda();
}

template<class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::VectorSingleStruct<Impl::QthreadTeamPolicyMember>& single_struct, const FunctorType& lambda, ValueType& val) {
  lambda(val);
}

template<class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::ThreadSingleStruct<Impl::QthreadTeamPolicyMember>& single_struct, const FunctorType& lambda, ValueType& val) {
  if(single_struct.team_member.team_rank()==0) {
    lambda(val);
  }
  single_struct.team_member.team_broadcast(val,0);
}

} // namespace Kokkos


#endif /* #define KOKKOS_QTHREAD_PARALLEL_HPP */

