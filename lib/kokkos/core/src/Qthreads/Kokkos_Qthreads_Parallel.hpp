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

#ifndef KOKKOS_QTHREADS_PARALLEL_HPP
#define KOKKOS_QTHREADS_PARALLEL_HPP

#include <Kokkos_Macros.hpp>
#if defined( KOKKOS_ENABLE_QTHREADS )

#include <vector>

#include <Kokkos_Parallel.hpp>

#include <impl/Kokkos_StaticAssert.hpp>
#include <impl/Kokkos_FunctorAdapter.hpp>

#include <Qthreads/Kokkos_QthreadsExec.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template< class FunctorType , class ... Traits >
class ParallelFor< FunctorType
                 , Kokkos::RangePolicy< Traits ... >
                 , Kokkos::Qthreads
                 >
{
private:

  typedef Kokkos::RangePolicy< Traits ... >  Policy ;

  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::member_type  Member ;
  typedef typename Policy::WorkRange    WorkRange ;

  const FunctorType  m_functor ;
  const Policy       m_policy ;

  template< class TagType >
  inline static
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_range( const FunctorType & functor , const Member ibeg , const Member iend )
    {
      for ( Member i = ibeg ; i < iend ; ++i ) {
        functor( i );
      }
    }

  template< class TagType >
  inline static
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec_range( const FunctorType & functor , const Member ibeg , const Member iend )
    {
      const TagType t{} ;
      for ( Member i = ibeg ; i < iend ; ++i ) {
        functor( t , i );
      }
    }

  // Function is called once by every concurrent thread.
  static void exec( QthreadsExec & exec , const void * arg )
  {
    const ParallelFor & self = * ((const ParallelFor *) arg );

    const WorkRange range( self.m_policy, exec.worker_rank(), exec.worker_size() );

    ParallelFor::template exec_range< WorkTag > ( self.m_functor , range.begin() , range.end() );

    // All threads wait for completion.
    exec.exec_all_barrier();
  }

public:

  inline
  void execute() const
    {
      Impl::QthreadsExec::exec_all( Qthreads::instance() , & ParallelFor::exec , this );

    }

  ParallelFor( const FunctorType & arg_functor
             , const Policy      & arg_policy
             )
    : m_functor( arg_functor )
    , m_policy(  arg_policy )
    { }
};

//----------------------------------------------------------------------------

template< class FunctorType , class ReducerType , class ... Traits >
class ParallelReduce< FunctorType
                    , Kokkos::RangePolicy< Traits ... >
                    , ReducerType
                    , Kokkos::Qthreads
                    >
{
private:

  typedef Kokkos::RangePolicy< Traits ... >  Policy ;

  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::WorkRange    WorkRange ;
  typedef typename Policy::member_type  Member ;

  typedef Kokkos::Impl::if_c< std::is_same<InvalidType, ReducerType>::value, FunctorType, ReducerType > ReducerConditional;
  typedef typename ReducerConditional::type ReducerTypeFwd;
  typedef typename Kokkos::Impl::if_c< std::is_same<InvalidType, ReducerType>::value, WorkTag, void >::type WorkTagFwd;

  // Static Assert WorkTag void if ReducerType not InvalidType

  typedef Kokkos::Impl::FunctorValueTraits< ReducerTypeFwd , WorkTagFwd > ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   ReducerTypeFwd , WorkTagFwd > ValueInit ;

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;

  const FunctorType   m_functor ;
  const Policy        m_policy ;
  const ReducerType   m_reducer ;
  const pointer_type  m_result_ptr ;

  template< class TagType >
  inline static
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_range( const FunctorType & functor
            , const Member ibeg , const Member iend
            , reference_type update )
    {
      for ( Member i = ibeg ; i < iend ; ++i ) {
        functor( i , update );
      }
    }

  template< class TagType >
  inline static
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec_range( const FunctorType & functor
            , const Member ibeg , const Member iend
            , reference_type update )
    {
      const TagType t{} ;
      for ( Member i = ibeg ; i < iend ; ++i ) {
        functor( t , i , update );
      }
    }

  static void exec( QthreadsExec & exec , const void * arg )
  {
    const ParallelReduce & self = * ((const ParallelReduce *) arg );

    const WorkRange range( self.m_policy, exec.worker_rank(), exec.worker_size() );

    ParallelReduce::template exec_range< WorkTag >(
      self.m_functor, range.begin(), range.end(),
      ValueInit::init( ReducerConditional::select(self.m_functor , self.m_reducer)
                     , exec.exec_all_reduce_value() ) );

    exec.template exec_all_reduce< FunctorType, ReducerType, WorkTag >( self.m_functor, self.m_reducer );
  }

public:

  inline
  void execute() const
    {
      QthreadsExec::resize_worker_scratch( ValueTraits::value_size( ReducerConditional::select(m_functor , m_reducer) ) , 0 );
      Impl::QthreadsExec::exec_all( Qthreads::instance() , & ParallelReduce::exec , this );

      const pointer_type data = (pointer_type) QthreadsExec::exec_all_reduce_result();

      Kokkos::Impl::FunctorFinal< ReducerTypeFwd , WorkTagFwd >::final( ReducerConditional::select(m_functor , m_reducer) , data );

      if ( m_result_ptr ) {
        const unsigned n = ValueTraits::value_count( ReducerConditional::select(m_functor , m_reducer) );
        for ( unsigned i = 0 ; i < n ; ++i ) { m_result_ptr[i] = data[i]; }
      }
    }

  template< class ViewType >
  ParallelReduce( const FunctorType  & arg_functor
                , const Policy       & arg_policy
                , const ViewType & arg_result_view
                , typename std::enable_if<Kokkos::is_view< ViewType >::value &&
                                          !Kokkos::is_reducer_type< ReducerType >::value
                                          , void*>::type = NULL)
    : m_functor( arg_functor )
    , m_policy( arg_policy )
    , m_reducer( InvalidType() )
    , m_result_ptr( arg_result_view.data() )
    { }

  ParallelReduce( const FunctorType & arg_functor
                , Policy       arg_policy
                , const ReducerType& reducer )
    : m_functor( arg_functor )
    , m_policy( arg_policy )
    , m_reducer( reducer )
    , m_result_ptr( reducer.result_view().data() )
    { }
};

//----------------------------------------------------------------------------

template< class FunctorType , class ... Properties >
class ParallelFor< FunctorType
                 , TeamPolicy< Properties ... >
                 , Kokkos::Qthreads >
{
private:

  typedef Kokkos::Impl::TeamPolicyInternal< Kokkos::Qthreads , Properties ... > Policy ;
  typedef typename Policy::member_type  Member ;
  typedef typename Policy::work_tag     WorkTag ;

  const FunctorType  m_functor ;
  const Policy       m_policy ;

  template< class TagType >
  inline static
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_team( const FunctorType & functor , Member member )
    {
      while ( member ) {
        functor( member );
        member.team_barrier();
        member.next_team();
      }
    }

  template< class TagType >
  inline static
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec_team( const FunctorType & functor , Member member )
    {
      const TagType t{} ;
      while ( member ) {
        functor( t , member );
        member.team_barrier();
        member.next_team();
      }
    }

  static void exec( QthreadsExec & exec , const void * arg )
  {
    const ParallelFor & self = * ((const ParallelFor *) arg );

    ParallelFor::template exec_team< WorkTag >
      ( self.m_functor , Member( exec , self.m_policy ) );

    exec.exec_all_barrier();
  }

public:

  inline
  void execute() const
    {
      QthreadsExec::resize_worker_scratch
        ( /* reduction   memory */ 0
        , /* team shared memory */ FunctorTeamShmemSize< FunctorType >::value( m_functor , m_policy.team_size() ) );
      Impl::QthreadsExec::exec_all( Qthreads::instance() , & ParallelFor::exec , this );
    }

  ParallelFor( const FunctorType & arg_functor ,
               const Policy      & arg_policy )
    : m_functor( arg_functor )
    , m_policy( arg_policy )
    { }
};

//----------------------------------------------------------------------------

template< class FunctorType , class ReducerType , class ... Properties >
class ParallelReduce< FunctorType
                    , TeamPolicy< Properties... >
                    , ReducerType
                    , Kokkos::Qthreads
                    >
{
private:

  typedef Kokkos::Impl::TeamPolicyInternal< Kokkos::Qthreads , Properties ... > Policy ;

  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::member_type  Member ;

  typedef Kokkos::Impl::if_c< std::is_same<InvalidType,ReducerType>::value, FunctorType, ReducerType> ReducerConditional;
  typedef typename ReducerConditional::type ReducerTypeFwd;
  typedef typename Kokkos::Impl::if_c< std::is_same<InvalidType, ReducerType>::value, WorkTag, void >::type WorkTagFwd;

  typedef Kokkos::Impl::FunctorValueTraits< ReducerTypeFwd , WorkTagFwd >  ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   ReducerTypeFwd , WorkTagFwd >  ValueInit ;

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;

  const FunctorType  m_functor ;
  const Policy       m_policy ;
  const ReducerType  m_reducer ;
  const pointer_type m_result_ptr ;

  template< class TagType >
  inline static
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_team( const FunctorType & functor , Member member , reference_type update )
    {
      while ( member ) {
        functor( member , update );
        member.team_barrier();
        member.next_team();
      }
    }

  template< class TagType >
  inline static
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec_team( const FunctorType & functor , Member member , reference_type update )
    {
      const TagType t{} ;
      while ( member ) {
        functor( t , member , update );
        member.team_barrier();
        member.next_team();
      }
    }

  static void exec( QthreadsExec & exec , const void * arg )
  {
    const ParallelReduce & self = * ((const ParallelReduce *) arg );

    ParallelReduce::template exec_team< WorkTag >
      ( self.m_functor
      , Member( exec , self.m_policy )
      , ValueInit::init( ReducerConditional::select( self.m_functor , self.m_reducer )
                       , exec.exec_all_reduce_value() ) );

    exec.template exec_all_reduce< FunctorType, ReducerType, WorkTag >( self.m_functor, self.m_reducer );
  }

public:

  inline
  void execute() const
    {
      QthreadsExec::resize_worker_scratch
        ( /* reduction   memory */ ValueTraits::value_size( ReducerConditional::select(m_functor , m_reducer) )
        , /* team shared memory */ FunctorTeamShmemSize< FunctorType >::value( m_functor , m_policy.team_size() ) );

      Impl::QthreadsExec::exec_all( Qthreads::instance() , & ParallelReduce::exec , this );

      const pointer_type data = (pointer_type) QthreadsExec::exec_all_reduce_result();

      Kokkos::Impl::FunctorFinal< ReducerTypeFwd , WorkTagFwd >::final( ReducerConditional::select(m_functor , m_reducer), data );

      if ( m_result_ptr ) {
        const unsigned n = ValueTraits::value_count( ReducerConditional::select(m_functor , m_reducer) );
        for ( unsigned i = 0 ; i < n ; ++i ) { m_result_ptr[i] = data[i]; }
      }
    }

  template< class ViewType >
  ParallelReduce( const FunctorType & arg_functor
                , const Policy      & arg_policy
                , const ViewType    & arg_result
                , typename std::enable_if<Kokkos::is_view< ViewType >::value &&
                                          !Kokkos::is_reducer_type< ReducerType >::value
                                          , void*>::type = NULL)
    : m_functor( arg_functor )
    , m_policy( arg_policy )
    , m_reducer( InvalidType() )
    , m_result_ptr( arg_result.ptr_on_device() )
    { }

  inline
  ParallelReduce( const FunctorType & arg_functor
                , Policy       arg_policy
                , const ReducerType& reducer )
  : m_functor( arg_functor )
  , m_policy( arg_policy )
  , m_reducer( reducer )
  , m_result_ptr( reducer.result_view().data() )
  { }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class FunctorType , class ... Traits >
class ParallelScan< FunctorType
                  , Kokkos::RangePolicy< Traits ... >
                  , Kokkos::Qthreads
                  >
{
private:

  typedef Kokkos::RangePolicy< Traits ... >  Policy ;

  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::WorkRange    WorkRange ;
  typedef typename Policy::member_type  Member ;

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
            , const Member ibeg , const Member iend
            , reference_type update , const bool final )
    {
      for ( Member i = ibeg ; i < iend ; ++i ) {
        functor( i , update , final );
      }
    }

  template< class TagType >
  inline static
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec_range( const FunctorType & functor
            , const Member ibeg , const Member iend
            , reference_type update , const bool final )
    {
      const TagType t{} ;
      for ( Member i = ibeg ; i < iend ; ++i ) {
        functor( t , i , update , final );
      }
    }

  static void exec( QthreadsExec & exec , const void * arg )
  {
    const ParallelScan & self = * ((const ParallelScan *) arg );

    const WorkRange range( self.m_policy , exec.worker_rank() , exec.worker_size() );

    // Initialize thread-local value
    reference_type update = ValueInit::init( self.m_functor , exec.exec_all_reduce_value() );

    ParallelScan::template exec_range< WorkTag >( self.m_functor, range.begin() , range.end() , update , false );

    exec.template exec_all_scan< FunctorType , typename Policy::work_tag >( self.m_functor );

    ParallelScan::template exec_range< WorkTag >( self.m_functor , range.begin() , range.end() , update , true );

    exec.exec_all_barrier();
  }

public:

  inline
  void execute() const
    {
      QthreadsExec::resize_worker_scratch( ValueTraits::value_size( m_functor ) , 0 );
      Impl::QthreadsExec::exec_all( Qthreads::instance() , & ParallelScan::exec , this );
    }

  ParallelScan( const FunctorType & arg_functor
              , const Policy      & arg_policy
              )
    : m_functor( arg_functor )
    , m_policy( arg_policy )
    {
    }
};

} // namespace Impl

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template< typename iType >
KOKKOS_INLINE_FUNCTION
Impl::TeamThreadRangeBoundariesStruct< iType, Impl::QthreadsTeamPolicyMember >
TeamThreadRange( const Impl::QthreadsTeamPolicyMember& thread, const iType& count )
{
  return Impl::TeamThreadRangeBoundariesStruct< iType, Impl::QthreadsTeamPolicyMember >( thread, count );
}

template< typename iType1, typename iType2 >
KOKKOS_INLINE_FUNCTION
Impl::TeamThreadRangeBoundariesStruct< typename std::common_type< iType1, iType2 >::type,
                                       Impl::QthreadsTeamPolicyMember >
TeamThreadRange( const Impl::QthreadsTeamPolicyMember& thread, const iType1 & begin, const iType2 & end )
{
  typedef typename std::common_type< iType1, iType2 >::type iType;
  return Impl::TeamThreadRangeBoundariesStruct< iType, Impl::QthreadsTeamPolicyMember >( thread, iType(begin), iType(end) );
}

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::QthreadsTeamPolicyMember >
  ThreadVectorRange(const Impl::QthreadsTeamPolicyMember& thread, const iType& count) {
  return Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::QthreadsTeamPolicyMember >(thread,count);
}

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::QthreadsTeamPolicyMember >
  ThreadVectorRange(const Impl::QthreadsTeamPolicyMember& thread, const iType& arg_begin, const iType& arg_end) {
  return Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::QthreadsTeamPolicyMember >(thread,arg_begin,arg_end);
}

KOKKOS_INLINE_FUNCTION
Impl::ThreadSingleStruct<Impl::QthreadsTeamPolicyMember> PerTeam(const Impl::QthreadsTeamPolicyMember& thread) {
  return Impl::ThreadSingleStruct<Impl::QthreadsTeamPolicyMember>(thread);
}

KOKKOS_INLINE_FUNCTION
Impl::VectorSingleStruct<Impl::QthreadsTeamPolicyMember> PerThread(const Impl::QthreadsTeamPolicyMember& thread) {
  return Impl::VectorSingleStruct<Impl::QthreadsTeamPolicyMember>(thread);
}

/** \brief  Inter-thread parallel_for. Executes lambda(iType i) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team.
 * This functionality requires C++11 support.*/
template<typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION
void parallel_for(const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::QthreadsTeamPolicyMember>& loop_boundaries, const Lambda& lambda) {
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment)
    lambda(i);
}

/** \brief  Inter-thread vector parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team and a summation of
 * val is performed and put into result. This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ValueType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::QthreadsTeamPolicyMember>& loop_boundaries,
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
void parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::QthreadsTeamPolicyMember>& loop_boundaries,
                     const Lambda & lambda, const JoinType& join, ValueType& init_result) {

  ValueType result = init_result;

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    ValueType tmp = ValueType();
    lambda(i,tmp);
    join(result,tmp);
  }

  init_result = loop_boundaries.thread.team_reduce(result,Impl::JoinLambdaAdapter<ValueType,JoinType>(join));
}

/** \brief  Intra-thread vector parallel_for. Executes lambda(iType i) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread.
 * This functionality requires C++11 support.*/
template<typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION
void parallel_for(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::QthreadsTeamPolicyMember >&
    loop_boundaries, const Lambda& lambda) {
  #ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
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
void parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::QthreadsTeamPolicyMember >&
      loop_boundaries, const Lambda & lambda, ValueType& result) {
  result = ValueType();
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
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
void parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::QthreadsTeamPolicyMember >&
      loop_boundaries, const Lambda & lambda, const JoinType& join, ValueType& init_result) {

  ValueType result = init_result;
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
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
void parallel_scan(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::QthreadsTeamPolicyMember >&
      loop_boundaries, const FunctorType & lambda) {

  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , void > ValueTraits ;
  typedef typename ValueTraits::value_type value_type ;

  value_type scan_val = value_type();

#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    lambda(i,scan_val,true);
  }
}

template<class FunctorType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::VectorSingleStruct<Impl::QthreadsTeamPolicyMember>& single_struct, const FunctorType& lambda) {
  lambda();
}

template<class FunctorType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::ThreadSingleStruct<Impl::QthreadsTeamPolicyMember>& single_struct, const FunctorType& lambda) {
  if(single_struct.team_member.team_rank()==0) lambda();
}

template<class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::VectorSingleStruct<Impl::QthreadsTeamPolicyMember>& single_struct, const FunctorType& lambda, ValueType& val) {
  lambda(val);
}

template<class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::ThreadSingleStruct<Impl::QthreadsTeamPolicyMember>& single_struct, const FunctorType& lambda, ValueType& val) {
  if(single_struct.team_member.team_rank()==0) {
    lambda(val);
  }
  single_struct.team_member.team_broadcast(val,0);
}

} // namespace Kokkos

#endif
#endif /* #define KOKKOS_QTHREADS_PARALLEL_HPP */

