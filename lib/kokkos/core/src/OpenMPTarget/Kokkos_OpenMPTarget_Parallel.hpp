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

#ifndef KOKKOS_OPENMPTARGET_PARALLEL_HPP
#define KOKKOS_OPENMPTARGET_PARALLEL_HPP

#include <omp.h>
#include <iostream>
#include <Kokkos_Parallel.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_Exec.hpp>
#include <impl/Kokkos_FunctorAdapter.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class ... Traits >
class ParallelFor< FunctorType
                 , Kokkos::RangePolicy< Traits ... >
                 , Kokkos::Experimental::OpenMPTarget 
                 >
{
private:

  typedef Kokkos::RangePolicy< Traits ...  > Policy ;
  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::WorkRange    WorkRange ;
  typedef typename Policy::member_type  Member ;

  const FunctorType m_functor ;
  const Policy      m_policy ;


public:

  inline void execute() const {
    execute_impl<WorkTag>();
  }

  template< class TagType >
  inline
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  execute_impl() const
    {
      OpenMPTargetExec::verify_is_process("Kokkos::Experimental::OpenMPTarget parallel_for");
      OpenMPTargetExec::verify_initialized("Kokkos::Experimental::OpenMPTarget parallel_for");
      const typename Policy::member_type begin = m_policy.begin();
      const typename Policy::member_type end = m_policy.end();
      
      #pragma omp target teams distribute parallel for map(to:this->m_functor)
      for(int i=begin; i<end; i++)
        m_functor(i);
    }


  template< class TagType >
  inline
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  execute_impl() const
    {
      OpenMPTargetExec::verify_is_process("Kokkos::Experimental::OpenMPTarget parallel_for");
      OpenMPTargetExec::verify_initialized("Kokkos::Experimental::OpenMPTarget parallel_for");
      const typename Policy::member_type begin = m_policy.begin();
      const typename Policy::member_type end = m_policy.end();

      #pragma omp target teams distribute parallel for num_threads(128) map(to:this->m_functor)
      for(int i=begin; i<end; i++)
        m_functor(TagType(),i);
    }

  inline
  ParallelFor( const FunctorType & arg_functor
             , Policy arg_policy )
    : m_functor( arg_functor )
    , m_policy(  arg_policy )
    {}
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template<class FunctorType, class PolicyType, class ReducerType, class PointerType, class ValueType, int FunctorHasJoin, int UseReducerType>
struct ParallelReduceSpecialize {
  static inline void execute(const FunctorType& f, const PolicyType& p , PointerType result_ptr) {
    printf("Error: Invalid Specialization %i %i\n",FunctorHasJoin,UseReducerType);
  }
};

template<class FunctorType, class ReducerType, class PointerType, class ValueType, class ... PolicyArgs>
struct ParallelReduceSpecialize<FunctorType, Kokkos::RangePolicy<PolicyArgs...>, ReducerType, PointerType, ValueType, 0,0> {
  typedef Kokkos::RangePolicy<PolicyArgs...> PolicyType;
  template< class TagType >
  inline static
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  execute_impl(const FunctorType& f, const PolicyType& p, PointerType result_ptr)
    {
      OpenMPTargetExec::verify_is_process("Kokkos::Experimental::OpenMPTarget parallel_for");
      OpenMPTargetExec::verify_initialized("Kokkos::Experimental::OpenMPTarget parallel_for");
      const typename PolicyType::member_type begin = p.begin();
      const typename PolicyType::member_type end = p.end();
      
      ValueType result = ValueType();
      #pragma omp target teams distribute parallel for num_teams(512) map(to:f) map(tofrom:result) reduction(+: result)
      for(int i=begin; i<end; i++)
        f(i,result);

      *result_ptr=result;
    }


  template< class TagType >
  inline static
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  execute_impl(const FunctorType& f, const PolicyType& p, PointerType result_ptr)
    {
      OpenMPTargetExec::verify_is_process("Kokkos::Experimental::OpenMPTarget parallel_for");
      OpenMPTargetExec::verify_initialized("Kokkos::Experimental::OpenMPTarget parallel_for");
      const typename PolicyType::member_type begin = p.begin();
      const typename PolicyType::member_type end = p.end();

      ValueType result = ValueType();
      #pragma omp target teams distribute parallel for num_teams(512) map(to:f) map(tofrom: result) reduction(+: result)
      for(int i=begin; i<end; i++)
        f(TagType(),i,result);
      
      *result_ptr=result;
    }


    inline static
    void execute(const FunctorType& f, const PolicyType& p, PointerType ptr) {
      execute_impl<typename PolicyType::work_tag>(f,p,ptr);
    }
};
/*
template<class FunctorType, class PolicyType, class ReducerType, class PointerType, class ValueType>
struct ParallelReduceSpecialize<FunctorType, PolicyType, ReducerType, PointerType, ValueType, 0,1> {

  #pragma omp declare reduction(custom: ValueType : ReducerType::join(omp_out, omp_in)) initializer ( ReducerType::init(omp_priv) )

  template< class TagType >
  inline static
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  execute_impl(const FunctorType& f, const PolicyType& p, PointerType result_ptr)
    {
      OpenMPTargetExec::verify_is_process("Kokkos::Experimental::OpenMPTarget parallel_for");
      OpenMPTargetExec::verify_initialized("Kokkos::Experimental::OpenMPTarget parallel_for");
      const typename PolicyType::member_type begin = p.begin();
      const typename PolicyType::member_type end = p.end();

      ValueType result = ValueType();
      #pragma omp target teams distribute parallel for num_teams(512) map(to:f) map(tofrom:result) reduction(custom: result)
      for(int i=begin; i<end; i++)
        f(i,result);

      *result_ptr=result;
    }


  template< class TagType >
  inline static
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  execute_impl(const FunctorType& f, const PolicyType& p, PointerType result_ptr)
    {
      OpenMPTargetExec::verify_is_process("Kokkos::Experimental::OpenMPTarget parallel_for");
      OpenMPTargetExec::verify_initialized("Kokkos::Experimental::OpenMPTarget parallel_for");
      const typename PolicyType::member_type begin = p.begin();
      const typename PolicyType::member_type end = p.end();

      ValueType result = ValueType();
      #pragma omp target teams distribute parallel for num_teams(512) map(to:f) map(tofrom: result) reduction(custom: result)
      for(int i=begin; i<end; i++)
        f(TagType(),i,result);

      *result_ptr=result;
    }


    inline static
    void execute(const FunctorType& f, const PolicyType& p, PointerType ptr) {
      execute_impl<typename PolicyType::work_tag>(f,p,ptr);
    }
};
*/

template< class FunctorType , class ReducerType, class ... Traits >
class ParallelReduce< FunctorType
                    , Kokkos::RangePolicy< Traits ...>
                    , ReducerType
                    , Kokkos::Experimental::OpenMPTarget
                    >
{
private:

  typedef Kokkos::RangePolicy< Traits ... > Policy ;

  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::WorkRange    WorkRange ;
  typedef typename Policy::member_type  Member ;

  typedef Kokkos::Impl::if_c< std::is_same<InvalidType,ReducerType>::value, FunctorType, ReducerType> ReducerConditional;
  typedef typename ReducerConditional::type ReducerTypeFwd;
  typedef typename Kokkos::Impl::if_c< std::is_same<InvalidType,ReducerType>::value, WorkTag, void>::type WorkTagFwd;

  // Static Assert WorkTag void if ReducerType not InvalidType

  typedef Kokkos::Impl::FunctorValueTraits< ReducerTypeFwd , WorkTagFwd > ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   ReducerTypeFwd , WorkTagFwd > ValueInit ;
  typedef Kokkos::Impl::FunctorValueJoin<   ReducerTypeFwd , WorkTagFwd > ValueJoin ;

  enum {HasJoin = ReduceFunctorHasJoin<FunctorType>::value };
  enum {UseReducer = is_reducer_type<ReducerType>::value };

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;
  
  typedef ParallelReduceSpecialize<FunctorType,Policy,ReducerType,pointer_type,typename ValueTraits::value_type,HasJoin,UseReducer> ParForSpecialize;

  const FunctorType   m_functor ;
  const Policy        m_policy ;
  const ReducerType   m_reducer ;
  const pointer_type  m_result_ptr ;

public: 
  inline void execute() const {
    ParForSpecialize::execute(m_functor,m_policy,m_result_ptr);    
  }

  template< class ViewType >
  inline
  ParallelReduce( const FunctorType & arg_functor
                , Policy       arg_policy
                , const ViewType    & arg_result_view
                , typename std::enable_if<
                           Kokkos::is_view< ViewType >::value &&
                           !Kokkos::is_reducer_type<ReducerType>::value
                  ,void*>::type = NULL)
    : m_functor( arg_functor )
    , m_policy(  arg_policy )
    , m_reducer( InvalidType() )
    , m_result_ptr(  arg_result_view.data() )
    {
      /*static_assert( std::is_same< typename ViewType::memory_space
                                      , Kokkos::HostSpace >::value
        , "Reduction result on Kokkos::Experimental::OpenMPTarget must be a Kokkos::View in HostSpace" );*/
    }

  inline
  ParallelReduce( const FunctorType & arg_functor
                , Policy       arg_policy
                , const ReducerType& reducer )
    : m_functor( arg_functor )
    , m_policy(  arg_policy )
    , m_reducer( reducer )
    , m_result_ptr(  reducer.result_view().data() )
    {
      /*static_assert( std::is_same< typename ViewType::memory_space
                                      , Kokkos::HostSpace >::value
        , "Reduction result on Kokkos::Experimental::OpenMPTarget must be a Kokkos::View in HostSpace" );*/
    }

};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class ... Traits >
class ParallelScan< FunctorType
                  , Kokkos::RangePolicy< Traits ... >
                  , Kokkos::Experimental::OpenMPTarget
                  >
{
private:

  typedef Kokkos::RangePolicy< Traits ... > Policy ;

  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::WorkRange    WorkRange ;
  typedef typename Policy::member_type  Member ;

  typedef Kokkos::Impl::FunctorValueTraits< FunctorType, WorkTag > ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   FunctorType, WorkTag > ValueInit ;
  typedef Kokkos::Impl::FunctorValueJoin<   FunctorType, WorkTag > ValueJoin ;
  typedef Kokkos::Impl::FunctorValueOps<    FunctorType, WorkTag > ValueOps ;

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;

  const FunctorType   m_functor ;
  const Policy        m_policy ;
/*
  template< class TagType >
  inline static
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_range( const FunctorType & functor
            , const Member ibeg , const Member iend
            , reference_type update , const bool final )
    {
      #ifdef KOKKOS_OPT_RANGE_AGGRESSIVE_VECTORIZATION
      #ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
      #pragma ivdep
      #endif
      #endif
      for ( Member iwork = ibeg ; iwork < iend ; ++iwork ) {
        functor( iwork , update , final );
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
      #ifdef KOKKOS_OPT_RANGE_AGGRESSIVE_VECTORIZATION
      #ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
      #pragma ivdep
      #endif
      #endif
      for ( Member iwork = ibeg ; iwork < iend ; ++iwork ) {
        functor( t , iwork , update , final );
      }
    }
*/
public:

  inline
  void execute() const
    {
/*      OpenMPTargetExec::verify_is_process("Kokkos::Experimental::OpenMPTarget parallel_scan");
      OpenMPTargetExec::verify_initialized("Kokkos::Experimental::OpenMPTarget parallel_scan");

      OpenMPTargetExec::resize_scratch( 2 * ValueTraits::value_size( m_functor ) , 0 );

#pragma omp parallel
      {
        OpenMPTargetExec & exec = * OpenMPTargetExec::get_thread_omp();
        const WorkRange range( m_policy, exec.pool_rank(), exec.pool_size() );
        const pointer_type ptr =
          pointer_type( exec.scratch_reduce() ) +
          ValueTraits::value_count( m_functor );
        ParallelScan::template exec_range< WorkTag >
          ( m_functor , range.begin() , range.end()
          , ValueInit::init( m_functor , ptr ) , false );
      }

      {
        const unsigned thread_count = OpenMPTargetExec::pool_size();
        const unsigned value_count  = ValueTraits::value_count( m_functor );

        pointer_type ptr_prev = 0 ;

        for ( unsigned rank_rev = thread_count ; rank_rev-- ; ) {

          pointer_type ptr = pointer_type( OpenMPTargetExec::pool_rev(rank_rev)->scratch_reduce() );

          if ( ptr_prev ) {
            for ( unsigned i = 0 ; i < value_count ; ++i ) { ptr[i] = ptr_prev[ i + value_count ] ; }
            ValueJoin::join( m_functor , ptr + value_count , ptr );
          }
          else {
            ValueInit::init( m_functor , ptr );
          }

          ptr_prev = ptr ;
        }
      }

#pragma omp parallel
      {
        OpenMPTargetExec & exec = * OpenMPTargetExec::get_thread_omp();
        const WorkRange range( m_policy, exec.pool_rank(), exec.pool_size() );
        const pointer_type ptr = pointer_type( exec.scratch_reduce() );
        ParallelScan::template exec_range< WorkTag >
          ( m_functor , range.begin() , range.end()
          , ValueOps::reference( ptr ) , true );
      }
*/
    }

  //----------------------------------------

  inline
  ParallelScan( const FunctorType & arg_functor
              , const Policy      & arg_policy )
    : m_functor( arg_functor )
    , m_policy(  arg_policy )
  {}

  //----------------------------------------
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class ... Properties >
class ParallelFor< FunctorType
                 , Kokkos::TeamPolicy< Properties ... >
                 , Kokkos::Experimental::OpenMPTarget
                 >
{
private:

  typedef Kokkos::Impl::TeamPolicyInternal< Kokkos::Experimental::OpenMPTarget, Properties ... > Policy ;
  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::member_type  Member ;

  const FunctorType  m_functor ;
  const Policy       m_policy ;
  const int          m_shmem_size ;

public:

  inline void execute() const {
    OpenMPTargetExec::verify_is_process("Kokkos::Experimental::OpenMPTarget parallel_for");
    OpenMPTargetExec::verify_initialized("Kokkos::Experimental::OpenMPTarget parallel_for");
    execute_impl<WorkTag>();
  }

private:
  template< class TagType >
  inline
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  execute_impl() const
    {
      OpenMPTargetExec::verify_is_process("Kokkos::Experimental::OpenMPTarget parallel_for");
      OpenMPTargetExec::verify_initialized("Kokkos::Experimental::OpenMPTarget parallel_for");
      const int league_size = m_policy.league_size();
      const int team_size = m_policy.team_size();
      const int vector_length = m_policy.vector_length();
      const int nteams = OpenMPTargetExec::MAX_ACTIVE_TEAMS<league_size?OpenMPTargetExec::MAX_ACTIVE_TEAMS:league_size;

      OpenMPTargetExec::resize_scratch(0,Policy::member_type::TEAM_REDUCE_SIZE,0,0);
      void* scratch_ptr = OpenMPTargetExec::get_scratch_ptr();

      #pragma omp target teams distribute parallel for num_teams(league_size) num_threads(team_size*vector_length) schedule(static,1) \
          map(to:this->m_functor,scratch_ptr) 
      for(int i=0 ; i<league_size*team_size*vector_length ; i++) {
        typename Policy::member_type team(i/(team_size*vector_length),league_size,team_size,vector_length, scratch_ptr, 0,0);
        m_functor(team);
      }
    }


  template< class TagType >
  inline
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  execute_impl() const
    {
      OpenMPTargetExec::verify_is_process("Kokkos::Experimental::OpenMPTarget parallel_for");
      OpenMPTargetExec::verify_initialized("Kokkos::Experimental::OpenMPTarget parallel_for");
      const int league_size = m_policy.league_size();
      const int team_size = m_policy.team_size();
      const int vector_length = m_policy.vector_length();
      const int nteams = OpenMPTargetExec::MAX_ACTIVE_TEAMS<league_size?OpenMPTargetExec::MAX_ACTIVE_TEAMS:league_size;

      OpenMPTargetExec::resize_scratch(0,Policy::member_type::TEAM_REDUCE_SIZE,0,0);
      void* scratch_ptr = OpenMPTargetExec::get_scratch_ptr();
      #pragma omp target teams distribute parallel for num_teams(league_size) num_threads(team_size*vector_length) schedule(static,1) \
         map(to:this->m_functor,scratch_ptr)
      for(int i=0 ; i<league_size ; i++) {
        typename Policy::member_type team(i/(team_size*vector_length),league_size,team_size,vector_length, scratch_ptr, 0,0);
        m_functor(TagType(), team);
      }
    }

public:

  inline
  ParallelFor( const FunctorType & arg_functor ,
               const Policy      & arg_policy )
    : m_functor( arg_functor )
    , m_policy(  arg_policy )
    , m_shmem_size( arg_policy.scratch_size(0) + arg_policy.scratch_size(1) + FunctorTeamShmemSize< FunctorType >::value( arg_functor , arg_policy.team_size() ) )
    {}
};

template<class FunctorType, class ReducerType, class PointerType, class ValueType, class ... PolicyArgs>
struct ParallelReduceSpecialize<FunctorType, TeamPolicyInternal<PolicyArgs...>, ReducerType, PointerType, ValueType, 0,0> {
  typedef TeamPolicyInternal<PolicyArgs...> PolicyType;

  template< class TagType >
  inline static
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  execute_impl(const FunctorType& f, const PolicyType& p, PointerType result_ptr)
    {
      OpenMPTargetExec::verify_is_process("Kokkos::Experimental::OpenMPTarget parallel_for");
      OpenMPTargetExec::verify_initialized("Kokkos::Experimental::OpenMPTarget parallel_for");
      
      const int league_size = p.league_size();
      const int team_size = p.team_size();
      const int vector_length = p.vector_length();
      const int nteams = OpenMPTargetExec::MAX_ACTIVE_TEAMS<league_size?OpenMPTargetExec::MAX_ACTIVE_TEAMS:league_size;
      
      OpenMPTargetExec::resize_scratch(0,PolicyType::member_type::TEAM_REDUCE_SIZE,0,0);
      void* scratch_ptr = OpenMPTargetExec::get_scratch_ptr(); 

      ValueType result = ValueType();
      #pragma omp target teams distribute parallel for num_teams(nteams) num_threads(team_size*vector_length) \
         map(to:f,scratch_ptr) map(tofrom:result) reduction(+: result) schedule(static,1)
      for(int i=0 ; i<league_size*team_size*vector_length ; i++) {
        typename PolicyType::member_type team(i/(team_size*vector_length),league_size,team_size,vector_length, scratch_ptr, 0,0);
        f(team,result);
        if(team.m_vector_lane!=0) result = 0;
      }

      *result_ptr=result;
    }


  template< class TagType >
  inline static
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  execute_impl(const FunctorType& f, const PolicyType& p, PointerType result_ptr)
    {
      OpenMPTargetExec::verify_is_process("Kokkos::Experimental::OpenMPTarget parallel_for");
      OpenMPTargetExec::verify_initialized("Kokkos::Experimental::OpenMPTarget parallel_for");

      const int league_size = p.league_size();
      const int team_size = p.team_size();
      const int vector_length = p.vector_length();
      const int nteams = OpenMPTargetExec::MAX_ACTIVE_TEAMS<league_size?OpenMPTargetExec::MAX_ACTIVE_TEAMS:league_size;

      OpenMPTargetExec::resize_scratch(0,PolicyType::member_type::TEAM_REDUCE_SIZE,0,0);
      void* scratch_ptr = OpenMPTargetExec::get_scratch_ptr();

      ValueType result = ValueType();
      #pragma omp target teams distribute parallel for num_teams(nteams) num_threads(team_size*vector_length) \
         map(to:f,scratch_ptr) map(tofrom:result) reduction(+: result) schedule(static,1)
      for(int i=0 ; i<league_size*team_size*vector_length ; i++) {
        typename PolicyType::member_type team(i/(team_size*vector_length),league_size,team_size,vector_length, scratch_ptr, 0,0);
        f(TagType(),team,result);
        if(team.vector_lane!=0) result = 0;
      }
      *result_ptr=result;
    }


    inline static
    void execute(const FunctorType& f, const PolicyType& p, PointerType ptr) {
      execute_impl<typename PolicyType::work_tag>(f,p,ptr);
    }
};


template< class FunctorType , class ReducerType, class ... Properties >
class ParallelReduce< FunctorType
                    , Kokkos::TeamPolicy< Properties ... >
                    , ReducerType
                    , Kokkos::Experimental::OpenMPTarget
                    >
{
private:

  typedef Kokkos::Impl::TeamPolicyInternal< Kokkos::Experimental::OpenMPTarget, Properties ... >         Policy ;

  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::member_type  Member ;

  typedef Kokkos::Impl::if_c< std::is_same<InvalidType,ReducerType>::value, FunctorType, ReducerType> ReducerConditional;
  typedef typename ReducerConditional::type ReducerTypeFwd;
  typedef typename Kokkos::Impl::if_c< std::is_same<InvalidType,ReducerType>::value, WorkTag, void>::type WorkTagFwd;

  typedef Kokkos::Impl::FunctorValueTraits< ReducerTypeFwd , WorkTagFwd >  ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   ReducerTypeFwd , WorkTagFwd >  ValueInit ;
  typedef Kokkos::Impl::FunctorValueJoin<   ReducerTypeFwd , WorkTagFwd >  ValueJoin ;

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;
  typedef typename ValueTraits::value_type      value_type ;

  enum {HasJoin = ReduceFunctorHasJoin<FunctorType>::value };
  enum {UseReducer = is_reducer_type<ReducerType>::value };

  typedef ParallelReduceSpecialize<FunctorType,Policy,ReducerType,pointer_type,typename ValueTraits::value_type,HasJoin,UseReducer> ParForSpecialize;

  const FunctorType  m_functor ;
  const Policy       m_policy ;
  const ReducerType  m_reducer ;
  const pointer_type m_result_ptr ;
  const int          m_shmem_size ;

public:

  inline
  void execute() const {
    ParForSpecialize::execute(m_functor,m_policy,m_result_ptr);   
  }

  template< class ViewType >
  inline
  ParallelReduce( const FunctorType  & arg_functor ,
                  const Policy       & arg_policy ,
                  const ViewType     & arg_result ,
                  typename std::enable_if<
                    Kokkos::is_view< ViewType >::value &&
                    !Kokkos::is_reducer_type<ReducerType>::value
                    ,void*>::type = NULL)
    : m_functor( arg_functor )
    , m_policy(  arg_policy )
    , m_reducer( InvalidType() )
    , m_result_ptr( arg_result.ptr_on_device() )
    , m_shmem_size( arg_policy.scratch_size(0) + arg_policy.scratch_size(1) + FunctorTeamShmemSize< FunctorType >::value( arg_functor , arg_policy.team_size() ) )
    {}

  inline
  ParallelReduce( const FunctorType & arg_functor
    , Policy       arg_policy
    , const ReducerType& reducer )
  : m_functor( arg_functor )
  , m_policy(  arg_policy )
  , m_reducer( reducer )
  , m_result_ptr(  reducer.result_view().data() )
  , m_shmem_size( arg_policy.scratch_size(0) + arg_policy.scratch_size(1) + FunctorTeamShmemSize< FunctorType >::value( arg_functor , arg_policy.team_size() ) )
  {
  /*static_assert( std::is_same< typename ViewType::memory_space
                          , Kokkos::HostSpace >::value
  , "Reduction result on Kokkos::Experimental::OpenMPTarget must be a Kokkos::View in HostSpace" );*/
  }

};

} // namespace Impl
} // namespace Kokkos


namespace Kokkos {
namespace Impl {

  template<typename iType>
  struct TeamThreadRangeBoundariesStruct<iType,OpenMPTargetExecTeamMember> {
    typedef iType index_type;
    const iType start;
    const iType end;
    const iType increment;

    inline
    TeamThreadRangeBoundariesStruct (const OpenMPTargetExecTeamMember& thread_, const iType& count):
      start( thread_.team_rank() ),
      end( count ),
      increment( thread_.team_size() )
    {}
    inline
    TeamThreadRangeBoundariesStruct (const OpenMPTargetExecTeamMember& thread_, const iType& begin_, const iType& end_):
      start( begin_+thread_.team_rank() ),
      end( end_ ),
      increment( thread_.team_size() )
    {}
  };

  template<typename iType>
  struct ThreadVectorRangeBoundariesStruct<iType,OpenMPTargetExecTeamMember> {
    typedef iType index_type;
    const index_type start;
    const index_type end;
    const index_type increment;

    inline
    ThreadVectorRangeBoundariesStruct (const OpenMPTargetExecTeamMember& thread_, const index_type& count):
      start( thread_.m_vector_lane ),
      end( count ),
      increment( thread_.m_vector_length )
    {}
    inline
    ThreadVectorRangeBoundariesStruct (const OpenMPTargetExecTeamMember& thread_, const index_type& begin_, const index_type& end_):
      start( begin_+thread_.m_vector_lane ),
      end( end_ ),
      increment( thread_.m_vector_length )
    {}
  };

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

  template<typename iType>
  KOKKOS_INLINE_FUNCTION
  Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::OpenMPTargetExecTeamMember>
    ThreadVectorRange(const Impl::OpenMPTargetExecTeamMember& thread, const iType& begin, const iType& end) {
    return Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::OpenMPTargetExecTeamMember>(thread,begin,end);
  }

}

}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* KOKKOS_OPENMPTARGET_PARALLEL_HPP */

