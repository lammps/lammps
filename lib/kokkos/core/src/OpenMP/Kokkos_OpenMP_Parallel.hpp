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

#ifndef KOKKOS_OPENMP_PARALLEL_HPP
#define KOKKOS_OPENMP_PARALLEL_HPP

#include <omp.h>
#include <iostream>
#include <OpenMP/Kokkos_OpenMPexec.hpp>
#include <impl/Kokkos_FunctorAdapter.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class ... Traits >
class ParallelFor< FunctorType
                 , Kokkos::RangePolicy< Traits ... >
                 , Kokkos::OpenMP 
                 >
{
private:

  typedef Kokkos::RangePolicy< Traits ...  > Policy ;
  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::WorkRange    WorkRange ;
  typedef typename Policy::member_type  Member ;

  const FunctorType m_functor ;
  const Policy      m_policy ;

  template< class TagType >
  inline static
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_range( const FunctorType & functor
            , const Member ibeg , const Member iend )
    {
      #ifdef KOKKOS_OPT_RANGE_AGGRESSIVE_VECTORIZATION
      #ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
      #pragma ivdep
      #endif
      #endif
      for ( Member iwork = ibeg ; iwork < iend ; ++iwork ) {
        functor( iwork );
      }
    }

  template< class TagType >
  inline static
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec_range( const FunctorType & functor
            , const Member ibeg , const Member iend )
    {
      const TagType t{} ;
      #ifdef KOKKOS_OPT_RANGE_AGGRESSIVE_VECTORIZATION
      #ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
      #pragma ivdep
      #endif
      #endif
      for ( Member iwork = ibeg ; iwork < iend ; ++iwork ) {
        functor( t , iwork );
      }
    }

public:

  inline void execute() const
    {
      enum { is_dynamic = std::is_same< typename Policy::schedule_type::type
                                      , Kokkos::Dynamic >::value };

      OpenMPexec::verify_is_process("Kokkos::OpenMP parallel_for");
      OpenMPexec::verify_initialized("Kokkos::OpenMP parallel_for");

#pragma omp parallel
      {
        HostThreadTeamData & data = *OpenMPexec::get_thread_data();

        data.set_work_partition( m_policy.end() - m_policy.begin()
                               , m_policy.chunk_size() );

        if ( is_dynamic ) {
          // Make sure work partition is set before stealing
          if ( data.pool_rendezvous() ) data.pool_rendezvous_release();
        }

        std::pair<int64_t,int64_t> range(0,0);

        do {

          range = is_dynamic ? data.get_work_stealing_chunk()
                             : data.get_work_partition();

          ParallelFor::template
            exec_range< WorkTag >( m_functor
                                 , range.first  + m_policy.begin()
                                 , range.second + m_policy.begin() );

        } while ( is_dynamic && 0 <= range.first );
      }
      // END #pragma omp parallel
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

template< class FunctorType , class ReducerType, class ... Traits >
class ParallelReduce< FunctorType
                    , Kokkos::RangePolicy< Traits ...>
                    , ReducerType
                    , Kokkos::OpenMP
                    >
{
private:

  typedef Kokkos::RangePolicy< Traits ... > Policy ;

  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::WorkRange    WorkRange ;
  typedef typename Policy::member_type  Member ;

  typedef FunctorAnalysis< FunctorPatternInterface::REDUCE , Policy , FunctorType > Analysis ;

  typedef Kokkos::Impl::if_c< std::is_same<InvalidType,ReducerType>::value, FunctorType, ReducerType> ReducerConditional;
  typedef typename ReducerConditional::type ReducerTypeFwd;

  // Static Assert WorkTag void if ReducerType not InvalidType

  typedef Kokkos::Impl::FunctorValueInit<   ReducerTypeFwd, WorkTag > ValueInit ;
  typedef Kokkos::Impl::FunctorValueJoin<   ReducerTypeFwd, WorkTag > ValueJoin ;

  typedef typename Analysis::pointer_type    pointer_type ;
  typedef typename Analysis::reference_type  reference_type ;

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
      #ifdef KOKKOS_OPT_RANGE_AGGRESSIVE_VECTORIZATION
      #ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
      #pragma ivdep
      #endif
      #endif
      for ( Member iwork = ibeg ; iwork < iend ; ++iwork ) {
        functor( iwork , update );
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
      #ifdef KOKKOS_OPT_RANGE_AGGRESSIVE_VECTORIZATION
      #ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
      #pragma ivdep
      #endif
      #endif
      for ( Member iwork = ibeg ; iwork < iend ; ++iwork ) {
        functor( t , iwork , update );
      }
    }

public:

  inline void execute() const
    {
      enum { is_dynamic = std::is_same< typename Policy::schedule_type::type
                                      , Kokkos::Dynamic >::value };

      OpenMPexec::verify_is_process("Kokkos::OpenMP parallel_for");
      OpenMPexec::verify_initialized("Kokkos::OpenMP parallel_for");

      const size_t pool_reduce_bytes =
        Analysis::value_size( ReducerConditional::select(m_functor, m_reducer));

      OpenMPexec::resize_thread_data( pool_reduce_bytes
                                    , 0 // team_reduce_bytes
                                    , 0 // team_shared_bytes
                                    , 0 // thread_local_bytes
                                    );

#pragma omp parallel
      {
        HostThreadTeamData & data = *OpenMPexec::get_thread_data();

        data.set_work_partition( m_policy.end() - m_policy.begin()
                               , m_policy.chunk_size() );

        if ( is_dynamic ) {
          // Make sure work partition is set before stealing
          if ( data.pool_rendezvous() ) data.pool_rendezvous_release();
        }

        reference_type update =
          ValueInit::init( ReducerConditional::select(m_functor , m_reducer)
                         , data.pool_reduce_local() );

        std::pair<int64_t,int64_t> range(0,0);

        do {

          range = is_dynamic ? data.get_work_stealing_chunk()
                             : data.get_work_partition();

          ParallelReduce::template
            exec_range< WorkTag >( m_functor
                                 , range.first  + m_policy.begin()
                                 , range.second + m_policy.begin()
                                 , update );

        } while ( is_dynamic && 0 <= range.first );
      }
// END #pragma omp parallel

      // Reduction:

      const pointer_type ptr = pointer_type( OpenMPexec::get_thread_data(0)->pool_reduce_local() );

      for ( int i = 1 ; i < OpenMPexec::pool_size() ; ++i ) {
        ValueJoin::join( ReducerConditional::select(m_functor , m_reducer)
                       , ptr
                       , OpenMPexec::get_thread_data(i)->pool_reduce_local() );
      }

      Kokkos::Impl::FunctorFinal<  ReducerTypeFwd , WorkTag >::final( ReducerConditional::select(m_functor , m_reducer) , ptr );

      if ( m_result_ptr ) {
        const int n = Analysis::value_count( ReducerConditional::select(m_functor , m_reducer) );

        for ( int j = 0 ; j < n ; ++j ) { m_result_ptr[j] = ptr[j] ; }
      }
    }

  //----------------------------------------

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
        , "Reduction result on Kokkos::OpenMP must be a Kokkos::View in HostSpace" );*/
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
        , "Reduction result on Kokkos::OpenMP must be a Kokkos::View in HostSpace" );*/
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
                  , Kokkos::OpenMP
                  >
{
private:

  typedef Kokkos::RangePolicy< Traits ... > Policy ;

  typedef FunctorAnalysis< FunctorPatternInterface::SCAN , Policy , FunctorType > Analysis ;

  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::WorkRange    WorkRange ;
  typedef typename Policy::member_type  Member ;

  typedef Kokkos::Impl::FunctorValueInit<   FunctorType, WorkTag > ValueInit ;
  typedef Kokkos::Impl::FunctorValueJoin<   FunctorType, WorkTag > ValueJoin ;
  typedef Kokkos::Impl::FunctorValueOps<    FunctorType, WorkTag > ValueOps ;

  typedef typename Analysis::pointer_type    pointer_type ;
  typedef typename Analysis::reference_type  reference_type ;

  const FunctorType   m_functor ;
  const Policy        m_policy ;

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

public:

  inline
  void execute() const
    {
      OpenMPexec::verify_is_process("Kokkos::OpenMP parallel_scan");
      OpenMPexec::verify_initialized("Kokkos::OpenMP parallel_scan");

      const int    value_count       = Analysis::value_count( m_functor );
      const size_t pool_reduce_bytes = 2 * Analysis::value_size( m_functor );

      OpenMPexec::resize_thread_data( pool_reduce_bytes
                                    , 0 // team_reduce_bytes
                                    , 0 // team_shared_bytes
                                    , 0 // thread_local_bytes
                                    );

#pragma omp parallel
      {
        HostThreadTeamData & data = *OpenMPexec::get_thread_data();

        const WorkRange range( m_policy, data.pool_rank(), data.pool_size() );

        reference_type update_sum =
          ValueInit::init( m_functor , data.pool_reduce_local() );

        ParallelScan::template exec_range< WorkTag >
          ( m_functor , range.begin() , range.end() , update_sum , false );

        if ( data.pool_rendezvous() ) {

          pointer_type ptr_prev = 0 ;

          const int n = data.pool_size();

          for ( int i = 0 ; i < n ; ++i ) {

            pointer_type ptr = (pointer_type)
              data.pool_member(i)->pool_reduce_local();

            if ( i ) {
              for ( int j = 0 ; j < value_count ; ++j ) {
                ptr[j+value_count] = ptr_prev[j+value_count] ;
              }
              ValueJoin::join( m_functor , ptr + value_count , ptr_prev );
            }
            else {
              ValueInit::init( m_functor , ptr + value_count );
            }

            ptr_prev = ptr ;
          }

          data.pool_rendezvous_release();
        }

        reference_type update_base =
          ValueOps::reference
            ( ((pointer_type)data.pool_reduce_local()) + value_count );

        ParallelScan::template exec_range< WorkTag >
          ( m_functor , range.begin() , range.end() , update_base , true );
      }
/* END #pragma omp parallel */

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
                 , Kokkos::OpenMP
                 >
{
private:

  enum { TEAM_REDUCE_SIZE = 512 };

  typedef Kokkos::Impl::TeamPolicyInternal< Kokkos::OpenMP, Properties ... > Policy ;
  typedef typename Policy::work_tag             WorkTag ;
  typedef typename Policy::schedule_type::type  SchedTag ;
  typedef typename Policy::member_type          Member ;

  const FunctorType  m_functor ;
  const Policy       m_policy ;
  const int          m_shmem_size ;

  template< class TagType >
  inline static
  typename std::enable_if< ( std::is_same< TagType , void >::value ) >::type
  exec_team( const FunctorType & functor
           , HostThreadTeamData & data
           , const int league_rank_begin
           , const int league_rank_end
           , const int league_size )
    {
      for ( int r = league_rank_begin ; r < league_rank_end ; ) {

        functor( Member( data, r , league_size ) );

        if ( ++r < league_rank_end ) {
          // Don't allow team members to lap one another
          // so that they don't overwrite shared memory.
          if ( data.team_rendezvous() ) { data.team_rendezvous_release(); }
        }
      }
    }


  template< class TagType >
  inline static
  typename std::enable_if< ( ! std::is_same< TagType , void >::value ) >::type
  exec_team( const FunctorType & functor
           , HostThreadTeamData & data
           , const int league_rank_begin
           , const int league_rank_end
           , const int league_size )
    {
      const TagType t{};

      for ( int r = league_rank_begin ; r < league_rank_end ; ) {

        functor( t , Member( data, r , league_size ) );

        if ( ++r < league_rank_end ) {
          // Don't allow team members to lap one another
          // so that they don't overwrite shared memory.
          if ( data.team_rendezvous() ) { data.team_rendezvous_release(); }
        }
      }
    }

public:

  inline
  void execute() const
    {
      enum { is_dynamic = std::is_same< SchedTag , Kokkos::Dynamic >::value };

      OpenMPexec::verify_is_process("Kokkos::OpenMP parallel_for");
      OpenMPexec::verify_initialized("Kokkos::OpenMP parallel_for");

      const size_t pool_reduce_size = 0 ; // Never shrinks
      const size_t team_reduce_size = TEAM_REDUCE_SIZE * m_policy.team_size();
      const size_t team_shared_size = m_shmem_size + m_policy.scratch_size(1);
      const size_t thread_local_size = 0 ; // Never shrinks

      OpenMPexec::resize_thread_data( pool_reduce_size
                                    , team_reduce_size
                                    , team_shared_size
                                    , thread_local_size );

#pragma omp parallel
      {
        HostThreadTeamData & data = *OpenMPexec::get_thread_data();

        const int active = data.organize_team( m_policy.team_size() );

        if ( active ) {
          data.set_work_partition( m_policy.league_size()
                                 , ( 0 < m_policy.chunk_size()
                                   ? m_policy.chunk_size()
                                   : m_policy.team_iter() ) );
        }

        if ( is_dynamic ) {
          // Must synchronize to make sure each team has set its
          // partition before begining the work stealing loop.
          if ( data.pool_rendezvous() ) data.pool_rendezvous_release();
        }

        if ( active ) {

          std::pair<int64_t,int64_t> range(0,0);

          do {

            range = is_dynamic ? data.get_work_stealing_chunk()
                               : data.get_work_partition();

            ParallelFor::template exec_team< WorkTag >
              ( m_functor , data
              , range.first , range.second , m_policy.league_size() );

          } while ( is_dynamic && 0 <= range.first );
        }

        data.disband_team();
      }
// END #pragma omp parallel
    }


  inline
  ParallelFor( const FunctorType & arg_functor ,
               const Policy      & arg_policy )
    : m_functor( arg_functor )
    , m_policy(  arg_policy )
    , m_shmem_size( arg_policy.scratch_size(0) +
                    arg_policy.scratch_size(1) +
                    FunctorTeamShmemSize< FunctorType >
                      ::value( arg_functor , arg_policy.team_size() ) )
    {}
};

//----------------------------------------------------------------------------

template< class FunctorType , class ReducerType, class ... Properties >
class ParallelReduce< FunctorType
                    , Kokkos::TeamPolicy< Properties ... >
                    , ReducerType
                    , Kokkos::OpenMP
                    >
{
private:

  enum { TEAM_REDUCE_SIZE = 512 };

  typedef Kokkos::Impl::TeamPolicyInternal< Kokkos::OpenMP, Properties ... >         Policy ;

  typedef FunctorAnalysis< FunctorPatternInterface::REDUCE , Policy , FunctorType > Analysis ;

  typedef typename Policy::work_tag             WorkTag ;
  typedef typename Policy::schedule_type::type  SchedTag ;
  typedef typename Policy::member_type          Member ;

  typedef Kokkos::Impl::if_c< std::is_same<InvalidType,ReducerType>::value
                            , FunctorType, ReducerType> ReducerConditional;

  typedef typename ReducerConditional::type ReducerTypeFwd;

  typedef Kokkos::Impl::FunctorValueInit<   ReducerTypeFwd , WorkTag >  ValueInit ;
  typedef Kokkos::Impl::FunctorValueJoin<   ReducerTypeFwd , WorkTag >  ValueJoin ;

  typedef typename Analysis::pointer_type    pointer_type ;
  typedef typename Analysis::reference_type  reference_type ;

  const FunctorType  m_functor ;
  const Policy       m_policy ;
  const ReducerType  m_reducer ;
  const pointer_type m_result_ptr ;
  const int          m_shmem_size ;

  template< class TagType >
  inline static
  typename std::enable_if< ( std::is_same< TagType , void >::value ) >::type
  exec_team( const FunctorType & functor
           , HostThreadTeamData & data
           , reference_type     & update
           , const int league_rank_begin
           , const int league_rank_end
           , const int league_size )
    {
      for ( int r = league_rank_begin ; r < league_rank_end ; ) {

        functor( Member( data, r , league_size ) , update );

        if ( ++r < league_rank_end ) {
          // Don't allow team members to lap one another
          // so that they don't overwrite shared memory.
          if ( data.team_rendezvous() ) { data.team_rendezvous_release(); }
        }
      }
    }


  template< class TagType >
  inline static
  typename std::enable_if< ( ! std::is_same< TagType , void >::value ) >::type
  exec_team( const FunctorType & functor
           , HostThreadTeamData & data
           , reference_type     & update
           , const int league_rank_begin
           , const int league_rank_end
           , const int league_size )
    {
      const TagType t{};

      for ( int r = league_rank_begin ; r < league_rank_end ; ) {

        functor( t , Member( data, r , league_size ) , update );

        if ( ++r < league_rank_end ) {
          // Don't allow team members to lap one another
          // so that they don't overwrite shared memory.
          if ( data.team_rendezvous() ) { data.team_rendezvous_release(); }
        }
      }
    }

public:

  inline
  void execute() const
    {
      enum { is_dynamic = std::is_same< SchedTag , Kokkos::Dynamic >::value };

      OpenMPexec::verify_is_process("Kokkos::OpenMP parallel_reduce");
      OpenMPexec::verify_initialized("Kokkos::OpenMP parallel_reduce");

      const size_t pool_reduce_size =
        Analysis::value_size( ReducerConditional::select(m_functor, m_reducer));

      const size_t team_reduce_size = TEAM_REDUCE_SIZE * m_policy.team_size();
      const size_t team_shared_size = m_shmem_size + m_policy.scratch_size(1);
      const size_t thread_local_size = 0 ; // Never shrinks

      OpenMPexec::resize_thread_data( pool_reduce_size
                                    , team_reduce_size
                                    , team_shared_size
                                    , thread_local_size );

#pragma omp parallel
      {
        HostThreadTeamData & data = *OpenMPexec::get_thread_data();

        const int active = data.organize_team( m_policy.team_size() );

        if ( active ) {
          data.set_work_partition( m_policy.league_size()
                                 , ( 0 < m_policy.chunk_size()
                                   ? m_policy.chunk_size()
                                   : m_policy.team_iter() ) );
        }

        if ( is_dynamic ) {
          // Must synchronize to make sure each team has set its
          // partition before begining the work stealing loop.
          if ( data.pool_rendezvous() ) data.pool_rendezvous_release();
        }

        if ( active ) {
          reference_type update =
            ValueInit::init( ReducerConditional::select(m_functor , m_reducer)
                           , data.pool_reduce_local() );

          std::pair<int64_t,int64_t> range(0,0);

          do {

            range = is_dynamic ? data.get_work_stealing_chunk()
                               : data.get_work_partition();

            ParallelReduce::template exec_team< WorkTag >
              ( m_functor , data , update
              , range.first , range.second , m_policy.league_size() );

          } while ( is_dynamic && 0 <= range.first );
        } else {
          ValueInit::init( ReducerConditional::select(m_functor , m_reducer)
                           , data.pool_reduce_local() );
        }

        data.disband_team();
      }
// END #pragma omp parallel

      // Reduction:

      const pointer_type ptr = pointer_type( OpenMPexec::get_thread_data(0)->pool_reduce_local() );

      for ( int i = 1 ; i < OpenMPexec::pool_size() ; ++i ) {
        ValueJoin::join( ReducerConditional::select(m_functor , m_reducer)
                       , ptr
                       , OpenMPexec::get_thread_data(i)->pool_reduce_local() );
      }

      Kokkos::Impl::FunctorFinal<  ReducerTypeFwd , WorkTag >::final( ReducerConditional::select(m_functor , m_reducer) , ptr );

      if ( m_result_ptr ) {
        const int n = Analysis::value_count( ReducerConditional::select(m_functor , m_reducer) );

        for ( int j = 0 ; j < n ; ++j ) { m_result_ptr[j] = ptr[j] ; }
      }
    }

  //----------------------------------------

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
    , m_shmem_size( arg_policy.scratch_size(0) +
                    arg_policy.scratch_size(1) +
                    FunctorTeamShmemSize< FunctorType >
                      ::value( arg_functor , arg_policy.team_size() ) )
    {}

  inline
  ParallelReduce( const FunctorType & arg_functor
    , Policy       arg_policy
    , const ReducerType& reducer )
  : m_functor( arg_functor )
  , m_policy(  arg_policy )
  , m_reducer( reducer )
  , m_result_ptr(  reducer.result_view().data() )
  , m_shmem_size( arg_policy.scratch_size(0) +
                  arg_policy.scratch_size(1) +
                  FunctorTeamShmemSize< FunctorType >
                    ::value( arg_functor , arg_policy.team_size() ) )
  {
  /*static_assert( std::is_same< typename ViewType::memory_space
                          , Kokkos::HostSpace >::value
  , "Reduction result on Kokkos::OpenMP must be a Kokkos::View in HostSpace" );*/
  }

};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* KOKKOS_OPENMP_PARALLEL_HPP */

