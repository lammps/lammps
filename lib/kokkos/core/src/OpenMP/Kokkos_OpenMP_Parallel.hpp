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

#include <Kokkos_Parallel.hpp>
#include <OpenMP/Kokkos_OpenMPexec.hpp>
#include <impl/Kokkos_FunctorAdapter.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelFor< FunctorType
                 , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::OpenMP >
                 >
{
private:

  typedef Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::OpenMP > Policy ;
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
      #ifdef KOKKOS_HAVE_PRAGMA_IVDEP
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
      #ifdef KOKKOS_HAVE_PRAGMA_IVDEP
      #pragma ivdep
      #endif
      #endif
      for ( Member iwork = ibeg ; iwork < iend ; ++iwork ) {
        functor( t , iwork );
      }
    }

public:

  inline
  void execute() const
    {
      OpenMPexec::verify_is_process("Kokkos::OpenMP parallel_for");
      OpenMPexec::verify_initialized("Kokkos::OpenMP parallel_for");

#pragma omp parallel
      {
        OpenMPexec & exec = * OpenMPexec::get_thread_omp();

        const WorkRange range( m_policy, exec.pool_rank(), exec.pool_size() );

        ParallelFor::template exec_range< WorkTag >( m_functor , range.begin() , range.end() );
      }
/* END #pragma omp parallel */
    }

  inline
  ParallelFor( const FunctorType & arg_functor
             , const Policy      & arg_policy )
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

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelReduce< FunctorType
                    , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::OpenMP >
                    >
{
private:

  typedef Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::OpenMP > Policy ;

  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::WorkRange    WorkRange ;
  typedef typename Policy::member_type  Member ;

  typedef Kokkos::Impl::FunctorValueTraits< FunctorType, WorkTag > ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   FunctorType, WorkTag > ValueInit ;
  typedef Kokkos::Impl::FunctorValueJoin<   FunctorType, WorkTag > ValueJoin ;

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;

  const FunctorType   m_functor ;
  const Policy        m_policy ;
  const pointer_type  m_result_ptr ;

  template< class TagType >
  inline static
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_range( const FunctorType & functor
            , const Member ibeg , const Member iend
            , reference_type update )
    {
      #ifdef KOKKOS_OPT_RANGE_AGGRESSIVE_VECTORIZATION
      #ifdef KOKKOS_HAVE_PRAGMA_IVDEP
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
      #ifdef KOKKOS_HAVE_PRAGMA_IVDEP
      #pragma ivdep
      #endif
      #endif
      for ( Member iwork = ibeg ; iwork < iend ; ++iwork ) {
        functor( t , iwork , update );
      }
    }

public:

  inline
  void execute() const
    {
      OpenMPexec::verify_is_process("Kokkos::OpenMP parallel_reduce");
      OpenMPexec::verify_initialized("Kokkos::OpenMP parallel_reduce");

      OpenMPexec::resize_scratch( ValueTraits::value_size( m_functor ) , 0 );

#pragma omp parallel
      {
        OpenMPexec & exec = * OpenMPexec::get_thread_omp();
        const WorkRange range( m_policy, exec.pool_rank(), exec.pool_size() );
        ParallelReduce::template exec_range< WorkTag >
          ( m_functor , range.begin() , range.end()
          , ValueInit::init( m_functor , exec.scratch_reduce() ) );
      }
/* END #pragma omp parallel */

      // Reduction:

      const pointer_type ptr = pointer_type( OpenMPexec::pool_rev(0)->scratch_reduce() );

      for ( int i = 1 ; i < OpenMPexec::pool_size() ; ++i ) {
        ValueJoin::join( m_functor , ptr , OpenMPexec::pool_rev(i)->scratch_reduce() );
      }

      Kokkos::Impl::FunctorFinal<  FunctorType , WorkTag >::final( m_functor , ptr );

      if ( m_result_ptr ) {
        const int n = ValueTraits::value_count( m_functor );

        for ( int j = 0 ; j < n ; ++j ) { m_result_ptr[j] = ptr[j] ; }
      }
    }

  //----------------------------------------

  template< class ViewType >
  inline
  ParallelReduce( const FunctorType & arg_functor
                , const Policy      & arg_policy
                , const ViewType    & arg_result_view )
    : m_functor( arg_functor )
    , m_policy(  arg_policy )
    , m_result_ptr(  arg_result_view.ptr_on_device() )
    {
      static_assert( Kokkos::is_view< ViewType >::value
        , "Reduction result on Kokkos::OpenMP must be a Kokkos::View" );

      static_assert( std::is_same< typename ViewType::memory_space
                                      , Kokkos::HostSpace >::value
        , "Reduction result on Kokkos::OpenMP must be a Kokkos::View in HostSpace" );
    }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelScan< FunctorType
                  , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::OpenMP >
                  >
{
private:

  typedef Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::OpenMP > Policy ;

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

  template< class TagType >
  inline static
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_range( const FunctorType & functor
            , const Member ibeg , const Member iend
            , reference_type update , const bool final )
    {
      #ifdef KOKKOS_OPT_RANGE_AGGRESSIVE_VECTORIZATION
      #ifdef KOKKOS_HAVE_PRAGMA_IVDEP
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
      #ifdef KOKKOS_HAVE_PRAGMA_IVDEP
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

      OpenMPexec::resize_scratch( 2 * ValueTraits::value_size( m_functor ) , 0 );

#pragma omp parallel
      {
        OpenMPexec & exec = * OpenMPexec::get_thread_omp();
        const WorkRange range( m_policy, exec.pool_rank(), exec.pool_size() );
        const pointer_type ptr =
          pointer_type( exec.scratch_reduce() ) +
          ValueTraits::value_count( m_functor );
        ParallelScan::template exec_range< WorkTag >
          ( m_functor , range.begin() , range.end()
          , ValueInit::init( m_functor , ptr ) , false );
      }
/* END #pragma omp parallel */

      {
        const unsigned thread_count = OpenMPexec::pool_size();
        const unsigned value_count  = ValueTraits::value_count( m_functor );

        pointer_type ptr_prev = 0 ;

        for ( unsigned rank_rev = thread_count ; rank_rev-- ; ) {

          pointer_type ptr = pointer_type( OpenMPexec::pool_rev(rank_rev)->scratch_reduce() );

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
        OpenMPexec & exec = * OpenMPexec::get_thread_omp();
        const WorkRange range( m_policy, exec.pool_rank(), exec.pool_size() );
        const pointer_type ptr = pointer_type( exec.scratch_reduce() );
        ParallelScan::template exec_range< WorkTag >
          ( m_functor , range.begin() , range.end()
          , ValueOps::reference( ptr ) , true );
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

template< class FunctorType , class Arg0 , class Arg1 >
class ParallelFor< FunctorType
                 , Kokkos::TeamPolicy< Arg0 , Arg1 , Kokkos::OpenMP >
                 >
{
private:

  typedef Kokkos::TeamPolicy< Arg0 , Arg1 , Kokkos::OpenMP > Policy ;
  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::member_type  Member ;

  const FunctorType  m_functor ;
  const Policy       m_policy ;
  const int          m_shmem_size ;

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

public:

  inline
  void execute() const
    {
      OpenMPexec::verify_is_process("Kokkos::OpenMP parallel_for");
      OpenMPexec::verify_initialized("Kokkos::OpenMP parallel_for");

      const size_t team_reduce_size = Policy::member_type::team_reduce_size();

      OpenMPexec::resize_scratch( 0 , team_reduce_size + m_shmem_size );

#pragma omp parallel
      {
        ParallelFor::template exec_team< WorkTag >
          ( m_functor
          , Member( * OpenMPexec::get_thread_omp(), m_policy, m_shmem_size) );
      }
/* END #pragma omp parallel */
    }

  inline
  ParallelFor( const FunctorType & arg_functor ,
               const Policy      & arg_policy )
    : m_functor( arg_functor )
    , m_policy(  arg_policy )
    , m_shmem_size( arg_policy.scratch_size() + FunctorTeamShmemSize< FunctorType >::value( arg_functor , arg_policy.team_size() ) )
    {}
};


template< class FunctorType , class Arg0 , class Arg1 >
class ParallelReduce< FunctorType
                    , Kokkos::TeamPolicy< Arg0 , Arg1 , Kokkos::OpenMP >
                    >
{
private:

  typedef Kokkos::TeamPolicy< Arg0 , Arg1 , Kokkos::OpenMP >         Policy ;

  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::member_type  Member ;

  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , WorkTag >  ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   FunctorType , WorkTag >  ValueInit ;
  typedef Kokkos::Impl::FunctorValueJoin<   FunctorType , WorkTag >  ValueJoin ;

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;

  const FunctorType  m_functor ;
  const Policy       m_policy ;
  const pointer_type m_result_ptr ;
  const int          m_shmem_size ;

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

public:

  inline
  void execute() const
    {
      OpenMPexec::verify_is_process("Kokkos::OpenMP parallel_reduce");

      const size_t team_reduce_size = Policy::member_type::team_reduce_size();

      OpenMPexec::resize_scratch( ValueTraits::value_size( m_functor ) , team_reduce_size + m_shmem_size );

#pragma omp parallel
      {
        OpenMPexec & exec = * OpenMPexec::get_thread_omp();

        ParallelReduce::template exec_team< WorkTag >
          ( m_functor
          , Member( exec , m_policy , m_shmem_size )
          , ValueInit::init( m_functor , exec.scratch_reduce() ) );
      }
/* END #pragma omp parallel */

      {
        const pointer_type ptr = pointer_type( OpenMPexec::pool_rev(0)->scratch_reduce() );

        int max_active_threads = OpenMPexec::pool_size();
        if( max_active_threads > m_policy.league_size()* m_policy.team_size() )
          max_active_threads = m_policy.league_size()* m_policy.team_size();

        for ( int i = 1 ; i < max_active_threads ; ++i ) {
          ValueJoin::join( m_functor , ptr , OpenMPexec::pool_rev(i)->scratch_reduce() );
        }

        Kokkos::Impl::FunctorFinal< FunctorType , WorkTag >::final( m_functor , ptr );

        if ( m_result_ptr ) {
          const int n = ValueTraits::value_count( m_functor );

          for ( int j = 0 ; j < n ; ++j ) { m_result_ptr[j] = ptr[j] ; }
        }
      }
    }

  template< class ViewType >
  inline
  ParallelReduce( const FunctorType  & arg_functor ,
                  const Policy       & arg_policy ,
                  const ViewType     & arg_result )
    : m_functor( arg_functor )
    , m_policy(  arg_policy )
    , m_result_ptr( arg_result.ptr_on_device() )
    , m_shmem_size( arg_policy.scratch_size() + FunctorTeamShmemSize< FunctorType >::value( arg_functor , arg_policy.team_size() ) )
    {}
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* KOKKOS_OPENMP_PARALLEL_HPP */

