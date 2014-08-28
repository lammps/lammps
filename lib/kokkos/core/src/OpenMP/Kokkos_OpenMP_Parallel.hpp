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

#ifndef KOKKOS_OPENMP_PARALLEL_HPP
#define KOKKOS_OPENMP_PARALLEL_HPP

#include <omp.h>

#include <Kokkos_Parallel.hpp>
#include <OpenMP/Kokkos_OpenMPexec.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , typename IntType , unsigned P >
class ParallelFor< FunctorType
                 , Kokkos::RangePolicy< Kokkos::OpenMP , void , IntType , P >
                 , Kokkos::OpenMP
                 >
{
public:
  typedef Kokkos::RangePolicy< Kokkos::OpenMP , void , IntType , P > Policy ;

  inline
  ParallelFor( const FunctorType & functor
             , const Policy      & policy )
  {
    OpenMPexec::verify_is_process("Kokkos::OpenMP parallel_for");
    OpenMPexec::verify_initialized("Kokkos::OpenMP parallel_for");

#pragma omp parallel
    {
      OpenMPexec & exec = * OpenMPexec::get_thread_omp();

      const Policy range( policy , exec.pool_rank() , exec.pool_size() );

      const typename Policy::member_type work_end = range.end();
      for ( typename Policy::member_type iwork = range.begin() ; iwork < work_end ; ++iwork ) {
        functor( iwork );
      }
    }
/* END #pragma omp parallel */
  }
};
} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , typename IntType , unsigned P >
class ParallelReduce< FunctorType
                    , Kokkos::RangePolicy< Kokkos::OpenMP , void , IntType , P >
                    , Kokkos::OpenMP
                    >
{
public:
  typedef ReduceAdapter< FunctorType >   Reduce ;
  typedef typename Reduce::pointer_type  pointer_type ;
  typedef Kokkos::RangePolicy< Kokkos::OpenMP , void , IntType , P > Policy ;

  template< class HostViewType >
  inline
  ParallelReduce( const FunctorType  & functor
                , const Policy       & policy
                , const HostViewType & result_view )
  {
    OpenMPexec::verify_is_process("Kokkos::OpenMP parallel_reduce");
    OpenMPexec::verify_initialized("Kokkos::OpenMP parallel_reduce");

    OpenMPexec::resize_scratch( Reduce::value_size( functor ) , 0 );

#pragma omp parallel
    {
      OpenMPexec & exec = * OpenMPexec::get_thread_omp();

      typename Reduce::reference_type update = Reduce::init( functor , exec.scratch_reduce() );

      const Policy range( policy , exec.pool_rank() , exec.pool_size() );

      const typename Policy::member_type work_end = range.end();
      for ( typename Policy::member_type iwork = range.begin() ; iwork < work_end ; ++iwork ) {
        functor( iwork , update );
      }
    }
/* END #pragma omp parallel */

    {
      const pointer_type ptr = pointer_type( OpenMPexec::pool_rev(0)->scratch_reduce() );

      for ( int i = 1 ; i < OpenMPexec::pool_size() ; ++i ) {
        Reduce::join( functor , ptr , OpenMPexec::pool_rev(i)->scratch_reduce() );
      }

      Reduce::final( functor , ptr );

      if ( result_view.ptr_on_device() ) {
        const int n = Reduce::value_count( functor );

        for ( int j = 0 ; j < n ; ++j ) { result_view.ptr_on_device()[j] = ptr[j] ; }
      }
    }
  }

  void wait() {}
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , typename IntType , unsigned P >
class ParallelScan< FunctorType
                  , Kokkos::RangePolicy< Kokkos::OpenMP , void , IntType , P >
                  , Kokkos::OpenMP
                  >
{
public:
  typedef ReduceAdapter< FunctorType >   Reduce ;
  typedef typename Reduce::pointer_type  pointer_type ;
  typedef Kokkos::RangePolicy< Kokkos::OpenMP , void , IntType , P > Policy ;

  inline
  ParallelScan( const FunctorType & functor , const Policy & policy )
  {
    OpenMPexec::verify_is_process("Kokkos::OpenMP parallel_scan");
    OpenMPexec::verify_initialized("Kokkos::OpenMP parallel_scan");

    OpenMPexec::resize_scratch( 2 * Reduce::value_size( functor ) , 0 );

#pragma omp parallel
    {
      OpenMPexec & exec = * OpenMPexec::get_thread_omp();

      const Policy range( policy , exec.pool_rank() , exec.pool_size() );

      typename Reduce::reference_type update =
        Reduce::init( functor ,
                      pointer_type( exec.scratch_reduce() ) + Reduce::value_count( functor ) );

      const typename Policy::member_type work_end = range.end();
      for ( typename Policy::member_type iwork = range.begin() ; iwork < work_end ; ++iwork ) {
        functor( iwork , update , false );
      }
    }
/* END #pragma omp parallel */

    {
      const unsigned thread_count = OpenMPexec::pool_size();
      const unsigned value_count  = Reduce::value_count( functor );

      pointer_type ptr_prev = 0 ;

      for ( unsigned rank_rev = thread_count ; rank_rev-- ; ) {

        pointer_type ptr = pointer_type( OpenMPexec::pool_rev(rank_rev)->scratch_reduce() );

        if ( ptr_prev ) {
          for ( unsigned i = 0 ; i < value_count ; ++i ) { ptr[i] = ptr_prev[ i + value_count ] ; }
          Reduce::join( functor , ptr + value_count , ptr );
        }
        else {
          Reduce::init( functor , ptr );
        }

        ptr_prev = ptr ;
      }
    }

#pragma omp parallel
    {
      OpenMPexec & exec = * OpenMPexec::get_thread_omp();

      const Policy range( policy , exec.pool_rank() , exec.pool_size() );

      typename Reduce::reference_type update =
        Reduce::reference( pointer_type( exec.scratch_reduce() ) );

      const typename Policy::member_type work_end = range.end();
      for ( typename Policy::member_type iwork = range.begin() ; iwork < work_end ; ++iwork ) {
        functor( iwork , update , true );
      }
    }
/* END #pragma omp parallel */

  }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType >
class ParallelFor< FunctorType , Kokkos::TeamPolicy< Kokkos::OpenMP , void > , Kokkos::OpenMP >
{
public:
  typedef Kokkos::TeamPolicy< Kokkos::OpenMP , void > Policy ;

  inline
  ParallelFor( const FunctorType & functor ,
               const Policy      & policy )
  {
    OpenMPexec::verify_is_process("Kokkos::OpenMP parallel_for");
    OpenMPexec::verify_initialized("Kokkos::OpenMP parallel_for");

    const size_t team_reduce_size = Policy::member_type::team_reduce_size();
    const size_t team_shmem_size  = FunctorTeamShmemSize< FunctorType >::value( functor , policy.team_size() );

    OpenMPexec::resize_scratch( 0 , team_reduce_size + team_shmem_size );

#pragma omp parallel
    {
      typename Policy::member_type member( * OpenMPexec::get_thread_omp() , policy , team_shmem_size );;

      for ( ; member.valid() ; member.next() ) {
        functor( member );
      }
    }
/* END #pragma omp parallel */
  }

  void wait() {}
};

template< class FunctorType >
class ParallelReduce< FunctorType , Kokkos::TeamPolicy< Kokkos::OpenMP , void > , Kokkos::OpenMP >
{
public:
  typedef Kokkos::TeamPolicy< Kokkos::OpenMP , void > Policy ;
  typedef ReduceAdapter< FunctorType >   Reduce ;
  typedef typename Reduce::pointer_type  pointer_type ;

  inline
  ParallelReduce( const FunctorType  & functor ,
                  const Policy       & policy )
  {
    OpenMPexec::verify_is_process("Kokkos::OpenMP parallel_reduce");

    const size_t team_reduce_size = Policy::member_type::team_reduce_size();
    const size_t team_shmem_size  = FunctorTeamShmemSize< FunctorType >::value( functor , policy.team_size() );

    OpenMPexec::resize_scratch( Reduce::value_size( functor ) , team_reduce_size + team_shmem_size );

#pragma omp parallel
    {
      OpenMPexec & exec = * OpenMPexec::get_thread_omp();

      typename Reduce::reference_type update = Reduce::init( functor , exec.scratch_reduce() );

      for ( typename Policy::member_type member( exec , policy , team_shmem_size ); member.valid() ; member.next() ) {
        functor( member , update );
      }
    }
/* END #pragma omp parallel */

    {
      const pointer_type ptr = pointer_type( OpenMPexec::pool_rev(0)->scratch_reduce() );

      for ( int i = 1 ; i < OpenMPexec::pool_size() ; ++i ) {
        Reduce::join( functor , ptr , OpenMPexec::pool_rev(i)->scratch_reduce() );
      }

      Reduce::final( functor , ptr );
    }
  }

  template< class ViewType >
  inline
  ParallelReduce( const FunctorType  & functor ,
                  const Policy       & policy ,
                  const ViewType     & result )
  {
    OpenMPexec::verify_is_process("Kokkos::OpenMP parallel_reduce");

    const size_t team_reduce_size = Policy::member_type::team_reduce_size();
    const size_t team_shmem_size  = FunctorTeamShmemSize< FunctorType >::value( functor , policy.team_size() );

    OpenMPexec::resize_scratch( Reduce::value_size( functor ) , team_reduce_size + team_shmem_size );

#pragma omp parallel
    {
      OpenMPexec & exec = * OpenMPexec::get_thread_omp();

      typename Reduce::reference_type update = Reduce::init( functor , exec.scratch_reduce() );

      for ( typename Policy::member_type member( exec , policy , team_shmem_size ); member.valid() ; member.next() ) {
        functor( member , update );
      }
    }
/* END #pragma omp parallel */

    {
      const pointer_type ptr = pointer_type( OpenMPexec::pool_rev(0)->scratch_reduce() );

      for ( int i = 1 ; i < OpenMPexec::pool_size() ; ++i ) {
        Reduce::join( functor , ptr , OpenMPexec::pool_rev(i)->scratch_reduce() );
      }

      Reduce::final( functor , ptr );

      const int n = Reduce::value_count( functor );

      for ( int j = 0 ; j < n ; ++j ) { result.ptr_on_device()[j] = ptr[j] ; }
    }
  }

  void wait() {}
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* KOKKOS_OPENMP_PARALLEL_HPP */

