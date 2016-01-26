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

#ifndef KOKKOS_CUDA_PARALLEL_HPP
#define KOKKOS_CUDA_PARALLEL_HPP

#include <iostream>
#include <algorithm>
#include <stdio.h>

#include <Kokkos_Macros.hpp>

/* only compile this file if CUDA is enabled for Kokkos */
#if defined( __CUDACC__ ) && defined( KOKKOS_HAVE_CUDA )

#include <utility>
#include <Kokkos_Parallel.hpp>

#include <Cuda/Kokkos_CudaExec.hpp>
#include <Cuda/Kokkos_Cuda_ReduceScan.hpp>
#include <Cuda/Kokkos_Cuda_Internal.hpp>
#include <Kokkos_Vectorization.hpp>

#ifdef KOKKOSP_ENABLE_PROFILING
#include <impl/Kokkos_Profiling_Interface.hpp>
#include <typeinfo>
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< typename Type >
struct CudaJoinFunctor {
  typedef Type value_type ;

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    volatile const value_type & input )
    { update += input ; }
};

class CudaTeamMember {
private:

  typedef Kokkos::Cuda                           execution_space ;
  typedef execution_space::scratch_memory_space  scratch_memory_space ;

  void                * m_team_reduce ;
  scratch_memory_space  m_team_shared ;
  int                   m_league_rank ;
  int                   m_league_size ;

public:

#if defined( __CUDA_ARCH__ )

  __device__ inline
  const execution_space::scratch_memory_space & team_shmem() const
    { return m_team_shared ; }

  __device__ inline int league_rank() const { return m_league_rank ; }
  __device__ inline int league_size() const { return m_league_size ; }
  __device__ inline int team_rank() const { return threadIdx.y ; }
  __device__ inline int team_size() const { return blockDim.y ; }

  __device__ inline void team_barrier() const { __syncthreads(); }

  template<class ValueType>
  __device__ inline void team_broadcast(ValueType& value, const int& thread_id) const {
    __shared__ ValueType sh_val;
    if(threadIdx.x == 0 && threadIdx.y == thread_id) {
      sh_val = val;
    }
    team_barrier();
    val = sh_val;
  }

#ifdef KOKKOS_HAVE_CXX11
  template< class ValueType, class JoinOp >
  __device__ inline
  typename JoinOp::value_type team_reduce( const ValueType & value
                                         , const JoinOp & op_in ) const
    {
      typedef JoinLambdaAdapter<ValueType,JoinOp> JoinOpFunctor ;
      const JoinOpFunctor op(op_in);
      ValueType * const base_data = (ValueType *) m_team_reduce ;
#else
  template< class JoinOp >
  __device__ inline
  typename JoinOp::value_type team_reduce( const typename JoinOp::value_type & value
                                         , const JoinOp & op ) const
    {
      typedef JoinOp JoinOpFunctor ;
      typename JoinOp::value_type * const base_data = (typename JoinOp::value_type *) m_team_reduce ;
#endif

      __syncthreads(); // Don't write in to shared data until all threads have entered this function

      if ( 0 == threadIdx.y ) { base_data[0] = 0 ; }

      base_data[ threadIdx.y ] = value ;

      Impl::cuda_intra_block_reduce_scan<false,JoinOpFunctor,void>( op , base_data );

      return base_data[ blockDim.y - 1 ];
    }

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
  __device__ inline Type team_scan( const Type & value , Type * const global_accum ) const
    {
      Type * const base_data = (Type *) m_team_reduce ;

      __syncthreads(); // Don't write in to shared data until all threads have entered this function

      if ( 0 == threadIdx.y ) { base_data[0] = 0 ; }

      base_data[ threadIdx.y + 1 ] = value ;

      Impl::cuda_intra_block_reduce_scan<true,Impl::CudaJoinFunctor<Type>,void>( Impl::CudaJoinFunctor<Type>() , base_data + 1 );

      if ( global_accum ) {
        if ( blockDim.y == threadIdx.y + 1 ) {
          base_data[ blockDim.y ] = atomic_fetch_add( global_accum , base_data[ blockDim.y ] );
        }
        __syncthreads(); // Wait for atomic
        base_data[ threadIdx.y ] += base_data[ blockDim.y ] ;
      }

      return base_data[ threadIdx.y ];
    }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
   *
   *  The highest rank thread can compute the reduction total as
   *    reduction_total = dev.team_scan( value ) + value ;
   */
  template< typename Type >
  __device__ inline Type team_scan( const Type & value ) const
    { return this->template team_scan<Type>( value , 0 ); }

  //----------------------------------------
  // Private for the driver

  __device__ inline
  CudaTeamMember( void * shared
                , const int shared_begin
                , const int shared_size
                , const int arg_league_rank
                , const int arg_league_size )
    : m_team_reduce( shared )
    , m_team_shared( ((char *)shared) + shared_begin , shared_size )
    , m_league_rank( arg_league_rank ) 
    , m_league_size( arg_league_size ) 
    {}

#else

  const execution_space::scratch_memory_space & team_shmem() const {return m_team_shared;}

  int league_rank() const {return 0;}
  int league_size() const {return 1;}
  int team_rank() const {return 0;}
  int team_size() const {return 1;}

  void team_barrier() const {}
  template<class ValueType>
  void team_broadcast(ValueType& value, const int& thread_id) const {}

  template< class JoinOp >
  typename JoinOp::value_type team_reduce( const typename JoinOp::value_type & value
                                         , const JoinOp & op ) const {return typename JoinOp::value_type();}

  template< typename Type >
  Type team_scan( const Type & value , Type * const global_accum ) const {return Type();}

  template< typename Type >
  Type team_scan( const Type & value ) const {return Type();}

  //----------------------------------------
  // Private for the driver

  CudaTeamMember( void * shared
                , const int shared_begin
                , const int shared_end
                , const int arg_league_rank
                , const int arg_league_size );

#endif /* #if ! defined( __CUDA_ARCH__ ) */

};

} // namespace Impl

template< class Arg0 , class Arg1 >
class TeamPolicy< Arg0 , Arg1 , Kokkos::Cuda >
{
private:

  enum { MAX_WARP = 8 };

  const int m_league_size ;
  const int m_team_size ;
  const int m_vector_length ;
  const size_t m_scratch_size ;

public:

  //! Tag this class as a kokkos execution policy
  typedef TeamPolicy     execution_policy ;

  //! Execution space of this execution policy
  typedef Kokkos::Cuda  execution_space ;

  typedef typename
    Impl::if_c< ! Impl::is_same< Kokkos::Cuda , Arg0 >::value , Arg0 , Arg1 >::type
      work_tag ;

  //----------------------------------------

  template< class FunctorType >
  inline static
  int team_size_max( const FunctorType & functor )
    {
      int n = MAX_WARP * Impl::CudaTraits::WarpSize ;

      for ( ; n ; n >>= 1 ) {
        const int shmem_size =
          /* for global reduce */ Impl::cuda_single_inter_block_reduce_scan_shmem<false,FunctorType,work_tag>( functor , n )
          /* for team   reduce */ + ( n + 2 ) * sizeof(double)
          /* for team   shared */ + Impl::FunctorTeamShmemSize< FunctorType >::value( functor , n );

        if ( shmem_size < Impl::CudaTraits::SharedMemoryCapacity ) break ;
      }

      return n ;
    }

  template< class FunctorType >
  static int team_size_recommended( const FunctorType & functor )
    { return team_size_max( functor ); }

  template< class FunctorType >
  static int team_size_recommended( const FunctorType & functor , const int vector_length)
    {
      int max = team_size_max( functor )/vector_length;
      if(max<1) max = 1;
      return max;
    }

  inline static
  int vector_length_max()
    { return Impl::CudaTraits::WarpSize; }

  //----------------------------------------

  inline int vector_length()   const { return m_vector_length ; }
  inline int team_size()   const { return m_team_size ; }
  inline int league_size() const { return m_league_size ; }
  inline size_t scratch_size() const { return m_scratch_size ; }

  /** \brief  Specify league size, request team size */
  TeamPolicy( execution_space &
            , int league_size_
            , int team_size_request
            , int vector_length_request = 1 )
    : m_league_size( league_size_ )
    , m_team_size( team_size_request )
    , m_vector_length( vector_length_request )
    , m_scratch_size ( 0 )
    {
      // Allow only power-of-two vector_length
      if ( ! Kokkos::Impl::is_integral_power_of_two( vector_length_request ) ) {
        Impl::throw_runtime_exception( "Requested non-power-of-two vector length for TeamPolicy.");
      }

      // Make sure league size is permissable
      if(league_size_ >= int(Impl::cuda_internal_maximum_grid_count()))
        Impl::throw_runtime_exception( "Requested too large league_size for TeamPolicy on Cuda execution space.");

      // Make sure total block size is permissable
      if ( m_team_size * m_vector_length > 1024 ) {
        Impl::throw_runtime_exception(std::string("Kokkos::TeamPolicy< Cuda > the team size is too large. Team size x vector length must be smaller than 1024."));
      }
    }

  /** \brief  Specify league size, request team size */
  TeamPolicy( execution_space &
            , int league_size_
            , const Kokkos::AUTO_t & /* team_size_request */
            , int vector_length_request = 1 )
    : m_league_size( league_size_ )
    , m_team_size( -1 )
    , m_vector_length( vector_length_request )
    , m_scratch_size ( 0 )
    {
      // Allow only power-of-two vector_length
      if ( ! Kokkos::Impl::is_integral_power_of_two( vector_length_request ) ) {
        Impl::throw_runtime_exception( "Requested non-power-of-two vector length for TeamPolicy.");
      }

      // Make sure league size is permissable
      if(league_size_ >= int(Impl::cuda_internal_maximum_grid_count()))
        Impl::throw_runtime_exception( "Requested too large league_size for TeamPolicy on Cuda execution space.");
    }

  TeamPolicy( int league_size_
            , int team_size_request
            , int vector_length_request = 1 )
    : m_league_size( league_size_ )
    , m_team_size( team_size_request )
    , m_vector_length ( vector_length_request )
    , m_scratch_size ( 0 )
    {
      // Allow only power-of-two vector_length
      if ( ! Kokkos::Impl::is_integral_power_of_two( vector_length_request ) ) {
        Impl::throw_runtime_exception( "Requested non-power-of-two vector length for TeamPolicy.");
      }

      // Make sure league size is permissable
      if(league_size_ >= int(Impl::cuda_internal_maximum_grid_count()))
        Impl::throw_runtime_exception( "Requested too large league_size for TeamPolicy on Cuda execution space.");

      // Make sure total block size is permissable
      if ( m_team_size * m_vector_length > 1024 ) {
        Impl::throw_runtime_exception(std::string("Kokkos::TeamPolicy< Cuda > the team size is too large. Team size x vector length must be smaller than 1024."));
      }
    }

  TeamPolicy( int league_size_
            , const Kokkos::AUTO_t & /* team_size_request */
            , int vector_length_request = 1 )
    : m_league_size( league_size_ )
    , m_team_size( -1 )
    , m_vector_length ( vector_length_request )
    , m_scratch_size ( 0 )
    {
      // Allow only power-of-two vector_length
      if ( ! Kokkos::Impl::is_integral_power_of_two( vector_length_request ) ) {
        Impl::throw_runtime_exception( "Requested non-power-of-two vector length for TeamPolicy.");
      }

      // Make sure league size is permissable
      if(league_size_ >= int(Impl::cuda_internal_maximum_grid_count()))
        Impl::throw_runtime_exception( "Requested too large league_size for TeamPolicy on Cuda execution space.");
    }

  template<class MemorySpace>
  TeamPolicy( int league_size_
            , int team_size_request
            , const Experimental::TeamScratchRequest<MemorySpace> & scratch_request )
    : m_league_size( league_size_ )
    , m_team_size( team_size_request )
    , m_vector_length( 1 )
    , m_scratch_size(scratch_request.total(team_size_request))
    {
      // Allow only power-of-two vector_length
      if ( ! Kokkos::Impl::is_integral_power_of_two( m_vector_length ) ) {
        Impl::throw_runtime_exception( "Requested non-power-of-two vector length for TeamPolicy.");
      }

      // Make sure league size is permissable
      if(league_size_ >= int(Impl::cuda_internal_maximum_grid_count()))
        Impl::throw_runtime_exception( "Requested too large league_size for TeamPolicy on Cuda execution space.");

      // Make sure total block size is permissable
      if ( m_team_size * m_vector_length > 1024 ) {
        Impl::throw_runtime_exception(std::string("Kokkos::TeamPolicy< Cuda > the team size is too large. Team size x vector length must be smaller than 1024."));
      }
    }

  template<class MemorySpace>
  TeamPolicy( int league_size_
            , const Kokkos::AUTO_t & /* team_size_request */
            , const Experimental::TeamScratchRequest<MemorySpace> & scratch_request )
    : m_league_size( league_size_ )
    , m_team_size( 256 )
    , m_vector_length ( 1 )
    , m_scratch_size(scratch_request.total(2356))
    {
      // Allow only power-of-two vector_length
      if ( ! Kokkos::Impl::is_integral_power_of_two( m_vector_length ) ) {
        Impl::throw_runtime_exception( "Requested non-power-of-two vector length for TeamPolicy.");
      }

      // Make sure league size is permissable
      if(league_size_ >= int(Impl::cuda_internal_maximum_grid_count()))
        Impl::throw_runtime_exception( "Requested too large league_size for TeamPolicy on Cuda execution space.");
    }

  typedef Kokkos::Impl::CudaTeamMember member_type ;
};

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelFor< FunctorType
                 , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Cuda >
                 >
{
private:

  typedef Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Cuda > Policy ;
  typedef typename Policy::member_type  Member ;
  typedef typename Policy::work_tag     WorkTag ;

  const FunctorType  m_functor ;
  const Policy       m_policy ;  

  ParallelFor() = delete ;
  ParallelFor & operator = ( const ParallelFor & ) = delete ;

  template< class TagType >
  inline __device__
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_range( const Member i ) const
    { m_functor( i ); }

  template< class TagType >
  inline __device__
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec_range( const Member i ) const
    { m_functor( TagType() , i ); }

public:

  typedef FunctorType functor_type ;

  inline
  __device__
  void operator()(void) const
    {
      const Member work_stride = blockDim.y * gridDim.x ;
      const Member work_end    = m_policy.end();

      for ( Member
              iwork =  m_policy.begin() + threadIdx.y + blockDim.y * blockIdx.x ;
              iwork <  work_end ;
              iwork += work_stride ) {
        this-> template exec_range< WorkTag >( iwork );
      }
    }

  inline
  void execute() const
    {
      const int nwork = m_policy.end() - m_policy.begin();
      const dim3 block(  1 , CudaTraits::WarpSize * cuda_internal_maximum_warp_count(), 1);
      const dim3 grid( std::min( ( nwork + block.y - 1 ) / block.y , cuda_internal_maximum_grid_count() ) , 1 , 1);

      CudaParallelLaunch< ParallelFor >( *this , grid , block , 0 );
    }

  ParallelFor( const FunctorType  & arg_functor ,
               const Policy       & arg_policy )
    : m_functor( arg_functor )
    , m_policy(  arg_policy )
    { }
};

template< class FunctorType , class Arg0 , class Arg1 >
class ParallelFor< FunctorType
                 , Kokkos::TeamPolicy< Arg0 , Arg1 , Kokkos::Cuda >
                 >
{
private:

  typedef Kokkos::TeamPolicy< Arg0 , Arg1 , Kokkos::Cuda >   Policy ;
  typedef typename Policy::member_type  Member ;
  typedef typename Policy::work_tag     WorkTag ;

public:

  typedef FunctorType      functor_type ;
  typedef Cuda::size_type  size_type ;

private:

  // Algorithmic constraints: blockDim.y is a power of two AND blockDim.y == blockDim.z == 1
  // shared memory utilization:
  //
  //  [ team   reduce space ]
  //  [ team   shared space ]
  //

  const FunctorType m_functor ;
  const size_type   m_league_size ;
  const size_type   m_team_size ;
  const size_type   m_vector_size ;
  const size_type   m_shmem_begin ;
  const size_type   m_shmem_size ;

  template< class TagType >
  __device__ inline
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_team( const Member & member ) const
    { m_functor( member ); }

  template< class TagType >
  __device__ inline
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec_team( const Member & member ) const
    { m_functor( TagType() , member ); }

public:

  __device__ inline
  void operator()(void) const
  {
    // Iterate this block through the league
    for ( int league_rank = blockIdx.x ; league_rank < m_league_size ; league_rank += gridDim.x ) {

      this-> template exec_team< WorkTag >(
        typename Policy::member_type( kokkos_impl_cuda_shared_memory<void>()
                                    , m_shmem_begin
                                    , m_shmem_size
                                    , league_rank
                                    , m_league_size ) );
    }
  }

  inline
  void execute() const
    {
      const int shmem_size_total = m_shmem_begin + m_shmem_size ;
      const dim3 grid( int(m_league_size) , 1 , 1 );
      const dim3 block( int(m_vector_size) , int(m_team_size) , 1 );

      CudaParallelLaunch< ParallelFor >( *this, grid, block, shmem_size_total ); // copy to device and execute

    }

  ParallelFor( const FunctorType  & arg_functor 
             , const Policy       & arg_policy 
             )
    : m_functor( arg_functor )
    , m_league_size( arg_policy.league_size() )
    , m_team_size( 0 <= arg_policy.team_size() ? arg_policy.team_size() :
        Kokkos::Impl::cuda_get_opt_block_size< ParallelFor >( arg_functor , arg_policy.vector_length(), arg_policy.scratch_size() ) / arg_policy.vector_length() )
    , m_vector_size( arg_policy.vector_length() )
    , m_shmem_begin( sizeof(double) * ( m_team_size + 2 ) )
    , m_shmem_size( arg_policy.scratch_size() + FunctorTeamShmemSize< FunctorType >::value( m_functor , m_team_size ) )
    {
      // Functor's reduce memory, team scan memory, and team shared memory depend upon team size.

      const int shmem_size_total = m_shmem_begin + m_shmem_size ;

      if ( CudaTraits::SharedMemoryCapacity < shmem_size_total ) {
        Kokkos::Impl::throw_runtime_exception(std::string("Kokkos::Impl::ParallelFor< Cuda > insufficient shared memory"));
      }

      if ( m_team_size >
           Kokkos::Impl::cuda_get_max_block_size< ParallelFor >
                 ( arg_functor , arg_policy.vector_length(), arg_policy.scratch_size() ) / arg_policy.vector_length()) {
        Kokkos::Impl::throw_runtime_exception(std::string("Kokkos::Impl::ParallelFor< Cuda > requested too large team size."));
      }
    }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelReduce< FunctorType
                    , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Cuda >
                    >
{
private:

  typedef Kokkos::RangePolicy<Arg0,Arg1,Arg2, Kokkos::Cuda >         Policy ;

  typedef typename Policy::WorkRange    WorkRange ;
  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::member_type  Member ;

  typedef Kokkos::Impl::FunctorValueTraits< FunctorType, WorkTag > ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   FunctorType, WorkTag > ValueInit ;

public:

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::value_type      value_type ;
  typedef typename ValueTraits::reference_type  reference_type ;
  typedef FunctorType                           functor_type ;
  typedef Cuda::size_type                       size_type ;

  // Algorithmic constraints: blockSize is a power of two AND blockDim.y == blockDim.z == 1

  const FunctorType   m_functor ;
  const Policy        m_policy ;
  const pointer_type  m_result_ptr ;
  size_type *         m_scratch_space ;
  size_type *         m_scratch_flags ;
  size_type *         m_unified_space ;

  template< class TagType >
  __device__ inline
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_range( const Member & i , reference_type update ) const
    { m_functor( i , update ); }

  template< class TagType >
  __device__ inline
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec_range( const Member & i , reference_type update ) const
    { m_functor( TagType() , i , update ); }

#if ! defined( KOKKOS_EXPERIMENTAL_CUDA_SHFL_REDUCTION )

  __device__ inline
  void operator()(void) const
  {
    const integral_nonzero_constant< size_type , ValueTraits::StaticValueSize / sizeof(size_type) >
      word_count( ValueTraits::value_size( m_functor ) / sizeof(size_type) );

    {
      reference_type value =
        ValueInit::init( m_functor , kokkos_impl_cuda_shared_memory<size_type>() + threadIdx.y * word_count.value );

      // Number of blocks is bounded so that the reduction can be limited to two passes.
      // Each thread block is given an approximately equal amount of work to perform.
      // Accumulate the values for this block.
      // The accumulation ordering does not match the final pass, but is arithmatically equivalent.

      const WorkRange range( m_policy , blockIdx.x , gridDim.x );

      for ( Member iwork = range.begin() + threadIdx.y , iwork_end = range.end() ;
            iwork < iwork_end ; iwork += blockDim.y ) {
        this-> template exec_range< WorkTag >( iwork , value );
      }
    }

    // Reduce with final value at blockDim.y - 1 location.
    if ( cuda_single_inter_block_reduce_scan<false,FunctorType,WorkTag>(
           m_functor , blockIdx.x , gridDim.x ,
           kokkos_impl_cuda_shared_memory<size_type>() , m_scratch_space , m_scratch_flags ) ) {

      // This is the final block with the final result at the final threads' location

      size_type * const shared = kokkos_impl_cuda_shared_memory<size_type>() + ( blockDim.y - 1 ) * word_count.value ;
      size_type * const global = m_unified_space ? m_unified_space : m_scratch_space ;

      if ( threadIdx.y == 0 ) {
        Kokkos::Impl::FunctorFinal< FunctorType , WorkTag >::final( m_functor , shared );
      }

      if ( CudaTraits::WarpSize < word_count.value ) { __syncthreads(); }

      for ( unsigned i = threadIdx.y ; i < word_count.value ; i += blockDim.y ) { global[i] = shared[i]; }
    }
  }

#else /* defined( KOKKOS_EXPERIMENTAL_CUDA_SHFL_REDUCTION ) */

  __device__ inline
   void operator()(void) const
   {

     value_type value = 0;

     // Number of blocks is bounded so that the reduction can be limited to two passes.
     // Each thread block is given an approximately equal amount of work to perform.
     // Accumulate the values for this block.
     // The accumulation ordering does not match the final pass, but is arithmatically equivalent.

     const Policy range( m_policy , blockIdx.x , gridDim.x );

     for ( Member iwork = range.begin() + threadIdx.y , iwork_end = range.end() ;
           iwork < iwork_end ; iwork += blockDim.y ) {
       this-> template exec_range< WorkTag >( iwork , value );
     }

     pointer_type const result = (pointer_type) (m_unified_space ? m_unified_space : m_scratch_space) ;
     int max_active_thread = range.end()-range.begin() < blockDim.y ? range.end() - range.begin():blockDim.y;
     max_active_thread = max_active_thread == 0?blockDim.y:max_active_thread;
     if(Impl::cuda_inter_block_reduction<FunctorType,Impl::JoinAdd<value_type> >
            (value,Impl::JoinAdd<value_type>(),m_scratch_space,result,m_scratch_flags,max_active_thread)) {
       const unsigned id = threadIdx.y*blockDim.x + threadIdx.x;
       if(id==0) {
         Kokkos::Impl::FunctorFinal< FunctorType , WorkTag >::final( m_functor , (void*) &value );
         *result = value;
       }
     }
   }

#endif

  // Determine block size constrained by shared memory:
  static inline
  unsigned local_block_size( const FunctorType & f )
    {
      unsigned n = CudaTraits::WarpSize * 8 ;
      while ( n && CudaTraits::SharedMemoryCapacity < cuda_single_inter_block_reduce_scan_shmem<false,FunctorType,WorkTag>( f , n ) ) { n >>= 1 ; }
      return n ;
    }

  inline
  void execute()
    {
      const int block_size = local_block_size( m_functor );

      m_scratch_space = cuda_internal_scratch_space( ValueTraits::value_size( m_functor ) * block_size /* block_size == max block_count */ );
      m_scratch_flags = cuda_internal_scratch_flags( sizeof(size_type) );
      m_unified_space = cuda_internal_scratch_unified( ValueTraits::value_size( m_functor ) );

      const int nwork = m_policy.end() - m_policy.begin();
      // REQUIRED ( 1 , N , 1 )
      const dim3 block( 1 , block_size , 1 );
      // Required grid.x <= block.y
      const dim3 grid( std::min( int(block.y) , int( ( nwork + block.y - 1 ) / block.y ) ) , 1 , 1 );

#ifdef KOKKOS_EXPERIMENTAL_CUDA_SHFL_REDUCTION
    const int shmem = 0;
#else
    const int shmem = cuda_single_inter_block_reduce_scan_shmem<false,FunctorType,WorkTag>( m_functor , block.y );
#endif

    CudaParallelLaunch< ParallelReduce >( *this, grid, block, shmem ); // copy to device and execute

    Cuda::fence();

    if ( m_result_ptr ) {
      if ( m_unified_space ) {
        const int count = ValueTraits::value_count( m_functor );
        for ( int i = 0 ; i < count ; ++i ) { m_result_ptr[i] = pointer_type(m_unified_space)[i] ; }
      }
      else {
        const int size = ValueTraits::value_size( m_functor );
        DeepCopy<HostSpace,CudaSpace>( m_result_ptr , m_scratch_space , size );
      }
    }
  }

  template< class HostViewType >
  ParallelReduce( const FunctorType  & arg_functor 
                , const Policy       & arg_policy 
                , const HostViewType & arg_result
                )
  : m_functor( arg_functor )
  , m_policy(  arg_policy )
  , m_result_ptr( arg_result.ptr_on_device() )
  , m_scratch_space( 0 )
  , m_scratch_flags( 0 )
  , m_unified_space( 0 )
  { }
};

//----------------------------------------------------------------------------

template< class FunctorType , class Arg0 , class Arg1 >
class ParallelReduce< FunctorType
                    , Kokkos::TeamPolicy< Arg0 , Arg1 , Kokkos::Cuda >
                    >
{
private:

  typedef Kokkos::TeamPolicy<Arg0,Arg1,Kokkos::Cuda>  Policy ;
  typedef typename Policy::member_type  Member ;
  typedef typename Policy::work_tag     WorkTag ;

  typedef Kokkos::Impl::FunctorValueTraits< FunctorType, WorkTag > ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   FunctorType, WorkTag > ValueInit ;

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;

public:

  typedef FunctorType      functor_type ;
  typedef Cuda::size_type  size_type ;

private:

  // Algorithmic constraints: blockDim.y is a power of two AND blockDim.y == blockDim.z == 1
  // shared memory utilization:
  //
  //  [ global reduce space ]
  //  [ team   reduce space ]
  //  [ team   shared space ]
  //

  const FunctorType   m_functor ;
  const pointer_type  m_result_ptr ;
  size_type *         m_scratch_space ;
  size_type *         m_scratch_flags ;
  size_type *         m_unified_space ;
  size_type           m_team_begin ;
  size_type           m_shmem_begin ;
  size_type           m_shmem_size ;
  const size_type     m_league_size ;
  const size_type     m_team_size ;

  template< class TagType >
  __device__ inline
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_team( const Member & member , reference_type update ) const
    { m_functor( member , update ); }

  template< class TagType >
  __device__ inline
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec_team( const Member & member , reference_type update ) const
    { m_functor( TagType() , member , update ); }

public:

  __device__ inline
  void operator()(void) const
  {
    const integral_nonzero_constant< size_type , ValueTraits::StaticValueSize / sizeof(size_type) >
      word_count( ValueTraits::value_size( m_functor ) / sizeof(size_type) );

    reference_type value =
      ValueInit::init( m_functor , kokkos_impl_cuda_shared_memory<size_type>() + threadIdx.y * word_count.value );

    // Iterate this block through the league
    for ( int league_rank = blockIdx.x ; league_rank < m_league_size ; league_rank += gridDim.x ) {
      this-> template exec_team< WorkTag >
        ( Member( kokkos_impl_cuda_shared_memory<char>() + m_team_begin
                                        , m_shmem_begin
                                        , m_shmem_size
                                        , league_rank
                                        , m_league_size )
        , value );
    }

    // Reduce with final value at blockDim.y - 1 location.
    if ( cuda_single_inter_block_reduce_scan<false,FunctorType,WorkTag>(
           m_functor , blockIdx.x , gridDim.x ,
           kokkos_impl_cuda_shared_memory<size_type>() , m_scratch_space , m_scratch_flags ) ) {

      // This is the final block with the final result at the final threads' location

      size_type * const shared = kokkos_impl_cuda_shared_memory<size_type>() + ( blockDim.y - 1 ) * word_count.value ;
      size_type * const global = m_unified_space ? m_unified_space : m_scratch_space ;

      if ( threadIdx.y == 0 ) {
        Kokkos::Impl::FunctorFinal< FunctorType , WorkTag >::final( m_functor , shared );
      }

      if ( CudaTraits::WarpSize < word_count.value ) { __syncthreads(); }

      for ( unsigned i = threadIdx.y ; i < word_count.value ; i += blockDim.y ) { global[i] = shared[i]; }
    }
  }

  inline
  void execute()
    {
      const int block_count = std::min( m_league_size , m_team_size );

      m_scratch_space = cuda_internal_scratch_space( ValueTraits::value_size( m_functor ) * block_count );
      m_scratch_flags = cuda_internal_scratch_flags( sizeof(size_type) );
      m_unified_space = cuda_internal_scratch_unified( ValueTraits::value_size( m_functor ) );

      // REQUIRED DIMENSIONS ( 1 , N , 1 )
      const dim3 block( 1 , m_team_size , 1 );
      const dim3 grid( std::min( int(m_league_size) , int(m_team_size) ) , 1 , 1 );
      const int shmem_size_total = m_team_begin + m_shmem_begin + m_shmem_size ;

      CudaParallelLaunch< ParallelReduce >( *this, grid, block, shmem_size_total ); // copy to device and execute

      Cuda::fence();

      if ( m_result_ptr ) {
        if ( m_unified_space ) {
          const int count = ValueTraits::value_count( m_functor );
          for ( int i = 0 ; i < count ; ++i ) { m_result_ptr[i] = pointer_type(m_unified_space)[i] ; }
        }
        else {
          const int size = ValueTraits::value_size( m_functor );
          DeepCopy<HostSpace,CudaSpace>( m_result_ptr, m_scratch_space, size );
        }
      }
    }

  template< class HostViewType >
  ParallelReduce( const FunctorType  & arg_functor 
                , const Policy       & arg_policy 
                , const HostViewType & arg_result
                )
  : m_functor( arg_functor )
  , m_result_ptr( arg_result.ptr_on_device() )
  , m_scratch_space( 0 )
  , m_scratch_flags( 0 )
  , m_unified_space( 0 )
  , m_team_begin( 0 )
  , m_shmem_begin( 0 )
  , m_shmem_size( 0 )
  , m_league_size( arg_policy.league_size() )
  , m_team_size( 0 <= arg_policy.team_size() ? arg_policy.team_size() :
      Kokkos::Impl::cuda_get_opt_block_size< ParallelReduce >( arg_functor , arg_policy.vector_length(), arg_policy.scratch_size() ) / arg_policy.vector_length() )
  {
    // Return Init value if the number of worksets is zero
    if( arg_policy.league_size() == 0) {
      ValueInit::init( m_functor , arg_result.ptr_on_device() );
      return ;
    }

    m_team_begin = cuda_single_inter_block_reduce_scan_shmem<false,FunctorType,WorkTag>( arg_functor , m_team_size );
    m_shmem_begin = sizeof(double) * ( m_team_size + 2 );
    m_shmem_size = arg_policy.scratch_size() + FunctorTeamShmemSize< FunctorType >::value( arg_functor , m_team_size );

    // The global parallel_reduce does not support vector_length other than 1 at the moment
    if( arg_policy.vector_length() > 1)
      Impl::throw_runtime_exception( "Kokkos::parallel_reduce with a TeamPolicy using a vector length of greater than 1 is not currently supported for CUDA.");

    if( m_team_size < 32)
      Impl::throw_runtime_exception( "Kokkos::parallel_reduce with a TeamPolicy using a team_size smaller than 32 is not currently supported with CUDA.");

    // Functor's reduce memory, team scan memory, and team shared memory depend upon team size.

    const int shmem_size_total = m_team_begin + m_shmem_begin + m_shmem_size ;

    if ( ! Kokkos::Impl::is_integral_power_of_two( m_team_size ) ||
         CudaTraits::SharedMemoryCapacity < shmem_size_total ) {
      Kokkos::Impl::throw_runtime_exception(std::string("Kokkos::Impl::ParallelReduce< Cuda > bad team size"));
    }

    if ( m_team_size >
         Kokkos::Impl::cuda_get_max_block_size< ParallelReduce >
               ( arg_functor , arg_policy.vector_length(), arg_policy.scratch_size() ) / arg_policy.vector_length()) {
      Kokkos::Impl::throw_runtime_exception(std::string("Kokkos::Impl::ParallelReduce< Cuda > requested too large team size."));
    }

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
                  , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Cuda >
                  >
{
private:

  typedef Kokkos::RangePolicy<Arg0,Arg1,Arg2,Kokkos::Cuda>  Policy ;
  typedef typename Policy::member_type  Member ;
  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::WorkRange    WorkRange ;

  typedef Kokkos::Impl::FunctorValueTraits< FunctorType, WorkTag > ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   FunctorType, WorkTag > ValueInit ;
  typedef Kokkos::Impl::FunctorValueOps<    FunctorType, WorkTag > ValueOps ;

public:

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;
  typedef FunctorType                           functor_type ;
  typedef Cuda::size_type                       size_type ;

private:

  // Algorithmic constraints:
  //  (a) blockDim.y is a power of two
  //  (b) blockDim.y == blockDim.z == 1
  //  (c) gridDim.x  <= blockDim.y * blockDim.y
  //  (d) gridDim.y  == gridDim.z == 1

  const FunctorType m_functor ;
  const Policy      m_policy ;
  size_type *       m_scratch_space ;
  size_type *       m_scratch_flags ;
  size_type         m_final ;

  template< class TagType >
  __device__ inline
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_range( const Member & i , reference_type update , const bool final ) const
    { m_functor( i , update , final ); }

  template< class TagType >
  __device__ inline
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec_range( const Member & i , reference_type update , const bool final ) const
    { m_functor( TagType() , i , update , final ); }

  //----------------------------------------

  __device__ inline
  void initial(void) const
  {
    const integral_nonzero_constant< size_type , ValueTraits::StaticValueSize / sizeof(size_type) >
      word_count( ValueTraits::value_size( m_functor ) / sizeof(size_type) );

    size_type * const shared_value = kokkos_impl_cuda_shared_memory<size_type>() + word_count.value * threadIdx.y ;

    ValueInit::init( m_functor , shared_value );

    // Number of blocks is bounded so that the reduction can be limited to two passes.
    // Each thread block is given an approximately equal amount of work to perform.
    // Accumulate the values for this block.
    // The accumulation ordering does not match the final pass, but is arithmatically equivalent.

    const WorkRange range( m_policy , blockIdx.x , gridDim.x );

    for ( Member iwork = range.begin() + threadIdx.y , iwork_end = range.end() ;
          iwork < iwork_end ; iwork += blockDim.y ) {
      this-> template exec_range< WorkTag >( iwork , ValueOps::reference( shared_value ) , false );
    }

    // Reduce and scan, writing out scan of blocks' totals and block-groups' totals.
    // Blocks' scan values are written to 'blockIdx.x' location.
    // Block-groups' scan values are at: i = ( j * blockDim.y - 1 ) for i < gridDim.x
    cuda_single_inter_block_reduce_scan<true,FunctorType,WorkTag>( m_functor , blockIdx.x , gridDim.x , kokkos_impl_cuda_shared_memory<size_type>() , m_scratch_space , m_scratch_flags );
  }

  //----------------------------------------

  __device__ inline
  void final(void) const
  {
    const integral_nonzero_constant< size_type , ValueTraits::StaticValueSize / sizeof(size_type) >
      word_count( ValueTraits::value_size( m_functor ) / sizeof(size_type) );

    // Use shared memory as an exclusive scan: { 0 , value[0] , value[1] , value[2] , ... }
    size_type * const shared_data   = kokkos_impl_cuda_shared_memory<size_type>();
    size_type * const shared_prefix = shared_data + word_count.value * threadIdx.y ;
    size_type * const shared_accum  = shared_data + word_count.value * ( blockDim.y + 1 );

    // Starting value for this thread block is the previous block's total.
    if ( blockIdx.x ) {
      size_type * const block_total = m_scratch_space + word_count.value * ( blockIdx.x - 1 );
      for ( unsigned i = threadIdx.y ; i < word_count.value ; ++i ) { shared_accum[i] = block_total[i] ; }
    }
    else if ( 0 == threadIdx.y ) {
      ValueInit::init( m_functor , shared_accum );
    }

    const WorkRange range( m_policy , blockIdx.x , gridDim.x );

    for ( typename Policy::member_type iwork_base = range.begin(); iwork_base < range.end() ; iwork_base += blockDim.y ) {

      const typename Policy::member_type iwork = iwork_base + threadIdx.y ;

      __syncthreads(); // Don't overwrite previous iteration values until they are used

      ValueInit::init( m_functor , shared_prefix + word_count.value );

      // Copy previous block's accumulation total into thread[0] prefix and inclusive scan value of this block
      for ( unsigned i = threadIdx.y ; i < word_count.value ; ++i ) {
        shared_data[i + word_count.value] = shared_data[i] = shared_accum[i] ;
      }

      if ( CudaTraits::WarpSize < word_count.value ) { __syncthreads(); } // Protect against large scan values.

      // Call functor to accumulate inclusive scan value for this work item
      if ( iwork < range.end() ) {
        this-> template exec_range< WorkTag >( iwork , ValueOps::reference( shared_prefix + word_count.value ) , false );
      }

      // Scan block values into locations shared_data[1..blockDim.y]
      cuda_intra_block_reduce_scan<true,FunctorType,WorkTag>( m_functor , ValueTraits::pointer_type(shared_data+word_count.value) );

      {
        size_type * const block_total = shared_data + word_count.value * blockDim.y ;
        for ( unsigned i = threadIdx.y ; i < word_count.value ; ++i ) { shared_accum[i] = block_total[i]; }
      }

      // Call functor with exclusive scan value
      if ( iwork < range.end() ) {
        this-> template exec_range< WorkTag >( iwork , ValueOps::reference( shared_prefix ) , true );
      }
    }
  }

public:

  //----------------------------------------

  __device__ inline
  void operator()(void) const
  {
    if ( ! m_final ) {
      initial();
    }
    else {
      final();
    }
  }

  // Determine block size constrained by shared memory:
  static inline
  unsigned local_block_size( const FunctorType & f )
    {
      // blockDim.y must be power of two = 128 (4 warps) or 256 (8 warps) or 512 (16 warps)
      // gridDim.x <= blockDim.y * blockDim.y
      //
      // 4 warps was 10% faster than 8 warps and 20% faster than 16 warps in unit testing

      unsigned n = CudaTraits::WarpSize * 4 ;
      while ( n && CudaTraits::SharedMemoryCapacity < cuda_single_inter_block_reduce_scan_shmem<false,FunctorType,WorkTag>( f , n ) ) { n >>= 1 ; }
      return n ;
    }

  inline
  void execute()
    {
      enum { GridMaxComputeCapability_2x = 0x0ffff };

      const int block_size = local_block_size( m_functor );

      const int grid_max =
        ( block_size * block_size ) < GridMaxComputeCapability_2x ?
        ( block_size * block_size ) : GridMaxComputeCapability_2x ;

      // At most 'max_grid' blocks:
      const int nwork    = m_policy.end() - m_policy.begin();
      const int max_grid = std::min( int(grid_max) , int(( nwork + block_size - 1 ) / block_size ));

      // How much work per block:
      const int work_per_block = ( nwork + max_grid - 1 ) / max_grid ;

      // How many block are really needed for this much work:
      const int grid_x = ( nwork + work_per_block - 1 ) / work_per_block ;

      m_scratch_space = cuda_internal_scratch_space( ValueTraits::value_size( m_functor ) * grid_x );
      m_scratch_flags = cuda_internal_scratch_flags( sizeof(size_type) * 1 );

      const dim3 grid( grid_x , 1 , 1 );
      const dim3 block( 1 , block_size , 1 ); // REQUIRED DIMENSIONS ( 1 , N , 1 )
      const int shmem = ValueTraits::value_size( m_functor ) * ( block_size + 2 );

      m_final = false ;
      CudaParallelLaunch< ParallelScan >( *this, grid, block, shmem ); // copy to device and execute

      m_final = true ;
      CudaParallelLaunch< ParallelScan >( *this, grid, block, shmem ); // copy to device and execute
    }

  ParallelScan( const FunctorType  & arg_functor ,
                const Policy       & arg_policy )
  : m_functor( arg_functor )
  , m_policy( arg_policy )
  , m_scratch_space( 0 )
  , m_scratch_flags( 0 )
  , m_final( false )
  { }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {
  template<typename iType>
  struct TeamThreadRangeBoundariesStruct<iType,CudaTeamMember> {
    typedef iType index_type;
    const iType start;
    const iType end;
    const iType increment;
    const CudaTeamMember& thread;

#ifdef __CUDA_ARCH__
    __device__ inline
    TeamThreadRangeBoundariesStruct (const CudaTeamMember& thread_, const iType& count):
      start( threadIdx.y ),
      end( count ),
      increment( blockDim.y ),
      thread(thread_)
    {}
    __device__ inline
    TeamThreadRangeBoundariesStruct (const CudaTeamMember& thread_, const iType& begin_, const iType& end_):
      start( begin_+threadIdx.y ),
      end( end_ ),
      increment( blockDim.y ),
      thread(thread_)
    {}
#else
    KOKKOS_INLINE_FUNCTION
    TeamThreadRangeBoundariesStruct (const CudaTeamMember& thread_, const iType& count):
      start( 0 ),
      end( count ),
      increment( 1 ),
      thread(thread_)
    {}
    KOKKOS_INLINE_FUNCTION
    TeamThreadRangeBoundariesStruct (const CudaTeamMember& thread_,  const iType& begin_, const iType& end_):
      start( begin_ ),
      end( end_ ),
      increment( 1 ),
      thread(thread_)
    {}
#endif
  };

  template<typename iType>
  struct ThreadVectorRangeBoundariesStruct<iType,CudaTeamMember> {
    typedef iType index_type;
    const iType start;
    const iType end;
    const iType increment;

#ifdef __CUDA_ARCH__
    __device__ inline
    ThreadVectorRangeBoundariesStruct (const CudaTeamMember& thread, const iType& count):
    start( threadIdx.x ),
    end( count ),
    increment( blockDim.x )
    {}
#else
    KOKKOS_INLINE_FUNCTION
    ThreadVectorRangeBoundariesStruct (const CudaTeamMember& thread_, const iType& count):
      start( 0 ),
      end( count ),
      increment( 1 )
    {}
#endif
    };

} // namespace Impl

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::TeamThreadRangeBoundariesStruct<iType,Impl::CudaTeamMember>
  TeamThreadRange(const Impl::CudaTeamMember& thread, const iType& count) {
  return Impl::TeamThreadRangeBoundariesStruct<iType,Impl::CudaTeamMember>(thread,count);
}

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::TeamThreadRangeBoundariesStruct<iType,Impl::CudaTeamMember>
  TeamThreadRange(const Impl::CudaTeamMember& thread, const iType& begin, const iType& end) {
  return Impl::TeamThreadRangeBoundariesStruct<iType,Impl::CudaTeamMember>(thread,begin,end);
}

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::CudaTeamMember >
  ThreadVectorRange(const Impl::CudaTeamMember& thread, const iType& count) {
  return Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::CudaTeamMember >(thread,count);
}

KOKKOS_INLINE_FUNCTION
Impl::ThreadSingleStruct<Impl::CudaTeamMember> PerTeam(const Impl::CudaTeamMember& thread) {
  return Impl::ThreadSingleStruct<Impl::CudaTeamMember>(thread);
}

KOKKOS_INLINE_FUNCTION
Impl::VectorSingleStruct<Impl::CudaTeamMember> PerThread(const Impl::CudaTeamMember& thread) {
  return Impl::VectorSingleStruct<Impl::CudaTeamMember>(thread);
}

} // namespace Kokkos

namespace Kokkos {

  /** \brief  Inter-thread parallel_for. Executes lambda(iType i) for each i=0..N-1.
   *
   * The range i=0..N-1 is mapped to all threads of the the calling thread team.
   * This functionality requires C++11 support.*/
template<typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION
void parallel_for(const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::CudaTeamMember>& loop_boundaries, const Lambda& lambda) {
  #ifdef __CUDA_ARCH__
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment)
    lambda(i);
  #endif
}

/** \brief  Inter-thread vector parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team and a summation of
 * val is performed and put into result. This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ValueType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::CudaTeamMember>& loop_boundaries,
                     const Lambda & lambda, ValueType& result) {

#ifdef __CUDA_ARCH__
  result = ValueType();

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    lambda(i,result);
  }

  Impl::cuda_intra_warp_reduction(result,[&] (ValueType& dst, const ValueType& src) { dst+=src; });
  Impl::cuda_inter_warp_reduction(result,[&] (ValueType& dst, const ValueType& src) { dst+=src; });

#endif
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
void parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::CudaTeamMember>& loop_boundaries,
                     const Lambda & lambda, const JoinType& join, ValueType& init_result) {

#ifdef __CUDA_ARCH__
  ValueType result = init_result;

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    lambda(i,result);
  }

  Impl::cuda_intra_warp_reduction(result, join );
  Impl::cuda_inter_warp_reduction(result, join );

  init_result = result;
#endif
}

} //namespace Kokkos

namespace Kokkos {
/** \brief  Intra-thread vector parallel_for. Executes lambda(iType i) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread.
 * This functionality requires C++11 support.*/
template<typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION
void parallel_for(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::CudaTeamMember >&
    loop_boundaries, const Lambda& lambda) {
#ifdef __CUDA_ARCH__
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment)
    lambda(i);
#endif
}

/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread and a summation of
 * val is performed and put into result. This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ValueType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::CudaTeamMember >&
      loop_boundaries, const Lambda & lambda, ValueType& result) {
#ifdef __CUDA_ARCH__
  ValueType val = ValueType();

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    lambda(i,val);
  }

  result = val;

  if (loop_boundaries.increment > 1)
    result += shfl_down(result, 1,loop_boundaries.increment);
  if (loop_boundaries.increment > 2)
    result += shfl_down(result, 2,loop_boundaries.increment);
  if (loop_boundaries.increment > 4)
    result += shfl_down(result, 4,loop_boundaries.increment);
  if (loop_boundaries.increment > 8)
    result += shfl_down(result, 8,loop_boundaries.increment);
  if (loop_boundaries.increment > 16)
    result += shfl_down(result, 16,loop_boundaries.increment);

  result = shfl(result,0,loop_boundaries.increment);
#endif
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
void parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::CudaTeamMember >&
      loop_boundaries, const Lambda & lambda, const JoinType& join, ValueType& init_result) {

#ifdef __CUDA_ARCH__
  ValueType result = init_result;

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    lambda(i,result);
  }

  if (loop_boundaries.increment > 1)
    join( result, shfl_down(result, 1,loop_boundaries.increment));
  if (loop_boundaries.increment > 2)
    join( result, shfl_down(result, 2,loop_boundaries.increment));
  if (loop_boundaries.increment > 4)
    join( result, shfl_down(result, 4,loop_boundaries.increment));
  if (loop_boundaries.increment > 8)
    join( result, shfl_down(result, 8,loop_boundaries.increment));
  if (loop_boundaries.increment > 16)
    join( result, shfl_down(result, 16,loop_boundaries.increment));

  init_result = shfl(result,0,loop_boundaries.increment);
#endif
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
void parallel_scan(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::CudaTeamMember >&
      loop_boundaries, const FunctorType & lambda) {

#ifdef __CUDA_ARCH__
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , void > ValueTraits ;
  typedef typename ValueTraits::value_type value_type ;

  value_type scan_val = value_type();
  const int VectorLength = blockDim.x;

  iType loop_bound = ((loop_boundaries.end+VectorLength-1)/VectorLength) * VectorLength;
  for(int _i = threadIdx.x; _i < loop_bound; _i += VectorLength) {
    value_type val = value_type();
    if(_i<loop_boundaries.end)
      lambda(_i , val , false);

    value_type tmp = val;
    value_type result_i;

    if(threadIdx.x%VectorLength == 0)
      result_i = tmp;
    if (VectorLength > 1) {
      const value_type tmp2 = shfl_up(tmp, 1,VectorLength);
      if(threadIdx.x > 0)
        tmp+=tmp2;
    }
    if(threadIdx.x%VectorLength == 1)
      result_i = tmp;
    if (VectorLength > 3) {
      const value_type tmp2 = shfl_up(tmp, 2,VectorLength);
      if(threadIdx.x > 1)
        tmp+=tmp2;
    }
    if ((threadIdx.x%VectorLength >= 2) &&
        (threadIdx.x%VectorLength < 4))
      result_i = tmp;
    if (VectorLength > 7) {
      const value_type tmp2 = shfl_up(tmp, 4,VectorLength);
      if(threadIdx.x > 3)
        tmp+=tmp2;
    }
    if ((threadIdx.x%VectorLength >= 4) &&
        (threadIdx.x%VectorLength < 8))
      result_i = tmp;
    if (VectorLength > 15) {
      const value_type tmp2 = shfl_up(tmp, 8,VectorLength);
      if(threadIdx.x > 7)
        tmp+=tmp2;
    }
    if ((threadIdx.x%VectorLength >= 8) &&
        (threadIdx.x%VectorLength < 16))
      result_i = tmp;
    if (VectorLength > 31) {
      const value_type tmp2 = shfl_up(tmp, 16,VectorLength);
      if(threadIdx.x > 15)
        tmp+=tmp2;
    }
    if (threadIdx.x%VectorLength >= 16)
      result_i = tmp;

    val = scan_val + result_i - val;
    scan_val += shfl(tmp,VectorLength-1,VectorLength);
    if(_i<loop_boundaries.end)
      lambda(_i , val , true);
  }
#endif
}

}

namespace Kokkos {

template<class FunctorType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::VectorSingleStruct<Impl::CudaTeamMember>& , const FunctorType& lambda) {
#ifdef __CUDA_ARCH__
  if(threadIdx.x == 0) lambda();
#endif
}

template<class FunctorType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::ThreadSingleStruct<Impl::CudaTeamMember>& , const FunctorType& lambda) {
#ifdef __CUDA_ARCH__
  if(threadIdx.x == 0 && threadIdx.y == 0) lambda();
#endif
}

template<class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::VectorSingleStruct<Impl::CudaTeamMember>& , const FunctorType& lambda, ValueType& val) {
#ifdef __CUDA_ARCH__
  if(threadIdx.x == 0) lambda(val);
  val = shfl(val,0,blockDim.x);
#endif
}

template<class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::ThreadSingleStruct<Impl::CudaTeamMember>& single_struct, const FunctorType& lambda, ValueType& val) {
#ifdef __CUDA_ARCH__
  if(threadIdx.x == 0 && threadIdx.y == 0) {
    lambda(val);
  }
  single_struct.team_member.team_broadcast(val,0);
#endif
}

}

namespace Kokkos {

namespace Impl {
  template< class FunctorType, class ExecPolicy, class ValueType , class Tag = typename ExecPolicy::work_tag>
  struct CudaFunctorAdapter {
    const FunctorType f;
    typedef ValueType value_type;
    CudaFunctorAdapter(const FunctorType& f_):f(f_) {}

    __device__ inline
    void operator() (typename ExecPolicy::work_tag, const typename ExecPolicy::member_type& i, ValueType& val) const {
      //Insert Static Assert with decltype on ValueType equals third argument type of FunctorType::operator()
      f(typename ExecPolicy::work_tag(), i,val);
    }
  };

  template< class FunctorType, class ExecPolicy, class ValueType >
  struct CudaFunctorAdapter<FunctorType,ExecPolicy,ValueType,void> {
    const FunctorType f;
    typedef ValueType value_type;
    CudaFunctorAdapter(const FunctorType& f_):f(f_) {}

    __device__ inline
    void operator() (const typename ExecPolicy::member_type& i, ValueType& val) const {
      //Insert Static Assert with decltype on ValueType equals second argument type of FunctorType::operator()
      f(i,val);
    }

  };

  template< class FunctorType, class Enable = void>
  struct ReduceFunctorHasInit {
    enum {value = false};
  };

  template< class FunctorType>
  struct ReduceFunctorHasInit<FunctorType, typename Impl::enable_if< 0 < sizeof( & FunctorType::init ) >::type > {
    enum {value = true};
  };

  template< class FunctorType, class Enable = void>
  struct ReduceFunctorHasJoin {
    enum {value = false};
  };

  template< class FunctorType>
  struct ReduceFunctorHasJoin<FunctorType, typename Impl::enable_if< 0 < sizeof( & FunctorType::join ) >::type > {
    enum {value = true};
  };

  template< class FunctorType, class Enable = void>
  struct ReduceFunctorHasFinal {
    enum {value = false};
  };

  template< class FunctorType>
  struct ReduceFunctorHasFinal<FunctorType, typename Impl::enable_if< 0 < sizeof( & FunctorType::final ) >::type > {
    enum {value = true};
  };

  template< class FunctorType, bool Enable =
      ( FunctorDeclaresValueType<FunctorType,void>::value) ||
      ( ReduceFunctorHasInit<FunctorType>::value  ) ||
      ( ReduceFunctorHasJoin<FunctorType>::value  ) ||
      ( ReduceFunctorHasFinal<FunctorType>::value )
      >
  struct IsNonTrivialReduceFunctor {
    enum {value = false};
  };

  template< class FunctorType>
  struct IsNonTrivialReduceFunctor<FunctorType, true> {
    enum {value = true};
  };

  template<class FunctorType, class ResultType, class Tag, bool Enable = IsNonTrivialReduceFunctor<FunctorType>::value >
  struct FunctorReferenceType {
    typedef ResultType& reference_type;
  };

  template<class FunctorType, class ResultType, class Tag>
  struct FunctorReferenceType<FunctorType, ResultType, Tag, true> {
    typedef typename Kokkos::Impl::FunctorValueTraits< FunctorType ,Tag >::reference_type reference_type;
  };

}

// general policy and view ouput
template< class ExecPolicy , class FunctorTypeIn , class ViewType >
inline
void parallel_reduce( const ExecPolicy  & policy
                    , const FunctorTypeIn & functor_in
                    , const ViewType    & result_view
                    , const std::string& str = "" 
                    , typename Impl::enable_if<
                      ( Kokkos::is_view<ViewType>::value && ! Impl::is_integral< ExecPolicy >::value &&
                        Impl::is_same<typename ExecPolicy::execution_space,Kokkos::Cuda>::value
                      )>::type * = 0 )
{
  enum {FunctorHasValueType = Impl::IsNonTrivialReduceFunctor<FunctorTypeIn>::value };
  typedef typename Kokkos::Impl::if_c<FunctorHasValueType, FunctorTypeIn, Impl::CudaFunctorAdapter<FunctorTypeIn,ExecPolicy,typename ViewType::value_type> >::type FunctorType;
  FunctorType functor = Impl::if_c<FunctorHasValueType,FunctorTypeIn,FunctorType>::select(functor_in,FunctorType(functor_in));

#ifdef KOKKOSP_ENABLE_PROFILING
    uint64_t kpID = 0;
     if(Kokkos::Experimental::profileLibraryLoaded()) {
	Kokkos::Experimental::beginParallelScan("" == str ? typeid(FunctorType).name() : str, 0, &kpID);
     }
#endif

  Kokkos::Impl::shared_allocation_tracking_claim_and_disable();
  Impl::ParallelReduce< FunctorType, ExecPolicy > closure( functor , policy , result_view );
  Kokkos::Impl::shared_allocation_tracking_release_and_enable();

  closure.execute();
    
#ifdef KOKKOSP_ENABLE_PROFILING
     if(Kokkos::Experimental::profileLibraryLoaded()) {
	Kokkos::Experimental::endParallelScan(kpID);
     }
#endif
}

// general policy and pod or array of pod output
template< class ExecPolicy , class FunctorTypeIn , class ResultType>
inline
void parallel_reduce( const ExecPolicy  & policy
                    , const FunctorTypeIn & functor_in
                    , ResultType& result_ref
                    , const std::string& str = "" 
                    , typename Impl::enable_if<
                      ( ! Kokkos::is_view<ResultType>::value &&
                        ! Impl::IsNonTrivialReduceFunctor<FunctorTypeIn>::value &&
                        ! Impl::is_integral< ExecPolicy >::value  &&
                          Impl::is_same<typename ExecPolicy::execution_space,Kokkos::Cuda>::value )>::type * = 0 )
{
  typedef typename Impl::CudaFunctorAdapter<FunctorTypeIn,ExecPolicy,ResultType> FunctorType;

  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , typename ExecPolicy::work_tag >  ValueTraits ;
  typedef Kokkos::Impl::FunctorValueOps<    FunctorType , typename ExecPolicy::work_tag >  ValueOps ;

  // Wrap the result output request in a view to inform the implementation
  // of the type and memory space.

  typedef typename Kokkos::Impl::if_c< (ValueTraits::StaticValueSize != 0)
                                     , typename ValueTraits::value_type
                                     , typename ValueTraits::pointer_type
                                     >::type value_type ;
  Kokkos::View< value_type
              , HostSpace
              , Kokkos::MemoryUnmanaged
              >
    result_view( ValueOps::pointer( result_ref )
               , 1
               );

#ifdef KOKKOSP_ENABLE_PROFILING
    uint64_t kpID = 0;
     if(Kokkos::Experimental::profileLibraryLoaded()) {
	Kokkos::Experimental::beginParallelScan("" == str ? typeid(FunctorType).name() : str, 0, &kpID);
     }
#endif
    
  Kokkos::Impl::shared_allocation_tracking_claim_and_disable();
  Impl::ParallelReduce< FunctorType, ExecPolicy > closure( FunctorType(functor_in) , policy , result_view );
  Kokkos::Impl::shared_allocation_tracking_release_and_enable();

  closure.execute();
    
#ifdef KOKKOSP_ENABLE_PROFILING
     if(Kokkos::Experimental::profileLibraryLoaded()) {
	Kokkos::Experimental::endParallelScan(kpID);
     }
#endif
}

// general policy and pod or array of pod output
template< class ExecPolicy , class FunctorType>
inline
void parallel_reduce( const ExecPolicy  & policy
                    , const FunctorType & functor
                    , typename Kokkos::Impl::FunctorValueTraits< FunctorType , typename ExecPolicy::work_tag >::reference_type result_ref
                    , const std::string& str = "" 
                    , typename Impl::enable_if<
                      (   Impl::IsNonTrivialReduceFunctor<FunctorType>::value &&
                        ! Impl::is_integral< ExecPolicy >::value  &&
                          Impl::is_same<typename ExecPolicy::execution_space,Kokkos::Cuda>::value )>::type * = 0 )
{
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , typename ExecPolicy::work_tag >  ValueTraits ;
  typedef Kokkos::Impl::FunctorValueOps<    FunctorType , typename ExecPolicy::work_tag >  ValueOps ;

  // Wrap the result output request in a view to inform the implementation
  // of the type and memory space.

  typedef typename Kokkos::Impl::if_c< (ValueTraits::StaticValueSize != 0)
                                     , typename ValueTraits::value_type
                                     , typename ValueTraits::pointer_type
                                     >::type value_type ;

  Kokkos::View< value_type
              , HostSpace
              , Kokkos::MemoryUnmanaged
              >
    result_view( ValueOps::pointer( result_ref )
               , ValueTraits::value_count( functor )
               );

#ifdef KOKKOSP_ENABLE_PROFILING
    uint64_t kpID = 0;
     if(Kokkos::Experimental::profileLibraryLoaded()) {
	Kokkos::Experimental::beginParallelScan("" == str ? typeid(FunctorType).name() : str, 0, &kpID);
     }
#endif
    
  Kokkos::Impl::shared_allocation_tracking_claim_and_disable();
  Impl::ParallelReduce< FunctorType, ExecPolicy > closure( functor , policy , result_view );
  Kokkos::Impl::shared_allocation_tracking_release_and_enable();

  closure.execute();
    
#ifdef KOKKOSP_ENABLE_PROFILING
     if(Kokkos::Experimental::profileLibraryLoaded()) {
	Kokkos::Experimental::endParallelScan(kpID);
     }
#endif
}

// integral range policy and view ouput
template< class FunctorTypeIn , class ViewType >
inline
void parallel_reduce( const size_t        work_count
                    , const FunctorTypeIn & functor_in
                    , const ViewType    & result_view
                    , const std::string& str = "" 
                    , typename Impl::enable_if<( Kokkos::is_view<ViewType>::value &&
                                                 Impl::is_same<
                          typename Impl::FunctorPolicyExecutionSpace< FunctorTypeIn , void >::execution_space,
                          Kokkos::Cuda>::value
                        )>::type * = 0 )
{
  enum {FunctorHasValueType = Impl::IsNonTrivialReduceFunctor<FunctorTypeIn>::value };
  typedef typename
    Impl::FunctorPolicyExecutionSpace< FunctorTypeIn , void >::execution_space
      execution_space ;

  typedef RangePolicy< execution_space > ExecPolicy ;

  typedef typename Kokkos::Impl::if_c<FunctorHasValueType, FunctorTypeIn, Impl::CudaFunctorAdapter<FunctorTypeIn,ExecPolicy,typename ViewType::value_type> >::type FunctorType;

  FunctorType functor = Impl::if_c<FunctorHasValueType,FunctorTypeIn,FunctorType>::select(functor_in,FunctorType(functor_in));

#ifdef KOKKOSP_ENABLE_PROFILING
    uint64_t kpID = 0;
     if(Kokkos::Experimental::profileLibraryLoaded()) {
	Kokkos::Experimental::beginParallelScan("" == str ? typeid(FunctorType).name() : str, 0, &kpID);
     }
#endif
    
  Kokkos::Impl::shared_allocation_tracking_claim_and_disable();
  Impl::ParallelReduce< FunctorType, ExecPolicy > closure( functor , ExecPolicy(0,work_count) , result_view );
  Kokkos::Impl::shared_allocation_tracking_release_and_enable();

  closure.execute();

#ifdef KOKKOSP_ENABLE_PROFILING
     if(Kokkos::Experimental::profileLibraryLoaded()) {
	Kokkos::Experimental::endParallelScan(kpID);
     }
#endif

}

// integral range policy and pod or array of pod output
template< class FunctorTypeIn , class ResultType>
inline
void parallel_reduce( const size_t        work_count
                    , const FunctorTypeIn & functor_in
                    , ResultType& result
                    , const std::string& str = "" 
                    , typename Impl::enable_if< ! Kokkos::is_view<ResultType>::value &&
                                                ! Impl::IsNonTrivialReduceFunctor<FunctorTypeIn>::value &&
                                                Impl::is_same<
                             typename Impl::FunctorPolicyExecutionSpace< FunctorTypeIn , void >::execution_space,
                             Kokkos::Cuda>::value >::type * = 0 )
{
  typedef typename
    Kokkos::Impl::FunctorPolicyExecutionSpace< FunctorTypeIn , void >::execution_space
      execution_space ;
  typedef Kokkos::RangePolicy< execution_space > ExecPolicy ;

  typedef Impl::CudaFunctorAdapter<FunctorTypeIn,ExecPolicy,ResultType> FunctorType;


  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , void >  ValueTraits ;
  typedef Kokkos::Impl::FunctorValueOps<    FunctorType , void >  ValueOps ;


  // Wrap the result output request in a view to inform the implementation
  // of the type and memory space.

  typedef typename Kokkos::Impl::if_c< (ValueTraits::StaticValueSize != 0)
                                     , typename ValueTraits::value_type
                                     , typename ValueTraits::pointer_type
                                     >::type value_type ;

  Kokkos::View< value_type
              , HostSpace
              , Kokkos::MemoryUnmanaged
              >
    result_view( ValueOps::pointer( result )
               , 1
               );

#ifdef KOKKOSP_ENABLE_PROFILING
    uint64_t kpID = 0;
     if(Kokkos::Experimental::profileLibraryLoaded()) {
	Kokkos::Experimental::beginParallelScan("" == str ? typeid(FunctorType).name() : str, 0, &kpID);
     }
#endif
    
  Kokkos::Impl::shared_allocation_tracking_claim_and_disable();
  Impl::ParallelReduce< FunctorType , ExecPolicy > closure( FunctorType(functor_in) , ExecPolicy(0,work_count) , result_view );
  Kokkos::Impl::shared_allocation_tracking_release_and_enable();

  closure.execute();
    
#ifdef KOKKOSP_ENABLE_PROFILING
     if(Kokkos::Experimental::profileLibraryLoaded()) {
	Kokkos::Experimental::endParallelScan(kpID);
     }
#endif
}

template< class FunctorType>
inline
void parallel_reduce( const size_t        work_count
                    , const FunctorType & functor
                    , typename Kokkos::Impl::FunctorValueTraits< FunctorType , void >::reference_type result
                    , const std::string& str = "" 
                    , typename Impl::enable_if< Impl::IsNonTrivialReduceFunctor<FunctorType>::value &&
                                                Impl::is_same<
                             typename Impl::FunctorPolicyExecutionSpace< FunctorType , void >::execution_space,
                             Kokkos::Cuda>::value >::type * = 0 )
{

  typedef typename
    Kokkos::Impl::FunctorPolicyExecutionSpace< FunctorType , void >::execution_space
      execution_space ;
  typedef Kokkos::RangePolicy< execution_space > ExecPolicy ;



  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , void >  ValueTraits ;
  typedef Kokkos::Impl::FunctorValueOps<    FunctorType , void >  ValueOps ;


  // Wrap the result output request in a view to inform the implementation
  // of the type and memory space.

  typedef typename Kokkos::Impl::if_c< (ValueTraits::StaticValueSize != 0)
                                     , typename ValueTraits::value_type
                                     , typename ValueTraits::pointer_type
                                     >::type value_type ;

  Kokkos::View< value_type
              , HostSpace
              , Kokkos::MemoryUnmanaged
              >
    result_view( ValueOps::pointer( result )
               , ValueTraits::value_count( functor )
               );

#ifdef KOKKOSP_ENABLE_PROFILING
    uint64_t kpID = 0;
     if(Kokkos::Experimental::profileLibraryLoaded()) {
	Kokkos::Experimental::beginParallelScan("" == str ? typeid(FunctorType).name() : str, 0, &kpID);
     }
#endif
    
  Kokkos::Impl::shared_allocation_tracking_claim_and_disable();
  Impl::ParallelReduce< FunctorType , ExecPolicy > closure( functor , ExecPolicy(0,work_count) , result_view );
  Kokkos::Impl::shared_allocation_tracking_release_and_enable();

  closure.execute();
    
#ifdef KOKKOSP_ENABLE_PROFILING
     if(Kokkos::Experimental::profileLibraryLoaded()) {
	Kokkos::Experimental::endParallelScan(kpID);
     }
#endif
}

#ifdef KOKKOS_HAVE_CUDA
template< class ExecPolicy , class FunctorType , class ResultType >
inline
void parallel_reduce( const std::string & str
                    , const ExecPolicy  & policy
                    , const FunctorType & functor
                    , ResultType * result)
{
  #if KOKKOS_ENABLE_DEBUG_PRINT_KERNEL_NAMES
  Kokkos::fence();
  std::cout << "KOKKOS_DEBUG Start parallel_reduce kernel: " << str << std::endl;
  #endif

  parallel_reduce(policy,functor,result,str);

  #if KOKKOS_ENABLE_DEBUG_PRINT_KERNEL_NAMES
  Kokkos::fence();
  std::cout << "KOKKOS_DEBUG End   parallel_reduce kernel: " << str << std::endl;
  #endif
  (void) str;
}

template< class ExecPolicy , class FunctorType , class ResultType >
inline
void parallel_reduce( const std::string & str
                    , const ExecPolicy  & policy
                    , const FunctorType & functor
                    , ResultType & result)
{
  #if KOKKOS_ENABLE_DEBUG_PRINT_KERNEL_NAMES
  Kokkos::fence();
  std::cout << "KOKKOS_DEBUG Start parallel_reduce kernel: " << str << std::endl;
  #endif

  parallel_reduce(policy,functor,result,str);

  #if KOKKOS_ENABLE_DEBUG_PRINT_KERNEL_NAMES
  Kokkos::fence();
  std::cout << "KOKKOS_DEBUG End   parallel_reduce kernel: " << str << std::endl;
  #endif
  (void) str;
}

template< class ExecPolicy , class FunctorType >
inline
void parallel_reduce( const std::string & str
                    , const ExecPolicy  & policy
                    , const FunctorType & functor)
{
  #if KOKKOS_ENABLE_DEBUG_PRINT_KERNEL_NAMES
  Kokkos::fence();
  std::cout << "KOKKOS_DEBUG Start parallel_reduce kernel: " << str << std::endl;
  #endif

  parallel_reduce(policy,functor,str);

  #if KOKKOS_ENABLE_DEBUG_PRINT_KERNEL_NAMES
  Kokkos::fence();
  std::cout << "KOKKOS_DEBUG End   parallel_reduce kernel: " << str << std::endl;
  #endif
  (void) str;
}
#endif
} // namespace Kokkos
#endif /* defined( __CUDACC__ ) */

#endif /* #ifndef KOKKOS_CUDA_PARALLEL_HPP */

