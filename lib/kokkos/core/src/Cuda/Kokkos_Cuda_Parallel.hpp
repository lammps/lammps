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

#ifndef KOKKOS_CUDA_PARALLEL_HPP
#define KOKKOS_CUDA_PARALLEL_HPP

#include <iostream>
#include <stdio.h>

#if defined( __CUDACC__ )

#include <utility>
#include <Kokkos_Parallel.hpp>

#include <Cuda/Kokkos_CudaExec.hpp>
#include <Cuda/Kokkos_Cuda_ReduceScan.hpp>
#include <Cuda/Kokkos_Cuda_Internal.hpp>
#include <Kokkos_Vectorization.hpp>
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


#ifdef KOKKOS_HAVE_CXX11
  template< class Operation >
  __device__ inline void vector_single(const Operation & op) const {
    if(threadIdx.x == 0)
      op();
  }

  template< class Operation, typename ValueType>
  __device__ inline void vector_single(const Operation & op, ValueType& bcast) const {
    if(threadIdx.x == 0)
      op();
    bcast = shfl(bcast,0,blockDim.x);
  }

#endif

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

#ifdef KOKKOS_HAVE_CXX11
  template< class Operation >
  void vector_single(const Operation & op) const {}

  template< class Operation , typename ValueType>
  void vector_single(const Operation & op, ValueType& val) const {}
#endif
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

  inline static
  int vector_length_max()
    { return Impl::CudaTraits::WarpSize; }

  //----------------------------------------

  inline int vector_length()   const { return m_vector_length ; }
  inline int team_size()   const { return m_team_size ; }
  inline int league_size() const { return m_league_size ; }

  /** \brief  Specify league size, request team size */
  TeamPolicy( execution_space & , int league_size , int team_size_request , int vector_length_request = 1 )
    : m_league_size( league_size )
    , m_team_size( team_size_request )
    , m_vector_length ( vector_length_request )
    {
      // Allow only power-of-two vector_length
      int check = 0;
      for(int k = 1; k < vector_length_max(); k*=2)
        if(k == vector_length_request)
          check = 1;
      if(!check)
        Impl::throw_runtime_exception( "Requested non-power-of-two vector length for TeamPolicy.");

      // Make sure league size is permissable
      if(league_size >= int(Impl::cuda_internal_maximum_grid_count()))
        Impl::throw_runtime_exception( "Requested too large league_size for TeamPolicy on Cuda execution space.");
    }

  TeamPolicy( int league_size , int team_size_request , int vector_length_request = 1 )
    : m_league_size( league_size )
    , m_team_size( team_size_request )
    , m_vector_length ( vector_length_request )
    {
      // Allow only power-of-two vector_length
      int check = 0;
      for(int k = 1; k < vector_length_max(); k*=2)
        if(k == vector_length_request)
          check = 1;
      if(!check)
        Impl::throw_runtime_exception( "Requested non-power-of-two vector length for TeamPolicy.");

      // Make sure league size is permissable
      if(league_size >= int(Impl::cuda_internal_maximum_grid_count()))
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
class ParallelFor< FunctorType , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Cuda > >
{
private:

  typedef Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Cuda > Policy ;

  const FunctorType  m_functor ;
  const Policy       m_policy ;  

  ParallelFor();
  ParallelFor & operator = ( const ParallelFor & );

  template< class Tag >
  inline static
  __device__
  void driver( const FunctorType & functor
             , typename Impl::enable_if< Impl::is_same< Tag , void >::value
               , typename Policy::member_type const & >::type iwork
             )
    { functor( iwork ); }

  template< class Tag >
  inline static
  __device__
  void driver( const FunctorType & functor
             , typename Impl::enable_if< ! Impl::is_same< Tag , void >::value
               , typename Policy::member_type const & >::type iwork
             )
    { functor( Tag() , iwork ); }

public:

  typedef FunctorType functor_type ;

  inline
  __device__
  void operator()(void) const
    {
      const typename Policy::member_type work_stride = blockDim.y * gridDim.x ;
      const typename Policy::member_type work_end    = m_policy.end();

      for ( typename Policy::member_type
              iwork =  m_policy.begin() + threadIdx.y + blockDim.y * blockIdx.x ;
              iwork <  work_end ;
              iwork += work_stride ) {
        ParallelFor::template driver< typename Policy::work_tag >( m_functor, iwork );
      }
    }

  ParallelFor( const FunctorType  & functor ,
               const Policy       & policy )
    : m_functor( functor )
    , m_policy(  policy )
    {
      const dim3 block(  1 , CudaTraits::WarpSize * cuda_internal_maximum_warp_count(), 1);
      const dim3 grid( std::min( ( int( policy.end() - policy.begin() ) + block.y - 1 ) / block.y
                               , cuda_internal_maximum_grid_count() )
                     , 1 , 1);

      CudaParallelLaunch< ParallelFor >( *this , grid , block , 0 );
    }
};

template< class FunctorType , class Arg0 , class Arg1 >
class ParallelFor< FunctorType , Kokkos::TeamPolicy< Arg0 , Arg1 , Kokkos::Cuda > >
{
private:

  typedef Kokkos::TeamPolicy< Arg0 , Arg1 , Kokkos::Cuda >   Policy ;

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
  size_type         m_shmem_begin ;
  size_type         m_shmem_size ;
  size_type         m_league_size ;

  template< class TagType >
  KOKKOS_FORCEINLINE_FUNCTION
  void driver( typename Impl::enable_if< Impl::is_same< TagType , void >::value ,
                 const typename Policy::member_type & >::type member ) const
    { m_functor( member ); }

  template< class TagType >
  KOKKOS_FORCEINLINE_FUNCTION
  void driver( typename Impl::enable_if< ! Impl::is_same< TagType , void >::value ,
                 const typename Policy::member_type & >::type  member ) const
    { m_functor( TagType() , member ); }

public:

  __device__ inline
  void operator()(void) const
  {
    // Iterate this block through the league
    for ( int league_rank = blockIdx.x ; league_rank < m_league_size ; league_rank += gridDim.x ) {

      ParallelFor::template driver< typename Policy::work_tag >(
        typename Policy::member_type( kokkos_impl_cuda_shared_memory<void>()
                                    , m_shmem_begin
                                    , m_shmem_size
                                    , league_rank
                                    , m_league_size ) );
    }
  }


  ParallelFor( const FunctorType  & functor 
             , const Policy       & policy 
             )
  : m_functor( functor )
  , m_shmem_begin( sizeof(double) * ( policy.team_size() + 2 ) )
  , m_shmem_size( FunctorTeamShmemSize< FunctorType >::value( functor , policy.team_size() ) )
  , m_league_size( policy.league_size() )
  {
    // Functor's reduce memory, team scan memory, and team shared memory depend upon team size.

    const int shmem_size_total = m_shmem_begin + m_shmem_size ;

    if ( CudaTraits::SharedMemoryCapacity < shmem_size_total ) {
      Kokkos::Impl::throw_runtime_exception(std::string("Kokkos::Impl::ParallelFor< Cuda > insufficient shared memory"));
    }

    const dim3 grid( int(policy.league_size()) , 1 , 1 );
    const dim3 block( policy.vector_length() , policy.team_size() , 1 );

    CudaParallelLaunch< ParallelFor >( *this, grid, block, shmem_size_total ); // copy to device and execute
  }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelReduce< FunctorType , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Cuda > >
{
private:

  typedef Kokkos::RangePolicy<Arg0,Arg1,Arg2, Kokkos::Cuda >         Policy ;
  typedef typename Policy::WorkRange                                 work_range ;
  typedef typename Policy::work_tag                                  work_tag ;
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , work_tag > ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   FunctorType , work_tag > ValueInit ;

public:

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::value_type      value_type ;
  typedef typename ValueTraits::reference_type  reference_type ;
  typedef FunctorType                           functor_type ;
  typedef Cuda::size_type                       size_type ;

  // Algorithmic constraints: blockSize is a power of two AND blockDim.y == blockDim.z == 1

  const FunctorType m_functor ;
  const Policy      m_policy ;
  size_type *       m_scratch_space ;
  size_type *       m_scratch_flags ;
  size_type *       m_unified_space ;

  // Determine block size constrained by shared memory:
  static inline
  unsigned local_block_size( const FunctorType & f )
    {
      unsigned n = CudaTraits::WarpSize * 8 ;
      while ( n && CudaTraits::SharedMemoryCapacity < cuda_single_inter_block_reduce_scan_shmem<false,FunctorType,work_tag>( f , n ) ) { n >>= 1 ; }
      return n ;
    }

  template< class Tag >
  inline static
  __device__
  void driver( const FunctorType & functor
             , typename Impl::enable_if< Impl::is_same< Tag , void >::value
               , typename Policy::member_type const & >::type iwork
             , reference_type value )
    { functor( iwork , value ); }

  template< class Tag >
  inline static
  __device__
  void driver( const FunctorType & functor
             , typename Impl::enable_if< ! Impl::is_same< Tag , void >::value
               , typename Policy::member_type const & >::type iwork
             , reference_type value )
    { functor( Tag() , iwork , value ); }

#ifndef KOKKOS_EXPERIMENTAL_CUDA_SHFL_REDUCTION
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

      const work_range range( m_policy , blockIdx.x , gridDim.x );

      for ( typename work_range::member_type iwork = range.begin() + threadIdx.y , iwork_end = range.end() ;
            iwork < iwork_end ; iwork += blockDim.y ) {
        ParallelReduce::template driver< work_tag >( m_functor , iwork , value );
      }
    }

    // Reduce with final value at blockDim.y - 1 location.
    if ( cuda_single_inter_block_reduce_scan<false,FunctorType,work_tag>(
           m_functor , blockIdx.x , gridDim.x ,
           kokkos_impl_cuda_shared_memory<size_type>() , m_scratch_space , m_scratch_flags ) ) {

      // This is the final block with the final result at the final threads' location

      size_type * const shared = kokkos_impl_cuda_shared_memory<size_type>() + ( blockDim.y - 1 ) * word_count.value ;
      size_type * const global = m_unified_space ? m_unified_space : m_scratch_space ;

      if ( threadIdx.y == 0 ) {
        Kokkos::Impl::FunctorFinal< FunctorType , work_tag >::final( m_functor , shared );
      }

      if ( CudaTraits::WarpSize < word_count.value ) { __syncthreads(); }

      for ( unsigned i = threadIdx.y ; i < word_count.value ; i += blockDim.y ) { global[i] = shared[i]; }
    }
  }
#else
  __device__ inline
   void operator()(void) const
   {

     value_type value = 0;

     // Number of blocks is bounded so that the reduction can be limited to two passes.
     // Each thread block is given an approximately equal amount of work to perform.
     // Accumulate the values for this block.
     // The accumulation ordering does not match the final pass, but is arithmatically equivalent.

     const Policy range( m_policy , blockIdx.x , gridDim.x );

     for ( typename Policy::member_type iwork = range.begin() + threadIdx.y , iwork_end = range.end() ;
           iwork < iwork_end ; iwork += blockDim.y ) {
       ParallelReduce::template driver< work_tag >( m_functor , iwork , value );
     }

     pointer_type const result = (pointer_type) (m_unified_space ? m_unified_space : m_scratch_space) ;
     int max_active_thread = range.end()-range.begin() < blockDim.y ? range.end() - range.begin():blockDim.y;
     max_active_thread = max_active_thread == 0?blockDim.y:max_active_thread;
     if(Impl::cuda_inter_block_reduction<FunctorType,Impl::JoinAdd<value_type> >
            (value,Impl::JoinAdd<value_type>(),m_scratch_space,result,m_scratch_flags,max_active_thread)) {
       const unsigned id = threadIdx.y*blockDim.x + threadIdx.x;
       if(id==0) {
         Kokkos::Impl::FunctorFinal< FunctorType , work_tag >::final( m_functor , (void*) &value );
         *result = value;
       }
     }
   }
#endif
  template< class HostViewType >
  ParallelReduce( const FunctorType  & functor 
                , const Policy       & policy 
                , const HostViewType & result
                )
  : m_functor( functor )
  , m_policy(  policy )
  , m_scratch_space( 0 )
  , m_scratch_flags( 0 )
  , m_unified_space( 0 )
  {
    const int block_size  = local_block_size( functor );
    const int block_count = std::min( int(block_size)
                                    , ( int(policy.end() - policy.begin()) + block_size - 1 ) / block_size
                                    );

    m_scratch_space = cuda_internal_scratch_space( ValueTraits::value_size( functor ) * block_count );
    m_scratch_flags = cuda_internal_scratch_flags( sizeof(size_type) );
    m_unified_space = cuda_internal_scratch_unified( ValueTraits::value_size( functor ) );

    const dim3 grid( block_count , 1 , 1 );
    const dim3 block( 1 , block_size , 1 ); // REQUIRED DIMENSIONS ( 1 , N , 1 )
#ifdef KOKKOS_EXPERIMENTAL_CUDA_SHFL_REDUCTION
    const int shmem = 0;
#else
    const int shmem = cuda_single_inter_block_reduce_scan_shmem<false,FunctorType,work_tag>( m_functor , block.y );
#endif

    CudaParallelLaunch< ParallelReduce >( *this, grid, block, shmem ); // copy to device and execute

    Cuda::fence();

    if ( result.ptr_on_device() ) {
      if ( m_unified_space ) {
        const int count = ValueTraits::value_count( m_functor );
        for ( int i = 0 ; i < count ; ++i ) { result.ptr_on_device()[i] = pointer_type(m_unified_space)[i] ; }
      }
      else {
        const int size = ValueTraits::value_size( m_functor );
        DeepCopy<HostSpace,CudaSpace>( result.ptr_on_device() , m_scratch_space , size );
      }
    }
  }
};

template< class FunctorType , class Arg0 , class Arg1 >
class ParallelReduce< FunctorType , Kokkos::TeamPolicy< Arg0 , Arg1 , Kokkos::Cuda > >
{
private:

  typedef Kokkos::TeamPolicy<Arg0,Arg1,Kokkos::Cuda>                  Policy ;
  typedef typename Policy::work_tag                                   work_tag ;
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , work_tag >  ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   FunctorType , work_tag >  ValueInit ;
  typedef typename ValueTraits::pointer_type                          pointer_type ;
  typedef typename ValueTraits::reference_type                        reference_type ;

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

  const FunctorType m_functor ;
  size_type *       m_scratch_space ;
  size_type *       m_scratch_flags ;
  size_type *       m_unified_space ;
  size_type         m_team_begin ;
  size_type         m_shmem_begin ;
  size_type         m_shmem_size ;
  size_type         m_league_size ;

  template< class TagType >
  KOKKOS_FORCEINLINE_FUNCTION
  void driver( typename Impl::enable_if< Impl::is_same< TagType , void >::value ,
                 const typename Policy::member_type & >::type  member 
             , reference_type update ) const
    { m_functor( member , update ); }

  template< class TagType >
  KOKKOS_FORCEINLINE_FUNCTION
  void driver( typename Impl::enable_if< ! Impl::is_same< TagType , void >::value ,
                 const typename Policy::member_type & >::type  member 
             , reference_type update ) const
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

      ParallelReduce::template driver< work_tag >
        ( typename Policy::member_type( kokkos_impl_cuda_shared_memory<char>() + m_team_begin
                                        , m_shmem_begin
                                        , m_shmem_size
                                        , league_rank
                                        , m_league_size )
        , value );
    }

    // Reduce with final value at blockDim.y - 1 location.
    if ( cuda_single_inter_block_reduce_scan<false,FunctorType,work_tag>(
           m_functor , blockIdx.x , gridDim.x ,
           kokkos_impl_cuda_shared_memory<size_type>() , m_scratch_space , m_scratch_flags ) ) {

      // This is the final block with the final result at the final threads' location

      size_type * const shared = kokkos_impl_cuda_shared_memory<size_type>() + ( blockDim.y - 1 ) * word_count.value ;
      size_type * const global = m_unified_space ? m_unified_space : m_scratch_space ;

      if ( threadIdx.y == 0 ) {
        Kokkos::Impl::FunctorFinal< FunctorType , work_tag >::final( m_functor , shared );
      }

      if ( CudaTraits::WarpSize < word_count.value ) { __syncthreads(); }

      for ( unsigned i = threadIdx.y ; i < word_count.value ; i += blockDim.y ) { global[i] = shared[i]; }
    }
  }


  template< class HostViewType >
  ParallelReduce( const FunctorType  & functor 
                , const Policy       & policy 
                , const HostViewType & result
                )
  : m_functor( functor )
  , m_scratch_space( 0 )
  , m_scratch_flags( 0 )
  , m_unified_space( 0 )
  , m_team_begin( cuda_single_inter_block_reduce_scan_shmem<false,FunctorType,work_tag>( functor , policy.team_size() ) )
  , m_shmem_begin( sizeof(double) * ( policy.team_size() + 2 ) )
  , m_shmem_size( FunctorTeamShmemSize< FunctorType >::value( functor , policy.team_size() ) )
  , m_league_size( policy.league_size() )
  {

    // The global parallel_reduce does not support vector_length other than 1 at the moment
    if(policy.vector_length() > 1)
      Impl::throw_runtime_exception( "Kokkos::parallel_reduce with a TeamPolicy using a vector length of greater than 1 is not currently supported for CUDA.");

    // Functor's reduce memory, team scan memory, and team shared memory depend upon team size.

    const int shmem_size_total = m_team_begin + m_shmem_begin + m_shmem_size ;
    const int not_power_of_two = 0 != ( policy.team_size() & ( policy.team_size() - 1 ) );

    if ( not_power_of_two ||  CudaTraits::SharedMemoryCapacity < shmem_size_total ) {
      Kokkos::Impl::throw_runtime_exception(std::string("Kokkos::Impl::ParallelReduce< Cuda > bad team size"));
    }

    const int block_count = std::min( policy.league_size() , policy.team_size() );

    m_scratch_space = cuda_internal_scratch_space( ValueTraits::value_size( functor ) * block_count );
    m_scratch_flags = cuda_internal_scratch_flags( sizeof(size_type) );
    m_unified_space = cuda_internal_scratch_unified( ValueTraits::value_size( functor ) );

    const dim3 grid( block_count , 1 , 1 );
    const dim3 block( 1 , policy.team_size() , 1 ); // REQUIRED DIMENSIONS ( 1 , N , 1 )

    CudaParallelLaunch< ParallelReduce >( *this, grid, block, shmem_size_total ); // copy to device and execute

    Cuda::fence();

    if ( result.ptr_on_device() ) {
      if ( m_unified_space ) {
        const int count = ValueTraits::value_count( m_functor );
        for ( int i = 0 ; i < count ; ++i ) { result.ptr_on_device()[i] = pointer_type(m_unified_space)[i] ; }
      }
      else {
        const int size = ValueTraits::value_size( m_functor );
        DeepCopy<HostSpace,CudaSpace>( result.ptr_on_device() , m_scratch_space , size );
      }
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
class ParallelScan< FunctorType , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Cuda > >
{
private:

  typedef Kokkos::RangePolicy<Arg0,Arg1,Arg2, Kokkos::Cuda >          Policy ;
  typedef typename Policy::WorkRange                                  work_range ;
  typedef typename Policy::work_tag                                   work_tag ;
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , work_tag >  ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   FunctorType , work_tag >  ValueInit ;
  typedef Kokkos::Impl::FunctorValueOps<    FunctorType , work_tag >  ValueOps ;

public:
  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;
  typedef FunctorType                           functor_type ;
  typedef Cuda::size_type                       size_type ;

  // Algorithmic constraints:
  //  (a) blockDim.y is a power of two
  //  (b) blockDim.y == blockDim.z == 1
  //  (c) gridDim.x  <= blockDim.y * blockDim.y
  //  (d) gridDim.y  == gridDim.z == 1

  // Determine block size constrained by shared memory:
  static inline
  unsigned local_block_size( const FunctorType & f )
    {
      // blockDim.y must be power of two = 128 (4 warps) or 256 (8 warps) or 512 (16 warps)
      // gridDim.x <= blockDim.y * blockDim.y
      //
      // 4 warps was 10% faster than 8 warps and 20% faster than 16 warps in unit testing

      unsigned n = CudaTraits::WarpSize * 4 ;
      while ( n && CudaTraits::SharedMemoryCapacity < cuda_single_inter_block_reduce_scan_shmem<false,FunctorType,work_tag>( f , n ) ) { n >>= 1 ; }
      return n ;
    }

  const FunctorType m_functor ;
  const Policy      m_policy ;
  size_type *       m_scratch_space ;
  size_type *       m_scratch_flags ;
        size_type   m_final ;
  
  template< class Tag >
  inline static
  __device__
  void driver( const FunctorType & functor
             , typename Impl::enable_if< Impl::is_same< Tag , void >::value
               , typename Policy::member_type const & >::type iwork
             , reference_type value 
             , const bool     final )
    { functor( iwork , value , final ); }

  template< class Tag >
  inline static
  __device__
  void driver( const FunctorType & functor
             , typename Impl::enable_if< ! Impl::is_same< Tag , void >::value
               , typename Policy::member_type const & >::type iwork
             , reference_type value
             , const bool     final )
    { functor( Tag() , iwork , value , final ); }

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

    const work_range range( m_policy , blockIdx.x , gridDim.x );

    for ( typename Policy::member_type iwork = range.begin() + threadIdx.y , iwork_end = range.end() ;
          iwork < iwork_end ; iwork += blockDim.y ) {
      ParallelScan::template driver< work_tag >
        ( m_functor , iwork , ValueOps::reference( shared_value ) , false );
    }

    // Reduce and scan, writing out scan of blocks' totals and block-groups' totals.
    // Blocks' scan values are written to 'blockIdx.x' location.
    // Block-groups' scan values are at: i = ( j * blockDim.y - 1 ) for i < gridDim.x
    cuda_single_inter_block_reduce_scan<true,FunctorType,work_tag>( m_functor , blockIdx.x , gridDim.x , kokkos_impl_cuda_shared_memory<size_type>() , m_scratch_space , m_scratch_flags );
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

    const work_range range( m_policy , blockIdx.x , gridDim.x );

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
        ParallelScan::template driver< work_tag >
          ( m_functor , iwork , ValueOps::reference( shared_prefix + word_count.value ) , false );
      }

      // Scan block values into locations shared_data[1..blockDim.y]
      cuda_intra_block_reduce_scan<true,FunctorType,work_tag>( m_functor , ValueTraits::pointer_type(shared_data+word_count.value) );

      {
        size_type * const block_total = shared_data + word_count.value * blockDim.y ;
        for ( unsigned i = threadIdx.y ; i < word_count.value ; ++i ) { shared_accum[i] = block_total[i]; }
      }

      // Call functor with exclusive scan value
      if ( iwork < range.end() ) {
        ParallelScan::template driver< work_tag >
          ( m_functor , iwork , ValueOps::reference( shared_prefix ) , true );
      }
    }
  }

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

  ParallelScan( const FunctorType  & functor ,
                const Policy       & policy )
  : m_functor( functor )
  , m_policy( policy )
  , m_scratch_space( 0 )
  , m_scratch_flags( 0 )
  , m_final( false )
  {
    enum { GridMaxComputeCapability_2x = 0x0ffff };

    const int block_size = local_block_size( functor );

    const int grid_max = ( block_size * block_size ) < GridMaxComputeCapability_2x ?
                         ( block_size * block_size ) : GridMaxComputeCapability_2x ;

    // At most 'max_grid' blocks:
    const int nwork    = policy.end() - policy.begin();
    const int max_grid = std::min( int(grid_max) , int(( nwork + block_size - 1 ) / block_size ));

    // How much work per block:
    const int work_per_block = ( nwork + max_grid - 1 ) / max_grid ;

    // How many block are really needed for this much work:
    const dim3 grid( ( nwork + work_per_block - 1 ) / work_per_block , 1 , 1 );
    const dim3 block( 1 , block_size , 1 ); // REQUIRED DIMENSIONS ( 1 , N , 1 )
    const int shmem = ValueTraits::value_size( functor ) * ( block_size + 2 );

    m_scratch_space = cuda_internal_scratch_space( ValueTraits::value_size( functor ) * grid.x );
    m_scratch_flags = cuda_internal_scratch_flags( sizeof(size_type) * 1 );

    m_final = false ;
    CudaParallelLaunch< ParallelScan >( *this, grid, block, shmem ); // copy to device and execute

    m_final = true ;
    CudaParallelLaunch< ParallelScan >( *this, grid, block, shmem ); // copy to device and execute
  }

  void wait() const { Cuda::fence(); }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifdef KOKKOS_HAVE_CXX11

namespace Kokkos {
namespace Impl {
  template<typename iType>
  struct TeamThreadLoopBoundariesStruct<iType,CudaTeamMember> {
    typedef iType index_type;
    const iType start;
    const iType end;
    const iType increment;
    const CudaTeamMember& thread;

#ifdef __CUDA_ARCH__
    __device__ inline
    TeamThreadLoopBoundariesStruct (const CudaTeamMember& thread_, const iType& count):
      start( threadIdx.y ),
      end( count ),
      increment( blockDim.y ),
      thread(thread_)
    {}
#else
    KOKKOS_INLINE_FUNCTION
    TeamThreadLoopBoundariesStruct (const CudaTeamMember& thread_, const iType& count):
      start( 0 ),
      end( count ),
      increment( 1 ),
      thread(thread_)
    {}
#endif
  };

  template<typename iType>
  struct ThreadVectorLoopBoundariesStruct<iType,CudaTeamMember> {
    typedef iType index_type;
    const iType start;
    const iType end;
    const iType increment;

#ifdef __CUDA_ARCH__
    __device__ inline
    ThreadVectorLoopBoundariesStruct (const CudaTeamMember& thread, const iType& count):
    start( threadIdx.x ),
    end( count ),
    increment( blockDim.x )
    {}
#else
    KOKKOS_INLINE_FUNCTION
    ThreadVectorLoopBoundariesStruct (const CudaTeamMember& thread_, const iType& count):
      start( 0 ),
      end( count ),
      increment( 1 )
    {}
#endif
    };

} // namespace Impl

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::TeamThreadLoopBoundariesStruct<iType,Impl::CudaTeamMember>
  TeamThreadLoop(const Impl::CudaTeamMember& thread, const iType& count) {
  return Impl::TeamThreadLoopBoundariesStruct<iType,Impl::CudaTeamMember>(thread,count);
}

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::ThreadVectorLoopBoundariesStruct<iType,Impl::CudaTeamMember >
  ThreadVectorLoop(Impl::CudaTeamMember thread, const iType count) {
  return Impl::ThreadVectorLoopBoundariesStruct<iType,Impl::CudaTeamMember >(thread,count);
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
void parallel_for(const Impl::TeamThreadLoopBoundariesStruct<iType,Impl::CudaTeamMember>& loop_boundaries, const Lambda& lambda) {
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment)
    lambda(i);
}

/** \brief  Inter-thread vector parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team and a summation of
 * val is performed and put into result. This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ValueType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::TeamThreadLoopBoundariesStruct<iType,Impl::CudaTeamMember>& loop_boundaries,
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
void parallel_reduce(const Impl::TeamThreadLoopBoundariesStruct<iType,Impl::CudaTeamMember>& loop_boundaries,
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
void parallel_for(const Impl::ThreadVectorLoopBoundariesStruct<iType,Impl::CudaTeamMember >&
    loop_boundaries, const Lambda& lambda) {

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment)
    lambda(i);
}

/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread and a summation of
 * val is performed and put into result. This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ValueType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::ThreadVectorLoopBoundariesStruct<iType,Impl::CudaTeamMember >&
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
void parallel_reduce(const Impl::ThreadVectorLoopBoundariesStruct<iType,Impl::CudaTeamMember >&
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
void parallel_scan(const Impl::ThreadVectorLoopBoundariesStruct<iType,Impl::CudaTeamMember >&
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

#endif // KOKKOS_HAVE_CXX11

namespace Kokkos {
template<int N>
struct Vectorization<Cuda,N> {
  typedef Kokkos::TeamPolicy< Cuda >         team_policy ;
  typedef typename team_policy::member_type  team_member ;
  enum {increment = N};

#ifdef __CUDA_ARCH__
  KOKKOS_FORCEINLINE_FUNCTION
  static int begin() { return threadIdx.y%N;}
#else
  KOKKOS_FORCEINLINE_FUNCTION
  static int begin() { return 0;}
#endif

  KOKKOS_FORCEINLINE_FUNCTION
  static int thread_rank(const team_member &dev) {
    return dev.team_rank()/increment;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  static int team_rank(const team_member &dev) {
    return dev.team_rank()/increment;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  static int team_size(const team_member &dev) {
    return dev.team_size()/increment;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  static int global_thread_rank(const team_member &dev) {
    return (dev.league_rank()*dev.team_size()+dev.team_rank())/increment;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  static bool is_lane_0(const team_member &dev) {
    return (dev.team_rank()%increment)==0;
  }

  template<class Scalar>
  KOKKOS_INLINE_FUNCTION
  static Scalar reduce(const Scalar& val) {
    #ifdef __CUDA_ARCH__
    __shared__ Scalar result[256];
    Scalar myresult;
    for(int k=0;k<blockDim.y;k+=256) {
      const int tid = threadIdx.y - k;
      if(tid > 0 && tid<256) {
        result[tid] = val;
        if ( (N > 1) && (tid%2==0) )
          result[tid] += result[tid+1];
        if ( (N > 2) && (tid%4==0) )
          result[tid] += result[tid+2];
        if ( (N > 4) && (tid%8==0) )
          result[tid] += result[tid+4];
        if ( (N > 8) && (tid%16==0) )
          result[tid] += result[tid+8];
        if ( (N > 16) && (tid%32==0) )
          result[tid] += result[tid+16];
        myresult = result[tid];
      }
      if(blockDim.y>256)
        __syncthreads();
    }
    return myresult;
    #else
    return val;
    #endif
  }

#ifdef __CUDA_ARCH__
  #if (__CUDA_ARCH__ >= 300)
  KOKKOS_INLINE_FUNCTION
  static int reduce(const int& val) {
    int result = val;
    if (N > 1)
      result += shfl_down(result, 1,N);
    if (N > 2)
      result += shfl_down(result, 2,N);
    if (N > 4)
      result += shfl_down(result, 4,N);
    if (N > 8)
      result += shfl_down(result, 8,N);
    if (N > 16)
      result += shfl_down(result, 16,N);
    return result;
  }

  KOKKOS_INLINE_FUNCTION
  static unsigned int reduce(const unsigned int& val) {
    unsigned int result = val;
    if (N > 1)
      result += shfl_down(result, 1,N);
    if (N > 2)
      result += shfl_down(result, 2,N);
    if (N > 4)
      result += shfl_down(result, 4,N);
    if (N > 8)
      result += shfl_down(result, 8,N);
    if (N > 16)
      result += shfl_down(result, 16,N);
    return result;
  }

  KOKKOS_INLINE_FUNCTION
  static long int reduce(const long int& val) {
    long int result = val;
    if (N > 1)
      result += shfl_down(result, 1,N);
    if (N > 2)
      result += shfl_down(result, 2,N);
    if (N > 4)
      result += shfl_down(result, 4,N);
    if (N > 8)
      result += shfl_down(result, 8,N);
    if (N > 16)
      result += shfl_down(result, 16,N);
    return result;
  }

  KOKKOS_INLINE_FUNCTION
  static unsigned long int reduce(const unsigned long int& val) {
    unsigned long int result = val;
    if (N > 1)
      result += shfl_down(result, 1,N);
    if (N > 2)
      result += shfl_down(result, 2,N);
    if (N > 4)
      result += shfl_down(result, 4,N);
    if (N > 8)
      result += shfl_down(result, 8,N);
    if (N > 16)
      result += shfl_down(result, 16,N);
    return result;
  }

  KOKKOS_INLINE_FUNCTION
  static float reduce(const float& val) {
    float result = val;
    if (N > 1)
      result += shfl_down(result, 1,N);
    if (N > 2)
      result += shfl_down(result, 2,N);
    if (N > 4)
      result += shfl_down(result, 4,N);
    if (N > 8)
      result += shfl_down(result, 8,N);
    if (N > 16)
      result += shfl_down(result, 16,N);
    return result;
  }

  KOKKOS_INLINE_FUNCTION
  static double reduce(const double& val) {
    double result = val;
    if (N > 1)
      result += shfl_down(result, 1,N);
    if (N > 2)
      result += shfl_down(result, 2,N);
    if (N > 4)
      result += shfl_down(result, 4,N);
    if (N > 8)
      result += shfl_down(result, 8,N);
    if (N > 16)
      result += shfl_down(result, 16,N);
    return result;
  }
  #endif
#endif

};
}

#endif /* defined( __CUDACC__ ) */

#endif /* #ifndef KOKKOS_CUDA_PARALLEL_HPP */

