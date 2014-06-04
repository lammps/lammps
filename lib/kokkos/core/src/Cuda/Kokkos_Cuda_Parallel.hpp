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

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class WorkSpec >
class ParallelFor< FunctorType , WorkSpec /* size_t */ , Cuda > {
private:

  const FunctorType     m_functor ;
  const Cuda::size_type m_work ;  

  ParallelFor();
  ParallelFor & operator = ( const ParallelFor & );

public:

  inline
  __device__
  void operator()(void) const
  {
    const Cuda::size_type work_stride = blockDim.x * gridDim.x ;

    for ( Cuda::size_type
            iwork = threadIdx.x + blockDim.x * blockIdx.x ;
            iwork < m_work ;
            iwork += work_stride ) {
      m_functor( iwork );
    }
  }

  ParallelFor( const FunctorType  & functor ,
               const size_t         work )
    : m_functor( functor )
    , m_work(    work )
    {
      const dim3 block( CudaTraits::WarpSize * cuda_internal_maximum_warp_count(), 1, 1);
      const dim3 grid( std::min( ( m_work + block.x - 1 ) / block.x , cuda_internal_maximum_grid_count() ) , 1 , 1 );

      CudaParallelLaunch< ParallelFor >( *this , grid , block , 0 );
    }
};

template< class FunctorType >
class ParallelFor< FunctorType , ParallelWorkRequest , Cuda > {
private:

  const FunctorType          m_functor ;
  const ParallelWorkRequest  m_work ;
  const int                  m_shmem ;

  ParallelFor();
  ParallelFor & operator = ( const ParallelFor & );

public:

  inline
  __device__
  void operator()(void) const
  {
    CudaExec exec( 0 , m_shmem );
    m_functor( Cuda( exec ) );
  }

  ParallelFor( const FunctorType         & functor ,
               const ParallelWorkRequest &  work )
    : m_functor( functor )
    , m_work( std::min( work.league_size , size_t(cuda_internal_maximum_grid_count()) ) ,
              std::min( work.team_size ,   size_t(CudaTraits::WarpSize * cuda_internal_maximum_warp_count()) ) )
    , m_shmem( FunctorShmemSize< FunctorType >::value( functor ) )
    {
      const dim3 grid(  m_work.league_size , 1 , 1 );
      const dim3 block( m_work.team_size , 1, 1 );

      CudaParallelLaunch< ParallelFor >( *this , grid , block , m_shmem );
    }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType >
class ParallelFor< FunctorType , CudaWorkConfig , Cuda > {
public:

  const FunctorType m_work_functor ;

  inline
  __device__
  void operator()(void) const
  {
    Cuda::size_type iwork = threadIdx.x + blockDim.x * (
                            threadIdx.y + blockDim.y * (
                            threadIdx.z + blockDim.z * (
                            blockIdx.x + gridDim.x * (
                            blockIdx.y + gridDim.y * (
                            blockIdx.z )))));

    m_work_functor( iwork );
  }

  ParallelFor( const FunctorType    & functor ,
               const CudaWorkConfig & work_config )
  : m_work_functor( functor )
  {
    const dim3 grid( work_config.grid[0] ,
                     work_config.grid[1] ,
                     work_config.grid[2] );

    const dim3 block( work_config.block[0] ,
                      work_config.block[1] ,
                      work_config.block[2] );

    CudaParallelLaunch< ParallelFor >( *this , grid , block , work_config.shared );
  }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class WorkSpec >
class ParallelReduce< FunctorType , WorkSpec , Cuda >
{
public:
  typedef ReduceAdapter< FunctorType >        Reduce ;
  typedef typename Reduce::pointer_type       pointer_type ;
  typedef typename Reduce::reference_type     reference_type ;
  typedef Cuda::size_type                     size_type ;

  // Algorithmic constraints:
  //  (a) blockSize is a power of two
  //  (b) blockDim.x == BlockSize == 1 << BlockSizeShift
  //  (c) blockDim.y == blockDim.z == 1

  enum { WarpCount      = 8 };
  enum { BlockSize      = CudaTraits::WarpSize << power_of_two< WarpCount >::value };
  enum { BlockSizeShift = power_of_two< BlockSize >::value };
  enum { BlockSizeMask  = BlockSize - 1 };

  enum { GridMaxComputeCapability_2x = 0x0ffff };
  enum { GridMax = BlockSize };

  const FunctorType m_functor ;
  size_type *       m_scratch_space ;
  size_type *       m_scratch_flags ;
  size_type *       m_unified_space ;
  pointer_type      m_host_pointer ;
  size_type         m_work ;
  size_type         m_work_per_block ;
  size_type         m_local_block_count ;
  size_type         m_global_block_begin ;
  size_type         m_global_block_count ;


  __device__ inline
  void operator()(void) const
  {
    extern __shared__ size_type shared_data[];

    const integral_nonzero_constant< size_type , Reduce::StaticValueSize / sizeof(size_type) >
      word_count( Reduce::value_size( m_functor ) / sizeof(size_type) );

    {
      reference_type value = Reduce::reference( shared_data + threadIdx.x * word_count.value );

      m_functor.init( value );

      // Number of blocks is bounded so that the reduction can be limited to two passes.
      // Each thread block is given an approximately equal amount of work to perform.
      // Accumulate the values for this block.
      // The accumulation ordering does not match the final pass, but is arithmatically equivalent.

      const size_type iwork_beg = blockIdx.x * m_work_per_block ;
      const size_type iwork_end = iwork_beg + m_work_per_block < m_work
                                ? iwork_beg + m_work_per_block : m_work ;

      for ( size_type iwork = threadIdx.x + iwork_beg ; iwork < iwork_end ; iwork += BlockSize ) {
        m_functor( iwork , value );
      }
    }

    // Reduce with final value at BlockSize - 1 location.
    if ( cuda_single_inter_block_reduce_scan<false,BlockSize>(
           m_functor , m_global_block_begin + blockIdx.x , m_global_block_count ,
           shared_data , m_scratch_space , m_scratch_flags ) ) {

      // This is the final block with the final result at the final threads' location

      size_type * const shared = shared_data + BlockSizeMask * word_count.value ;
      size_type * const global = m_unified_space ? m_unified_space : m_scratch_space ;

      if ( threadIdx.x == 0 ) { Reduce::final( m_functor , shared ); }

      if ( CudaTraits::WarpSize < word_count.value ) { __syncthreads(); }

      for ( unsigned i = threadIdx.x ; i < word_count.value ; i += BlockSize ) { global[i] = shared[i]; }
    }
  }

  ParallelReduce( const FunctorType  & functor ,
                  const size_t         nwork ,
                  const pointer_type   result = 0 ,
                  const bool execute_immediately = true )
  : m_functor( functor )
  , m_scratch_space( 0 )
  , m_scratch_flags( 0 )
  , m_unified_space( 0 )
  , m_host_pointer( result )
  , m_work( nwork )
  , m_work_per_block( 0 )
  , m_local_block_count( 0 )
  , m_global_block_begin( 0 )
  , m_global_block_count( 0 )
  {
    // At most 'max_grid' blocks:
    const int max_grid = std::min( int(GridMax) , int(( nwork + BlockSizeMask ) / BlockSize ));

    // How much work per block:
    m_work_per_block = ( nwork + max_grid - 1 ) / max_grid ;

    // How many block are really needed for this much work:
    m_local_block_count  = ( nwork + m_work_per_block - 1 ) / m_work_per_block ;
    m_global_block_count = m_local_block_count ;

    m_scratch_space = cuda_internal_scratch_space( Reduce::value_size( functor ) * m_local_block_count );
    m_scratch_flags = cuda_internal_scratch_flags( sizeof(size_type) );
    m_unified_space = cuda_internal_scratch_unified( Reduce::value_size( functor ) );

    if ( execute_immediately ) { execute(); }
  }

  inline
  void execute() const
  {
    const dim3 grid( m_local_block_count , 1 , 1 );
    const dim3 block( BlockSize , 1 , 1 );
    const int shmem = cuda_single_inter_block_reduce_scan_shmem<false,BlockSize>( m_functor );

    CudaParallelLaunch< ParallelReduce >( *this, grid, block, shmem ); // copy to device and execute
  }

  void wait() const
  {
    Cuda::fence();

    if ( m_host_pointer ) {
      if ( m_unified_space ) {
        const int count = Reduce::value_count( m_functor );
        for ( int i = 0 ; i < count ; ++i ) { m_host_pointer[i] = pointer_type(m_unified_space)[i] ; }
      }
      else {
        const int size = Reduce::value_size( m_functor );
        DeepCopy<HostSpace,CudaSpace>( m_host_pointer , m_scratch_space , size );
      }
    }
  }
};


template< class FunctorType >
class ParallelReduce< FunctorType , ParallelWorkRequest , Cuda >
{
public:
  typedef ReduceAdapter< FunctorType >        Reduce ;
  typedef typename Reduce::pointer_type       pointer_type ;
  typedef typename Reduce::reference_type     reference_type ;
  typedef Cuda::size_type                     size_type ;

  // Algorithmic constraints:
  //  (a) blockSize is a power of two
  //  (b) blockDim.x == BlockSize == 1 << BlockSizeShift
  //  (b) blockDim.y == blockDim.z == 1

  enum { WarpCount      = 8 };
  enum { BlockSize      = CudaTraits::WarpSize << power_of_two< WarpCount >::value };
  enum { BlockSizeShift = power_of_two< BlockSize >::value };
  enum { BlockSizeMask  = BlockSize - 1 };

  enum { GridMaxComputeCapability_2x = 0x0ffff };
  enum { GridMax = BlockSize };

  const FunctorType m_functor ;
  size_type *       m_scratch_space ;
  size_type *       m_scratch_flags ;
  size_type *       m_unified_space ;
  pointer_type      m_host_pointer ;
  size_type         m_shmem_begin ;
  size_type         m_shmem_end ;
  size_type         m_local_block_count ;
  size_type         m_global_block_begin ;
  size_type         m_global_block_count ;

  __device__ inline
  void operator()(void) const
  {
    extern __shared__ size_type shared_data[];

    const integral_nonzero_constant< size_type , Reduce::StaticValueSize / sizeof(size_type) >
      word_count( Reduce::value_size( m_functor ) / sizeof(size_type) );

    {
      reference_type value = Reduce::reference( shared_data + threadIdx.x * word_count.value );

      m_functor.init( value );

      CudaExec exec( m_shmem_begin , m_shmem_end );

      m_functor( Cuda( exec ) , value );
    }

    // Reduce with final value at BlockSize - 1 location.
    if ( cuda_single_inter_block_reduce_scan<false,BlockSize>(
           m_functor , m_global_block_begin + blockIdx.x , m_global_block_count ,
           shared_data , m_scratch_space , m_scratch_flags ) ) {

      // This is the final block with the final result at the final threads' location

      size_type * const shared = shared_data + BlockSizeMask * word_count.value ;
      size_type * const global = m_unified_space ? m_unified_space : m_scratch_space ;

      if ( threadIdx.x == 0 ) { Reduce::final( m_functor , shared ); }

      if ( CudaTraits::WarpSize < word_count.value ) { __syncthreads(); }

      for ( unsigned i = threadIdx.x ; i < word_count.value ; i += BlockSize ) { global[i] = shared[i]; }
    }
  }


  ParallelReduce( const FunctorType         & functor ,
                  const ParallelWorkRequest & work ,
                  const pointer_type          result = 0 ,
                  const bool execute_immediately = true )
  : m_functor( functor )
  , m_scratch_space( 0 )
  , m_scratch_flags( 0 )
  , m_unified_space( 0 )
  , m_host_pointer( result )
  , m_shmem_begin( cuda_single_inter_block_reduce_scan_shmem<false,BlockSize>( functor ) )
  , m_shmem_end(   cuda_single_inter_block_reduce_scan_shmem<false,BlockSize>( functor )
                   + FunctorShmemSize< FunctorType >::value( functor ) )
  , m_local_block_count( 0 )
  , m_global_block_begin( 0 )
  , m_global_block_count( 0 )
  {
    m_local_block_count  = std::min( int(GridMax) , int(work.league_size) );
    m_global_block_count = std::min( int(GridMax) , int(work.league_size) );
    m_scratch_space = cuda_internal_scratch_space( Reduce::value_size( functor ) * m_local_block_count );
    m_scratch_flags = cuda_internal_scratch_flags( sizeof(size_type) );
    m_unified_space = cuda_internal_scratch_unified( Reduce::value_size( functor ) );

    if ( execute_immediately ) { execute(); }
  }

  inline
  void execute() const
  {
    const dim3 grid( m_local_block_count , 1 , 1 );
    const dim3 block( BlockSize , 1 , 1 );

    CudaParallelLaunch< ParallelReduce >( *this, grid, block, m_shmem_end ); // copy to device and execute
  }

  void wait() const
  {
    Cuda::fence();

    if ( m_host_pointer ) {
      if ( m_unified_space ) {
        const int count = Reduce::value_count( m_functor );
        for ( int i = 0 ; i < count ; ++i ) { m_host_pointer[i] = pointer_type(m_unified_space)[i] ; }
      }
      else {
        const int size = Reduce::value_size( m_functor );
        DeepCopy<HostSpace,CudaSpace>( m_host_pointer , m_scratch_space , size );
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

template< class Functor >
class MultiFunctorParallelReduceMember ;

template<>
class MultiFunctorParallelReduceMember<void>
{
private:

  MultiFunctorParallelReduceMember( const MultiFunctorParallelReduceMember & );
  MultiFunctorParallelReduceMember & operator = ( const MultiFunctorParallelReduceMember & );

protected:

  MultiFunctorParallelReduceMember() {}

public:

  virtual unsigned block_count() const = 0 ;

  virtual ~MultiFunctorParallelReduceMember() {}

  virtual void execute( void * const host_pointer ,
                        const unsigned global_block_begin ,
                        const unsigned global_block_count ) = 0 ;

  virtual void wait() const = 0 ;
};

template< class Functor >
class MultiFunctorParallelReduceMember : public MultiFunctorParallelReduceMember<void> {
public:
  ParallelReduce< Functor , size_t , Cuda >  m_functor ;

  MultiFunctorParallelReduceMember( const Functor & f , size_t nwork )
    : MultiFunctorParallelReduceMember<void>()
    , m_functor( f , nwork , 0 , false )
    {}

  virtual unsigned block_count() const { return m_functor.m_local_block_count ; }

  virtual void execute( void * const host_pointer ,
                        const unsigned global_block_begin ,
                        const unsigned global_block_count )
  {
    m_functor.m_host_pointer = typename ReduceAdapter< Functor >::pointer_type(host_pointer);
    m_functor.m_global_block_begin = global_block_begin ;
    m_functor.m_global_block_count = global_block_count ;
    m_functor.execute();
  }

  virtual void wait() const { m_functor.wait(); }
};

} // namespace Impl
} // namespace Kokkos

namespace Kokkos {

template<>
class MultiFunctorParallelReduce< Cuda >
{
private:

  typedef std::vector< Impl::MultiFunctorParallelReduceMember<void> * > MemberVector ;

  MemberVector m_functors ;

public:

  MultiFunctorParallelReduce()
    : m_functors()
    {}

  ~MultiFunctorParallelReduce()
  {
    while ( ! m_functors.empty() ) {
      delete m_functors.back();
      m_functors.pop_back();
    }
  }

  template< class FunctorType >
  void push_back( const size_t work_count , const FunctorType & f )
  {
    m_functors.push_back( new Impl::MultiFunctorParallelReduceMember<FunctorType>( f , work_count ) );
  }

  void execute( void * host_pointer )
  {
    typename MemberVector::iterator m ;

    Cuda::size_type block_count = 0 ;

    for ( m = m_functors.begin() ; m != m_functors.end() ; ++m ) {
      block_count += (*m)->block_count();
    }

    Cuda::size_type block_offset = 0 ;

    for ( m = m_functors.begin() ; m != m_functors.end() ; ++m ) {
      (*m)->execute( host_pointer , block_offset , block_count );
      block_offset += (*m)->block_count();
    }
  }

  void wait() const
  {
    if ( ! m_functors.empty() ) { (m_functors.back())->wait(); }
  }
};

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class WorkSpec >
class ParallelScan< FunctorType , WorkSpec , Cuda >
{
public:
  typedef ReduceAdapter< FunctorType >        Reduce ;
  typedef typename Reduce::pointer_type       pointer_type ;
  typedef typename Reduce::reference_type     reference_type ;
  typedef Cuda::size_type                     size_type ;

  // Algorithmic constraints:
  //  (a) blockSize is a power of two
  //  (b) blockDim.x == BlockSize == 1 << BlockSizeShift
  //  (b) blockDim.y == blockDim.z == 1
  //  (c) gridDim.x  <= blockDim.x * blockDim.x
  //  (d) gridDim.y  == gridDim.z == 1

  // blockDim.x must be power of two = 128 (4 warps) or 256 (8 warps) or 512 (16 warps)
  // gridDim.x <= blockDim.x * blockDim.x
  //
  // 4 warps was 10% faster than 8 warps and 20% faster than 16 warps in unit testing

  enum { WarpCount      = 4 };
  enum { BlockSize      = CudaTraits::WarpSize << power_of_two< WarpCount >::value };
  enum { BlockSizeShift = power_of_two< BlockSize >::value };
  enum { BlockSizeMask  = BlockSize - 1 };

  enum { GridMaxComputeCapability_2x = 0x0ffff };
  enum { GridMax = ( BlockSize * BlockSize ) < GridMaxComputeCapability_2x
                 ? ( BlockSize * BlockSize ) : GridMaxComputeCapability_2x };

  const FunctorType m_functor ;
  size_type *       m_scratch_space ;
  size_type *       m_scratch_flags ;
  const size_type   m_work ;
        size_type   m_work_per_block ;
        size_type   m_final ;
  
  //----------------------------------------

  __device__ inline
  void initial(void) const
  {
    extern __shared__ size_type shared_data[];

    const integral_nonzero_constant< size_type , Reduce::StaticValueSize / sizeof(size_type) >
      word_count( Reduce::value_size( m_functor ) / sizeof(size_type) );

    size_type * const shared_value = shared_data + word_count.value * threadIdx.x ;

    m_functor.init( Reduce::reference( shared_value ) );

    // Number of blocks is bounded so that the reduction can be limited to two passes.
    // Each thread block is given an approximately equal amount of work to perform.
    // Accumulate the values for this block.
    // The accumulation ordering does not match the final pass, but is arithmatically equivalent.

    const size_type iwork_beg = blockIdx.x * m_work_per_block ;
    const size_type iwork_end = iwork_beg + m_work_per_block < m_work 
                              ? iwork_beg + m_work_per_block : m_work ;

    for ( size_type iwork = threadIdx.x + iwork_beg ; iwork < iwork_end ; iwork += BlockSize ) {
      m_functor( iwork , Reduce::reference( shared_value ) , false );
    }

    // Reduce and scan, writing out scan of blocks' totals and block-groups' totals.
    // Blocks' scan values are written to 'blockIdx.x' location.
    // Block-groups' scan values are at: i = ( j * BlockSize - 1 ) for i < gridDim.x
    cuda_single_inter_block_reduce_scan<true,BlockSize>( m_functor , blockIdx.x , gridDim.x , shared_data , m_scratch_space , m_scratch_flags );
  }

  //----------------------------------------

  __device__ inline
  void final(void) const
  {
    extern __shared__ size_type shared_data[];

    const integral_nonzero_constant< size_type , Reduce::StaticValueSize / sizeof(size_type) >
      word_count( Reduce::value_size( m_functor ) / sizeof(size_type) );

    // Use shared memory as an exclusive scan: { 0 , value[0] , value[1] , value[2] , ... }
    size_type * const shared_prefix = shared_data + word_count.value * threadIdx.x ;
    size_type * const shared_accum  = shared_data + word_count.value * ( BlockSize + 1 );

    // Starting value for this thread block is the previous block's total.
    if ( blockIdx.x ) {
      size_type * const block_total = m_scratch_space + word_count.value * ( blockIdx.x - 1 );
      for ( unsigned i = threadIdx.x ; i < word_count.value ; ++i ) { shared_accum[i] = block_total[i] ; }
    }
    else if ( 0 == threadIdx.x ) {
      m_functor.init( Reduce::reference( shared_accum ) );
    }

          unsigned iwork_beg = blockIdx.x * m_work_per_block ;
    const unsigned iwork_end = iwork_beg + m_work_per_block ;

    for ( ; iwork_beg < iwork_end ; iwork_beg += BlockSize ) {

      const unsigned iwork = threadIdx.x + iwork_beg ;

      __syncthreads(); // Don't overwrite previous iteration values until they are used

      m_functor.init( Reduce::reference( shared_prefix + word_count.value ) );

      // Copy previous block's accumulation total into thread[0] prefix and inclusive scan value of this block
      for ( unsigned i = threadIdx.x ; i < word_count.value ; ++i ) {
        shared_data[i + word_count.value] = shared_data[i] = shared_accum[i] ;
      }

      if ( CudaTraits::WarpSize < word_count.value ) { __syncthreads(); } // Protect against large scan values.

      // Call functor to accumulate inclusive scan value for this work item
      if ( iwork < m_work ) { m_functor( iwork , Reduce::reference( shared_prefix + word_count.value ) , false ); }

      // Scan block values into locations shared_data[1..BlockSize]
      cuda_intra_block_reduce_scan<true>( m_functor , Reduce::pointer_type(shared_data+word_count.value) );

      {
        size_type * const block_total = shared_data + word_count.value * blockDim.x ;
        for ( unsigned i = threadIdx.x ; i < word_count.value ; ++i ) { shared_accum[i] = block_total[i]; }
      }

      // Call functor with exclusive scan value
      if ( iwork < m_work ) { m_functor( iwork , Reduce::reference( shared_prefix ) , true ); }
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
                const size_t         nwork )
  : m_functor( functor )
  , m_scratch_space( 0 )
  , m_scratch_flags( 0 )
  , m_work( nwork )
  , m_work_per_block( 0 )
  , m_final( false )
  {
    // At most 'max_grid' blocks:
    const int max_grid = std::min( int(GridMax) , int(( nwork + BlockSizeMask ) / BlockSize ));

    // How much work per block:
    m_work_per_block = ( nwork + max_grid - 1 ) / max_grid ;

    // How many block are really needed for this much work:
    const dim3 grid( ( nwork + m_work_per_block - 1 ) / m_work_per_block , 1 , 1 );
    const dim3 block( BlockSize , 1 , 1 );
    const int shmem = Reduce::value_size( functor ) * ( BlockSize + 2 );

    m_scratch_space = cuda_internal_scratch_space( Reduce::value_size( functor ) * grid.x );
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

#if defined( __CUDA_ARCH__ )

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

} // namespace Impl
} // namespace Kokkos

namespace Kokkos {

template< typename TypeLocal , typename TypeGlobal >
__device__ inline TypeGlobal Cuda::team_scan( const TypeLocal & value , TypeGlobal * const global_accum )
{
  enum { BlockSizeMax = 512 };

  __shared__ TypeGlobal base_data[ BlockSizeMax + 1 ];

  __syncthreads(); // Don't write in to shared data until all threads have entered this function

  if ( 0 == threadIdx.x ) { base_data[0] = 0 ; }

  base_data[ threadIdx.x + 1 ] = value ;

  Impl::cuda_intra_block_reduce_scan<true>( Impl::CudaJoinFunctor<TypeGlobal>() , base_data + 1 );

  if ( global_accum ) {
    if ( blockDim.x == threadIdx.x + 1 ) {
      base_data[ blockDim.x ] = atomic_fetch_add( global_accum , base_data[ blockDim.x ] );
    }
    __syncthreads(); // Wait for atomic
    base_data[ threadIdx.x ] += base_data[ blockDim.x ] ;
  }

  return base_data[ threadIdx.x ];
}

template< typename Type >
__device__ inline Type Cuda::team_scan( const Type & value )
{ return team_scan( value , (Type*) 0 ); }

} // namespace Kokkos

#else /* ! defined( __CUDA_ARCH__ ) */

namespace Kokkos {

template< typename Type > inline Type Cuda::team_scan( const Type & ) { return 0 ; }

template< typename TypeLocal , typename TypeGlobal >
inline TypeGlobal Cuda::team_scan( const TypeLocal & , TypeGlobal * const ) { return 0 ; }

} // namespace Kokkos

#endif /* ! defined( __CUDA_ARCH__ ) */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* defined( __CUDACC__ ) */

#endif /* #ifndef KOKKOS_CUDA_PARALLEL_HPP */

