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

#ifndef KOKKOS_CUDA_REDUCESCAN_HPP
#define KOKKOS_CUDA_REDUCESCAN_HPP

#if defined( __CUDACC__ )

#include <utility>

#include <Kokkos_Parallel.hpp>
#include <impl/Kokkos_Error.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
// See section B.17 of Cuda C Programming Guide Version 3.2
// for discussion of
//   __launch_bounds__(maxThreadsPerBlock,minBlocksPerMultiprocessor)
// function qualifier which could be used to improve performance.
//----------------------------------------------------------------------------
// Maximize shared memory and minimize L1 cache:
//   cudaFuncSetCacheConfig(MyKernel, cudaFuncCachePreferShared );
// For 2.0 capability: 48 KB shared and 16 KB L1
//----------------------------------------------------------------------------
// Must have consistent '__shared__' statement across all device kernels.
// Since there may be more than one kernel in a file then have to make this
// a simple array of words.
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/*
 *  Algorithmic constraints:
 *   (a) blockDim.x is a power of two
 *   (b) blockDim.x <= 512
 *   (c) blockDim.y == blockDim.z == 1
 */
template< bool DoScan , class FunctorType >
__device__
void cuda_intra_block_reduce_scan( const FunctorType & functor ,
                                   const typename ReduceAdapter< FunctorType >::pointer_type base_data )
{
  typedef ReduceAdapter< FunctorType >   Reduce ;
  typedef typename Reduce::pointer_type  pointer_type ;

  const unsigned value_count   = Reduce::value_count( functor );
  const unsigned BlockSizeMask = blockDim.x - 1 ;

  // Must have power of two thread count

  if ( BlockSizeMask & blockDim.x ) { cuda_abort("Cuda::cuda_intra_block_scan requires power-of-two blockDim"); }

#define BLOCK_REDUCE_STEP( R , TD , S )  \
  if ( ! ( R & ((1<<(S+1))-1) ) ) \
    { functor.join( Reduce::reference(TD) , Reduce::reference(TD - (value_count<<S))); }

#define BLOCK_SCAN_STEP( TD , N , S )  \
  if ( N == (1<<S) ) \
    { functor.join( Reduce::reference(TD) , Reduce::reference(TD - (value_count<<S))); }

  const unsigned     rtid_intra = threadIdx.x ^ BlockSizeMask ;
  const pointer_type tdata_intra = base_data + value_count * threadIdx.x ;

  { // Intra-warp reduction:
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,0)
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,1)
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,2)
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,3)
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,4)
  }

  __syncthreads(); // Wait for all warps to reduce

  { // Inter-warp reduce-scan by a single warp to avoid extra synchronizations
    const unsigned rtid_inter = ( threadIdx.x ^ BlockSizeMask ) << CudaTraits::WarpIndexShift ;

    if ( rtid_inter < blockDim.x ) {

      const pointer_type tdata_inter = base_data + value_count * ( rtid_inter ^ BlockSizeMask );

      if ( (1<<5) < BlockSizeMask ) {                        BLOCK_REDUCE_STEP(rtid_inter,tdata_inter,5) }
      if ( (1<<6) < BlockSizeMask ) { __threadfence_block(); BLOCK_REDUCE_STEP(rtid_inter,tdata_inter,6) }
      if ( (1<<7) < BlockSizeMask ) { __threadfence_block(); BLOCK_REDUCE_STEP(rtid_inter,tdata_inter,7) }
      if ( (1<<8) < BlockSizeMask ) { __threadfence_block(); BLOCK_REDUCE_STEP(rtid_inter,tdata_inter,8) }

      if ( DoScan ) {

        int n = ( rtid_inter &  32 ) ?  32 : (
                ( rtid_inter &  64 ) ?  64 : (
                ( rtid_inter & 128 ) ? 128 : (
                ( rtid_inter & 256 ) ? 256 : 0 )));

        if ( ! ( rtid_inter + n < blockDim.x ) ) n = 0 ;

        BLOCK_SCAN_STEP(tdata_inter,n,8)
        BLOCK_SCAN_STEP(tdata_inter,n,7)
        BLOCK_SCAN_STEP(tdata_inter,n,6)
        BLOCK_SCAN_STEP(tdata_inter,n,5)
      }
    }
  }

  __syncthreads(); // Wait for inter-warp reduce-scan to complete

  if ( DoScan ) {
    int n = ( rtid_intra &  1 ) ?  1 : (
            ( rtid_intra &  2 ) ?  2 : (
            ( rtid_intra &  4 ) ?  4 : (
            ( rtid_intra &  8 ) ?  8 : (
            ( rtid_intra & 16 ) ? 16 : 0 ))));

    if ( ! ( rtid_intra + n < blockDim.x ) ) n = 0 ;

    BLOCK_SCAN_STEP(tdata_intra,n,4) __threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,3) __threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,2) __threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,1) __threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,0)
  }

#undef BLOCK_SCAN_STEP
#undef BLOCK_REDUCE_STEP
}

//----------------------------------------------------------------------------
/**\brief  Input value-per-thread starting at 'shared_data'.
 *         Reduction value at last thread's location.
 *
 *  If 'DoScan' then write blocks' scan values and block-groups' scan values.
 *
 *  Global reduce result is in the last threads' 'shared_data' location.
 */
template< bool DoScan , unsigned ArgBlockSize , class FunctorType >
__device__
bool cuda_single_inter_block_reduce_scan( const FunctorType     & functor ,
                                          const Cuda::size_type   block_id ,
                                          const Cuda::size_type   block_count ,
                                          Cuda::size_type * const shared_data ,
                                          Cuda::size_type * const global_data ,
                                          Cuda::size_type * const global_flags )
{
  typedef Cuda::size_type                  size_type ;
  typedef ReduceAdapter< FunctorType >     Reduce ;
  typedef typename Reduce::pointer_type    pointer_type ;
  typedef typename Reduce::reference_type  reference_type ;

  enum { BlockSize      = ArgBlockSize };
  enum { BlockSizeMask  = BlockSize - 1 };
  enum { BlockSizeShift = power_of_two< BlockSize >::value };

  const integral_nonzero_constant< size_type , Reduce::StaticValueSize / sizeof(size_type) >
    word_count( Reduce::value_size( functor ) / sizeof(size_type) );

  // Must have power of two thread count
  if ( BlockSize != blockDim.x ) { cuda_abort("Cuda::cuda_inter_block_scan wrong blockDim.x"); }

  // Reduce the accumulation for the entire block.
  cuda_intra_block_reduce_scan<false>( functor , pointer_type(shared_data) );

  {
    // Write accumulation total to global scratch space.
    // Accumulation total is the last thread's data.
    size_type * const shared = shared_data + word_count.value * BlockSizeMask ;
    size_type * const global = global_data + word_count.value * block_id ;

    for ( size_type i = threadIdx.x ; i < word_count.value ; i += BlockSize ) { global[i] = shared[i] ; }
  }

  // Contributing blocks note that their contribution has been completed via an atomic-increment flag
  // If this block is not the last block to contribute to this group then the block is done.
  const bool is_last_block =
    ! __syncthreads_or( threadIdx.x ? 0 : ( 1 + atomicInc( global_flags , block_count - 1 ) < block_count ) );

  if ( is_last_block ) {

    const size_type b = ( long(block_count) * long(threadIdx.x) ) >> BlockSizeShift ;
    const size_type e = ( long(block_count) * long( threadIdx.x + 1 ) ) >> BlockSizeShift ;

    {
      reference_type shared_value = Reduce::reference( shared_data + word_count.value * threadIdx.x );

      functor.init( shared_value );

      for ( size_type i = b ; i < e ; ++i ) {
        functor.join( shared_value , Reduce::reference( global_data + word_count.value * i ) );
      }
    }

    cuda_intra_block_reduce_scan<DoScan>( functor , pointer_type(shared_data) );

    if ( DoScan ) {

      size_type * const shared_value = shared_data + word_count.value * ( threadIdx.x ? threadIdx.x - 1 : BlockSize );

      if ( ! threadIdx.x ) { functor.init( Reduce::reference( shared_value ) ); }

      // Join previous inclusive scan value to each member
      for ( size_type i = b ; i < e ; ++i ) {
        size_type * const global_value = global_data + word_count.value * i ;
        functor.join( Reduce::reference( shared_value ) , Reduce::reference( global_value ) );
        Reduce::copy( functor , global_value , shared_value );
      }
    }
  }

  return is_last_block ;
}

template< bool DoScan , unsigned ArgBlockSize , class FunctorType >
inline
unsigned cuda_single_inter_block_reduce_scan_shmem( const FunctorType & functor )
{
  return ( ArgBlockSize + 2 ) * ReduceAdapter< FunctorType >::value_size( functor );
}

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( __CUDACC__ ) */
#endif /* KOKKOS_CUDA_REDUCESCAN_HPP */

