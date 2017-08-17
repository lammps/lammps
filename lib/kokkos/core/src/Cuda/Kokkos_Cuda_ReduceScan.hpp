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

#ifndef KOKKOS_CUDA_REDUCESCAN_HPP
#define KOKKOS_CUDA_REDUCESCAN_HPP

#include <Kokkos_Macros.hpp>
#if defined( __CUDACC__ ) && defined( KOKKOS_ENABLE_CUDA )

#include <utility>

#include <Kokkos_Parallel.hpp>
#include <impl/Kokkos_FunctorAdapter.hpp>
#include <impl/Kokkos_Error.hpp>
#include <Cuda/Kokkos_Cuda_Vectorization.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
// Shuffle operations require input to be a register (stack) variable

template< typename T >
__device__ inline
void cuda_shfl( T & out , T const & in , int lane ,
  typename std::enable_if< sizeof(int) == sizeof(T) , int >::type width )
{
  *reinterpret_cast<int*>(&out) =
    __shfl( *reinterpret_cast<int const *>(&in) , lane , width );
}

template< typename T >
__device__ inline
void cuda_shfl( T & out , T const & in , int lane ,
  typename std::enable_if
    < ( sizeof(int) < sizeof(T) ) && ( 0 == ( sizeof(T) % sizeof(int) ) )
    , int >::type width )
{
  enum : int { N = sizeof(T) / sizeof(int) };

  for ( int i = 0 ; i < N ; ++i ) {
    reinterpret_cast<int*>(&out)[i] =
      __shfl( reinterpret_cast<int const *>(&in)[i] , lane , width );
  }
}

//----------------------------------------------------------------------------

template< typename T >
__device__ inline
void cuda_shfl_down( T & out , T const & in , int delta ,
  typename std::enable_if< sizeof(int) == sizeof(T) , int >::type width )
{
  *reinterpret_cast<int*>(&out) =
    __shfl_down( *reinterpret_cast<int const *>(&in) , delta , width );
}

template< typename T >
__device__ inline
void cuda_shfl_down( T & out , T const & in , int delta ,
  typename std::enable_if
    < ( sizeof(int) < sizeof(T) ) && ( 0 == ( sizeof(T) % sizeof(int) ) )
    , int >::type width )
{
  enum : int { N = sizeof(T) / sizeof(int) };

  for ( int i = 0 ; i < N ; ++i ) {
    reinterpret_cast<int*>(&out)[i] =
      __shfl_down( reinterpret_cast<int const *>(&in)[i] , delta , width );
  }
}

//----------------------------------------------------------------------------

template< typename T >
__device__ inline
void cuda_shfl_up( T & out , T const & in , int delta ,
  typename std::enable_if< sizeof(int) == sizeof(T) , int >::type width )
{
  *reinterpret_cast<int*>(&out) =
    __shfl_up( *reinterpret_cast<int const *>(&in) , delta , width );
}

template< typename T >
__device__ inline
void cuda_shfl_up( T & out , T const & in , int delta ,
  typename std::enable_if
    < ( sizeof(int) < sizeof(T) ) && ( 0 == ( sizeof(T) % sizeof(int) ) )
    , int >::type width )
{
  enum : int { N = sizeof(T) / sizeof(int) };

  for ( int i = 0 ; i < N ; ++i ) {
    reinterpret_cast<int*>(&out)[i] =
      __shfl_up( reinterpret_cast<int const *>(&in)[i] , delta , width );
  }
}

//----------------------------------------------------------------------------
/*
 *  Algorithmic constraints:
 *   (a) threads with same threadIdx.y have same value
 *   (b) blockDim.x == power of two
 *   (c) blockDim.z == 1
 */

template< class ValueType , class JoinOp>
__device__
inline void cuda_intra_warp_reduction( ValueType& result,
                                       const JoinOp& join,
                                       const uint32_t max_active_thread = blockDim.y) {

  unsigned int shift = 1;

  //Reduce over values from threads with different threadIdx.y
  while(blockDim.x * shift < 32 ) {
    const ValueType tmp = shfl_down(result, blockDim.x*shift,32u);
    //Only join if upper thread is active (this allows non power of two for blockDim.y
    if(threadIdx.y + shift < max_active_thread)
      join(result , tmp);
    shift*=2;
  }

  result = shfl(result,0,32);
}

template< class ValueType , class JoinOp>
__device__
inline void cuda_inter_warp_reduction( ValueType& value,
                                       const JoinOp& join,
                                       const int max_active_thread = blockDim.y) {

  #define STEP_WIDTH 4
  // Depending on the ValueType _shared__ memory must be aligned up to 8byte boundaries
  // The reason not to use ValueType directly is that for types with constructors it
  // could lead to race conditions
  __shared__ double sh_result[(sizeof(ValueType)+7)/8*STEP_WIDTH];
  ValueType* result = (ValueType*) & sh_result;
  const int step = 32 / blockDim.x;
  int shift = STEP_WIDTH;
  const int id = threadIdx.y%step==0?threadIdx.y/step:65000;
  if(id < STEP_WIDTH ) {
    result[id] = value;
  }
  __syncthreads();
  while (shift<=max_active_thread/step) {
    if(shift<=id && shift+STEP_WIDTH>id && threadIdx.x==0) {
      join(result[id%STEP_WIDTH],value);
    }
    __syncthreads();
    shift+=STEP_WIDTH;
  }


  value = result[0];
  for(int i = 1; (i*step<max_active_thread) && i<STEP_WIDTH; i++)
    join(value,result[i]);
}

template< class ValueType , class JoinOp>
__device__
inline void cuda_intra_block_reduction( ValueType& value,
                                        const JoinOp& join,
                                        const int max_active_thread = blockDim.y) {
  cuda_intra_warp_reduction(value,join,max_active_thread);
  cuda_inter_warp_reduction(value,join,max_active_thread);
}

template< class FunctorType , class JoinOp , class ArgTag = void >
__device__
bool cuda_inter_block_reduction( typename FunctorValueTraits< FunctorType , ArgTag >::reference_type  value,
                                 typename FunctorValueTraits< FunctorType , ArgTag >::reference_type  neutral,
                                 const JoinOp& join,
                                 Cuda::size_type * const m_scratch_space,
                                 typename FunctorValueTraits< FunctorType , ArgTag >::pointer_type const result,
                                 Cuda::size_type * const m_scratch_flags,
                                 const int max_active_thread = blockDim.y) {
#ifdef __CUDA_ARCH__
  typedef typename FunctorValueTraits< FunctorType , ArgTag >::pointer_type pointer_type;
  typedef typename FunctorValueTraits< FunctorType , ArgTag >::value_type value_type;

  //Do the intra-block reduction with shfl operations and static shared memory
  cuda_intra_block_reduction(value,join,max_active_thread);

  const int id = threadIdx.y*blockDim.x + threadIdx.x;

  //One thread in the block writes block result to global scratch_memory
  if(id == 0 ) {
    pointer_type global = ((pointer_type) m_scratch_space) + blockIdx.x;
    *global = value;
  }

  //One warp of last block performs inter block reduction through loading the block values from global scratch_memory
  bool last_block = false;

  __syncthreads();
  if ( id < 32 ) {
    Cuda::size_type count;

    //Figure out whether this is the last block
    if(id == 0)
      count = Kokkos::atomic_fetch_add(m_scratch_flags,1);
    count = Kokkos::shfl(count,0,32);

    //Last block does the inter block reduction
    if( count == gridDim.x - 1) {
      //set flag back to zero
      if(id == 0)
        *m_scratch_flags = 0;
      last_block = true;
      value = neutral;

      pointer_type const volatile global = (pointer_type) m_scratch_space ;

      //Reduce all global values with splitting work over threads in one warp
      const int step_size = blockDim.x*blockDim.y < 32 ? blockDim.x*blockDim.y : 32;
      for(int i=id; i<(int)gridDim.x; i+=step_size) {
        value_type tmp = global[i];
        join(value, tmp);
      }

      //Perform shfl reductions within the warp only join if contribution is valid (allows gridDim.x non power of two and <32)
      if (int(blockDim.x*blockDim.y) > 1) {
        value_type tmp = Kokkos::shfl_down(value, 1,32);
        if( id + 1 < int(gridDim.x) )
          join(value, tmp);
      }
      int active = __ballot(1);
      if (int(blockDim.x*blockDim.y) > 2) {
        value_type tmp = Kokkos::shfl_down(value, 2,32);
        if( id + 2 < int(gridDim.x) )
          join(value, tmp);
      }
      active += __ballot(1);
      if (int(blockDim.x*blockDim.y) > 4) {
        value_type tmp = Kokkos::shfl_down(value, 4,32);
        if( id + 4 < int(gridDim.x) )
          join(value, tmp);
      }
      active += __ballot(1);
      if (int(blockDim.x*blockDim.y) > 8) {
        value_type tmp = Kokkos::shfl_down(value, 8,32);
        if( id + 8 < int(gridDim.x) )
          join(value, tmp);
      }
      active += __ballot(1);
      if (int(blockDim.x*blockDim.y) > 16) {
        value_type tmp = Kokkos::shfl_down(value, 16,32);
        if( id + 16 < int(gridDim.x) )
          join(value, tmp);
      }
      active += __ballot(1);
    }
  }
  //The last block has in its thread=0 the global reduction value through "value"
  return last_block;
#else
  return true;
#endif
}

template< class ReducerType >
__device__ inline
typename std::enable_if< Kokkos::is_reducer<ReducerType>::value >::type
cuda_intra_warp_reduction( const ReducerType& reducer,
                           const uint32_t max_active_thread = blockDim.y) {

  typedef typename ReducerType::value_type ValueType;

  unsigned int shift = 1;

  ValueType result = reducer.reference();
  //Reduce over values from threads with different threadIdx.y
  while(blockDim.x * shift < 32 ) {
    const ValueType tmp = shfl_down(result, blockDim.x*shift,32u);
    //Only join if upper thread is active (this allows non power of two for blockDim.y
    if(threadIdx.y + shift < max_active_thread)
      reducer.join(result , tmp);
    shift*=2;
  }

  result = shfl(result,0,32);
  reducer.reference() = result;
}

template< class ReducerType >
__device__ inline
typename std::enable_if< Kokkos::is_reducer<ReducerType>::value >::type
cuda_inter_warp_reduction( const ReducerType& reducer,
                           const int max_active_thread = blockDim.y) {

  typedef typename ReducerType::value_type ValueType;

  #define STEP_WIDTH 4
  // Depending on the ValueType _shared__ memory must be aligned up to 8byte boundaries
  // The reason not to use ValueType directly is that for types with constructors it
  // could lead to race conditions
  __shared__ double sh_result[(sizeof(ValueType)+7)/8*STEP_WIDTH];
  ValueType* result = (ValueType*) & sh_result;
  ValueType value = reducer.reference();
  const int step = 32 / blockDim.x;
  int shift = STEP_WIDTH;
  const int id = threadIdx.y%step==0?threadIdx.y/step:65000;
  if(id < STEP_WIDTH ) {
    result[id] = value;
  }
  __syncthreads();
  while (shift<=max_active_thread/step) {
    if(shift<=id && shift+STEP_WIDTH>id && threadIdx.x==0) {
      reducer.join(result[id%STEP_WIDTH],value);
    }
    __syncthreads();
    shift+=STEP_WIDTH;
  }


  value = result[0];
  for(int i = 1; (i*step<max_active_thread) && i<STEP_WIDTH; i++)
    reducer.join(value,result[i]);

  reducer.reference() = value;
}

template< class ReducerType >
__device__ inline
typename std::enable_if< Kokkos::is_reducer<ReducerType>::value >::type
cuda_intra_block_reduction( const ReducerType& reducer,
                            const int max_active_thread = blockDim.y) {
  cuda_intra_warp_reduction(reducer,max_active_thread);
  cuda_inter_warp_reduction(reducer,max_active_thread);
}

template< class ReducerType>
__device__ inline
typename std::enable_if< Kokkos::is_reducer<ReducerType>::value , bool >::type
cuda_inter_block_reduction( const ReducerType& reducer,
                            Cuda::size_type * const m_scratch_space,
                            Cuda::size_type * const m_scratch_flags,
                            const int max_active_thread = blockDim.y) {
#ifdef __CUDA_ARCH__
  typedef typename ReducerType::value_type* pointer_type;
  typedef typename ReducerType::value_type value_type;

  //Do the intra-block reduction with shfl operations and static shared memory
  cuda_intra_block_reduction(reducer,max_active_thread);

  value_type value = reducer.reference();

  const int id = threadIdx.y*blockDim.x + threadIdx.x;

  //One thread in the block writes block result to global scratch_memory
  if(id == 0 ) {
    pointer_type global = ((pointer_type) m_scratch_space) + blockIdx.x;
    *global = value;
  }

  //One warp of last block performs inter block reduction through loading the block values from global scratch_memory
  bool last_block = false;

  __syncthreads();
  if ( id < 32 ) {
    Cuda::size_type count;

    //Figure out whether this is the last block
    if(id == 0)
      count = Kokkos::atomic_fetch_add(m_scratch_flags,1);
    count = Kokkos::shfl(count,0,32);

    //Last block does the inter block reduction
    if( count == gridDim.x - 1) {
      //set flag back to zero
      if(id == 0)
        *m_scratch_flags = 0;
      last_block = true;
      reducer.init(value);

      pointer_type const volatile global = (pointer_type) m_scratch_space ;

      //Reduce all global values with splitting work over threads in one warp
      const int step_size = blockDim.x*blockDim.y < 32 ? blockDim.x*blockDim.y : 32;
      for(int i=id; i<(int)gridDim.x; i+=step_size) {
        value_type tmp = global[i];
        reducer.join(value, tmp);
      }

      //Perform shfl reductions within the warp only join if contribution is valid (allows gridDim.x non power of two and <32)
      if (int(blockDim.x*blockDim.y) > 1) {
        value_type tmp = Kokkos::shfl_down(value, 1,32);
        if( id + 1 < int(gridDim.x) )
          reducer.join(value, tmp);
      }
      int active = __ballot(1);
      if (int(blockDim.x*blockDim.y) > 2) {
        value_type tmp = Kokkos::shfl_down(value, 2,32);
        if( id + 2 < int(gridDim.x) )
          reducer.join(value, tmp);
      }
      active += __ballot(1);
      if (int(blockDim.x*blockDim.y) > 4) {
        value_type tmp = Kokkos::shfl_down(value, 4,32);
        if( id + 4 < int(gridDim.x) )
          reducer.join(value, tmp);
      }
      active += __ballot(1);
      if (int(blockDim.x*blockDim.y) > 8) {
        value_type tmp = Kokkos::shfl_down(value, 8,32);
        if( id + 8 < int(gridDim.x) )
          reducer.join(value, tmp);
      }
      active += __ballot(1);
      if (int(blockDim.x*blockDim.y) > 16) {
        value_type tmp = Kokkos::shfl_down(value, 16,32);
        if( id + 16 < int(gridDim.x) )
          reducer.join(value, tmp);
      }
      active += __ballot(1);
    }
  }

  //The last block has in its thread=0 the global reduction value through "value"
  return last_block;
#else
  return true;
#endif
}

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
//----------------------------------------------------------------------------
/*
 *  Algorithmic constraints:
 *   (a) blockDim.y is a power of two
 *   (b) blockDim.y <= 512
 *   (c) blockDim.x == blockDim.z == 1
 */

template< bool DoScan , class FunctorType , class ArgTag >
__device__
void cuda_intra_block_reduce_scan( const FunctorType & functor ,
                                   const typename FunctorValueTraits< FunctorType , ArgTag >::pointer_type base_data )
{
  typedef FunctorValueTraits< FunctorType , ArgTag >  ValueTraits ;
  typedef FunctorValueJoin<   FunctorType , ArgTag >  ValueJoin ;

  typedef typename ValueTraits::pointer_type  pointer_type ;

  const unsigned value_count   = ValueTraits::value_count( functor );
  const unsigned BlockSizeMask = blockDim.y - 1 ;

  // Must have power of two thread count

  if ( BlockSizeMask & blockDim.y ) { Kokkos::abort("Cuda::cuda_intra_block_scan requires power-of-two blockDim"); }

#define BLOCK_REDUCE_STEP( R , TD , S )  \
  if ( ! ( R & ((1<<(S+1))-1) ) ) { ValueJoin::join( functor , TD , (TD - (value_count<<S)) ); }

#define BLOCK_SCAN_STEP( TD , N , S )  \
  if ( N == (1<<S) ) { ValueJoin::join( functor , TD , (TD - (value_count<<S))); }

  const unsigned     rtid_intra = threadIdx.y ^ BlockSizeMask ;
  const pointer_type tdata_intra = base_data + value_count * threadIdx.y ;

  { // Intra-warp reduction:
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,0)
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,1)
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,2)
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,3)
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,4)
  }

  __syncthreads(); // Wait for all warps to reduce

  { // Inter-warp reduce-scan by a single warp to avoid extra synchronizations
    const unsigned rtid_inter = ( threadIdx.y ^ BlockSizeMask ) << CudaTraits::WarpIndexShift ;

    if ( rtid_inter < blockDim.y ) {

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

        if ( ! ( rtid_inter + n < blockDim.y ) ) n = 0 ;

        __threadfence_block(); BLOCK_SCAN_STEP(tdata_inter,n,8)
        __threadfence_block(); BLOCK_SCAN_STEP(tdata_inter,n,7)
        __threadfence_block(); BLOCK_SCAN_STEP(tdata_inter,n,6)
        __threadfence_block(); BLOCK_SCAN_STEP(tdata_inter,n,5)
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

    if ( ! ( rtid_intra + n < blockDim.y ) ) n = 0 ;
    #ifdef KOKKOS_IMPL_CUDA_CLANG_WORKAROUND
    BLOCK_SCAN_STEP(tdata_intra,n,4) __syncthreads();//__threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,3) __syncthreads();//__threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,2) __syncthreads();//__threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,1) __syncthreads();//__threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,0) __syncthreads();
    #else
    BLOCK_SCAN_STEP(tdata_intra,n,4) __threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,3) __threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,2) __threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,1) __threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,0) __threadfence_block();
    #endif
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
template< bool DoScan , class FunctorType , class ArgTag >
__device__
bool cuda_single_inter_block_reduce_scan( const FunctorType     & functor ,
                                          const Cuda::size_type   block_id ,
                                          const Cuda::size_type   block_count ,
                                          Cuda::size_type * const shared_data ,
                                          Cuda::size_type * const global_data ,
                                          Cuda::size_type * const global_flags )
{
  typedef Cuda::size_type                  size_type ;
  typedef FunctorValueTraits< FunctorType , ArgTag >  ValueTraits ;
  typedef FunctorValueJoin<   FunctorType , ArgTag >  ValueJoin ;
  typedef FunctorValueInit<   FunctorType , ArgTag >  ValueInit ;
  typedef FunctorValueOps<    FunctorType , ArgTag >  ValueOps ;

  typedef typename ValueTraits::pointer_type    pointer_type ;
  //typedef typename ValueTraits::reference_type  reference_type ;

  // '__ffs' = position of the least significant bit set to 1.
  // 'blockDim.y' is guaranteed to be a power of two so this
  // is the integral shift value that can replace an integral divide.
  const unsigned BlockSizeShift = __ffs( blockDim.y ) - 1 ;
  const unsigned BlockSizeMask  = blockDim.y - 1 ;

  // Must have power of two thread count
  if ( BlockSizeMask & blockDim.y ) { Kokkos::abort("Cuda::cuda_single_inter_block_reduce_scan requires power-of-two blockDim"); }

  const integral_nonzero_constant< size_type , ValueTraits::StaticValueSize / sizeof(size_type) >
    word_count( ValueTraits::value_size( functor ) / sizeof(size_type) );

  // Reduce the accumulation for the entire block.
  cuda_intra_block_reduce_scan<false,FunctorType,ArgTag>( functor , pointer_type(shared_data) );

  {
    // Write accumulation total to global scratch space.
    // Accumulation total is the last thread's data.
    size_type * const shared = shared_data + word_count.value * BlockSizeMask ;
    size_type * const global = global_data + word_count.value * block_id ;

//#if (__CUDA_ARCH__ < 500)
    for ( int i = int(threadIdx.y) ; i < int(word_count.value) ; i += int(blockDim.y) ) { global[i] = shared[i] ; }
//#else
//    for ( size_type i = 0 ; i < word_count.value ; i += 1 ) { global[i] = shared[i] ; }
//#endif

  }

  // Contributing blocks note that their contribution has been completed via an atomic-increment flag
  // If this block is not the last block to contribute to this group then the block is done.
  const bool is_last_block =
    ! __syncthreads_or( threadIdx.y ? 0 : ( 1 + atomicInc( global_flags , block_count - 1 ) < block_count ) );

  if ( is_last_block ) {

    const size_type b = ( long(block_count) * long(threadIdx.y) ) >> BlockSizeShift ;
    const size_type e = ( long(block_count) * long( threadIdx.y + 1 ) ) >> BlockSizeShift ;

    {
      void * const shared_ptr = shared_data + word_count.value * threadIdx.y ;
      /* reference_type shared_value = */ ValueInit::init( functor , shared_ptr );

      for ( size_type i = b ; i < e ; ++i ) {
        ValueJoin::join( functor , shared_ptr , global_data + word_count.value * i );
      }
    }

    cuda_intra_block_reduce_scan<DoScan,FunctorType,ArgTag>( functor , pointer_type(shared_data) );

    if ( DoScan ) {

      size_type * const shared_value = shared_data + word_count.value * ( threadIdx.y ? threadIdx.y - 1 : blockDim.y );

      if ( ! threadIdx.y ) { ValueInit::init( functor , shared_value ); }

      // Join previous inclusive scan value to each member
      for ( size_type i = b ; i < e ; ++i ) {
        size_type * const global_value = global_data + word_count.value * i ;
        ValueJoin::join( functor , shared_value , global_value );
        ValueOps ::copy( functor , global_value , shared_value );
      }
    }
  }

  return is_last_block ;
}

// Size in bytes required for inter block reduce or scan
template< bool DoScan , class FunctorType , class ArgTag >
inline
unsigned cuda_single_inter_block_reduce_scan_shmem( const FunctorType & functor , const unsigned BlockSize )
{
  return ( BlockSize + 2 ) * Impl::FunctorValueTraits< FunctorType , ArgTag >::value_size( functor );
}

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( __CUDACC__ ) */
#endif /* KOKKOS_CUDA_REDUCESCAN_HPP */

