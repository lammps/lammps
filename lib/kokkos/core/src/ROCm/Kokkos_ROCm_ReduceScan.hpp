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

#ifndef KOKKOS_ROCM_REDUCESCAN_HPP
#define KOKKOS_ROCM_REDUCESCAN_HPP

#include <Kokkos_Macros.hpp>

/* only compile this file if ROCM is enabled for Kokkos */
#if defined( __HCC__ ) && defined( KOKKOS_ENABLE_ROCM )

//#include <utility>

#include <Kokkos_Parallel.hpp>
#include <impl/Kokkos_FunctorAdapter.hpp>
#include <impl/Kokkos_Error.hpp>
#include <ROCm/Kokkos_ROCm_Vectorization.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//#if __KALMAR_ACCELERATOR__ == 1
KOKKOS_INLINE_FUNCTION
void __syncthreads() [[hc]]
{
   amp_barrier(CLK_LOCAL_MEM_FENCE);
}

#define LT0 ((threadIdx_x+threadIdx_y+threadIdx_z)?0:1)


// returns non-zero if and only if predicate is non-zero for all threads
// note that syncthreads_or uses the first 64 bits of dynamic group memory.
// this reserved memory must be accounted for everwhere 
// that get_dynamic_group_segment_base_pointer is called.
KOKKOS_INLINE_FUNCTION
uint64_t __syncthreads_or(uint64_t  pred) 
{
  uint64_t *shared_var = (uint64_t *)hc::get_dynamic_group_segment_base_pointer();
  if(LT0) *shared_var = 0;
  amp_barrier(CLK_LOCAL_MEM_FENCE);
#if __KALMAR_ACCELERATOR__ == 1
  if (pred) hc::atomic_or_uint64(shared_var,1);
#endif
  amp_barrier(CLK_LOCAL_MEM_FENCE);
  return (*shared_var);
}

KOKKOS_INLINE_FUNCTION
void __threadfence() 
{
   amp_barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
}

KOKKOS_INLINE_FUNCTION
void __threadfence_block() 
{
   amp_barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
}
//#endif
struct ROCm_atomic_CAS {
    template<class OP>
    KOKKOS_INLINE_FUNCTION
    unsigned long operator () (volatile unsigned long * dest, OP &&op){
       unsigned long read,compare,val;
       compare = *dest;
       read = compare;
       do {
         compare = read;
         val = op(compare);
#if __KALMAR_ACCELERATOR__ == 1
         hc::atomic_compare_exchange((uint64_t *)dest,&read,val);
#endif
       } while (read != compare);
       return val;
    }
};

  template<class OP>
  KOKKOS_INLINE_FUNCTION
  unsigned long atomic_cas_op (volatile unsigned long * dest, OP &&op) {
    ROCm_atomic_CAS cas_op;
    return cas_op(dest, std::forward<OP>(op));
  }

  KOKKOS_INLINE_FUNCTION
  unsigned long atomicInc (volatile unsigned long * dest, const unsigned long& val) {
    return atomic_cas_op(dest, [=](unsigned long old){return ((old>=val)?0:(old+1));});
  }


//----------------------------------------------------------------------------

template< typename T >
KOKKOS_INLINE_FUNCTION
void rocm_shfl( T & out , T const & in , int lane ,
  typename std::enable_if< sizeof(int) == sizeof(T) , int >::type width )
{
  *reinterpret_cast<int*>(&out) =
    __shfl( *reinterpret_cast<int const *>(&in) , lane , width );
}

template< typename T >
KOKKOS_INLINE_FUNCTION
void rocm_shfl( T & out , T const & in , int lane ,
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
KOKKOS_INLINE_FUNCTION
void rocm_shfl_down( T & out , T const & in , int delta ,
  typename std::enable_if< sizeof(int) == sizeof(T) , int >::type width )
{
  *reinterpret_cast<int*>(&out) =
    __shfl_down( *reinterpret_cast<int const *>(&in) , delta , width );
}

template< typename T >
KOKKOS_INLINE_FUNCTION
void rocm_shfl_down( T & out , T const & in , int delta ,
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
KOKKOS_INLINE_FUNCTION
void rocm_shfl_up( T & out , T const & in , int delta ,
  typename std::enable_if< sizeof(int) == sizeof(T) , int >::type width )
{
  *reinterpret_cast<int*>(&out) =
    __shfl_up( *reinterpret_cast<int const *>(&in) , delta , width );
}

template< typename T >
KOKKOS_INLINE_FUNCTION
void rocm_shfl_up( T & out , T const & in , int delta ,
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
#if 0
//----------------------------------------------------------------------------
/** \brief  Reduce within a workgroup over team.vector_length(), the "vector" dimension.
 *
 *  This will be called within a nested, intra-team parallel operation.
 *  Use shuffle operations to avoid conflicts with shared memory usage.
 *
 *  Requires:
 *    team.vector_length() is power of 2
 *    team.vector_length() <= 32 (one workgroup)
 *
 *  Cannot use "butterfly" pattern because floating point
 *  addition is non-associative.  Therefore, must broadcast
 *  the final result.
 */
template< class Reducer >
KOKKOS_INLINE_FUNCTION
void rocm_intra_workgroup_vector_reduce( Reducer const & reducer )
{
  static_assert(
    std::is_reference< typename Reducer::reference_type >::value , "" );

  if ( 1 < team.vector_length() ) {

    typename Reducer::value_type tmp ;

    for ( int i = team.vector_length() ; ( i >>= 1 ) ; ) {

      rocm_shfl_down( tmp , reducer.reference() , i , team.vector_length() );

      if ( team.vector_rank() < i ) { reducer.join( reducer.data() , & tmp ); }
    }

    // Broadcast from root "lane" to all other "lanes"

    rocm_shfl( reducer.reference() , reducer.reference() , 0 , team.vector_length() );
  }
}

/** \brief  Inclusive scan over team.vector_length(), the "vector" dimension.
 *
 *  This will be called within a nested, intra-team parallel operation.
 *  Use shuffle operations to avoid conflicts with shared memory usage.
 *
 *  Algorithm is concurrent bottom-up reductions in triangular pattern
 *  where each ROCM thread is the root of a reduction tree from the
 *  zeroth ROCM thread to itself.
 *
 *  Requires:
 *    team.vector_length() is power of 2
 *    team.vector_length() <= 32 (one workgroup)
 */
template< typename ValueType >
KOKKOS_INLINE_FUNCTION
void rocm_intra_workgroup_vector_inclusive_scan( ValueType & local )
{
  ValueType tmp ;

  // Bottom up:
  //   [t] += [t-1] if t >= 1
  //   [t] += [t-2] if t >= 2
  //   [t] += [t-4] if t >= 4
  // ...

  for ( int i = 1 ; i < team.vector_length() ; i <<= 1 ) {

    rocm_shfl_up( tmp , local , i , team.vector_length() );

    if ( i <= team.vector_rank() ) { local += tmp ; }
  }
}
#endif

//----------------------------------------------------------------------------
/*
 *  Algorithmic constraints:
 *   (a) threads with same team.team_rank() have same value
 *   (b) team.vector_length() == power of two
 *   (c) blockDim.z == 1
 */

template< class ValueType , class JoinOp>
KOKKOS_INLINE_FUNCTION
void rocm_intra_workgroup_reduction( const ROCmTeamMember& team, 
                                       ValueType& result,
                                       const JoinOp& join) {

  unsigned int shift = 1;
  int max_active_thread = team.team_size();

  //Reduce over values from threads with different team.team_rank()
  while(team.vector_length() * shift < 32 ) {
    const ValueType tmp = shfl_down(result, team.vector_length()*shift,32u);
    //Only join if upper thread is active (this allows non power of two for team.team_size()
    if(team.team_rank() + shift < max_active_thread)
      join(result , tmp);
    shift*=2;
  }

  result = shfl(result,0,32);
}

template< class ValueType , class JoinOp>
KOKKOS_INLINE_FUNCTION
void rocm_inter_workgroup_reduction( const ROCmTeamMember& team,
                                       ValueType& value,
                                       const JoinOp& join) {

  #define STEP_WIDTH 4
  
  tile_static ValueType sh_result[256];
  int max_active_thread = team.team_size();
  ValueType* result = (ValueType*) & sh_result;
  const unsigned step = 256 / team.vector_length();
  unsigned shift = STEP_WIDTH;
  const int id = team.team_rank()%step==0?team.team_rank()/step:65000;
  if(id < STEP_WIDTH ) {
    result[id] = value;
  }
  team.team_barrier();

  while (shift<=max_active_thread/step) {
    if(shift<=id && shift+STEP_WIDTH>id && team.vector_rank()==0) {
      join(result[id%STEP_WIDTH],value);
    }
    team.team_barrier();
    shift+=STEP_WIDTH;
  }


  value = result[0];
  for(int i = 1; (i*step<max_active_thread) && i<STEP_WIDTH; i++)
    join(value,result[i]);
}

#if 0
template< class ValueType , class JoinOp>
KOKKOS_INLINE_FUNCTION
void rocm_intra_block_reduction( ROCmTeamMember& team,
                                        ValueType& value,
                                        const JoinOp& join,
                                        const int max_active_thread) {
  rocm_intra_workgroup_reduction(team,value,join,max_active_thread);
  rocm_inter_workgroup_reduction(team,value,join,max_active_thread);
}

template< class FunctorType , class JoinOp , class ArgTag = void >
KOKKOS_INLINE_FUNCTION
bool rocm_inter_block_reduction( ROCmTeamMember& team,
                                 typename FunctorValueTraits< FunctorType , ArgTag >::reference_type  value,
                                 typename FunctorValueTraits< FunctorType , ArgTag >::reference_type  neutral,
                                 const JoinOp& join,
                                 ROCm::size_type * const m_scratch_space,
                                 typename FunctorValueTraits< FunctorType , ArgTag >::pointer_type const result,
                                 ROCm::size_type * const m_scratch_flags,
                                 const int max_active_thread) {
#ifdef __ROCM_ARCH__
  typedef typename FunctorValueTraits< FunctorType , ArgTag >::pointer_type pointer_type;
  typedef typename FunctorValueTraits< FunctorType , ArgTag >::value_type value_type;

  //Do the intra-block reduction with shfl operations and static shared memory
  rocm_intra_block_reduction(value,join,max_active_thread);

  const unsigned id = team.team_rank()*team.vector_length() + team.vector_rank();

  //One thread in the block writes block result to global scratch_memory
  if(id == 0 ) {
    pointer_type global = ((pointer_type) m_scratch_space) + blockIdx.x;
    *global = value;
  }

  //One workgroup of last block performs inter block reduction through loading the block values from global scratch_memory
  bool last_block = false;

  team.team_barrier();
  if ( id < 32 ) {
    ROCm::size_type count;

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

      //Reduce all global values with splitting work over threads in one workgroup
      const int step_size = team.vector_length()*team.team_size() < 32 ? team.vector_length()*team.team_size() : 32;
      for(int i=id; i<gridDim.x; i+=step_size) {
        value_type tmp = global[i];
        join(value, tmp);
      }

      //Perform shfl reductions within the workgroup only join if contribution is valid (allows gridDim.x non power of two and <32)
      if (team.vector_length()*team.team_size() > 1) {
        value_type tmp = Kokkos::shfl_down(value, 1,32);
        if( id + 1 < gridDim.x )
          join(value, tmp);
      }
      if (team.vector_length()*team.team_size() > 2) {
        value_type tmp = Kokkos::shfl_down(value, 2,32);
        if( id + 2 < gridDim.x )
          join(value, tmp);
      }
      if (team.vector_length()*team.team_size() > 4) {
        value_type tmp = Kokkos::shfl_down(value, 4,32);
        if( id + 4 < gridDim.x )
          join(value, tmp);
      }
      if (team.vector_length()*team.team_size() > 8) {
        value_type tmp = Kokkos::shfl_down(value, 8,32);
        if( id + 8 < gridDim.x )
          join(value, tmp);
      }
      if (team.vector_length()*team.team_size() > 16) {
        value_type tmp = Kokkos::shfl_down(value, 16,32);
        if( id + 16 < gridDim.x )
          join(value, tmp);
      }
    }
  }

  //The last block has in its thread=0 the global reduction value through "value"
  return last_block;
#else
  return true;
#endif
}
#endif

//----------------------------------------------------------------------------
/*
 *  Algorithmic constraints:
 *   (a) team.team_size() is a power of two
 *   (b) team.team_size() <= 512
 *   (c) team.vector_length() == blockDim.z == 1
 */

template< bool DoScan , class FunctorType , class ArgTag >
KOKKOS_INLINE_FUNCTION
void rocm_intra_block_reduce_scan( const FunctorType & functor ,
                                   const typename FunctorValueTraits< FunctorType , ArgTag >::pointer_type base_data )
{
  typedef FunctorValueTraits< FunctorType , ArgTag >  ValueTraits ;
  typedef FunctorValueJoin<   FunctorType , ArgTag >  ValueJoin ;

  typedef typename ValueTraits::pointer_type  pointer_type ;

  const unsigned value_count   = ValueTraits::value_count( functor );
  const unsigned BlockSizeMask = blockDim_y  - 1 ;

  // Must have power of two thread count

  if ( BlockSizeMask & blockDim_y ) { Kokkos::abort("ROCm::rocm_intra_block_scan requires power-of-two blockDim"); }

#define BLOCK_REDUCE_STEP( R , TD , S )  \
  if ( ! (( R & ((1<<(S+1))-1) )|(blockDim_y<(1<<(S+1)))) ) { ValueJoin::join( functor , TD , (TD - (value_count<<S)) ); }

#define BLOCK_SCAN_STEP( TD , N , S )  \
  if ( N == (1<<S) ) { ValueJoin::join( functor , TD , (TD - (value_count<<S))); }
#define KOKKOS_IMPL_ROCM_SYNCWF __threadfence_block()

  const unsigned     rtid_intra = threadIdx_y ^ BlockSizeMask ;
  const pointer_type tdata_intra = base_data + value_count * threadIdx_y ;

  { // Intra-workgroup reduction: min blocksize of 64
    KOKKOS_IMPL_ROCM_SYNCWF;
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,0)
    KOKKOS_IMPL_ROCM_SYNCWF;
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,1)
    KOKKOS_IMPL_ROCM_SYNCWF;
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,2)
    KOKKOS_IMPL_ROCM_SYNCWF;
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,3)
    KOKKOS_IMPL_ROCM_SYNCWF;
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,4)
    KOKKOS_IMPL_ROCM_SYNCWF;
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,5)
    KOKKOS_IMPL_ROCM_SYNCWF;
  }

  __syncthreads(); // Wait for all workgroups to reduce

  { // Inter-workgroup reduce-scan by a single workgroup to avoid extra synchronizations
    if(threadIdx_y < value_count) {
      for(int i=blockDim_y-65; i>0; i-= 64)
        ValueJoin::join( functor , base_data + (blockDim_y-1)*value_count + threadIdx_y ,  base_data + i*value_count + threadIdx_y );
    }
    __syncthreads();
#if 0
    const unsigned rtid_inter = ( threadIdx_y ^ BlockSizeMask ) << ROCmTraits::WavefrontIndexShift ;

    if ( rtid_inter < blockDim_y ) {


      const pointer_type tdata_inter = base_data + value_count * ( rtid_inter ^ BlockSizeMask );
//
// remove these comments
// for rocm, we start with a block size of 64, so the 5 step is already done.
// The remaining steps are only done if block size is > 64, so we leave them
// in place until we tune blocksize for performance, then remove the ones 
// that will never be used.
//      if ( (1<<6) < BlockSizeMask ) { __threadfence_block(); BLOCK_REDUCE_STEP(rtid_inter,tdata_inter,6) }
//      if ( (1<<7) < BlockSizeMask ) { __threadfence_block(); BLOCK_REDUCE_STEP(rtid_inter,tdata_inter,7) }
//      if ( (1<<8) < BlockSizeMask ) { __threadfence_block(); BLOCK_REDUCE_STEP(rtid_inter,tdata_inter,8) }
//      if ( (1<<9) < BlockSizeMask ) { __threadfence_block(); BLOCK_REDUCE_STEP(rtid_inter,tdata_inter,9) }


      if ( DoScan ) {

        int n = ( rtid_inter &  64 ) ?  64 : (
                ( rtid_inter & 128 ) ? 128 : (
                ( rtid_inter & 256 ) ? 256 : 0 ));

        if ( ! ( rtid_inter + n < blockDim_y ) ) n = 0 ;

        __threadfence_block(); BLOCK_SCAN_STEP(tdata_inter,n,8)
        __threadfence_block(); BLOCK_SCAN_STEP(tdata_inter,n,7)
        __threadfence_block(); BLOCK_SCAN_STEP(tdata_inter,n,6)
//        __threadfence_block(); BLOCK_SCAN_STEP(tdata_inter,n,5)
      }
    }
#endif
  }

  __syncthreads(); // Wait for inter-workgroup reduce-scan to complete

  if ( DoScan ) {
    int n = ( rtid_intra &  1 ) ?  1 : (
            ( rtid_intra &  2 ) ?  2 : (
            ( rtid_intra &  4 ) ?  4 : (
            ( rtid_intra &  8 ) ?  8 : (
            ( rtid_intra & 16 ) ? 16 : (
            ( rtid_intra & 32 ) ? 32 : 0 )))));

    if ( ! ( rtid_intra + n < blockDim_y ) ) n = 0 ;

//    BLOCK_SCAN_STEP(tdata_intra,n,5) __threadfence_block();
//    BLOCK_SCAN_STEP(tdata_intra,n,4) __threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,3) __threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,2) __threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,1) __threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,0) __threadfence_block();
  }

#undef BLOCK_SCAN_STEP
#undef BLOCK_REDUCE_STEP
#undef KOKKOS_IMPL_ROCM_SYNCWF
}

//----------------------------------------------------------------------------
/**\brief  Input value-per-thread starting at 'shared_data'.
 *         Reduction value at last thread's location.
 *
 *  If 'DoScan' then write blocks' scan values and block-groups' scan values.
 *
 *  Global reduce result is in the last threads' 'shared_data' location.
 */
using ROCM  = Kokkos::Experimental::ROCm ;

template< bool DoScan , class FunctorType , class ArgTag >
KOKKOS_INLINE_FUNCTION
bool rocm_single_inter_block_reduce_scan( const FunctorType     & functor ,
                                          const ROCM::size_type   block_id ,
                                          const ROCM::size_type   block_count ,
                                          typename FunctorValueTraits<FunctorType, ArgTag>::value_type * const shared_data ,
                                          typename FunctorValueTraits<FunctorType, ArgTag>::value_type * const global_data ,
                                          ROCM::size_type * const global_flags )
{
  typedef ROCM::size_type                  size_type ;
  typedef FunctorValueTraits< FunctorType , ArgTag >  ValueTraits ;
  typedef FunctorValueJoin<   FunctorType , ArgTag >  ValueJoin ;
  typedef FunctorValueInit<   FunctorType , ArgTag >  ValueInit ;
  typedef FunctorValueOps<    FunctorType , ArgTag >  ValueOps ;

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;
  typedef typename ValueTraits::value_type      value_type ;

  // '__ffs' = position of the least significant bit set to 1.
  // blockDim_y is guaranteed to be a power of two so this
  // is the integral shift value that can replace an integral divide.
  //  const unsigned long BlockSizeShift = __ffs( blockDim_y ) - 1 ;
  const unsigned long BlockSizeShift = __lastbit_u32_u32( blockDim_y )  ;
  const unsigned long BlockSizeMask  = blockDim_y - 1 ;

  // Must have power of two thread count
  if ( BlockSizeMask & blockDim_y ) { Kokkos::abort("ROCm::rocm_single_inter_block_reduce_scan requires power-of-two blockDim"); }

  const integral_nonzero_constant< size_type , ValueTraits::StaticValueSize / sizeof(value_type) >
    word_count( ValueTraits::value_size( functor )/ sizeof(value_type) );

  // Reduce the accumulation for the entire block.
  rocm_intra_block_reduce_scan<false,FunctorType,ArgTag>( functor , pointer_type(shared_data) );

  {
    // Write accumulation total to global scratch space.
    // Accumulation total is the last thread's data.
    value_type * const shared = shared_data +  
                                   word_count.value * BlockSizeMask ;
    value_type * const global = global_data + word_count.value * block_id ;

    for ( int i = int(threadIdx_y) ; i < word_count.value ; i += blockDim_y ) { global[i] = shared[i] ; }
  }

  // Contributing blocks note that their contribution has been completed via an atomic-increment flag
  // If this block is not the last block to contribute to this group then the block is done.
    
  const bool is_last_block =
    !  __syncthreads_or( threadIdx_y ? 0 : ( 1 + atomicInc( global_flags , block_count - 1 ) < block_count ) );
  if ( is_last_block ) {

    const size_type b = ( long(block_count) * long(threadIdx_y )) >> BlockSizeShift ;
    const size_type e = ( long(block_count) * long(threadIdx_y + 1 ) ) >> BlockSizeShift ;

    {
      value_type * const shared_ptr = shared_data + word_count.value * threadIdx_y ;
      ValueInit::init( functor , shared_ptr );


      for ( size_type i = b ; i < e ; ++i ) {
        ValueJoin::join( functor , shared_ptr , global_data + word_count.value * i );
      }
    }
    rocm_intra_block_reduce_scan<DoScan,FunctorType,ArgTag>( functor , pointer_type(shared_data) );

    if ( DoScan ) {
      value_type * const shared_value = shared_data + word_count.value * ( threadIdx_y ? threadIdx_y - 1 : blockDim_y );

      if ( ! threadIdx_y ) { ValueInit::init( functor , shared_value ); }

      // Join previous inclusive scan value to each member
      for ( size_type i = b ; i < e ; ++i ) {
        value_type * const global_value = global_data + word_count.value * i ;
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
unsigned rocm_single_inter_block_reduce_scan_shmem( const FunctorType & functor , const unsigned BlockSize )
{
  return ( BlockSize + 2 ) * Impl::FunctorValueTraits< FunctorType , ArgTag >::value_size( functor );
}

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( __ROCMCC__ ) */
#endif /* KOKKOS_ROCM_REDUCESCAN_HPP */

