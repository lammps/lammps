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

#ifndef KOKKOS_CUDA_INTERNAL_HPP
#define KOKKOS_CUDA_INTERNAL_HPP

namespace Kokkos {
namespace Impl {

void cuda_internal_error_throw( cudaError e , const char * name, const char * file = NULL, const int line = 0 );

void cuda_device_synchronize();

inline
void cuda_internal_safe_call( cudaError e , const char * name, const char * file = NULL, const int line = 0)
{
  if ( cudaSuccess != e ) { cuda_internal_error_throw( e , name, file, line ); }
}

template<class DriverType>
int cuda_get_max_block_size(const typename DriverType::functor_type & f) {
#if ( CUDA_VERSION < 6050 )
  return 256;
#else
  bool Large = ( CudaTraits::ConstantMemoryUseThreshold < sizeof(DriverType) );

  int numBlocks;
  if(Large) {
    int blockSize=32;
    int sharedmem = FunctorTeamShmemSize< typename DriverType::functor_type  >::value( f , blockSize );
    cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        &numBlocks,
        cuda_parallel_launch_constant_memory<DriverType>,
        blockSize,
        sharedmem);

    while (blockSize<1024 && numBlocks>0) {
      blockSize*=2;
      sharedmem = FunctorTeamShmemSize< typename DriverType::functor_type  >::value( f , blockSize );

      cudaOccupancyMaxActiveBlocksPerMultiprocessor(
          &numBlocks,
          cuda_parallel_launch_constant_memory<DriverType>,
          blockSize,
          sharedmem);
    }
    if(numBlocks>0) return blockSize;
    else return blockSize/2;
  } else {
    int blockSize=32;
    int sharedmem = FunctorTeamShmemSize< typename DriverType::functor_type  >::value( f , blockSize );
    cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        &numBlocks,
        cuda_parallel_launch_local_memory<DriverType>,
        blockSize,
        sharedmem);

    while (blockSize<1024 && numBlocks>0) {
      blockSize*=2;
      sharedmem = FunctorTeamShmemSize< typename DriverType::functor_type  >::value( f , blockSize );

      cudaOccupancyMaxActiveBlocksPerMultiprocessor(
          &numBlocks,
          cuda_parallel_launch_local_memory<DriverType>,
          blockSize,
          sharedmem);
    }
    if(numBlocks>0) return blockSize;
    else return blockSize/2;
  }
#endif
}

template<class DriverType>
int cuda_get_opt_block_size(const typename DriverType::functor_type & f) {
#if ( CUDA_VERSION < 6050 )
  return 256;
#else
  bool Large = ( CudaTraits::ConstantMemoryUseThreshold < sizeof(DriverType) );

  int blockSize=16;
  int numBlocks;
  int sharedmem;
  int maxOccupancy=0;
  int bestBlockSize=0;

  if(Large) {
    while(blockSize<1024) {
      blockSize*=2;

      //calculate the occupancy with that optBlockSize and check whether its larger than the largest one found so far
      sharedmem = FunctorTeamShmemSize< typename DriverType::functor_type  >::value( f , blockSize );
      cudaOccupancyMaxActiveBlocksPerMultiprocessor(
              &numBlocks,
              cuda_parallel_launch_constant_memory<DriverType>,
              blockSize,
              sharedmem);
      if(maxOccupancy < numBlocks*blockSize) {
        maxOccupancy = numBlocks*blockSize;
        bestBlockSize = blockSize;
      }
    }
  } else {
    while(blockSize<1024) {
      blockSize*=2;
      sharedmem = FunctorTeamShmemSize< typename DriverType::functor_type  >::value( f , blockSize );

      cudaOccupancyMaxActiveBlocksPerMultiprocessor(
              &numBlocks,
              cuda_parallel_launch_local_memory<DriverType>,
              blockSize,
              sharedmem);

      if(maxOccupancy < numBlocks*blockSize) {
        maxOccupancy = numBlocks*blockSize;
        bestBlockSize = blockSize;
      }
    }
  }
  return bestBlockSize;
#endif
}

}
}

#define CUDA_SAFE_CALL( call )  \
	Kokkos::Impl::cuda_internal_safe_call( call , #call, __FILE__, __LINE__ )

#endif /* #ifndef KOKKOS_CUDA_INTERNAL_HPP */

