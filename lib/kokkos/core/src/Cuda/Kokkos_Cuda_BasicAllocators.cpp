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

#include <Kokkos_Macros.hpp>

/* only compile this file if CUDA is enabled for Kokkos */
#ifdef KOKKOS_HAVE_CUDA

#include <impl/Kokkos_Error.hpp>
#include <Cuda/Kokkos_Cuda_BasicAllocators.hpp>
#include <Cuda/Kokkos_Cuda_Error.hpp>

#include <sstream>

namespace Kokkos { namespace Impl {


/*--------------------------------------------------------------------------*/
TextureAttribute::TextureAttribute(  void * const alloc_ptr
                                   , size_t alloc_size
                                   , cudaChannelFormatDesc const & desc
                                  )
  : m_tex_obj(0)
{
  cuda_device_synchronize();

  struct cudaResourceDesc resDesc ;
  struct cudaTextureDesc  texDesc ;

  memset( & resDesc , 0 , sizeof(resDesc) );
  memset( & texDesc , 0 , sizeof(texDesc) );

  resDesc.resType                = cudaResourceTypeLinear ;
  resDesc.res.linear.desc        = desc ;
  resDesc.res.linear.sizeInBytes = alloc_size ;
  resDesc.res.linear.devPtr      = alloc_ptr ;

  CUDA_SAFE_CALL( cudaCreateTextureObject( & m_tex_obj , & resDesc, & texDesc, NULL) );

  cuda_device_synchronize();
}


TextureAttribute::~TextureAttribute()
{
  if (m_tex_obj) {
    cudaDestroyTextureObject( m_tex_obj );
  }
}

/*--------------------------------------------------------------------------*/

void * CudaMallocAllocator::allocate( size_t size )
{
  void * ptr = NULL;

  CUDA_SAFE_CALL( cudaMalloc( &ptr, size ) );

  return ptr;
}

void CudaMallocAllocator::deallocate( void * ptr, size_t /*size*/ )
{
  try {
    CUDA_SAFE_CALL( cudaFree( ptr ) );
  } catch(...) {}
}

void * CudaMallocAllocator::reallocate(void * old_ptr, size_t old_size, size_t new_size)
{
  void * ptr = old_ptr;
  if (old_size != new_size) {
    ptr = allocate( new_size );
    size_t copy_size = old_size < new_size ? old_size : new_size;

    CUDA_SAFE_CALL( cudaMemcpy( ptr , old_ptr , copy_size , cudaMemcpyDefault ) );

    deallocate( old_ptr, old_size );
  }
  return ptr;
}

/*--------------------------------------------------------------------------*/

void * CudaUVMAllocator::allocate( size_t size )
{
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION )
  void * ptr = NULL;
  CUDA_SAFE_CALL( cudaMallocManaged( &ptr, size, cudaMemAttachGlobal ) );
  return ptr;
#else
  throw_runtime_exception( "CUDA VERSION does not support UVM" );
  return NULL;
#endif
}

void CudaUVMAllocator::deallocate( void * ptr, size_t /*size*/ )
{
  try {
    CUDA_SAFE_CALL( cudaFree( ptr ) );
  } catch(...) {}
}

void * CudaUVMAllocator::reallocate(void * old_ptr, size_t old_size, size_t new_size)
{
  void * ptr = old_ptr;
  if (old_size != new_size) {
    ptr = allocate( new_size );
    size_t copy_size = old_size < new_size ? old_size : new_size;

    CUDA_SAFE_CALL( cudaMemcpy( ptr , old_ptr , copy_size , cudaMemcpyDefault ) );

    deallocate( old_ptr, old_size );
  }
  return ptr;
}

/*--------------------------------------------------------------------------*/

void * CudaHostAllocator::allocate( size_t size )
{
  void * ptr = NULL;
  CUDA_SAFE_CALL( cudaHostAlloc( &ptr , size , cudaHostAllocDefault ) );
  return ptr;
}

void CudaHostAllocator::deallocate( void * ptr, size_t /*size*/ )
{
  try {
    CUDA_SAFE_CALL( cudaFreeHost( ptr ) );
  } catch(...) {}
}

void * CudaHostAllocator::reallocate(void * old_ptr, size_t old_size, size_t new_size)
{
  void * ptr = old_ptr;
  if (old_size != new_size) {
    ptr = allocate( new_size );
    size_t copy_size = old_size < new_size ? old_size : new_size;

    CUDA_SAFE_CALL( cudaMemcpy( ptr , old_ptr , copy_size , cudaMemcpyHostToHost ) );

    deallocate( old_ptr, old_size );
  }
  return ptr;
}

/*--------------------------------------------------------------------------*/

}} // namespace Kokkos::Impl

#endif //KOKKOS_HAVE_CUDA
