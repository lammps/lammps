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
#ifndef KOKKOS_ROCMEXEC_HPP
#define KOKKOS_ROCMEXEC_HPP

#include <algorithm>
#include <typeinfo>

#if defined(__HCC_ACCELERATOR__)
#define printf(...)
#endif

namespace Kokkos {
namespace Impl {

struct ROCmTraits {
// TODO: determine if needed
  enum { WavefrontSize       = 64  /* 64  */ };
  enum { WorkgroupSize       = 256 /* 256  */ };
  enum { WavefrontIndexMask  = 0x003f  /* Mask for wavefrontindex */ };
  enum { WavefrontIndexShift = 6   /* WavefrontSize == 1 << WavefrontShift */ };

  enum { SharedMemoryBanks    = 64      /* GCN */ };
  enum { SharedMemoryCapacity = 0x10000 /* 64k shared / 16k L1 Cache */ };
  enum { SharedMemoryUsage    = 0x04000 /* 64k shared / 16k L1 Cache */ };

  enum { UpperBoundExtentCount    = 4294967295 /* Hard upper bound */ };
#if 0
  KOKKOS_INLINE_FUNCTION static
  ROCmSpace::size_type wavefront_count( ROCmSpace::size_type i )
    { return ( i +  WavefrontIndexMask ) >>  WavefrontIndexShift ; }

  KOKKOS_INLINE_FUNCTION static
  ROCmSpace::size_type wavefront_align( ROCmSpace::size_type i )
    {
      enum { Mask = ~ROCmSpace::size_type(  WavefrontIndexMask ) };
      return ( i +  WavefrontIndexMask ) & Mask ;
    }
#endif
};
size_t rocm_internal_cu_count();
size_t rocm_internal_maximum_workgroup_count();

size_t * rocm_internal_scratch_flags( const size_t size );
size_t * rocm_internal_scratch_space( const size_t size );

// This pointer is the start of dynamic shared memory (LDS).
// Dynamic is at the end of LDS and it's size must be specified
// in a tile_block specification at kernel launch time.
template< typename T >
KOKKOS_INLINE_FUNCTION
T * kokkos_impl_rocm_shared_memory()
//{ return (T*) hc::get_group_segment_base_pointer() ; }
{ return (T*) hc::get_dynamic_group_segment_base_pointer() ; }


}
} // namespace Kokkos
#define ROCM_SPACE_ATOMIC_MASK      0x1FFFF
#define ROCM_SPACE_ATOMIC_XOR_MASK  0x15A39
//int rocm_space_atomic_locks[ROCM_SPACE_ATOMIC_MASK+1];
extern int
   *rocm_space_atomic_locks;

namespace Kokkos {
namespace Impl {
  void init_lock_arrays_rocm_space();

  void* rocm_resize_scratch_space(size_t bytes, bool force_shrink = false);

// TODO: determine if needed
KOKKOS_INLINE_FUNCTION
bool lock_address_rocm_space(void* ptr) {
#if 0
return(Kokkos::Impl::lock_address_host_space(ptr));
#else
  size_t offset = size_t(ptr);
  offset = offset >> 2;
  offset = offset & ROCM_SPACE_ATOMIC_MASK;
  return (0 == hc::atomic_compare_exchange(&rocm_space_atomic_locks[offset],0,1));
#endif
}
KOKKOS_INLINE_FUNCTION
void unlock_address_rocm_space(void* ptr) {
#if 0
Kokkos::Impl::unlock_address_host_space(ptr) ;
#else
  size_t offset = size_t(ptr);
  offset = offset >> 2;
  offset = offset & ROCM_SPACE_ATOMIC_MASK;
  hc::atomic_exchange( &rocm_space_atomic_locks[ offset ], 0);
#endif
}

}
} // namespace Kokkos

namespace Kokkos {
namespace Impl {
//extern 
//KOKKOS_INLINE_FUNCTION
//void init_lock_arrays_rocm_space(); 


}
} // namespace Kokkos
//#if defined(__HCC_ACCELERATOR__)
namespace Kokkos {
namespace Impl {
/*
template< class DriverType>
__global__
static void rocm_parallel_launch_constant_memory()
{
  const DriverType & driver =
    *((const DriverType *) kokkos_impl_rocm_constant_memory_buffer );

  driver();
}

template< class DriverType, unsigned int maxTperB, unsigned int minBperSM >
__global__
__launch_bounds__(maxTperB, minBperSM)
static void rocm_parallel_launch_constant_memory()
{
  const DriverType & driver =
    *((const DriverType *) kokkos_impl_rocm_constant_memory_buffer );

  driver();
}

template< class DriverType>
__global__
static void rocm_parallel_launch_local_memory( const DriverType driver )
{
  driver();
}

template< class DriverType, unsigned int maxTperB, unsigned int minBperSM >
__global__
__launch_bounds__(maxTperB, minBperSM)
static void rocm_parallel_launch_local_memory( const DriverType driver )
{
  driver();
}
*/
template < class DriverType
         , class LaunchBounds = Kokkos::LaunchBounds<> >
struct ROCmParallelLaunch ;

template < class DriverType
         , unsigned int MaxThreadsPerBlock
         , unsigned int MinBlocksPerSM >
struct ROCmParallelLaunch< DriverType
                         , Kokkos::LaunchBounds< MaxThreadsPerBlock
                                               , MinBlocksPerSM >>
{
  inline
  ROCmParallelLaunch( const DriverType & driver
                    , const dim3       & grid
                    , const dim3       & block
                    , const int      shmem )
  {
    if ( grid.x && ( block.x * block.y * block.z ) ) {
      if ( ROCmTraits::SharedMemoryCapacity < shmem ) {
        Kokkos::Impl::throw_runtime_exception( std::string("ROCmParallelLaunch FAILED: shared memory request is too large") );
      }
      DriverType * rocm_memory_buffer = (DriverType *)
                                     rocm_device_allocate(sizeof(DriverType));
      // Copy functor to constant memory on the device
      Kokkos::Impl::DeepCopy<HostSpace,Kokkos::Experimental::ROCmSpace>
              ( rocm_memory_buffer , (void *)&driver , sizeof(DriverType) );

//      KOKKOS_ENSURE_ROCM_LOCK_ARRAYS_ON_DEVICE();

      // Invoke the driver function on the device
      auto ext = hc::extent<3>(grid.z,grid.y,grid.x);
      size_t bx = (grid.x > block.x)? block.x : grid.x;
      size_t by = (grid.y > block.y)? block.y : grid.y;
      size_t bz = (grid.z > block.z)? block.z : grid.z;

      hc::parallel_for_each(ext.tile_with_dynamic(bz,by,bx,shmem), [=](const hc::index<3> & idx) [[hc]]
 
      { rocm_memory_buffer->operator()();
      }).wait();
      rocm_device_free(rocm_memory_buffer);

//#if defined( KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK )
//      ROCM_SAFE_CALL( rocmGetLastError() );
//      Kokkos::ROCm().fence();
//#endif
    }
  }
};
template < class DriverType >
struct ROCmParallelLaunch< DriverType
                         , Kokkos::LaunchBounds<>>
{
  inline
  ROCmParallelLaunch( const DriverType & driver
                    , const dim3       & grid
                    , const dim3       & block
                    , const int          shmem )
  {
    if ( grid.x && ( block.x * block.y * block.z ) ) {
      if ( ROCmTraits::SharedMemoryCapacity < shmem ) {
        Kokkos::Impl::throw_runtime_exception( std::string("ROCmParallelLaunch FAILED: shared memory request is too large") );
      }

      DriverType * rocm_memory_buffer = (DriverType *)
                                     rocm_device_allocate(sizeof(DriverType));
      // Copy functor to constant memory on the device
      Kokkos::Impl::DeepCopy<HostSpace,Kokkos::Experimental::ROCmSpace>
              ( rocm_memory_buffer , (void *)&driver , sizeof(DriverType) );

//      KOKKOS_ENSURE_ROCM_LOCK_ARRAYS_ON_DEVICE();
      // Invoke the driver function on the device
      auto ext = hc::extent<3>(grid.z,grid.y,grid.x);
      size_t bx = (grid.x > block.x)? block.x : grid.x;
      size_t by = (grid.y > block.y)? block.y : grid.y;
      size_t bz = (grid.z > block.z)? block.z : grid.z;
      hc::parallel_for_each(ext.tile_with_dynamic(bz,by,bx,shmem), [=](const hc::index<3> & idx) [[hc]]
 
 
      { rocm_memory_buffer->operator()();
      }).wait();
      rocm_device_free(rocm_memory_buffer);
    }
  }
};
} // namespace Impl
} // namespace Kokkos


#endif /* #ifndef KOKKOS_ROCMEXEC_HPP */
