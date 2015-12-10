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

#ifndef KOKKOS_CORE_HPP
#define KOKKOS_CORE_HPP

//----------------------------------------------------------------------------
// Include the execution space header files for the enabled execution spaces.

#include <Kokkos_Core_fwd.hpp>

#if defined( KOKKOS_HAVE_SERIAL )
#include <Kokkos_Serial.hpp>
#endif

#if defined( KOKKOS_HAVE_OPENMP )
#include <Kokkos_OpenMP.hpp>
#endif

#if defined( KOKKOS_HAVE_PTHREAD )
#include <Kokkos_Threads.hpp>
#endif

#if defined( KOKKOS_HAVE_CUDA )
#include <Kokkos_Cuda.hpp>
#endif

#include <Kokkos_Pair.hpp>
#include <Kokkos_Array.hpp>
#include <Kokkos_View.hpp>
#include <Kokkos_Vectorization.hpp>
#include <Kokkos_Atomic.hpp>
#include <Kokkos_hwloc.hpp>

#ifdef KOKKOS_HAVE_CXX11
#include <Kokkos_Complex.hpp>
#endif


//----------------------------------------------------------------------------

namespace Kokkos {

struct InitArguments {
  int num_threads;
  int num_numa;
  int device_id;

  InitArguments() {
    num_threads = -1;
    num_numa = -1;
    device_id = -1;
  }
};

void initialize(int& narg, char* arg[]);

void initialize(const InitArguments& args = InitArguments());

/** \brief  Finalize the spaces that were initialized via Kokkos::initialize */
void finalize();

/** \brief  Finalize all known execution spaces */
void finalize_all();

void fence();

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {

/* Allocate memory from a memory space.
 * The allocation is tracked in Kokkos memory tracking system, so
 * leaked memory can be identified.
 */
template< class Space = typename Kokkos::DefaultExecutionSpace::memory_space >
inline
void * kokkos_malloc( const std::string & arg_alloc_label
                    , const size_t arg_alloc_size )
{
  typedef typename Space::memory_space MemorySpace ;
  return Impl::SharedAllocationRecord< MemorySpace >::
    allocate_tracked( MemorySpace() , arg_alloc_label , arg_alloc_size );
}

template< class Space = typename Kokkos::DefaultExecutionSpace::memory_space >
inline
void * kokkos_malloc( const size_t arg_alloc_size )
{
  typedef typename Space::memory_space MemorySpace ;
  return Impl::SharedAllocationRecord< MemorySpace >::
    allocate_tracked( MemorySpace() , "no-label" , arg_alloc_size );
}

template< class Space = typename Kokkos::DefaultExecutionSpace::memory_space >
inline
void kokkos_free( void * arg_alloc )
{
  typedef typename Space::memory_space MemorySpace ;
  return Impl::SharedAllocationRecord< MemorySpace >::
    deallocate_tracked( arg_alloc );
}

template< class Space = typename Kokkos::DefaultExecutionSpace::memory_space >
inline
void * kokkos_realloc( void * arg_alloc , const size_t arg_alloc_size )
{
  typedef typename Space::memory_space MemorySpace ;
  return Impl::SharedAllocationRecord< MemorySpace >::
    reallocate_tracked( arg_alloc , arg_alloc_size );
}

} // namespace Experimental
} // namespace Kokkos

#if defined( KOKKOS_USING_EXPERIMENTAL_VIEW )

namespace Kokkos {

using Kokkos::Experimental::kokkos_malloc ;
using Kokkos::Experimental::kokkos_realloc ;
using Kokkos::Experimental::kokkos_free ;

}

#else

namespace Kokkos {

namespace Impl {
// should only by used by kokkos_malloc and kokkos_free
struct MallocHelper
{
  static void increment_ref_count( AllocationTracker const & tracker )
  {
    tracker.increment_ref_count();
  }

  static void decrement_ref_count( AllocationTracker const & tracker )
  {
    tracker.decrement_ref_count();
  }
};
} // namespace Impl

/* Allocate memory from a memory space.
 * The allocation is tracked in Kokkos memory tracking system, so
 * leaked memory can be identified.
 */
template< class Arg = DefaultExecutionSpace>
void* kokkos_malloc(const std::string label, size_t count) {
  if(count == 0) return NULL;
  typedef typename Arg::memory_space MemorySpace;
  Impl::AllocationTracker tracker = MemorySpace::allocate_and_track(label,count);;
  Impl::MallocHelper::increment_ref_count( tracker );
  return tracker.alloc_ptr();
}

template< class Arg = DefaultExecutionSpace>
void* kokkos_malloc(const size_t& count) {
  return kokkos_malloc<Arg>("DefaultLabel",count);
}


/* Free memory from a memory space.
 */
template< class Arg = DefaultExecutionSpace>
void kokkos_free(const void* ptr) {
  typedef typename Arg::memory_space MemorySpace;
  typedef typename MemorySpace::allocator allocator;
  Impl::AllocationTracker tracker = Impl::AllocationTracker::find<allocator>(ptr);
  if (tracker.is_valid()) {
    Impl::MallocHelper::decrement_ref_count( tracker );
  }
}


template< class Arg = DefaultExecutionSpace>
void* kokkos_realloc(const void* old_ptr, size_t size) {
  if(old_ptr == NULL)
    return kokkos_malloc<Arg>(size);

  typedef typename Arg::memory_space MemorySpace;
  typedef typename MemorySpace::allocator allocator;
  Impl::AllocationTracker tracker = Impl::AllocationTracker::find<allocator>(old_ptr);

  tracker.reallocate(size);

  return tracker.alloc_ptr();
}

} // namespace Kokkos

#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif

