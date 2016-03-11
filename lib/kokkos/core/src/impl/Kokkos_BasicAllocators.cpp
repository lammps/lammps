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

#include <Kokkos_HostSpace.hpp>

#if ! defined( KOKKOS_USING_EXPERIMENTAL_VIEW )

#include <impl/Kokkos_BasicAllocators.hpp>
#include <impl/Kokkos_Error.hpp>


#include <stdint.h>    // uintptr_t
#include <cstdlib>     // for malloc, realloc, and free
#include <cstring>     // for memcpy

#if defined(KOKKOS_POSIX_MEMALIGN_AVAILABLE)
#include <sys/mman.h>  // for mmap, munmap, MAP_ANON, etc
#include <unistd.h>    // for sysconf, _SC_PAGE_SIZE, _SC_PHYS_PAGES
#endif

#include <sstream>

namespace Kokkos { namespace Impl {

/*--------------------------------------------------------------------------*/

void* MallocAllocator::allocate( size_t size )
{
  void * ptr = NULL;
  if (size) {
    ptr = malloc(size);

    if (!ptr)
    {
      std::ostringstream msg ;
      msg << name() << ": allocate(" << size << ") FAILED";
      throw_runtime_exception( msg.str() );
    }
  }
  return ptr;
}

void MallocAllocator::deallocate( void * ptr, size_t /*size*/ )
{
  if (ptr) {
    free(ptr);
  }
}

void * MallocAllocator::reallocate(void * old_ptr, size_t /*old_size*/, size_t new_size)
{
  void * ptr = realloc(old_ptr, new_size);

  if (new_size > 0u && ptr == NULL) {
    throw_runtime_exception("Error: Malloc Allocator could not reallocate memory");
  }
  return ptr;
}

/*--------------------------------------------------------------------------*/

namespace {

void * raw_aligned_allocate( size_t size, size_t alignment )
{
  void * ptr = NULL;
  if ( size ) {
#if defined( __INTEL_COMPILER ) && !defined ( KOKKOS_HAVE_CUDA )
    ptr = _mm_malloc( size , alignment );

#elif defined(KOKKOS_POSIX_MEMALIGN_AVAILABLE)

    posix_memalign( & ptr, alignment , size );

#else
    // Over-allocate to and round up to guarantee proper alignment.
    size_t size_padded = size + alignment + sizeof(void *);
    void * alloc_ptr = malloc( size_padded );

    if (alloc_ptr) {
      uintptr_t address = reinterpret_cast<uintptr_t>(alloc_ptr);
      // offset enough to record the alloc_ptr
      address += sizeof(void *);
      uintptr_t rem = address % alignment;
      uintptr_t offset = rem ? (alignment - rem) : 0u;
      address += offset;
      ptr = reinterpret_cast<void *>(address);
      // record the alloc'd pointer
      address -= sizeof(void *);
      *reinterpret_cast<void **>(address) = alloc_ptr;
    }
#endif
  }
  return ptr;
}

void raw_aligned_deallocate( void * ptr, size_t /*size*/ )
{
  if ( ptr ) {
#if defined( __INTEL_COMPILER ) && !defined ( KOKKOS_HAVE_CUDA )
    _mm_free( ptr );

#elif defined(KOKKOS_POSIX_MEMALIGN_AVAILABLE)
    free( ptr );
#else
    // get the alloc'd pointer
    void * alloc_ptr = *(reinterpret_cast<void **>(ptr) -1);
    free( alloc_ptr );
#endif
  }

}

}

void* AlignedAllocator::allocate( size_t size )
{
  void * ptr = 0 ;

  if ( size ) {
    ptr = raw_aligned_allocate(size, MEMORY_ALIGNMENT);

    if (!ptr)
    {
      std::ostringstream msg ;
      msg << name() << ": allocate(" << size << ") FAILED";
      throw_runtime_exception( msg.str() );
    }
  }
  return ptr;
}

void AlignedAllocator::deallocate( void * ptr, size_t size )
{
  raw_aligned_deallocate( ptr, size);
}

void * AlignedAllocator::reallocate(void * old_ptr, size_t old_size, size_t new_size)
{
  void * ptr = old_ptr;;

  if (old_size < new_size) {
    ptr = allocate( new_size );

    memcpy(ptr, old_ptr, old_size );

    deallocate( old_ptr, old_size );
  }

  return ptr;
}

/*--------------------------------------------------------------------------*/

// mmap flags for private anonymous memory allocation
#if defined( MAP_ANONYMOUS ) && defined( MAP_PRIVATE )
  #define MMAP_FLAGS (MAP_PRIVATE | MAP_ANONYMOUS)
#elif defined( MAP_ANON) && defined( MAP_PRIVATE )
  #define MMAP_FLAGS (MAP_PRIVATE | MAP_ANON)
#else
  #define NO_MMAP
#endif

// huge page tables
#if !defined( NO_MMAP )
  #if defined( MAP_HUGETLB )
    #define MMAP_FLAGS_HUGE (MMAP_FLAGS | MAP_HUGETLB )
  #elif defined( MMAP_FLAGS )
    #define MMAP_FLAGS_HUGE MMAP_FLAGS
  #endif
  // threshold to use huge pages
  #define MMAP_USE_HUGE_PAGES (1u << 27)
#endif

// read write access to private memory
#if !defined( NO_MMAP )
  #define MMAP_PROTECTION (PROT_READ | PROT_WRITE)
#endif


void* PageAlignedAllocator::allocate( size_t size )
{
  void *ptr = NULL;
  if (size) {
#if !defined NO_MMAP
    if ( size < MMAP_USE_HUGE_PAGES ) {
      ptr = mmap( NULL, size, MMAP_PROTECTION, MMAP_FLAGS, -1 /*file descriptor*/, 0 /*offset*/);
    } else {
      ptr = mmap( NULL, size, MMAP_PROTECTION, MMAP_FLAGS_HUGE, -1 /*file descriptor*/, 0 /*offset*/);
    }
    if (ptr == MAP_FAILED) {
      ptr = NULL;
    }
#else
    static const size_t page_size = 4096; // TODO: read in from sysconf( _SC_PAGE_SIZE )

    ptr = raw_aligned_allocate( size, page_size);
#endif
    if (!ptr)
    {
      std::ostringstream msg ;
      msg << name() << ": allocate(" << size << ") FAILED";
      throw_runtime_exception( msg.str() );
    }
  }
  return ptr;
}

void PageAlignedAllocator::deallocate( void * ptr, size_t size )
{
#if !defined( NO_MMAP )
  munmap(ptr, size);
#else
  raw_aligned_deallocate(ptr, size);
#endif
}

void * PageAlignedAllocator::reallocate(void * old_ptr, size_t old_size, size_t new_size)
{
  void * ptr = NULL;
#if defined( NO_MMAP ) || defined( __APPLE__ ) || defined( __CYGWIN__ )

  if (old_size != new_size) {
    ptr = allocate( new_size );

    memcpy(ptr, old_ptr, (old_size < new_size ? old_size : new_size) );

    deallocate( old_ptr, old_size );
  }
  else {
    ptr = old_ptr;
  }
#else
  ptr = mremap( old_ptr, old_size, new_size, MREMAP_MAYMOVE );

  if (ptr == MAP_FAILED) {
    throw_runtime_exception("Error: Page Aligned Allocator could not reallocate memory");
  }
#endif

  return ptr;
}

}} // namespace Kokkos::Impl

#endif /* #if ! defined( KOKKOS_USING_EXPERIMENTAL_VIEW ) */

