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

#ifndef KOKKOS_MEMORYPOOL_CPP
#define KOKKOS_MEMORYPOOL_CPP

// How should errors be handled?  In general, production code should return a
// value indicating failure so the user can decide how the error is handled.
// While experimental, code can abort instead.  If KOKKOS_MEMPOOLLIST_PRINTERR
// is defined, the code will abort with an error message.  Otherwise, the code
// will return with a value indicating failure when possible, or do nothing
// instead.
//#define KOKKOS_MEMPOOLLIST_PRINTERR

//#define KOKKOS_MEMPOOLLIST_PRINT_INFO

//----------------------------------------------------------------------------

#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA )

/* This '.cpp' is being included by the header file
 * to inline these functions for Cuda.
 *
 *  Prefer to implement these functions in a separate
 *  compilation unit.  However, the 'nvcc' linker
 *  has an internal error when attempting separate compilation
 *  (--relocatable-device-code=true)
 *  of Kokkos unit tests.
 */

#define KOKKOS_MEMPOOLLIST_INLINE inline

#else

/*  This '.cpp' file is being separately compiled for the Host */

#include <Kokkos_MemoryPool.hpp>
#include <Kokkos_Atomic.hpp>

#define KOKKOS_MEMPOOLLIST_INLINE /* */

#endif

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {

#if defined(KOKKOS_MEMPOOLLIST_PRINT_INFO) && defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
long MemPoolList::m_count = 0;
#endif

KOKKOS_FUNCTION
KOKKOS_MEMPOOLLIST_INLINE
uint64_t
MemPoolList::acquire_lock( volatile uint64_t * freelist ) const
{
  uint64_t old_head;
  bool locked = false;

  while ( !locked ) {
    old_head = *freelist;

    if ( old_head != FREELIST_LOCK_HEAD ) {
      // In the initial look at the head, the freelist wasn't locked.
      // Attempt to lock the head of list.  If the list was changed (including
      // being locked) between the initial look and now, head will be different
      // than old_head.  This means the lock can't proceed and has to be
      // tried again.
      uint64_t head =
        atomic_compare_exchange( freelist, old_head, uint64_t(FREELIST_LOCK_HEAD) );

      if ( head == old_head ) locked = true;
    }
  }

  return old_head;
}

KOKKOS_FUNCTION
KOKKOS_MEMPOOLLIST_INLINE
void
MemPoolList::release_lock( volatile uint64_t * freelist, uint64_t new_head ) const
{
  // This function is only intended to be called if acquire_lock() has already
  // been called to acquire a lock on freelist.  Thus, we know that the value
  // pointed to by freelist is FREELIST_LOCK_HEAD.
#ifdef KOKKOS_MEMPOOLLIST_PRINTERR
  uint64_t head =
#endif
    atomic_compare_exchange( freelist, uint64_t(FREELIST_LOCK_HEAD), new_head );

#ifdef KOKKOS_MEMPOOLLIST_PRINTERR
  if ( head != FREELIST_LOCK_HEAD ) {
    // We shouldn't get here, but this check is here for sanity.
    printf( "\n** MemoryPool::allocate() UNLOCK_ERROR(0x%llx) **\n",
            reinterpret_cast<uint64_t>( freelist ) );
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    fflush( stdout );
#endif
    Kokkos::abort( "" );
  }
#endif
}

KOKKOS_FUNCTION
KOKKOS_MEMPOOLLIST_INLINE
void *
MemPoolList::refill_freelist( size_t l_exp ) const
{
  void * p = 0;
  volatile uint64_t * l_exp_freelist = m_freelist + l_exp;

  // The l_exp freelist is empty. Grab a lock on the freelist.
  uint64_t l_exp_old_head = acquire_lock( l_exp_freelist );

  uint64_t l_exp_old_head_off = get_head_offset( l_exp_old_head );

  if ( l_exp_old_head_off != FREELIST_END ) {
    // Another thread put some more entries on the freelist between when
    // this thread saw it empty and acquired the lock.  Just return an entry.
    uint64_t l_exp_old_head_tag = get_head_tag( l_exp_old_head );
    uint64_t new_head_tag = increment_tag( l_exp_old_head_tag );
    uint64_t new_head_off = *reinterpret_cast<uint64_t *>( m_data + l_exp_old_head_off );
    uint64_t new_head = create_head( new_head_off, new_head_tag );

    // Release the lock, replacing the head with the next entry on the list.
    release_lock( l_exp_freelist, new_head );

    // Set the chunk to return.
    p = m_data + l_exp_old_head_off;
  }
  else {
    // The l_exp freelist is empty.

    size_t l = l_exp + 1;
    bool done = false;

    while ( !done ) {
      // Find the next freelist that is either locked or not empty.  A locked
      // freelist will probably have memory available when the lock is
      // released.
      while ( m_chunk_size[l] > 0 &&
              get_head_offset( m_freelist[l] ) == FREELIST_END ) ++l;

      if ( m_chunk_size[l] == 0 ) {
        // We got to the end of the list of freelists without finding any
        // available memory which means the pool is empty.  Release the lock
        // on the l_exp freelist.
        release_lock( l_exp_freelist, l_exp_old_head );

        // Exit out of the loop.
        done = true;
      }
      else {
        volatile uint64_t * l_freelist = m_freelist + l;

        // Grab a lock on the l freelist.
        uint64_t l_old_head = acquire_lock( l_freelist );
        uint64_t l_old_head_off = get_head_offset( l_old_head );

        if ( l_old_head_off != FREELIST_END ) {
          // The l freelist has chunks.  Grab one to divide.

          // Create a new head for the l_freelist by using the second entry
          // in the list and incrementing the current tag.
          uint64_t l_old_head_tag = get_head_tag( l_old_head );
          uint64_t new_head_tag = increment_tag( l_old_head_tag );
          uint64_t new_head_off =
            *reinterpret_cast<volatile uint64_t *>( m_data + l_old_head_off );
          uint64_t new_head = create_head( new_head_off, new_head_tag );

          // Release the lock on the l freelist.
          release_lock( l_freelist, new_head );

          // Subdivide the chunk into smaller chunks.  The first chunk will
          // be returned to satisfy the allocaiton request.  The remainder
          // of the chunks will be inserted onto the appropriate freelist.
          size_t num_chunks = m_chunk_size[l] / m_chunk_size[l_exp];

          // Link the chunks following the first chunk to form a list.
          uint64_t lp_head = l_old_head_off + m_chunk_size[l_exp];
          uint64_t lp_tail = l_old_head_off + (num_chunks - 1) * m_chunk_size[l_exp];

          for ( uint64_t offset = lp_head; offset < lp_tail;
                offset += m_chunk_size[l_exp] )
          {
            *reinterpret_cast<uint64_t *>( m_data + offset ) =
              offset + m_chunk_size[l_exp];
          }

          // Set the tail entry to be the end of the list.
          *reinterpret_cast<volatile uint64_t *>( m_data + lp_tail ) = FREELIST_END;

          memory_fence();

          // Create a new head for the l_exp_freelist.
          new_head = create_head( lp_head, get_head_tag( l_exp_old_head ) );

          // This thread already has the lock on the l_exp freelist, so just
          // release the lock placing the divided memory on the list.
          release_lock( l_exp_freelist, new_head );

          // Set the chunk to return.
          p = m_data + l_old_head_off;
          done = true;
        }
        else {
          // Release the lock on the l freelist.  Put the old head back on.
          release_lock( l_freelist, l_old_head );
        }
      }
    }
  }

  return p;
}

KOKKOS_FUNCTION
KOKKOS_MEMPOOLLIST_INLINE
void *
MemPoolList::allocate( size_t alloc_size ) const
{
  void * p = 0;

  // Find the first freelist whose chunk size is big enough for allocation.
  size_t l_exp = 0;
  while ( m_chunk_size[l_exp] > 0 && alloc_size > m_chunk_size[l_exp] ) ++l_exp;

#ifdef KOKKOS_MEMPOOLLIST_PRINTERR
  if ( m_chunk_size[l_exp] == 0 ) {
    Kokkos::abort( "\n** MemoryPool::allocate() REQUESTED_SIZE_TOO_LARGE **\n" );
  }
#endif

  // Do a fast fail test for an empty list.  This checks for l_exp and all
  // higher freelists being empty.
  size_t l = l_exp;
  while ( m_chunk_size[l] > 0 &&
          get_head_offset( m_freelist[l] ) == FREELIST_END ) ++l;

  if ( m_chunk_size[l] != 0 ) {
    // Try to grab a chunk from the l_exp list.
    volatile uint64_t * l_exp_freelist = m_freelist + l_exp;

    bool done = false;

    while ( !done ) {
      uint64_t old_head = *l_exp_freelist;
      uint64_t old_head_off = get_head_offset( old_head );

      if ( old_head_off == FREELIST_END ) {
        // The list is empty.  Try to refill it and grab a chunk.
        p = refill_freelist(l_exp);

        done = true;
      }
      else if ( old_head_off != FREELIST_LOCK ) {
        // The freelist wasn't empty or locked, so try to pop off the head.
        uint64_t old_head_tag = get_head_tag( old_head );

        // Increment the tag by 1, wrapping around to 0 after 2^32-1.
        uint64_t new_head_tag = increment_tag( old_head_tag );
        uint64_t new_head_off = *reinterpret_cast<uint64_t *>( m_data + old_head_off );
        uint64_t new_head = create_head( new_head_off, new_head_tag );

        // Attempt to pull off the head of the list and put the next entry in
        // its place.  If the list was changed
        // (including being locked) between the initial look and now, head will
        // be different than old_head.  This means the insert can't proceed and
        // has to be tried again.
        uint64_t head = atomic_compare_exchange( l_exp_freelist, old_head, new_head );

        if ( head == old_head ) {
          done = true;
          p = m_data + old_head_off;
        }
      }
    }
  }

#ifdef KOKKOS_MEMPOOLLIST_PRINT_INFO
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
  long val = p == 0 ?
             *reinterpret_cast<volatile long *>( &m_count ) :
             Kokkos::atomic_fetch_add( &m_count, 1 );

  printf( "  allocate(): %6ld   size: %6lu    l: %2lu  %2lu   0x%llx\n", val,
          alloc_size, l_exp, l, reinterpret_cast<uint64_t>( p ) );
  fflush( stdout );
#else
  printf( "  allocate()   size: %6lu    l: %2lu  %2lu   0x%lx\n", alloc_size,
          l_exp, l, reinterpret_cast<uint64_t>( p ) );
#endif
#endif

#ifdef KOKKOS_MEMPOOLLIST_PRINTERR
  if ( p == 0 ) {
    printf( "** MemoryPool::allocate() NO_CHUNKS_BIG_ENOUGH **\n" );
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    fflush( stdout );
#endif
  }
#endif

  return p;
}

KOKKOS_FUNCTION
KOKKOS_MEMPOOLLIST_INLINE
void
MemPoolList::deallocate( void * alloc_ptr, size_t alloc_size ) const
{
#ifdef KOKKOS_MEMPOOLLIST_PRINTERR
  // Verify that the pointer is controlled by this pool.
  {
    char * ap = static_cast<char *>( alloc_ptr );

    if ( ap < m_data || ap + alloc_size > m_data + m_data_size ) {
      printf( "\n** MemoryPool::deallocate() ADDRESS_OUT_OF_RANGE(0x%llx) **\n",
              reinterpret_cast<uint64_t>( alloc_ptr ) );
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
      fflush( stdout );
#endif
      Kokkos::abort( "" );
    }
  }
#endif

  // Determine which freelist to place deallocated memory on.
  size_t l = 0;
  while ( m_chunk_size[l] > 0 && alloc_size > m_chunk_size[l] ) ++l;

#ifdef KOKKOS_MEMPOOLLIST_PRINTERR
  if ( m_chunk_size[l] == 0 ) {
    printf( "\n** MemoryPool::deallocate() CHUNK_TOO_LARGE(%lu) **\n", alloc_size );
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    fflush( stdout );
#endif
    Kokkos::abort( "" );
  }
#endif

  uint64_t offset = static_cast<char *>( alloc_ptr ) - m_data;

  // Insert a single chunk at the head of the freelist.
  volatile uint64_t * freelist = m_freelist + l;

  bool inserted = false;

  while ( !inserted ) {
    uint64_t old_head = *freelist;

    if ( old_head != FREELIST_LOCK_HEAD ) {
      // In the initial look at the head, the freelist wasn't locked.

      uint64_t old_head_off = get_head_offset(old_head);
      uint64_t old_head_tag = get_head_tag(old_head);
      uint64_t new_head = create_head( offset, old_head_tag );

      // Proactively point the new head to the old head assuming a successful
      // insertion into the list.
      *reinterpret_cast<volatile uint64_t *>( alloc_ptr ) = old_head_off;

      memory_fence();

      // Attempt to insert at head of list.  If the list was changed
      // (including being locked) between the initial look and now, head will
      // be different than old_head.  This means the insert can't proceed and
      // has to be tried again.
      uint64_t head = atomic_compare_exchange( freelist, old_head, new_head );

      if ( head == old_head ) inserted = true;
    }
  }

#ifdef KOKKOS_MEMPOOLLIST_PRINT_INFO
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
  long val = Kokkos::atomic_fetch_add( &m_count, -1 ) - 1;
  printf( "deallocate(): %6ld   size: %6lu    l: %2lu       0x%llx\n", val,
          alloc_size, l, reinterpret_cast<uint64_t>( alloc_ptr ) );
  fflush( stdout );
#else
  printf( "deallocate()   size: %6lu    l: %2lu       0x%lx\n", alloc_size, l,
          reinterpret_cast<uint64_t>( alloc_ptr ) );
#endif
#endif
}


} // namespace Impl
} // namespace Experimental
} // namespace Kokkos

#undef KOKKOS_MEMPOOLLIST_INLINE

#ifdef KOKKOS_MEMPOOLLIST_PRINTERR
#undef KOKKOS_MEMPOOLLIST_PRINTERR
#endif

#ifdef KOKKOS_MEMPOOLLIST_PRINT_INFO
#undef KOKKOS_MEMPOOLLIST_PRINT_INFO
#endif

#endif /* #ifndef KOKKOS_MEMORYPOOL_CPP */
