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


#ifndef KOKKOS_UNITTEST_MEMPOOL_HPP
#define KOKKOS_UNITTEST_MEMPOOL_HPP

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <algorithm>

#include <impl/Kokkos_Timer.hpp>

//#define TESTMEMORYPOOL_PRINT
//#define TESTMEMORYPOOL_PRINT_STATUS

#define STRIDE 1
#ifdef KOKKOS_ENABLE_CUDA
#define STRIDE_ALLOC 32
#else
#define STRIDE_ALLOC 1
#endif

namespace TestMemoryPool {

struct pointer_obj {
  uint64_t *  ptr;

  KOKKOS_INLINE_FUNCTION
  pointer_obj() : ptr( 0 ) {}
};

struct pointer_obj2 {
  void *  ptr;
  size_t  size;

  KOKKOS_INLINE_FUNCTION
  pointer_obj2() : ptr( 0 ), size( 0 ) {}
};

template < typename PointerView, typename Allocator >
struct allocate_memory {
  typedef typename PointerView::execution_space  execution_space;
  typedef typename execution_space::size_type    size_type;

  PointerView  m_pointers;
  size_t       m_chunk_size;
  Allocator    m_mempool;

  allocate_memory( PointerView & ptrs, size_t num_ptrs,
                   size_t cs, Allocator & m )
    : m_pointers( ptrs ), m_chunk_size( cs ), m_mempool( m )
  {
    // Initialize the view with the out degree of each vertex.
    Kokkos::parallel_for( num_ptrs * STRIDE_ALLOC, *this );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i ) const
  {
    if ( i % STRIDE_ALLOC == 0 ) {
      m_pointers[i / STRIDE_ALLOC].ptr =
        static_cast< uint64_t * >( m_mempool.allocate( m_chunk_size ) );
    }
  }
};

template < typename PointerView >
struct count_invalid_memory {
  typedef typename PointerView::execution_space  execution_space;
  typedef typename execution_space::size_type    size_type;
  typedef uint64_t                               value_type;

  PointerView  m_pointers;
  uint64_t &   m_result;

  count_invalid_memory( PointerView & ptrs, size_t num_ptrs, uint64_t & res )
    : m_pointers( ptrs ), m_result( res )
  {
    // Initialize the view with the out degree of each vertex.
    Kokkos::parallel_reduce( num_ptrs * STRIDE, *this, m_result );
  }

  KOKKOS_INLINE_FUNCTION
  void init( value_type & v ) const
  { v = 0; }

  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type & dst, volatile value_type const & src ) const
  { dst += src; }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i, value_type & r ) const
  {
    if ( i % STRIDE == 0 ) {
      r += ( m_pointers[i / STRIDE].ptr == 0 );
    }
  }
};

template < typename PointerView >
struct fill_memory {
  typedef typename PointerView::execution_space  execution_space;
  typedef typename execution_space::size_type    size_type;

  PointerView m_pointers;

  fill_memory( PointerView & ptrs, size_t num_ptrs ) : m_pointers( ptrs )
  {
    // Initialize the view with the out degree of each vertex.
    Kokkos::parallel_for( num_ptrs * STRIDE, *this );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i ) const
  {
    if ( i % STRIDE == 0 ) {
      *m_pointers[i / STRIDE].ptr = i / STRIDE;
    }
  }
};

template < typename PointerView >
struct sum_memory {
  typedef typename PointerView::execution_space  execution_space;
  typedef typename execution_space::size_type    size_type;
  typedef uint64_t                               value_type;

  PointerView  m_pointers;
  uint64_t &   m_result;

  sum_memory( PointerView & ptrs, size_t num_ptrs, uint64_t & res )
    : m_pointers( ptrs ), m_result( res )
  {
    // Initialize the view with the out degree of each vertex.
    Kokkos::parallel_reduce( num_ptrs * STRIDE, *this, m_result );
  }

  KOKKOS_INLINE_FUNCTION
  void init( value_type & v ) const
  { v = 0; }

  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type & dst, volatile value_type const & src ) const
  { dst += src; }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i, value_type & r ) const
  {
    if ( i % STRIDE == 0 ) {
      r += *m_pointers[i / STRIDE].ptr;
    }
  }
};

template < typename PointerView, typename Allocator >
struct deallocate_memory {
  typedef typename PointerView::execution_space  execution_space;
  typedef typename execution_space::size_type    size_type;

  PointerView  m_pointers;
  size_t       m_chunk_size;
  Allocator    m_mempool;

  deallocate_memory( PointerView & ptrs, size_t num_ptrs,
                     size_t cs, Allocator & m )
    : m_pointers( ptrs ), m_chunk_size( cs ), m_mempool( m )
  {
    // Initialize the view with the out degree of each vertex.
    Kokkos::parallel_for( num_ptrs * STRIDE, *this );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i ) const
  {
    if ( i % STRIDE == 0 ) {
      m_mempool.deallocate( m_pointers[i / STRIDE].ptr, m_chunk_size );
    }
  }
};

template < typename WorkView, typename PointerView, typename ScalarView,
           typename Allocator >
struct allocate_deallocate_memory {
  typedef typename WorkView::execution_space   execution_space;
  typedef typename execution_space::size_type  size_type;

  WorkView     m_work;
  PointerView  m_pointers;
  ScalarView   m_ptrs_front;
  ScalarView   m_ptrs_back;
  Allocator    m_mempool;

  allocate_deallocate_memory( WorkView & w, size_t work_size, PointerView & p,
                              ScalarView pf, ScalarView pb, Allocator & m )
    : m_work( w ), m_pointers( p ), m_ptrs_front( pf ), m_ptrs_back( pb ),
      m_mempool( m )
  {
    // Initialize the view with the out degree of each vertex.
    Kokkos::parallel_for( work_size * STRIDE_ALLOC, *this );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i ) const
  {
    if ( i % STRIDE_ALLOC == 0 ) {
      unsigned my_work = m_work[i / STRIDE_ALLOC];

      if ( ( my_work & 1 ) == 0 ) {
        // Allocation.
        size_t pos = Kokkos::atomic_fetch_add( &m_ptrs_back(), 1 );
        size_t alloc_size = my_work >> 1;
        m_pointers[pos].ptr = m_mempool.allocate( alloc_size );
        m_pointers[pos].size = alloc_size;
      }
      else {
        // Deallocation.
        size_t pos = Kokkos::atomic_fetch_add( &m_ptrs_front(), 1 );
        m_mempool.deallocate( m_pointers[pos].ptr, m_pointers[pos].size );
      }
    }
  }
};

#define PRECISION 6
#define SHIFTW 24
#define SHIFTW2 12

template < typename F >
void print_results( const std::string & text, F elapsed_time )
{
  std::cout << std::setw( SHIFTW ) << text << std::setw( SHIFTW2 )
            << std::fixed << std::setprecision( PRECISION ) << elapsed_time
            << std::endl;
}

template < typename F, typename T >
void print_results( const std::string & text, unsigned long long width,
                    F elapsed_time, T result )
{
  std::cout << std::setw( SHIFTW ) << text << std::setw( SHIFTW2 )
            << std::fixed << std::setprecision( PRECISION ) << elapsed_time
            << "     " << std::setw( width ) << result << std::endl;
}

template < typename F >
void print_results( const std::string & text, unsigned long long width,
                    F elapsed_time, const std::string & result )
{
  std::cout << std::setw( SHIFTW ) << text << std::setw( SHIFTW2 )
            << std::fixed << std::setprecision( PRECISION ) << elapsed_time
            << "     " << std::setw( width ) << result << std::endl;
}

// This test slams allocation and deallocation in a worse than real-world usage
// scenario to see how bad the thread-safety really is by having a loop where
// all threads allocate and a subsequent loop where all threads deallocate.
// All of the allocation requests are for equal-sized chunks that are the base
// chunk size of the memory pool.  It also tests initialization of the memory
// pool and breaking large chunks into smaller chunks to fulfill allocation
// requests.  It verifies that MemoryPool(), allocate(), and deallocate() work
// correctly.
template < class Device >
bool test_mempool( size_t chunk_size, size_t total_size )
{
  typedef typename Device::execution_space                 execution_space;
  typedef typename Device::memory_space                    memory_space;
  typedef Device                                           device_type;
  typedef Kokkos::View< pointer_obj *, device_type >       pointer_view;
  typedef Kokkos::Experimental::MemoryPool< device_type >  pool_memory_space;

  uint64_t result = 0;
  size_t num_chunks = total_size / chunk_size;
  bool return_val = true;

  pointer_view pointers( "pointers", num_chunks );

#ifdef TESTMEMORYPOOL_PRINT
  std::cout << "*** test_mempool() ***" << std::endl
            << std::setw( SHIFTW ) << "chunk_size: " << std::setw( 12 )
            << chunk_size << std::endl
            << std::setw( SHIFTW ) << "total_size: " << std::setw( 12 )
            << total_size << std::endl
            << std::setw( SHIFTW ) << "num_chunks: " << std::setw( 12 )
            << num_chunks << std::endl;

  double elapsed_time = 0;
  Kokkos::Timer timer;
#endif

  pool_memory_space mempool( memory_space(), total_size * 1.2, 20 );

#ifdef TESTMEMORYPOOL_PRINT
  execution_space::fence();
  elapsed_time = timer.seconds();
  print_results( "initialize mempool: ", elapsed_time );
#ifdef TESTMEMORYPOOL_PRINT_STATUS
  mempool.print_status();
#endif
  timer.reset();
#endif

  {
    allocate_memory< pointer_view, pool_memory_space >
      am( pointers, num_chunks, chunk_size, mempool );
  }

#ifdef TESTMEMORYPOOL_PRINT
  execution_space::fence();
  elapsed_time = timer.seconds();
  print_results( "allocate chunks: ", elapsed_time );
#ifdef TESTMEMORYPOOL_PRINT_STATUS
  mempool.print_status();
#endif
  timer.reset();
#endif

  {
    count_invalid_memory< pointer_view > sm( pointers, num_chunks, result );
  }

#ifdef TESTMEMORYPOOL_PRINT
  execution_space::fence();
  elapsed_time = timer.seconds();
  print_results( "invalid chunks: ", 16, elapsed_time, result );
  timer.reset();
#endif

  {
    fill_memory< pointer_view > fm( pointers, num_chunks );
  }

#ifdef TESTMEMORYPOOL_PRINT
  execution_space::fence();
  elapsed_time = timer.seconds();
  print_results( "fill chunks: ", elapsed_time );
  timer.reset();
#endif

  {
    sum_memory< pointer_view > sm( pointers, num_chunks, result );
  }

  execution_space::fence();

#ifdef TESTMEMORYPOOL_PRINT
  elapsed_time = timer.seconds();
  print_results( "sum chunks: ", 16, elapsed_time, result );
#endif

  if ( result != ( num_chunks * ( num_chunks - 1 ) ) / 2 ) {
    std::cerr << "Invalid sum value in memory." << std::endl;
    return_val = false;
  }

#ifdef TESTMEMORYPOOL_PRINT
  timer.reset();
#endif

  {
    deallocate_memory< pointer_view, pool_memory_space >
      dm( pointers, num_chunks, chunk_size, mempool );
  }

#ifdef TESTMEMORYPOOL_PRINT
  execution_space::fence();
  elapsed_time = timer.seconds();
  print_results( "deallocate chunks: ", elapsed_time );
#ifdef TESTMEMORYPOOL_PRINT_STATUS
  mempool.print_status();
#endif
  timer.reset();
#endif

  {
    allocate_memory< pointer_view, pool_memory_space >
      am( pointers, num_chunks, chunk_size, mempool );
  }

#ifdef TESTMEMORYPOOL_PRINT
  execution_space::fence();
  elapsed_time = timer.seconds();
  print_results( "allocate chunks: ", elapsed_time );
#ifdef TESTMEMORYPOOL_PRINT_STATUS
  mempool.print_status();
#endif
  timer.reset();
#endif

  {
    count_invalid_memory< pointer_view > sm( pointers, num_chunks, result );
  }

#ifdef TESTMEMORYPOOL_PRINT
  execution_space::fence();
  elapsed_time = timer.seconds();
  print_results( "invalid chunks: ", 16, elapsed_time, result );
  timer.reset();
#endif

  {
    fill_memory< pointer_view > fm( pointers, num_chunks );
  }

#ifdef TESTMEMORYPOOL_PRINT
  execution_space::fence();
  elapsed_time = timer.seconds();
  print_results( "fill chunks: ", elapsed_time );
  timer.reset();
#endif

  {
    sum_memory< pointer_view > sm( pointers, num_chunks, result );
  }

  execution_space::fence();

#ifdef TESTMEMORYPOOL_PRINT
  elapsed_time = timer.seconds();
  print_results( "sum chunks: ", 16, elapsed_time, result );
#endif

  if ( result != ( num_chunks * ( num_chunks - 1 ) ) / 2 ) {
    std::cerr << "Invalid sum value in memory." << std::endl;
    return_val = false;
  }

#ifdef TESTMEMORYPOOL_PRINT
  timer.reset();
#endif

  {
    deallocate_memory< pointer_view, pool_memory_space >
      dm( pointers, num_chunks, chunk_size, mempool );
  }

#ifdef TESTMEMORYPOOL_PRINT
  execution_space::fence();
  elapsed_time = timer.seconds();
  print_results( "deallocate chunks: ", elapsed_time );
#ifdef TESTMEMORYPOOL_PRINT_STATUS
  mempool.print_status();
#endif
#endif

  return return_val;
}

template < typename T >
T smallest_power2_ge( T val )
{
  // Find the most significant nonzero bit.
  int first_nonzero_bit = Kokkos::Impl::bit_scan_reverse( val );

  // If val is an integral power of 2, ceil( log2( val ) ) is equal to the
  // most significant nonzero bit.  Otherwise, you need to add 1.
  int lg2_size = first_nonzero_bit +
                 !Kokkos::Impl::is_integral_power_of_two( val );

  return T( 1 ) << T( lg2_size );
}

// This test makes allocation requests for multiple sizes and interleaves
// allocation and deallocation.
//
// There are 3 phases.  The first phase does only allocations to build up a
// working state for the allocator.  The second phase interleaves allocations
// and deletions.  The third phase does only deallocations to undo all the
// allocations from the first phase.  By building first to a working state,
// allocations and deallocations can happen in any order for the second phase.
// Each phase performs on multiple chunk sizes.
template < class Device >
void test_mempool2( unsigned base_chunk_size, size_t num_chunk_sizes,
                    size_t phase1_size, size_t phase2_size )
{
#ifdef TESTMEMORYPOOL_PRINT
  typedef typename Device::execution_space                 execution_space;
#endif
  typedef typename Device::memory_space                    memory_space;
  typedef Device                                           device_type;
  typedef Kokkos::View< unsigned *, device_type >          work_view;
  typedef Kokkos::View< size_t, device_type >              scalar_view;
  typedef Kokkos::View< pointer_obj2 *, device_type >      pointer_view;
  typedef Kokkos::Experimental::MemoryPool< device_type >  pool_memory_space;

  enum {
    MIN_CHUNK_SIZE      = 64,
    MIN_BASE_CHUNK_SIZE = MIN_CHUNK_SIZE / 2 + 1
  };

  // Make sure the base chunk size is at least MIN_BASE_CHUNK_SIZE bytes, so
  // all the different chunk sizes translate to different block sizes for the
  // allocator.
  if ( base_chunk_size < MIN_BASE_CHUNK_SIZE ) {
    base_chunk_size = MIN_BASE_CHUNK_SIZE;
  }

  // Get the smallest power of 2 >= the base chunk size.  The size must be
  // >= MIN_CHUNK_SIZE, though.
  unsigned ceil_base_chunk_size = smallest_power2_ge( base_chunk_size );
  if ( ceil_base_chunk_size < MIN_CHUNK_SIZE ) {
    ceil_base_chunk_size = MIN_CHUNK_SIZE;
  }

  // Make sure the phase 1 size is multiples of num_chunk_sizes.
  phase1_size = ( ( phase1_size + num_chunk_sizes - 1 ) / num_chunk_sizes ) *
                num_chunk_sizes;

  // Make sure the phase 2 size is multiples of ( 2 * num_chunk_sizes ).
  phase2_size =
    ( ( phase2_size + 2 * num_chunk_sizes - 1 ) / ( 2 * num_chunk_sizes ) ) *
    2 * num_chunk_sizes;

  // The phase2 size must be <= twice the phase1 size so that deallocations
  // can't happen before allocations.
  if ( phase2_size > 2 * phase1_size ) phase2_size = 2 * phase1_size;

  size_t phase3_size = phase1_size;
  size_t half_phase2_size = phase2_size / 2;

  // Each entry in the work views has the following format.  The least
  // significant bit indicates allocation (0) vs. deallocation (1).  For
  // allocation, the other bits indicate the desired allocation size.

  // Initialize the phase 1 work view with an equal number of allocations for
  // each chunk size.
  work_view phase1_work( "Phase 1 Work", phase1_size );
  typename work_view::HostMirror host_phase1_work =
    create_mirror_view( phase1_work );

  size_t inner_size = phase1_size / num_chunk_sizes;
  unsigned chunk_size = base_chunk_size;

  for ( size_t i = 0; i < num_chunk_sizes; ++i ) {
    for ( size_t j = 0; j < inner_size; ++j ) {
      host_phase1_work[i * inner_size + j] = chunk_size << 1;
    }

    chunk_size *= 2;
  }

  std::random_shuffle( host_phase1_work.ptr_on_device(),
                       host_phase1_work.ptr_on_device() + phase1_size );

  deep_copy( phase1_work, host_phase1_work );

  // Initialize the phase 2 work view with half allocations and half
  // deallocations with an equal number of allocations for each chunk size.
  work_view phase2_work( "Phase 2 Work", phase2_size );
  typename work_view::HostMirror host_phase2_work =
    create_mirror_view( phase2_work );

  inner_size = half_phase2_size / num_chunk_sizes;
  chunk_size = base_chunk_size;

  for ( size_t i = 0; i < num_chunk_sizes; ++i ) {
    for ( size_t j = 0; j < inner_size; ++j ) {
      host_phase2_work[i * inner_size + j] = chunk_size << 1;
    }

    chunk_size *= 2;
  }

  for ( size_t i = half_phase2_size; i < phase2_size; ++i ) {
    host_phase2_work[i] = 1;
  }

  std::random_shuffle( host_phase2_work.ptr_on_device(),
                       host_phase2_work.ptr_on_device() + phase2_size );

  deep_copy( phase2_work, host_phase2_work );

  // Initialize the phase 3 work view with all deallocations.
  work_view phase3_work( "Phase 3 Work", phase3_size );
  typename work_view::HostMirror host_phase3_work =
    create_mirror_view( phase3_work );

  inner_size = phase3_size / num_chunk_sizes;

  for ( size_t i = 0; i < phase3_size; ++i ) host_phase3_work[i] = 1;

  deep_copy( phase3_work, host_phase3_work );

  // Calculate the amount of memory needed for the allocator.  We need to know
  // the number of superblocks required for each chunk size and use that to
  // calculate the amount of memory for each chunk size.
  size_t lg_sb_size = 18;
  size_t sb_size = 1 << lg_sb_size;
  size_t total_size = 0;
  size_t allocs_per_size = phase1_size / num_chunk_sizes +
                           half_phase2_size / num_chunk_sizes;

  chunk_size = ceil_base_chunk_size;
  for ( size_t i = 0; i < num_chunk_sizes; ++i ) {
    size_t my_size = allocs_per_size * chunk_size;
    total_size += ( my_size + sb_size - 1 ) / sb_size * sb_size;
    chunk_size *= 2;
  }

  // Declare the queue to hold the records for allocated memory.  An allocation
  // adds a record to the back of the queue, and a deallocation removes a
  // record from the front of the queue.
  size_t num_allocations = phase1_size + half_phase2_size;
  scalar_view ptrs_front( "Pointers front" );
  scalar_view ptrs_back( "Pointers back" );

  pointer_view pointers( "pointers", num_allocations );

#ifdef TESTMEMORYPOOL_PRINT
  printf( "\n*** test_mempool2() ***\n" );
  printf( "       num_chunk_sizes: %12zu\n", num_chunk_sizes );
  printf( "       base_chunk_size: %12u\n", base_chunk_size );
  printf( "  ceil_base_chunk_size: %12u\n", ceil_base_chunk_size );
  printf( "           phase1_size: %12zu\n", phase1_size );
  printf( "           phase2_size: %12zu\n", phase2_size );
  printf( "           phase3_size: %12zu\n", phase3_size );
  printf( "       allocs_per_size: %12zu\n", allocs_per_size );
  printf( "       num_allocations: %12zu\n", num_allocations );
  printf( "            total_size: %12zu\n", total_size );
  fflush( stdout );

  double elapsed_time = 0;
  Kokkos::Timer timer;
#endif

  pool_memory_space mempool( memory_space(), total_size * 1.2, lg_sb_size );

#ifdef TESTMEMORYPOOL_PRINT
  execution_space::fence();
  elapsed_time = timer.seconds();
  print_results( "initialize mempool: ", elapsed_time );

#ifdef TESTMEMORYPOOL_PRINT_STATUS
  mempool.print_status();
#endif

  timer.reset();
#endif

  {
    allocate_deallocate_memory< work_view, pointer_view, scalar_view,
                                pool_memory_space >
      adm( phase1_work, phase1_size, pointers, ptrs_front, ptrs_back, mempool );
  }

#ifdef TESTMEMORYPOOL_PRINT
  execution_space::fence();
  elapsed_time = timer.seconds();
  print_results( "phase1: ", elapsed_time );

#ifdef TESTMEMORYPOOL_PRINT_STATUS
  mempool.print_status();
#endif

  timer.reset();
#endif

  {
    allocate_deallocate_memory< work_view, pointer_view, scalar_view,
                                pool_memory_space >
      adm( phase2_work, phase2_size, pointers, ptrs_front, ptrs_back, mempool );
  }

#ifdef TESTMEMORYPOOL_PRINT
  execution_space::fence();
  elapsed_time = timer.seconds();
  print_results( "phase2: ", elapsed_time );

#ifdef TESTMEMORYPOOL_PRINT_STATUS
  mempool.print_status();
#endif

  timer.reset();
#endif

  {
    allocate_deallocate_memory< work_view, pointer_view, scalar_view,
                                pool_memory_space >
      adm( phase3_work, phase3_size, pointers, ptrs_front, ptrs_back, mempool );
  }

#ifdef TESTMEMORYPOOL_PRINT
  execution_space::fence();
  elapsed_time = timer.seconds();
  print_results( "phase3: ", elapsed_time );
#ifdef TESTMEMORYPOOL_PRINT_STATUS
  mempool.print_status();
#endif
#endif
}

// Tests for correct behavior when the allocator is out of memory.
template < class Device >
void test_memory_exhaustion()
{
#ifdef TESTMEMORYPOOL_PRINT
  typedef typename Device::execution_space                 execution_space;
#endif
  typedef typename Device::memory_space                    memory_space;
  typedef Device                                           device_type;
  typedef Kokkos::View< pointer_obj *, device_type >       pointer_view;
  typedef Kokkos::Experimental::MemoryPool< device_type >  pool_memory_space;

  // The allocator will have a single superblock, and allocations will all be
  // of the same chunk size.  The allocation loop will attempt to allocate
  // twice the number of chunks as are available in the allocator.  The
  // deallocation loop will only free the successfully allocated chunks.

  size_t chunk_size = 128;
  size_t num_chunks = 128;
  size_t half_num_chunks = num_chunks / 2;
  size_t superblock_size = chunk_size * half_num_chunks;
  size_t lg_superblock_size =
    Kokkos::Impl::integral_power_of_two( superblock_size );

  pointer_view pointers( "pointers", num_chunks );

#ifdef TESTMEMORYPOOL_PRINT
  std::cout << "\n*** test_memory_exhaustion() ***" << std::endl;

  double elapsed_time = 0;
  Kokkos::Timer timer;
#endif

  pool_memory_space mempool( memory_space(), superblock_size,
                             lg_superblock_size );

#ifdef TESTMEMORYPOOL_PRINT
  execution_space::fence();
  elapsed_time = timer.seconds();
  print_results( "initialize mempool: ", elapsed_time );
#ifdef TESTMEMORYPOOL_PRINT_STATUS
  mempool.print_status();
#endif
  timer.reset();
#endif

  {
    allocate_memory< pointer_view, pool_memory_space >
      am( pointers, num_chunks, chunk_size, mempool );
  }

#ifdef TESTMEMORYPOOL_PRINT
  execution_space::fence();
  elapsed_time = timer.seconds();
  print_results( "allocate chunks: ", elapsed_time );
#ifdef TESTMEMORYPOOL_PRINT_STATUS
  mempool.print_status();
#endif
  timer.reset();
#endif

  {
    // In parallel, the allocations that succeeded were not put contiguously
    // into the pointers View.  The whole View can still be looped over and
    // have deallocate called because deallocate will just do nothing for NULL
    // pointers.
    deallocate_memory< pointer_view, pool_memory_space >
      dm( pointers, num_chunks, chunk_size, mempool );
  }

#ifdef TESTMEMORYPOOL_PRINT
  execution_space::fence();
  elapsed_time = timer.seconds();
  print_results( "deallocate chunks: ", elapsed_time );
#ifdef TESTMEMORYPOOL_PRINT_STATUS
  mempool.print_status();
#endif
#endif
}

}

#undef TESTMEMORYPOOL_PRINT
#undef TESTMEMORYPOOL_PRINT_STATUS
#undef STRIDE
#undef STRIDE_ALLOC

#endif
