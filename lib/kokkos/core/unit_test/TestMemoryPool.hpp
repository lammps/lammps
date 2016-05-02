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

#include <impl/Kokkos_Timer.hpp>

//#define TESTMEMORYPOOL_PRINT
//#define TESTMEMORYPOOL_PRINT_STATUS

namespace TestMemoryPool {

struct pointer_obj {
  uint64_t * ptr;
};

template < typename PointerView, typename MemorySpace >
struct allocate_memory {
  typedef typename PointerView::execution_space  execution_space;
  typedef typename execution_space::size_type    size_type;

  enum { STRIDE = 32 };

  PointerView m_pointers;
  size_t m_num_ptrs;
  size_t m_chunk_size;
  MemorySpace m_space;

  allocate_memory( PointerView & ptrs, size_t nptrs,
                   size_t cs, MemorySpace & sp )
    : m_pointers( ptrs ), m_num_ptrs( nptrs ),
      m_chunk_size( cs ), m_space( sp )
  {
    // Initialize the view with the out degree of each vertex.
    Kokkos::parallel_for( m_num_ptrs * STRIDE , *this );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i ) const
  {
    if ( i % STRIDE == 0 ) {
      m_pointers[i / STRIDE].ptr =
        static_cast< uint64_t * >( m_space.allocate( m_chunk_size ) );
    }
  }
};

template < typename PointerView >
struct fill_memory {
  typedef typename PointerView::execution_space  execution_space;
  typedef typename execution_space::size_type    size_type;

  enum { STRIDE = 32 };

  PointerView m_pointers;
  size_t m_num_ptrs;

  fill_memory( PointerView & ptrs, size_t nptrs )
    : m_pointers( ptrs ), m_num_ptrs( nptrs )
  {
    // Initialize the view with the out degree of each vertex.
    Kokkos::parallel_for( m_num_ptrs * STRIDE , *this );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i ) const
  {
    if ( i % STRIDE == 0 ) {
      *m_pointers[i / STRIDE].ptr = i / STRIDE ;
    }
  }
};

template < typename PointerView >
struct sum_memory {
  typedef typename PointerView::execution_space  execution_space;
  typedef typename execution_space::size_type    size_type;
  typedef uint64_t                               value_type;

  enum { STRIDE = 32 };

  PointerView m_pointers;
  size_t m_num_ptrs;
  uint64_t & result;

  sum_memory( PointerView & ptrs, size_t nptrs, uint64_t & res )
    : m_pointers( ptrs ), m_num_ptrs( nptrs ), result( res )
  {
    // Initialize the view with the out degree of each vertex.
    Kokkos::parallel_reduce( m_num_ptrs * STRIDE , *this, result );
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

template < typename PointerView, typename MemorySpace >
struct deallocate_memory {
  typedef typename PointerView::execution_space  execution_space;
  typedef typename execution_space::size_type    size_type;

  enum { STRIDE = 32 };

  PointerView m_pointers;
  size_t m_num_ptrs;
  size_t m_chunk_size;
  MemorySpace m_space;

  deallocate_memory( PointerView & ptrs, size_t nptrs,
                     size_t cs, MemorySpace & sp )
    : m_pointers( ptrs ), m_num_ptrs( nptrs ), m_chunk_size( cs ), m_space( sp )
  {
    // Initialize the view with the out degree of each vertex.
    Kokkos::parallel_for( m_num_ptrs * STRIDE , *this );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i ) const
  {
    if ( i % STRIDE == 0 ) {
      m_space.deallocate( m_pointers[i / STRIDE].ptr, m_chunk_size );
    }
  }
};

template < typename ExecutionSpace, typename MemorySpace >
struct allocate_deallocate_memory {
  typedef ExecutionSpace                       execution_space;
  typedef typename execution_space::size_type  size_type;

  enum { STRIDE = 32 };

  size_t m_num_max_chunks;
  size_t m_max_chunk_size;
  size_t m_min_chunk_size;
  size_t m_chunk_spacing;
  MemorySpace m_space;

  allocate_deallocate_memory( size_t nmc, size_t max_cs,
                              size_t min_cs, size_t cs, MemorySpace & sp )
    : m_num_max_chunks( nmc ), m_max_chunk_size( max_cs ),
      m_min_chunk_size( min_cs ), m_chunk_spacing( cs ), m_space( sp )
  {
    Kokkos::parallel_for( m_num_max_chunks * STRIDE, *this );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i ) const
  {
    if ( i % STRIDE == 0 ) {
      for ( size_t j = m_max_chunk_size; j >= m_min_chunk_size; j /= m_chunk_spacing ) {
        for ( size_t k = 0; k < 10; ++k ) {
          void * mem = m_space.allocate( j );
          m_space.deallocate( mem, j );
        }
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
template < class ExecSpace, class MemorySpace = typename ExecSpace::memory_space >
bool test_mempool( size_t chunk_size, size_t total_size )
{
  typedef Kokkos::View< pointer_obj *, ExecSpace >       pointer_view;
  typedef Kokkos::Experimental::MemoryPool< MemorySpace , ExecSpace >
   pool_memory_space;

  uint64_t result;
  size_t num_chunks = total_size / chunk_size;
  bool return_val = true;

  pointer_view pointers( "pointers", num_chunks );

#ifdef TESTMEMORYPOOL_PRINT
  std::cout << std::setw( SHIFTW ) << "chunk_size: " << std::setw( 12 )
            << chunk_size << std::endl
            << std::setw( SHIFTW ) << "total_size: " << std::setw( 12 )
            << total_size << std::endl
            << std::setw( SHIFTW ) << "num_chunks: " << std::setw( 12 )
            << num_chunks << std::endl;

  double elapsed_time = 0;
  Kokkos::Impl::Timer timer;
#endif

  pool_memory_space m_space( MemorySpace(), chunk_size, total_size );

#ifdef TESTMEMORYPOOL_PRINT
  ExecSpace::fence();
  elapsed_time = timer.seconds();
  print_results( "initialize mempool: ", elapsed_time );
#ifdef TESTMEMORYPOOL_PRINT_STATUS
  m_space.print_status();
#endif
  timer.reset();
#endif

  // Tests:
  //   test for correct behvior when out of memory
  //   test for correct behvior when interleaving allocate() and deallocate()

  {
    allocate_memory< pointer_view, pool_memory_space >
      am( pointers, num_chunks, chunk_size, m_space );
  }

#ifdef TESTMEMORYPOOL_PRINT
  ExecSpace::fence();
  elapsed_time = timer.seconds();
  print_results( "allocate chunks: ", elapsed_time );
#ifdef TESTMEMORYPOOL_PRINT_STATUS
  m_space.print_status();
#endif
  timer.reset();
#endif

  {
    fill_memory< pointer_view > fm( pointers, num_chunks );
  }

#ifdef TESTMEMORYPOOL_PRINT
  ExecSpace::fence();
  elapsed_time = timer.seconds();
  print_results( "fill chunks: ", elapsed_time );
  timer.reset();
#endif

  {
    sum_memory< pointer_view > sm( pointers, num_chunks, result );
  }

#ifdef TESTMEMORYPOOL_PRINT
  ExecSpace::fence();
  elapsed_time = timer.seconds();
  print_results( "sum chunks: ", 10, elapsed_time, result );
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
      dm( pointers, num_chunks, chunk_size, m_space );
  }

#ifdef TESTMEMORYPOOL_PRINT
  ExecSpace::fence();
  elapsed_time = timer.seconds();
  print_results( "deallocate chunks: ", elapsed_time );
#ifdef TESTMEMORYPOOL_PRINT_STATUS
  m_space.print_status();
#endif
  timer.reset();
#endif

  {
    allocate_memory< pointer_view, pool_memory_space >
      am( pointers, num_chunks, chunk_size, m_space );
  }

#ifdef TESTMEMORYPOOL_PRINT
  ExecSpace::fence();
  elapsed_time = timer.seconds();
  print_results( "allocate chunks: ", elapsed_time );
#ifdef TESTMEMORYPOOL_PRINT_STATUS
  m_space.print_status();
#endif
  timer.reset();
#endif

  {
    fill_memory< pointer_view > fm( pointers, num_chunks );
  }

#ifdef TESTMEMORYPOOL_PRINT
  ExecSpace::fence();
  elapsed_time = timer.seconds();
  print_results( "fill chunks: ", elapsed_time );
  timer.reset();
#endif

  {
    sum_memory< pointer_view > sm( pointers, num_chunks, result );
  }

#ifdef TESTMEMORYPOOL_PRINT
  ExecSpace::fence();
  elapsed_time = timer.seconds();
  print_results( "sum chunks: ", 10, elapsed_time, result );
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
      dm( pointers, num_chunks, chunk_size, m_space );
  }

#ifdef TESTMEMORYPOOL_PRINT
  ExecSpace::fence();
  elapsed_time = timer.seconds();
  print_results( "deallocate chunks: ", elapsed_time );
#ifdef TESTMEMORYPOOL_PRINT_STATUS
  m_space.print_status();
#endif
#endif

  return return_val;
}

// This test makes allocation requests for multiple sizes and interleaves
// allocation and deallocation.
template < class ExecSpace, class MemorySpace = typename ExecSpace::memory_space >
void test_mempool2( size_t chunk_size, size_t total_size )
{
  typedef Kokkos::Experimental::MemoryPool< MemorySpace , ExecSpace >
    pool_memory_space;

  size_t num_chunk_sizes = 4;
  size_t chunk_spacing = 4;

#ifdef TESTMEMORYPOOL_PRINT
  double elapsed_time = 0;
  Kokkos::Impl::Timer timer;
#endif

  pool_memory_space m_space( MemorySpace(), chunk_size, total_size,
                             num_chunk_sizes, chunk_spacing );

#ifdef TESTMEMORYPOOL_PRINT
  ExecSpace::fence();
  elapsed_time = timer.seconds();
  print_results( "initialize mempool: ", elapsed_time );
#endif

  chunk_size = m_space.get_min_chunk_size();
  total_size = m_space.get_mem_size();

  // Get the chunk size for the largest possible chunk.
  //   max_chunk_size =
  //     chunk_size * (MEMPOOL_CHUNK_SPACING ^ (MEMPOOL_NUM_CHUNK_SIZES - 1))
  size_t max_chunk_size = chunk_size;
  for (size_t i = 1; i < num_chunk_sizes; ++i) {
    max_chunk_size *= chunk_spacing;
  }

  size_t num_max_chunks = total_size / ( max_chunk_size * num_chunk_sizes );

#ifdef TESTMEMORYPOOL_PRINT_STATUS
  m_space.print_status();
#endif

#ifdef TESTMEMORYPOOL_PRINT
  timer.reset();
#endif

  {
    allocate_deallocate_memory< ExecSpace, pool_memory_space >
      am( num_max_chunks, max_chunk_size, chunk_size, chunk_spacing, m_space );
  }

#ifdef TESTMEMORYPOOL_PRINT
  ExecSpace::fence();
  elapsed_time = timer.seconds();
  print_results( "allocate / deallocate: ", elapsed_time );
#ifdef TESTMEMORYPOOL_PRINT_STATUS
  m_space.print_status();
#endif
#endif
}

}

#ifdef TESTMEMORYPOOL_PRINT
#undef TESTMEMORYPOOL_PRINT
#endif

#ifdef TESTMEMORYPOOL_PRINT_STATUS
#undef TESTMEMORYPOOL_PRINT_STATUS
#endif

#endif
