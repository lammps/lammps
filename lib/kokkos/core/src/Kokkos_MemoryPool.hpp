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

#ifndef KOKKOS_MEMORYPOOL_HPP
#define KOKKOS_MEMORYPOOL_HPP

#include <vector>

#include <Kokkos_Core_fwd.hpp>
#include <impl/Kokkos_Error.hpp>
#include <impl/KokkosExp_SharedAlloc.hpp>
#include <Kokkos_ExecPolicy.hpp>
#include <Kokkos_Atomic.hpp>

// How should errors be handled?  In general, production code should return a
// value indicating failure so the user can decide how the error is handled.
// While experimental, code can abort instead.  If KOKKOS_MEMPOOL_PRINTERR is
// defined, the code will abort with an error message.  Otherwise, the code will
// return with a value indicating failure when possible, or do nothing instead.
//#define KOKKOS_MEMPOOL_PRINTERR

//#define KOKKOS_MEMPOOL_PRINT_INFO

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {

template < typename Space , typename ExecSpace = typename Space::execution_space >
class MemoryPool;

namespace Impl {

#ifdef KOKKOS_MEMPOOL_PRINT_INFO
template < typename MemPool >
struct print_mempool {
  size_t      m_num_chunk_sizes;
  size_t *    m_chunk_size;
  uint64_t *  m_freelist;
  char *      m_data;

  print_mempool( size_t ncs, size_t * cs, uint64_t * f, char * d )
    : m_num_chunk_sizes(ncs), m_chunk_size(cs), m_freelist(f), m_data(d)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator()( size_t i ) const
  {
    if ( i == 0 ) {
      printf( "*** ON DEVICE ***\n");
      printf( "m_chunk_size: 0x%llx\n", reinterpret_cast<uint64_t>( m_chunk_size ) );
      printf( "  m_freelist: 0x%llx\n", reinterpret_cast<uint64_t>( m_freelist ) );
      printf( "      m_data: 0x%llx\n", reinterpret_cast<uint64_t>( m_data ) );
      for ( size_t l = 0; l < m_num_chunk_sizes; ++l ) {
        printf( "%2lu    freelist: %10llu    chunk_size: %6lu\n",
               l, get_head_offset( m_freelist[l] ), m_chunk_size[l] );
      }
      printf( "                              chunk_size: %6lu\n\n",
              m_chunk_size[m_num_chunk_sizes] );
    }
  }

  // This is only redefined here to avoid having to pass a MemPoolList object
  // to the class.
  KOKKOS_INLINE_FUNCTION
  uint64_t get_head_offset(uint64_t head) const
  { return ( head >> MemPool::TAGBITS ) << MemPool::LG_MIN_CHUNKSIZE; }
};
#endif

template < typename MemPool >
struct initialize_mempool {
  char *  m_data;
  size_t  m_chunk_size;
  size_t  m_last_chunk;
  size_t  m_base_offset;

  initialize_mempool( char * d, size_t cs, size_t lc, size_t bo )
    : m_data(d), m_chunk_size(cs), m_last_chunk(lc), m_base_offset(bo)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator()( size_t i ) const
  {
    uint64_t * lp =
      reinterpret_cast<uint64_t *>( m_data + m_base_offset + i * m_chunk_size );

    // All entries in the list point to the next entry except the last which
    // uses a reserved value to indicate the end of the list.  The offset from
    // the base pointer is stored in increments of the minimum chunk size.
    *lp = i < m_last_chunk ?
          m_base_offset + (i + 1) * m_chunk_size :
          MemPool::FREELIST_END;
  }
};

class MemPoolList {
private:

  typedef Impl::SharedAllocationTracker  Tracker;

  template < typename , typename > friend class Kokkos::Experimental::MemoryPool;
  template < typename > friend struct initialize_mempool;
#ifdef KOKKOS_MEMPOOL_PRINT_INFO
  template < typename > friend struct print_mempool;
#endif

  // Define some constants.
  enum {
    // The head of a freelist is a 64 bit unsigned interger.  We divide it
    // into 2 pieces.  The upper (64-TAGBITS) bits is the offset from the base
    // data pointer of the allocator in increments of the minimum chunk size.
    // The lower TAGBITS bits is the tag used to prevent ABA problems.  The
    // largest two values that fit in the offset portion are reserved to
    // represent the end of the freelist and that the freelist is locked.
    //
    // Using 32 bits for both the tag and offset and with a minimum chunk size
    // of 128 bytes, the offset can address 549755813632 bytes (app. 512 GB)
    // of memory.  This should be more than enough to address the whole address
    // space of a GPU or MIC for the foreseeable future.
    TAGBITS            = 32,
    MIN_CHUNKSIZE      = 128,

    TAGBITS_MASK       = ( uint64_t( 1 ) << TAGBITS ) - 1,
    LG_MIN_CHUNKSIZE   = Kokkos::Impl::integral_power_of_two(MIN_CHUNKSIZE),

    // The largest two values of the offset are reserved to indicate the end of a
    // freelist (2^TAGBITS - 2) and that the freelist is locked (2^TAGBITS - 1).
    // They are shifted so they can be compared directly to the result of
    // get_head_offset().
    FREELIST_END       = uint64_t( TAGBITS_MASK - 1 ) << LG_MIN_CHUNKSIZE,
    FREELIST_LOCK      = uint64_t( TAGBITS_MASK ) << LG_MIN_CHUNKSIZE,

    // This is the head value for a locked freelist.  It uses the lock value for
    // the offset and 0 for the tagbits.
    FREELIST_LOCK_HEAD = uint64_t( TAGBITS_MASK ) << TAGBITS
  };

  Tracker   m_track;

  // These three variables are pointers into device memory.
  size_t *    m_chunk_size; // Array of chunk sizes of freelists.
  uint64_t *  m_freelist;   // Array of freelist heads.
  char *      m_data;       // Beginning memory location used for chunks.

  size_t      m_data_size;
  size_t      m_chunk_spacing;

#if defined(KOKKOS_MEMPOOL_PRINT_INFO) && defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
  static long m_count;
#endif

  ~MemPoolList() = default;
  MemPoolList() = default;
  MemPoolList( MemPoolList && ) = default;
  MemPoolList( const MemPoolList & ) = default;
  MemPoolList & operator = ( MemPoolList && ) = default;
  MemPoolList & operator = ( const MemPoolList & ) = default;

  template < typename MemorySpace, typename ExecutionSpace >
  inline
  MemPoolList( const MemorySpace & memspace, const ExecutionSpace &,
               size_t arg_base_chunk_size, size_t arg_total_size,
               size_t num_chunk_sizes, size_t chunk_spacing )
    : m_track(), m_chunk_size(0), m_freelist(0), m_data(0), m_data_size(0),
      m_chunk_spacing(chunk_spacing)
  {
    static_assert( sizeof(size_t) <= sizeof(void*), "" );

    typedef Impl::SharedAllocationRecord< MemorySpace, void >  SharedRecord;
    typedef Kokkos::RangePolicy< ExecutionSpace >              Range;

    size_t base_chunk_size = arg_base_chunk_size;

    // The base chunk size must be at least MIN_CHUNKSIZE bytes as this is the
    // cache-line size for NVIDA GPUs.
    if ( base_chunk_size < MIN_CHUNKSIZE ) {

#ifdef KOKKOS_MEMPOOL_PRINT_INFO
      printf( "** Chunk size must be at least %u bytes.  Setting to %u. **\n",
              MIN_CHUNKSIZE, MIN_CHUNKSIZE);
      fflush( stdout );
#endif

      base_chunk_size = MIN_CHUNKSIZE;
    }

    // The base chunk size must also be a multiple of MIN_CHUNKSIZE bytes for
    // correct memory alignment of the chunks.  If it isn't a multiple of
    // MIN_CHUNKSIZE, set it to the smallest multiple of MIN_CHUNKSIZE
    // greater than the given chunk size.
    if ( base_chunk_size % MIN_CHUNKSIZE != 0 ) {
      size_t old_chunk_size = base_chunk_size;
      base_chunk_size = ( ( old_chunk_size + MIN_CHUNKSIZE - 1 ) / MIN_CHUNKSIZE ) *
                        MIN_CHUNKSIZE;

#ifdef KOKKOS_MEMPOOL_PRINT_INFO
      printf( "** Chunk size must be a multiple of %u bytes.  Given: %lu  Using: %lu. **\n",
              MIN_CHUNKSIZE, old_chunk_size, base_chunk_size);
      fflush( stdout );
#endif

    }

    // Force total_size to be a multiple of base_chunk_size.
    // Preserve the number of chunks originally requested.
    size_t total_size = base_chunk_size *
      ( ( arg_total_size + arg_base_chunk_size - 1 ) / arg_base_chunk_size );

    m_data_size = total_size;

    // Get the chunk size for the largest possible chunk.
    //   max_chunk_size =
    //     base_chunk_size * (m_chunk_spacing ^ (num_chunk_sizes - 1))
    size_t max_chunk_size = base_chunk_size;
    for (size_t i = 1; i < num_chunk_sizes; ++i) {
      max_chunk_size *= m_chunk_spacing;
    }

    // We want each chunk size to use total_size / num_chunk_sizes memory.  If
    // the total size of the pool is not enough to accomodate this, keep making
    // the next lower chunk size the max_chunk_size until it is.
    while ( max_chunk_size > total_size / num_chunk_sizes ) {
      max_chunk_size /= m_chunk_spacing;
      --num_chunk_sizes;
    }

    // We put a header at the beginnig of the device memory and use extra
    // chunks to store the header.  The header contains:
    //   size_t     chunk_size[num_chunk_sizes+1]
    //   uint64_t  freelist[num_chunk_sizes]

    // Calculate the size of the header where the size is rounded up to the
    // smallest multiple of base_chunk_size >= the needed size.  The size of the
    // chunk size array is calculated using sizeof(void*) to guarantee alignment
    // for the freelist array.  This assumes sizeof(size_t) <= sizeof(void*).
    size_t header_bytes = ( 2 * num_chunk_sizes + 1 ) * sizeof(void*);
    size_t header_size =
      ( header_bytes + base_chunk_size - 1 ) / base_chunk_size * base_chunk_size;

    // Allocate the memory including the header.
    size_t alloc_size = total_size + header_size;

#ifdef KOKKOS_MEMPOOL_PRINT_INFO
      printf( "** Allocating total %ld bytes\n", long(alloc_size));
      fflush( stdout );
#endif

    SharedRecord * rec =
      SharedRecord::allocate( memspace, "mempool", alloc_size );

#ifdef KOKKOS_MEMPOOL_PRINT_INFO
      printf( "** Allocated total %ld bytes at 0x%lx\n",
              long(alloc_size), long(rec->data()) );
      fflush( stdout );
#endif

    m_track.assign_allocated_record_to_uninitialized( rec );

    {
      // Get the pointers into the allocated memory.
      char * mem = reinterpret_cast<char *>( rec->data() );
      m_chunk_size = reinterpret_cast<size_t *>( mem );
      m_freelist = reinterpret_cast<uint64_t *>(
                   mem + ( num_chunk_sizes + 1 ) * sizeof(void*) );
      m_data = mem + header_size;

#ifdef KOKKOS_MEMPOOL_PRINT_INFO
      printf( "** Partitioning allocation 0x%lx : m_chunk_size[0x%lx] m_freelist[0x%lx] m_data[0x%lx]\n",
              (unsigned long) mem, (unsigned long) m_chunk_size,
              (unsigned long) m_freelist, (unsigned long) m_data );
      fflush( stdout );
#endif
    }

    // Initialize the chunk sizes array.  Create num_chunk_sizes different
    // chunk sizes where each successive chunk size is
    // m_chunk_spacing * previous chunk size.  The last entry in the array is
    // 0 and is used for a stopping condition.
    m_chunk_size[0] = base_chunk_size;
    for ( size_t i = 1; i < num_chunk_sizes; ++i ) {
      m_chunk_size[i] = m_chunk_size[i - 1] * m_chunk_spacing;
    }
    m_chunk_size[num_chunk_sizes] = 0;

    std::vector<size_t> num_chunks(num_chunk_sizes);

    // Set the starting point in memory and get the number of chunks for each
    // freelist.  Start with the largest chunk size to ensure usage of all the
    // memory.  If there is leftover memory for a chunk size, it will be used
    // by a smaller chunk size.
    size_t used_memory = 0;
    for ( size_t i = num_chunk_sizes; i > 0; --i ) {
      // Set the starting position in the memory for the current chunk sizes's
      // freelist and initialize the tag to 0.
      m_freelist[i - 1] = create_head( used_memory, 0UL );

      size_t mem_avail =
        total_size - (i - 1) * ( total_size / num_chunk_sizes ) - used_memory;

      // Set the number of chunks for the current chunk sizes's freelist.
      num_chunks[i - 1] = mem_avail / m_chunk_size[i - 1];

      used_memory += num_chunks[i - 1] * m_chunk_size[i - 1];
    }

#ifdef KOKKOS_MEMPOOL_PRINT_INFO
    printf( "\n" );
    printf( "*** ON HOST ***\n");
    printf( "m_chunk_size: 0x%llx\n", reinterpret_cast<uint64_t>( m_chunk_size ) );
    printf( "  m_freelist: 0x%llx\n", reinterpret_cast<uint64_t>( m_freelist ) );
    printf( "      m_data: 0x%llx\n", reinterpret_cast<uint64_t>( m_data ) );
    for ( size_t i = 0; i < num_chunk_sizes; ++i ) {
      printf( "%2lu    freelist: %10llu    chunk_size: %6lu    num_chunks: %8lu\n",
              i, get_head_offset( m_freelist[i] ), m_chunk_size[i], num_chunks[i] );
    }
    printf( "                              chunk_size: %6lu\n\n",
            m_chunk_size[num_chunk_sizes] );
    fflush( stdout );
#endif

#ifdef KOKKOS_MEMPOOL_PRINTERR
    if ( used_memory != total_size ) {
      printf( "\n** MemoryPool::MemoryPool() USED_MEMORY(%lu) != TOTAL_SIZE(%lu) **\n",
              used_memory, total_size );
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
      fflush( stdout );
#endif
      Kokkos::abort( "" );
    }
#endif

    // Create the chunks for each freelist.
    for ( size_t i = 0; i < num_chunk_sizes; ++i ) {
      // Initialize the next pointers to point to the next chunk for all but the
      // last chunk which uses a reserved value to indicate the end of the list.
      initialize_mempool<MemPoolList> im( m_data, m_chunk_size[i], num_chunks[i] - 1,
                                          get_head_offset( m_freelist[i] ) );

      Kokkos::Impl::ParallelFor< initialize_mempool<MemPoolList>, Range >
        closure( im, Range( 0, num_chunks[i] ) );

      closure.execute();

      ExecutionSpace::fence();
    }

#ifdef KOKKOS_MEMPOOL_PRINT_INFO
    print_mempool<MemPoolList> pm( num_chunk_sizes, m_chunk_size, m_freelist, m_data );

    Kokkos::Impl::ParallelFor< print_mempool<MemPoolList>, Range >
      closure( pm, Range( 0, 10 ) );

    closure.execute();

    ExecutionSpace::fence();
#endif
  }

  /// \brief Releases a lock on a freelist.
  KOKKOS_FUNCTION
  uint64_t acquire_lock( volatile uint64_t * freelist ) const;

  /// \brief Releases a lock on a freelist.
  KOKKOS_FUNCTION
  void release_lock( volatile uint64_t * freelist, uint64_t new_head ) const;

  /// \brief Tries to refill a freelist using a chunk from another freelist.
  KOKKOS_FUNCTION
  void * refill_freelist( size_t l_exp ) const;

  /// \brief Claim chunks of untracked memory from the pool.
  KOKKOS_FUNCTION
  void * allocate( size_t alloc_size ) const;

  /// \brief Release claimed memory back into the pool.
  KOKKOS_FUNCTION
  void deallocate( void * alloc_ptr, size_t alloc_size ) const;

  // \brief Pulls the offset from a freelist head.
  KOKKOS_INLINE_FUNCTION
  uint64_t get_head_offset(uint64_t head) const
  { return ( head >> TAGBITS ) << LG_MIN_CHUNKSIZE; }

  // \brief Pulls the tag from a freelist head.
  KOKKOS_INLINE_FUNCTION
  uint64_t get_head_tag(uint64_t head) const { return head & TAGBITS_MASK; }
  // \brief Creates a freelist head from a offset and tag.
  KOKKOS_INLINE_FUNCTION
  uint64_t create_head(uint64_t offset, uint64_t tag) const
  { return ( ( offset >> LG_MIN_CHUNKSIZE ) << TAGBITS ) | tag; }

  // \brief Increments a tag.
  KOKKOS_INLINE_FUNCTION
  uint64_t increment_tag(uint64_t tag) const { return ( tag + 1 ) & TAGBITS_MASK; }

  /// \brief Tests if the memory pool is empty.
  KOKKOS_INLINE_FUNCTION
  bool is_empty() const
  {
    size_t l = 0;
    while ( m_chunk_size[l] > 0 &&
            get_head_offset( m_freelist[l] ) == FREELIST_END )
    {
      ++l;
    }

    return m_chunk_size[l] == 0;
  }

  // The following functions are used for debugging.
  void print_status() const
  {
    for ( size_t l = 0; m_chunk_size[l] > 0; ++l ) {
      size_t count = 0;
      uint64_t chunk = get_head_offset( m_freelist[l] );

      while ( chunk != FREELIST_END ) {
        ++count;
        chunk = *reinterpret_cast<uint64_t *>( m_data + chunk );
      }

      printf( "chunk_size: %6lu    num_chunks: %8lu\n", m_chunk_size[l], count );
      fflush(stdout);
    }
  }

  KOKKOS_INLINE_FUNCTION
  size_t get_min_chunk_size() const { return m_chunk_size[0]; }

  size_t get_mem_size() const { return m_data_size; }
};

} // namespace Impl
} // namespace Experimental
} // namespace Kokkos

//----------------------------------------------------------------------------
/*  Prefer to implement these functions in a separate
 *  compilation unit.  For CUDA this requires nvcc command
 *  --relocatable-device-code=true
 *  When this command is set then the macro
 *  KOKKOS_CUDA_USE_RELOCATABLE_DEVICE_CODE
 *  is also set.
 */
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA ) && \
    ! defined( KOKKOS_CUDA_USE_RELOCATABLE_DEVICE_CODE )

#include <impl/Kokkos_MemoryPool_Inline.hpp>

#endif

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {

/// \class MemoryPool
/// \brief Memory management for pool of same-sized chunks of memory.
///
/// MemoryPool is a memory space that can be on host or device.  It provides a
/// pool memory allocator for fast allocation of same-sized chunks of memory.
/// The memory is only accessible on the host / device this allocator is
/// associated with.
template < typename Space , typename ExecSpace >
class MemoryPool {
private:

  Impl::MemPoolList  m_memory;

  typedef ExecSpace                     execution_space;
  typedef typename Space::memory_space  backend_memory_space;

#if defined( KOKKOS_HAVE_CUDA )

  // Current implementation requires CudaUVM memory space
  // for Cuda memory pool.

  static_assert(
    ! std::is_same< typename Space::memory_space , Kokkos::CudaSpace >::value ,
    "Kokkos::MemoryPool currently cannot use Kokkos::CudaSpace, you must use Kokkos::CudaUVMSpace" );

#endif

public:

  //! Tag this class as a kokkos memory space
  typedef MemoryPool  memory_space;

  //------------------------------------

  MemoryPool() = default;
  MemoryPool( MemoryPool && rhs ) = default;
  MemoryPool( const MemoryPool & rhs ) = default;
  MemoryPool & operator = ( MemoryPool && ) = default;
  MemoryPool & operator = ( const MemoryPool & ) = default;
  ~MemoryPool() = default;

  /// \brief Allocate memory pool
  /// \param memspace         From where to allocate the pool.
  /// \param base_chunk_size  Hand out memory in chunks of this size.
  /// \param total_size       Total size of the pool.
  MemoryPool( const backend_memory_space & memspace,
              size_t base_chunk_size, size_t total_size,
              size_t num_chunk_sizes = 4, size_t chunk_spacing = 4 )
    : m_memory( memspace, execution_space(), base_chunk_size, total_size,
                num_chunk_sizes, chunk_spacing )
  {}

  /// \brief Claim chunks of untracked memory from the pool.
  /// Can only be called from device.
  KOKKOS_INLINE_FUNCTION
  void * allocate( const size_t alloc_size ) const
  { return m_memory.allocate( alloc_size ); }

  /// \brief Release claimed memory back into the pool
  /// Can only be called from device.
  KOKKOS_INLINE_FUNCTION
  void deallocate( void * const alloc_ptr, const size_t alloc_size ) const
  { m_memory.deallocate( alloc_ptr, alloc_size ); }

  /// \brief Is out of memory at this instant
  KOKKOS_INLINE_FUNCTION
  bool is_empty() const { return m_memory.is_empty(); }

  /// \brief Minimum chunk size allocatable.
  KOKKOS_INLINE_FUNCTION
  size_t get_min_chunk_size() const { return m_memory.get_min_chunk_size(); }

  // The following functions are used for debugging.
  void print_status() const { m_memory.print_status(); }
  size_t get_mem_size() const { return m_memory.get_mem_size(); }
};

} // namespace Experimental
} // namespace Kokkos

#ifdef KOKKOS_MEMPOOL_PRINTERR
#undef KOKKOS_MEMPOOL_PRINTERR
#endif

#ifdef KOKKOS_MEMPOOL_PRINT_INFO
#undef KOKKOS_MEMPOOL_PRINT_INFO
#endif

#endif /* #define KOKKOS_MEMORYPOOL_HPP */
