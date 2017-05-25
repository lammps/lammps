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

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_Atomic.hpp>
#include <impl/Kokkos_BitOps.hpp>
#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_SharedAlloc.hpp>

#include <limits>
#include <algorithm>
#include <chrono>

// How should errors be handled?  In general, production code should return a
// value indicating failure so the user can decide how the error is handled.
// While experimental, code can abort instead.  If KOKKOS_ENABLE_MEMPOOL_PRINTERR is
// defined, the code will abort with an error message.  Otherwise, the code will
// return with a value indicating failure when possible, or do nothing instead.
//#define KOKKOS_ENABLE_MEMPOOL_PRINTERR

//#define KOKKOS_ENABLE_MEMPOOL_PRINT_INFO
//#define KOKKOS_ENABLE_MEMPOOL_PRINT_CONSTRUCTOR_INFO
//#define KOKKOS_ENABLE_MEMPOOL_PRINT_BLOCKSIZE_INFO
//#define KOKKOS_ENABLE_MEMPOOL_PRINT_SUPERBLOCK_INFO
//#define KOKKOS_ENABLE_MEMPOOL_PRINT_ACTIVE_SUPERBLOCKS
//#define KOKKOS_ENABLE_MEMPOOL_PRINT_PAGE_INFO
//#define KOKKOS_ENABLE_MEMPOOL_PRINT_INDIVIDUAL_PAGE_INFO

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {

namespace MempoolImpl {

template < typename T, typename ExecutionSpace >
struct initialize_array {
  typedef ExecutionSpace                      execution_space;
  typedef typename ExecutionSpace::size_type  size_type;

  T *  m_data;
  T    m_value;

  initialize_array( T * d, size_t size, T v ) : m_data( d ), m_value( v )
  {
    Kokkos::parallel_for( size, *this );

    execution_space::fence();
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i ) const { m_data[i] = m_value; }
};

template <typename Bitset>
struct bitset_count
{
  typedef typename Bitset::execution_space     execution_space;
  typedef typename execution_space::size_type  size_type;
  typedef typename Bitset::size_type           value_type;
  typedef typename Bitset::word_type           word_type;

  word_type *   m_words;
  value_type &  m_result;

  bitset_count( word_type * w, value_type num_words, value_type & r )
    : m_words( w ), m_result( r )
  {
    parallel_reduce( num_words, *this, m_result );
  }

  KOKKOS_INLINE_FUNCTION
  void init( value_type & v ) const
  { v = 0; }

  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type & dst, volatile value_type const & src ) const
  { dst += src; }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i, value_type & count ) const
  {
    count += Kokkos::Impl::bit_count( m_words[i] );
  }
};

template < typename Device >
class Bitset {
public:
  typedef typename Device::execution_space  execution_space;
  typedef typename Device::memory_space     memory_space;
  typedef unsigned                          word_type;
  typedef unsigned                          size_type;

  typedef Kokkos::Impl::DeepCopy< memory_space, Kokkos::HostSpace > raw_deep_copy;

  // Define some constants.
  enum {
    // Size of bitset word.  Should be 32.
    WORD_SIZE    = sizeof(word_type) * CHAR_BIT,
    LG_WORD_SIZE = Kokkos::Impl::integral_power_of_two( WORD_SIZE ),
    WORD_MASK    = WORD_SIZE - 1
  };

private:
  word_type *  m_words;
  size_type    m_size;
  size_type    m_num_words;
  word_type    m_last_word_mask;

public:
  ~Bitset() = default;
  Bitset() = default;
  Bitset( Bitset && ) = default;
  Bitset( const Bitset & ) = default;
  Bitset & operator = ( Bitset && ) = default;
  Bitset & operator = ( const Bitset & ) = default;

  void init( void * w, size_type s )
  {
    // Assumption: The size of the memory pointed to by w is a multiple of
    //             sizeof(word_type).

    m_words = reinterpret_cast<word_type*>( w );
    m_size = s;
    m_num_words = ( s + WORD_SIZE - 1 ) >> LG_WORD_SIZE;
    m_last_word_mask = m_size & WORD_MASK ? ( word_type(1) << ( m_size & WORD_MASK ) ) - 1 : 0;

    reset();
  }

  size_type size() const { return m_size; }

  size_type count() const
  {
    size_type val = 0;
    bitset_count< Bitset > bc( m_words, m_num_words, val );
    return val;
  }

  void set()
  {
    // Set all the bits.
    initialize_array< word_type, execution_space > ia( m_words, m_num_words, ~word_type(0) );

    if ( m_last_word_mask ) {
      // Clear the unused bits in the last block.
      raw_deep_copy( m_words + ( m_num_words - 1 ), &m_last_word_mask, sizeof(word_type) );
    }
  }

  void reset()
  {
    initialize_array< word_type, execution_space > ia( m_words, m_num_words, word_type(0) );
  }

  KOKKOS_FORCEINLINE_FUNCTION
  bool test( size_type i ) const
  {
    size_type word_pos = i >> LG_WORD_SIZE;
    word_type word = volatile_load( &m_words[ word_pos ] );
    word_type mask = word_type(1) << ( i & WORD_MASK );

    return word & mask;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  bool set( size_type i ) const
  {
    size_type word_pos = i >> LG_WORD_SIZE;
    word_type mask = word_type(1) << ( i & WORD_MASK );

    return !( atomic_fetch_or( &m_words[ word_pos ], mask ) & mask );
  }

  KOKKOS_FORCEINLINE_FUNCTION
  bool reset( size_type i ) const
  {
    size_type word_pos = i >> LG_WORD_SIZE;
    word_type mask = word_type(1) << ( i & WORD_MASK );

    return atomic_fetch_and( &m_words[ word_pos ], ~mask ) & mask;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  Kokkos::pair< bool, word_type >
  fetch_word_set( size_type i ) const
  {
    size_type word_pos = i >> LG_WORD_SIZE;
    word_type mask = word_type(1) << ( i & WORD_MASK );

    Kokkos::pair<bool, word_type> result;
    result.second = atomic_fetch_or( &m_words[ word_pos ], mask );
    result.first = !( result.second & mask );

    return result;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  Kokkos::pair< bool, word_type >
  fetch_word_reset( size_type i ) const
  {
    size_type word_pos = i >> LG_WORD_SIZE;
    word_type mask = word_type(1) << ( i & WORD_MASK );

    Kokkos::pair<bool, word_type> result;
    result.second = atomic_fetch_and( &m_words[ word_pos ], ~mask );
    result.first = result.second & mask;

    return result;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  Kokkos::pair< bool, word_type >
  set_any_in_word( size_type & pos ) const
  {
    size_type word_pos = pos >> LG_WORD_SIZE;
    word_type word = volatile_load( &m_words[ word_pos ] );

    // Loop until there are no more unset bits in the word.
    while ( ~word ) {
      // Find the first unset bit in the word.
      size_type bit = Kokkos::Impl::bit_scan_forward( ~word );

      // Try to set the bit.
      word_type mask = word_type(1) << bit;
      word = atomic_fetch_or( &m_words[ word_pos ], mask );

      if ( !( word & mask ) ) {
        // Successfully set the bit.
        pos = ( word_pos << LG_WORD_SIZE ) + bit;

        return Kokkos::pair<bool, word_type>( true, word );
      }
    }

    // Didn't find a free bit in this word.
    return Kokkos::pair<bool, word_type>( false, word_type(0) );
  }

  KOKKOS_FORCEINLINE_FUNCTION
  Kokkos::pair< bool, word_type >
  set_any_in_word( size_type & pos, word_type word_mask ) const
  {
    size_type word_pos = pos >> LG_WORD_SIZE;
    word_type word = volatile_load( &m_words[ word_pos ] );
    word = ( ~word ) & word_mask;

    // Loop until there are no more unset bits in the word.
    while ( word ) {
      // Find the first unset bit in the word.
      size_type bit = Kokkos::Impl::bit_scan_forward( word );

      // Try to set the bit.
      word_type mask = word_type(1) << bit;
      word = atomic_fetch_or( &m_words[ word_pos ], mask );

      if ( !( word & mask ) ) {
        // Successfully set the bit.
        pos = ( word_pos << LG_WORD_SIZE ) + bit;

        return Kokkos::pair<bool, word_type>( true, word );
      }

      word = ( ~word ) & word_mask;
    }

    // Didn't find a free bit in this word.
    return Kokkos::pair<bool, word_type>( false, word_type(0) );
  }

  KOKKOS_FORCEINLINE_FUNCTION
  Kokkos::pair< bool, word_type >
  reset_any_in_word( size_type & pos ) const
  {
    size_type word_pos = pos >> LG_WORD_SIZE;
    word_type word = volatile_load( &m_words[ word_pos ] );

    // Loop until there are no more set bits in the word.
    while ( word ) {
      // Find the first unset bit in the word.
      size_type bit = Kokkos::Impl::bit_scan_forward( word );

      // Try to reset the bit.
      word_type mask = word_type(1) << bit;
      word = atomic_fetch_and( &m_words[ word_pos ], ~mask );

      if ( word & mask ) {
        // Successfully reset the bit.
        pos = ( word_pos << LG_WORD_SIZE ) + bit;

        return Kokkos::pair<bool, word_type>( true, word );
      }
    }

    // Didn't find a free bit in this word.
    return Kokkos::pair<bool, word_type>( false, word_type(0) );
  }

  KOKKOS_FORCEINLINE_FUNCTION
  Kokkos::pair< bool, word_type >
  reset_any_in_word( size_type & pos, word_type word_mask ) const
  {
    size_type word_pos = pos >> LG_WORD_SIZE;
    word_type word = volatile_load( &m_words[ word_pos ] );
    word = word & word_mask;

    // Loop until there are no more set bits in the word.
    while ( word ) {
      // Find the first unset bit in the word.
      size_type bit = Kokkos::Impl::bit_scan_forward( word );

      // Try to reset the bit.
      word_type mask = word_type(1) << bit;
      word = atomic_fetch_and( &m_words[ word_pos ], ~mask );

      if ( word & mask ) {
        // Successfully reset the bit.
        pos = ( word_pos << LG_WORD_SIZE ) + bit;

        return Kokkos::pair<bool, word_type>( true, word );
      }

      word = word & word_mask;
    }

    // Didn't find a free bit in this word.
    return Kokkos::pair<bool, word_type>( false, word_type(0) );
  }
};

template < typename UInt32View, typename BSHeaderView, typename SBHeaderView,
           typename MempoolBitset >
struct create_histogram {
  typedef typename UInt32View::execution_space  execution_space;
  typedef typename execution_space::size_type   size_type;
  typedef Kokkos::pair< double, uint32_t >      value_type;

  size_t         m_start;
  UInt32View     m_page_histogram;
  BSHeaderView   m_blocksize_info;
  SBHeaderView   m_sb_header;
  MempoolBitset  m_sb_blocks;
  size_t         m_lg_max_sb_blocks;
  uint32_t       m_lg_min_block_size;
  uint32_t       m_blocks_per_page;
  value_type &   m_result;

  create_histogram( size_t start, size_t end, UInt32View ph, BSHeaderView bsi,
                    SBHeaderView sbh, MempoolBitset sbb, size_t lmsb,
                    uint32_t lmbs, uint32_t bpp, value_type & r )
    : m_start( start ), m_page_histogram( ph ), m_blocksize_info( bsi ),
      m_sb_header( sbh ), m_sb_blocks( sbb ), m_lg_max_sb_blocks( lmsb ),
      m_lg_min_block_size( lmbs ), m_blocks_per_page( bpp ), m_result( r )
  {
    Kokkos::parallel_reduce( end - start, *this, m_result );

    execution_space::fence();
  }

  KOKKOS_INLINE_FUNCTION
  void init( value_type & v ) const
  {
    v.first  = 0.0;
    v.second = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type & dst, volatile value_type const & src ) const
  {
    dst.first += src.first;
    dst.second += src.second;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i, value_type & r ) const
  {
    size_type i2 = i + m_start;

    uint32_t lg_block_size = m_sb_header(i2).m_lg_block_size;

    // A superblock only has a block size of 0 when it is empty.
    if ( lg_block_size != 0 ) {
      uint32_t block_size_id = lg_block_size - m_lg_min_block_size;
      uint32_t blocks_per_sb = m_blocksize_info[block_size_id].m_blocks_per_sb;
      uint32_t pages_per_sb = m_blocksize_info[block_size_id].m_pages_per_sb;

      uint32_t total_allocated_blocks = 0;

      for ( uint32_t j = 0; j < pages_per_sb; ++j ) {
        unsigned start_pos = ( i2 << m_lg_max_sb_blocks ) + j * m_blocks_per_page;
        unsigned end_pos = start_pos + m_blocks_per_page;
        uint32_t page_allocated_blocks = 0;

        for ( unsigned k = start_pos; k < end_pos; ++k ) {
          page_allocated_blocks += m_sb_blocks.test( k );
        }

        total_allocated_blocks += page_allocated_blocks;

        atomic_increment( &m_page_histogram(page_allocated_blocks) );
      }

      r.first += double(total_allocated_blocks) / blocks_per_sb;
      r.second += blocks_per_sb;
    }
  }
};

#ifdef KOKKOS_ENABLE_MEMPOOL_PRINT_SUPERBLOCK_INFO
template < typename UInt32View, typename SBHeaderView, typename MempoolBitset >
struct count_allocated_blocks {
  typedef typename UInt32View::execution_space  execution_space;
  typedef typename execution_space::size_type   size_type;

  UInt32View     m_num_allocated_blocks;
  SBHeaderView   m_sb_header;
  MempoolBitset  m_sb_blocks;
  size_t         m_sb_size;
  size_t         m_lg_max_sb_blocks;

  count_allocated_blocks( size_t num_sb, UInt32View nab, SBHeaderView sbh,
                          MempoolBitset sbb, size_t sbs, size_t lmsb )
    : m_num_allocated_blocks( nab ), m_sb_header( sbh ),
      m_sb_blocks( sbb ), m_sb_size( sbs ), m_lg_max_sb_blocks( lmsb )
  {
    Kokkos::parallel_for( num_sb, *this );

    execution_space::fence();
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i ) const
  {
    uint32_t lg_block_size = m_sb_header(i).m_lg_block_size;

    // A superblock only has a block size of 0 when it is empty.
    if ( lg_block_size != 0 ) {
      // Count the allocated blocks in the superblock.
      uint32_t blocks_per_sb = lg_block_size > 0 ? m_sb_size >> lg_block_size : 0;
      unsigned start_pos = i << m_lg_max_sb_blocks;
      unsigned end_pos = start_pos + blocks_per_sb;
      uint32_t count = 0;

      for ( unsigned j = start_pos; j < end_pos; ++j ) {
        count += m_sb_blocks.test( j );
      }

      m_num_allocated_blocks(i) = count;
    }
  }
};
#endif

}

/// \class MemoryPool
/// \brief Bitset based memory manager for pools of same-sized chunks of memory.
/// \tparam Device Kokkos device that gives the execution and memory space the
///                allocator will be used in.
///
/// MemoryPool is a memory space that can be on host or device.  It provides a
/// pool memory allocator for fast allocation of same-sized chunks of memory.
/// The memory is only accessible on the host / device this allocator is
/// associated with.
///
/// This allocator is based on ideas from the following GPU allocators:
///   Halloc (https://github.com/canonizer/halloc).
///   ScatterAlloc (https://github.com/ComputationalRadiationPhysics/scatteralloc)
template < typename Device >
class MemoryPool {
private:
  // The allocator uses superblocks.  A superblock is divided into pages, and a
  // page is divided into blocks.  A block is the chunk of memory that is given
  // out by the allocator.  A page always has a number of blocks equal to the
  // size of the word used by the bitset.  Thus, the pagesize can vary between
  // superblocks as it is based on the block size of the superblock.  The
  // allocator supports all powers of 2 from MIN_BLOCK_SIZE to the size of a
  // superblock as block sizes.

  // Superblocks are divided into 4 categories:
  //   1. empty    - is completely empty; there are no active allocations
  //   2. partfull - partially full; there are some active allocations
  //   3. full     - full enough with active allocations that new allocations
  //                 will likely fail
  //   4. active   - is currently the active superblock for a block size
  //
  // An inactive superblock is one that is empty, partfull, or full.
  //
  // New allocations occur only from an active superblock.  If a superblock is
  // made inactive after an allocation request is made to it but before the
  // allocation request is fulfilled, the allocation will still be attempted
  // from that superblock.  Deallocations can  occur to partfull, full, or
  // active superblocks.  Superblocks move between categories as allocations
  // and deallocations happen.  Superblocks all start empty.
  //
  // Here are the possible moves between categories:
  //   empty    -> active    During allocation, there is no active superblock
  //                         or the active superblock is full.
  //   active   -> full      During allocation, the full threshold of the
  //                         superblock is reached when increasing the fill
  //                         level.
  //   full     -> partfull  During deallocation, the full threshold of the
  //                         superblock is crossed when decreasing the fill
  //                         level.
  //   partfull -> empty     Deallocation of the last allocated block of an
  //                         inactive superblock.
  //   partfull -> active    During allocation, the active superblock is full.
  //
  // When a new active superblock is needed, partfull superblocks of the same
  // block size are chosen over empty superblocks.
  //
  // The empty and partfull superblocks are tracked using bitsets that represent
  // the superblocks in those repsective categories.  Empty superblocks use a
  // single bitset, while partfull superblocks use a bitset per block size
  // (contained sequentially in a single bitset).  Active superblocks are
  // tracked by the active superblocks array.  Full superblocks aren't tracked
  // at all.

  typedef typename Device::execution_space    execution_space;
  typedef typename Device::memory_space       backend_memory_space;
  typedef Device                              device_type;
  typedef MempoolImpl::Bitset< device_type >  MempoolBitset;

  // Define some constants.
  enum {
    MIN_BLOCK_SIZE     = 64,
    LG_MIN_BLOCK_SIZE  = Kokkos::Impl::integral_power_of_two( MIN_BLOCK_SIZE ),
    MAX_BLOCK_SIZES    = 31 - LG_MIN_BLOCK_SIZE + 1,

    // Size of bitset word.
    BLOCKS_PER_PAGE    = MempoolBitset::WORD_SIZE,
    LG_BLOCKS_PER_PAGE = MempoolBitset::LG_WORD_SIZE,

    INVALID_SUPERBLOCK = ~uint32_t(0),
    SUPERBLOCK_LOCK    = ~uint32_t(0) - 1,

    MAX_TRIES          = 32             // Cap on the number of pages searched
                                        // before an allocation returns empty.
  };

public:
  // Stores information about each superblock.
  struct SuperblockHeader {
    uint32_t  m_full_pages;
    uint32_t  m_empty_pages;
    uint32_t  m_lg_block_size;
    uint32_t  m_is_active;

    KOKKOS_FUNCTION
    SuperblockHeader() :
      m_full_pages(0), m_empty_pages(0), m_lg_block_size(0), m_is_active(false) {}
  };

  // Stores information about each block size.
  struct BlockSizeHeader {
    uint32_t  m_blocks_per_sb;
    uint32_t  m_pages_per_sb;
    uint32_t  m_sb_full_level;
    uint32_t  m_page_full_level;

    KOKKOS_FUNCTION
    BlockSizeHeader() :
      m_blocks_per_sb(0), m_pages_per_sb(0), m_sb_full_level(0), m_page_full_level(0) {}
  };

private:
  typedef Kokkos::Impl::SharedAllocationTracker    Tracker;
  typedef View< uint32_t *, device_type >          UInt32View;
  typedef View< SuperblockHeader *, device_type >  SBHeaderView;

  // The letters 'sb' used in any variable name mean superblock.

  size_t           m_lg_sb_size;        // Log2 of superblock size.
  size_t           m_sb_size;           // Superblock size.
  size_t           m_lg_max_sb_blocks;  // Log2 of the number of blocks of the
                                        // minimum block size in a superblock.
  size_t           m_num_sb;            // Number of superblocks.
  size_t           m_ceil_num_sb;       // Number of superblocks rounded up to the smallest
                                        // multiple of the bitset word size.  Used by
                                        // bitsets representing superblock categories to
                                        // ensure different block sizes never share a word
                                        // in the bitset.
  size_t           m_num_block_size;    // Number of block sizes supported.
  size_t           m_data_size;         // Amount of memory available to the allocator.
  size_t           m_sb_blocks_size;    // Amount of memory for free / empty blocks bitset.
  size_t           m_empty_sb_size;     // Amount of memory for empty superblocks bitset.
  size_t           m_partfull_sb_size;  // Amount of memory for partfull superblocks bitset.
  size_t           m_total_size;        // Total amount of memory allocated.
  char *           m_data;              // Beginning device memory location used for
                                        // superblocks.
  UInt32View       m_active;            // Active superblocks IDs.
  SBHeaderView     m_sb_header;         // Header info for superblocks.
  MempoolBitset    m_sb_blocks;         // Bitsets representing free / allocated status
                                        // of blocks in superblocks.
  MempoolBitset    m_empty_sb;          // Bitset representing empty superblocks.
  MempoolBitset    m_partfull_sb;       // Bitsets representing partially full superblocks.
  Tracker          m_track;             // Tracker for superblock memory.
  BlockSizeHeader  m_blocksize_info[MAX_BLOCK_SIZES];  // Header info for block sizes.

  // There were several methods tried for storing the block size header info: in a View,
  // in a View of const data, and in a RandomAccess View.  All of these were slower than
  // storing it in a static array that is a member variable to the class.  In the latter
  // case, the block size info gets copied into the constant memory on the GPU along with
  // the class when it is copied there for exeucting a parallel loop.  Instead of storing
  // the values, computing the values every time they were needed was also tried.  This
  // method was slightly slower than storing them in the static array.

public:
  //! Tag this class as a kokkos memory space
  typedef MemoryPool  memory_space;

  ~MemoryPool() = default;
  MemoryPool() = default;
  MemoryPool( MemoryPool && ) = default;
  MemoryPool( const MemoryPool & ) = default;
  MemoryPool & operator = ( MemoryPool && ) = default;
  MemoryPool & operator = ( const MemoryPool & ) = default;

  /// \brief Initializes the memory pool.
  /// \param memspace The memory space from which the memory pool will allocate memory.
  /// \param total_size The requested memory amount controlled by the allocator.  The
  ///                   actual amount is rounded up to the smallest multiple of the
  ///                   superblock size >= the requested size.
  /// \param log2_superblock_size Log2 of the size of superblocks used by the allocator.
  ///                             In most use cases, the default value should work.
  inline
  MemoryPool( const backend_memory_space & memspace,
              size_t total_size, size_t log2_superblock_size = 20 )
    : m_lg_sb_size( log2_superblock_size ),
      m_sb_size( size_t(1) << m_lg_sb_size ),
      m_lg_max_sb_blocks( m_lg_sb_size - LG_MIN_BLOCK_SIZE ),
      m_num_sb( ( total_size + m_sb_size - 1 ) >> m_lg_sb_size ),
      m_ceil_num_sb( ( ( m_num_sb + BLOCKS_PER_PAGE - 1 ) >> LG_BLOCKS_PER_PAGE ) <<
                     LG_BLOCKS_PER_PAGE ),
      m_num_block_size( m_lg_sb_size - LG_MIN_BLOCK_SIZE + 1 ),
      m_data_size( m_num_sb * m_sb_size ),
      m_sb_blocks_size( ( m_num_sb << m_lg_max_sb_blocks ) / CHAR_BIT ),
      m_empty_sb_size( m_ceil_num_sb / CHAR_BIT ),
      m_partfull_sb_size( m_ceil_num_sb * m_num_block_size / CHAR_BIT ),
      m_total_size( m_data_size +  m_sb_blocks_size + m_empty_sb_size + m_partfull_sb_size ),
      m_data(0),
      m_active( "Active superblocks" ),
      m_sb_header( "Superblock headers" ),
      m_track()
  {
    // Assumption.  The minimum block size must be a power of 2.
    static_assert( Kokkos::Impl::is_integral_power_of_two( MIN_BLOCK_SIZE ), "" );

    // Assumption.  Require a superblock be large enough so it takes at least 1
    // whole bitset word to represent it using the minimum blocksize.
    if ( m_sb_size < MIN_BLOCK_SIZE * BLOCKS_PER_PAGE ) {
      printf( "\n** MemoryPool::MemoryPool() Superblock size must be >= %u **\n",
              MIN_BLOCK_SIZE * BLOCKS_PER_PAGE );
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
      fflush( stdout );
#endif
      Kokkos::abort( "" );
    }

    // Assumption.  A superblock's size can be at most 2^31.  Verify this.
    if ( m_lg_sb_size > 31 ) {
      printf( "\n** MemoryPool::MemoryPool() Superblock size must be < %u **\n",
              ( uint32_t(1) << 31 ) );
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
      fflush( stdout );
#endif
      Kokkos::abort( "" );
    }

    // Assumption.  The Bitset only uses unsigned for size types which limits
    // the amount of memory the allocator can manage.  Verify the memory size
    // is below this limit.
    if ( m_data_size > size_t(MIN_BLOCK_SIZE) * std::numeric_limits<unsigned>::max() ) {
      printf( "\n** MemoryPool::MemoryPool() Allocator can only manage %lu bytes of memory; requested %lu **\n",
              size_t(MIN_BLOCK_SIZE) * std::numeric_limits<unsigned>::max(), total_size );
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
      fflush( stdout );
#endif
      Kokkos::abort( "" );
    }

    // Allocate memory for Views.  This is done here instead of at construction
    // so that the runtime checks can be performed before allocating memory.
    resize( m_active, m_num_block_size );
    resize( m_sb_header, m_num_sb );

    // Allocate superblock memory.
    typedef Kokkos::Impl::SharedAllocationRecord< backend_memory_space, void >  SharedRecord;
    SharedRecord * rec =
      SharedRecord::allocate( memspace, "mempool", m_total_size );

    m_track.assign_allocated_record_to_uninitialized( rec );
    m_data = reinterpret_cast<char *>( rec->data() );

    // Set and initialize the free / empty block bitset memory.
    m_sb_blocks.init( m_data + m_data_size, m_num_sb << m_lg_max_sb_blocks );

    // Set and initialize the empty superblock block bitset memory.
    m_empty_sb.init( m_data + m_data_size + m_sb_blocks_size, m_num_sb );

    // Start with all superblocks in the empty category.
    m_empty_sb.set();

    // Set and initialize the partfull superblock block bitset memory.
    m_partfull_sb.init( m_data + m_data_size + m_sb_blocks_size + m_empty_sb_size,
                        m_ceil_num_sb * m_num_block_size );

    // Initialize all active superblocks to be invalid.
    typename UInt32View::HostMirror host_active = create_mirror_view( m_active );
    for ( size_t i = 0; i < m_num_block_size; ++i ) host_active(i) = INVALID_SUPERBLOCK;
    deep_copy( m_active, host_active );

    // A superblock is considered full when this percentage of its pages are full.
    const double superblock_full_fraction = .8;

    // A page is considered full when this percentage of its blocks are full.
    const double page_full_fraction = .875;

    // Initialize the blocksize info.
    for ( size_t i = 0; i < m_num_block_size; ++i ) {
      uint32_t lg_block_size = i + LG_MIN_BLOCK_SIZE;
      uint32_t blocks_per_sb = m_sb_size >> lg_block_size;
      uint32_t pages_per_sb = ( blocks_per_sb + BLOCKS_PER_PAGE - 1 ) >> LG_BLOCKS_PER_PAGE;

      m_blocksize_info[i].m_blocks_per_sb = blocks_per_sb;
      m_blocksize_info[i].m_pages_per_sb = pages_per_sb;

      // Set the full level for the superblock.
      m_blocksize_info[i].m_sb_full_level =
        static_cast<uint32_t>( pages_per_sb * superblock_full_fraction );

      if ( m_blocksize_info[i].m_sb_full_level == 0 ) {
        m_blocksize_info[i].m_sb_full_level = 1;
      }

      // Set the full level for the page.
      uint32_t blocks_per_page =
        blocks_per_sb < BLOCKS_PER_PAGE ? blocks_per_sb : BLOCKS_PER_PAGE;

      m_blocksize_info[i].m_page_full_level =
        static_cast<uint32_t>( blocks_per_page * page_full_fraction );

      if ( m_blocksize_info[i].m_page_full_level == 0 ) {
        m_blocksize_info[i].m_page_full_level = 1;
      }
    }

#ifdef KOKKOS_ENABLE_MEMPOOL_PRINT_CONSTRUCTOR_INFO
    printf( "\n" );
    printf( "      m_lg_sb_size: %12lu\n", m_lg_sb_size );
    printf( "         m_sb_size: %12lu\n", m_sb_size );
    printf( "   m_max_sb_blocks: %12lu\n", size_t(1) << m_lg_max_sb_blocks );
    printf( "m_lg_max_sb_blocks: %12lu\n", m_lg_max_sb_blocks );
    printf( "          m_num_sb: %12lu\n", m_num_sb );
    printf( "     m_ceil_num_sb: %12lu\n", m_ceil_num_sb );
    printf( "  m_num_block_size: %12lu\n", m_num_block_size );
    printf( "        data bytes: %12lu\n", m_data_size );
    printf( "   sb_blocks bytes: %12lu\n", m_sb_blocks_size );
    printf( "    empty_sb bytes: %12lu\n", m_empty_sb_size );
    printf( " partfull_sb bytes: %12lu\n", m_partfull_sb_size );
    printf( "       total bytes: %12lu\n", m_total_size );
    printf( "   m_empty_sb size: %12u\n", m_empty_sb.size() );
    printf( "m_partfull_sb size: %12u\n", m_partfull_sb.size() );
    printf( "\n" );
    fflush( stdout );
#endif

#ifdef KOKKOS_ENABLE_MEMPOOL_PRINT_BLOCKSIZE_INFO
    // Print the blocksize info for all the block sizes.
    printf( "SIZE    BLOCKS_PER_SB    PAGES_PER_SB    SB_FULL_LEVEL    PAGE_FULL_LEVEL\n" );
    for ( size_t i = 0; i < m_num_block_size; ++i ) {
      printf( "%4zu    %13u    %12u    %13u    %15u\n", i + LG_MIN_BLOCK_SIZE,
              m_blocksize_info[i].m_blocks_per_sb, m_blocksize_info[i].m_pages_per_sb,
              m_blocksize_info[i].m_sb_full_level, m_blocksize_info[i].m_page_full_level );
    }
    printf( "\n" );
#endif
  }

  /// \brief  The actual block size allocated given alloc_size.
  KOKKOS_INLINE_FUNCTION
  size_t allocate_block_size( const size_t alloc_size ) const
  { return size_t(1) << ( get_block_size_index( alloc_size ) + LG_MIN_BLOCK_SIZE ); }

  /// \brief Allocate a chunk of memory.
  /// \param alloc_size Size of the requested allocated in number of bytes.
  ///
  /// The function returns a void pointer to a memory location on success and
  /// NULL on failure.
  KOKKOS_FUNCTION
  void * allocate( size_t alloc_size ) const
  {
    void * p = 0;

    // Only support allocations up to the superblock size.  Just return 0
    // (failed allocation) for any size above this.
    if ( alloc_size <= m_sb_size )
    {
      int block_size_id = get_block_size_index( alloc_size );
      uint32_t blocks_per_sb = m_blocksize_info[block_size_id].m_blocks_per_sb;
      uint32_t pages_per_sb = m_blocksize_info[block_size_id].m_pages_per_sb;

#ifdef KOKKOS_IMPL_CUDA_CLANG_WORKAROUND
      // Without this test it looks like pages_per_sb might come back wrong.
      if ( pages_per_sb == 0 ) return NULL;
#endif

      unsigned word_size = blocks_per_sb > 32 ? 32 : blocks_per_sb;
      unsigned word_mask = ( uint64_t(1) << word_size ) - 1;

      // Instead of forcing an atomic read to guarantee the updated value,
      // reading the old value is actually beneficial because more threads will
      // attempt allocations on the old active superblock instead of waiting on
      // the new active superblock.  This will help hide the latency of
      // switching the active superblock.
      uint32_t sb_id = volatile_load( &m_active(block_size_id) );

      // If the active is locked, keep reading it atomically until the lock is
      // released.
      while ( sb_id == SUPERBLOCK_LOCK ) {
        sb_id = atomic_fetch_or( &m_active(block_size_id), uint32_t(0) );
      }

      load_fence();

      bool allocation_done = false;

      while ( !allocation_done ) {
        bool need_new_sb = false;

        if ( sb_id != INVALID_SUPERBLOCK ) {
          // Use the value from the clock register as the hash value.
          uint64_t hash_val = get_clock_register();

          // Get the starting position for this superblock's bits in the bitset.
          uint32_t pos_base = sb_id << m_lg_max_sb_blocks;

          // Mod the hash value to choose a page in the superblock.  The
          // initial block searched is the first block of that page.
          uint32_t pos_rel = uint32_t( hash_val & ( pages_per_sb - 1 ) ) << LG_BLOCKS_PER_PAGE;

          // Get the absolute starting position for this superblock's bits in the bitset.
          uint32_t pos = pos_base + pos_rel;

          // Keep track of the number of pages searched.  Pages in the superblock are
          // searched linearly from the starting page.  All pages in the superblock are
          // searched until either a location is found, or it is proven empty.
          uint32_t pages_searched = 0;

          bool search_done = false;

          while ( !search_done ) {
            bool success = false;
            unsigned prev_val = 0;

            Kokkos::tie( success, prev_val ) = m_sb_blocks.set_any_in_word( pos, word_mask );

            if ( !success ) {
              if ( ++pages_searched >= pages_per_sb ) {
                // Searched all the pages in this superblock.  Look for a new superblock.
                //
                // The previous method tried limiting the number of pages searched, but
                // that caused a huge performance issue in CUDA where the outer loop
                // executed massive numbers of times.  Threads weren't able to find a
                // free location when the superblock wasn't full and were able to execute
                // the outer loop many times before the superblock was switched for a new
                // one.  Switching to an exhaustive search eliminated this possiblity and
                // didn't slow anything down for the tests.
                need_new_sb = true;
                search_done = true;
              }
              else {
                // Move to the next page making sure the new search position
                // doesn't go past this superblock's bits.
                pos += BLOCKS_PER_PAGE;
                pos = ( pos < pos_base + blocks_per_sb ) ? pos : pos_base;
              }
            }
            else {
              // Reserved a memory location to allocate.
              memory_fence();

              search_done = true;
              allocation_done = true;

              uint32_t lg_block_size = block_size_id + LG_MIN_BLOCK_SIZE;

              p = m_data + ( size_t(sb_id) << m_lg_sb_size ) +
                  ( ( pos - pos_base ) << lg_block_size );

              uint32_t used_bits = Kokkos::Impl::bit_count( prev_val );

              if ( used_bits == 0 ) {
                // This page was empty.  Decrement the number of empty pages for
                // the superblock.
                atomic_decrement( &m_sb_header(sb_id).m_empty_pages );
              }
              else if ( used_bits == m_blocksize_info[block_size_id].m_page_full_level - 1 )
              {
                // This page is full.  Increment the number of full pages for
                // the superblock.
                uint32_t full_pages = atomic_fetch_add( &m_sb_header(sb_id).m_full_pages, 1 );

                // This allocation made the superblock full, so a new one needs to be found.
                if ( full_pages == m_blocksize_info[block_size_id].m_sb_full_level - 1 ) {
                  need_new_sb = true;
                }
              }
            }
          }
        }
        else {
          // This is the first allocation for this block size.  A superblock needs
          // to be set as the active one.  If this point is reached any other time,
          // it is an error.
          need_new_sb = true;
        }

        if ( need_new_sb ) {
          uint32_t new_sb_id = find_superblock( block_size_id, sb_id );

          if ( new_sb_id == sb_id ) {
            allocation_done = true;
#ifdef KOKKOS_ENABLE_MEMPOOL_PRINT_INFO
            printf( "** No superblocks available. **\n" );
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
            fflush( stdout );
#endif
#endif
          }
          else {
            sb_id = new_sb_id;
          }
        }
      }
    }
#ifdef KOKKOS_ENABLE_MEMPOOL_PRINT_INFO
    else {
      printf( "** Requested allocation size (%zu) larger than superblock size (%lu). **\n",
              alloc_size, m_sb_size );
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
      fflush( stdout );
#endif
    }
#endif

    return p;
  }

  /// \brief Release allocated memory back to the pool.
  /// \param alloc_ptr Pointer to chunk of memory previously allocated by
  ///                  the allocator.
  /// \param alloc_size Size of the allocated memory in number of bytes.
  KOKKOS_FUNCTION
  void deallocate( void * alloc_ptr, size_t alloc_size ) const
  {
    char * ap = static_cast<char *>( alloc_ptr );

    // Only deallocate memory controlled by this pool.
    if ( ap >= m_data && ap + alloc_size <= m_data + m_data_size ) {
      // Get the superblock for the address.  This can be calculated by math on
      // the address since the superblocks are stored contiguously in one memory
      // chunk.
      uint32_t sb_id = ( ap - m_data ) >> m_lg_sb_size;

      // Get the starting position for this superblock's bits in the bitset.
      uint32_t pos_base = sb_id << m_lg_max_sb_blocks;

      // Get the relative position for this memory location's bit in the bitset.
      uint32_t offset = ( ap - m_data ) - ( size_t(sb_id) << m_lg_sb_size );
      uint32_t lg_block_size = m_sb_header(sb_id).m_lg_block_size;
      uint32_t block_size_id = lg_block_size - LG_MIN_BLOCK_SIZE;
      uint32_t pos_rel = offset >> lg_block_size;

      bool success = false;
      unsigned prev_val = 0;

      memory_fence();

      Kokkos::tie( success, prev_val ) = m_sb_blocks.fetch_word_reset( pos_base + pos_rel );

      // If the memory location was previously deallocated, do nothing.
      if ( success ) {
        uint32_t page_fill_level = Kokkos::Impl::bit_count( prev_val );

        if ( page_fill_level == 1 ) {
          // This page is now empty.  Increment the number of empty pages for the
          // superblock.
          uint32_t empty_pages = atomic_fetch_add( &m_sb_header(sb_id).m_empty_pages, 1 );

          if ( !volatile_load( &m_sb_header(sb_id).m_is_active ) &&
               empty_pages == m_blocksize_info[block_size_id].m_pages_per_sb - 1 )
          {
            // This deallocation caused the superblock to be empty.  Change the
            // superblock category from partially full to empty.
            unsigned pos = block_size_id * m_ceil_num_sb + sb_id;

            if ( m_partfull_sb.reset( pos ) ) {
              // Reset the empty pages and block size for the superblock.
              volatile_store( &m_sb_header(sb_id).m_empty_pages, uint32_t(0) );
              volatile_store( &m_sb_header(sb_id).m_lg_block_size, uint32_t(0) );

              store_fence();

              m_empty_sb.set( sb_id );
            }
          }
        }
        else if ( page_fill_level == m_blocksize_info[block_size_id].m_page_full_level ) {
          // This page is no longer full.  Decrement the number of full pages for
          // the superblock.
          uint32_t full_pages = atomic_fetch_sub( &m_sb_header(sb_id).m_full_pages, 1 );

          if ( !volatile_load( &m_sb_header(sb_id).m_is_active ) &&
               full_pages == m_blocksize_info[block_size_id].m_sb_full_level )
          {
            // This deallocation caused the number of full pages to decrease below
            // the full threshold.  Change the superblock category from full to
            // partially full.
            unsigned pos = block_size_id * m_ceil_num_sb + sb_id;
            m_partfull_sb.set( pos );
          }
        }
      }
    }
#ifdef KOKKOS_ENABLE_MEMPOOL_PRINTERR
    else {
      printf( "\n** MemoryPool::deallocate() ADDRESS_OUT_OF_RANGE(0x%llx) **\n",
              reinterpret_cast<uint64_t>( alloc_ptr ) );
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
      fflush( stdout );
#endif
    }
#endif
  }

  /// \brief Tests if the memory pool has no more memory available to allocate.
  KOKKOS_INLINE_FUNCTION
  bool is_empty() const
  {
    // The allocator is empty if all superblocks are full.  A superblock is
    // full if it has >= 80% of its pages allocated.

    // Look at all the superblocks.  If one is not full, then the allocator
    // isn't empty.
    for ( size_t i = 0; i < m_num_sb; ++i ) {
      uint32_t lg_block_size = m_sb_header(i).m_lg_block_size;

      // A superblock only has a block size of 0 when it is empty.
      if ( lg_block_size == 0 ) return false;

      uint32_t block_size_id = lg_block_size - LG_MIN_BLOCK_SIZE;
      uint32_t full_pages = volatile_load( &m_sb_header(i).m_full_pages );

      if ( full_pages < m_blocksize_info[block_size_id].m_sb_full_level ) return false;
    }

    // All the superblocks were full.  The allocator is empty.
    return true;
  }

  // The following functions are used for debugging.
  void print_status() const
  {
    printf( "\n" );

#ifdef KOKKOS_ENABLE_MEMPOOL_PRINT_SUPERBLOCK_INFO
    typename SBHeaderView::HostMirror host_sb_header = create_mirror_view( m_sb_header );
    deep_copy( host_sb_header, m_sb_header );

    UInt32View num_allocated_blocks( "Allocated Blocks", m_num_sb );

    // Count the number of allocated blocks per superblock.
    {
      MempoolImpl::count_allocated_blocks< UInt32View, SBHeaderView, MempoolBitset >
        mch( m_num_sb, num_allocated_blocks, m_sb_header,
             m_sb_blocks, m_sb_size, m_lg_max_sb_blocks );
    }

    typename UInt32View::HostMirror host_num_allocated_blocks =
      create_mirror_view( num_allocated_blocks );
    deep_copy( host_num_allocated_blocks, num_allocated_blocks );

    // Print header info of all superblocks.
    printf( "SB_ID    SIZE    ACTIVE    EMPTY_PAGES    FULL_PAGES    USED_BLOCKS\n" );
    for ( size_t i = 0; i < m_num_sb; ++i ) {
      printf( "%5zu    %4u    %6d    %11u    %10u     %10u\n", i,
              host_sb_header(i).m_lg_block_size, host_sb_header(i).m_is_active,
              host_sb_header(i).m_empty_pages, host_sb_header(i).m_full_pages,
              host_num_allocated_blocks(i) );
    }

    printf( "\n" );
#endif

    UInt32View page_histogram( "Page Histogram", 33 );

    // Get a View version of the blocksize info.
    typedef View< BlockSizeHeader *, device_type >  BSHeaderView;
    BSHeaderView blocksize_info( "BlockSize Headers", MAX_BLOCK_SIZES );

    Kokkos::Impl::DeepCopy< backend_memory_space, Kokkos::HostSpace >
      dc( blocksize_info.ptr_on_device(), m_blocksize_info,
          sizeof(BlockSizeHeader) * m_num_block_size );

    Kokkos::pair< double, uint32_t > result = Kokkos::pair< double, uint32_t >( 0.0, 0 );

    // Create the page histogram.
    {
      MempoolImpl::create_histogram< UInt32View, BSHeaderView, SBHeaderView, MempoolBitset >
        mch( 0, m_num_sb, page_histogram, blocksize_info, m_sb_header, m_sb_blocks,
             m_lg_max_sb_blocks, LG_MIN_BLOCK_SIZE, BLOCKS_PER_PAGE, result );
    }

    typename UInt32View::HostMirror host_page_histogram = create_mirror_view( page_histogram );
    deep_copy( host_page_histogram, page_histogram );

    // Find the used and total pages and blocks.
    uint32_t used_pages = 0;
    uint32_t used_blocks = 0;
    for ( uint32_t i = 1; i < 33; ++i ) {
      used_pages += host_page_histogram(i);
      used_blocks += i * host_page_histogram(i);
    }
    uint32_t total_pages = used_pages + host_page_histogram(0);

    unsigned num_empty_sb = m_empty_sb.count();
    unsigned num_non_empty_sb = m_num_sb - num_empty_sb;
    unsigned num_partfull_sb = m_partfull_sb.count();

    uint32_t total_blocks = result.second;
    double ave_sb_full = num_non_empty_sb == 0 ? 0.0 : result.first / num_non_empty_sb;
    double percent_used_sb = double( m_num_sb - num_empty_sb ) / m_num_sb;
    double percent_used_pages = total_pages == 0 ? 0.0 : double(used_pages) / total_pages;
    double percent_used_blocks = total_blocks == 0 ? 0.0 : double(used_blocks) / total_blocks;

    // Count active superblocks.
    typename UInt32View::HostMirror host_active = create_mirror_view( m_active );
    deep_copy( host_active, m_active );

    unsigned num_active_sb = 0;
    for ( size_t i = 0; i < m_num_block_size; ++i ) {
      num_active_sb += host_active(i) != INVALID_SUPERBLOCK;
    }

#ifdef KOKKOS_ENABLE_MEMPOOL_PRINT_ACTIVE_SUPERBLOCKS
    // Print active superblocks.
    printf( "BS_ID      SB_ID\n" );
    for ( size_t i = 0; i < m_num_block_size; ++i ) {
      uint32_t sb_id = host_active(i);

      if ( sb_id == INVALID_SUPERBLOCK ) {
        printf( "%5zu          I\n", i );
      }
      else if ( sb_id == SUPERBLOCK_LOCK ) {
        printf( "%5zu          L\n", i );
      }
      else {
        printf( "%5zu    %7u\n", i, sb_id );
      }
    }
    printf( "\n" );
    fflush( stdout );
#endif

#ifdef KOKKOS_ENABLE_MEMPOOL_PRINT_PAGE_INFO
    // Print the summary page histogram.
    printf( "USED_BLOCKS    PAGE_COUNT\n" );
    for ( uint32_t i = 0; i < 33; ++i ) {
      printf( "%10u    %10u\n", i, host_page_histogram[i] );
    }
    printf( "\n" );
#endif

#ifdef KOKKOS_ENABLE_MEMPOOL_PRINT_INDIVIDUAL_PAGE_INFO
    // Print the page histogram for a few individual superblocks.
//    const uint32_t num_sb_id = 2;
//    uint32_t sb_id[num_sb_id] = { 0, 10 };
    const uint32_t num_sb_id = 1;
    uint32_t sb_id[num_sb_id] = { 0 };

    for ( uint32_t i = 0; i < num_sb_id; ++i ) {
      deep_copy( page_histogram, 0 );

      {
        MempoolImpl::create_histogram< UInt32View, BSHeaderView, SBHeaderView, MempoolBitset >
          mch( sb_id[i], sb_id[i] + 1, page_histogram, blocksize_info, m_sb_header,
               m_sb_blocks, m_lg_max_sb_blocks, LG_MIN_BLOCK_SIZE, BLOCKS_PER_PAGE, result );
      }

      deep_copy( host_page_histogram, page_histogram );

      printf( "SB_ID    USED_BLOCKS    PAGE_COUNT\n" );
      for ( uint32_t j = 0; j < 33; ++j ) {
        printf( "%5u    %10u    %10u\n", sb_id[i], j, host_page_histogram[j] );
      }
      printf( "\n" );
    }

/*
    // Print the blocks used for each page of a few individual superblocks.
    for ( uint32_t i = 0; i < num_sb_id; ++i ) {
      uint32_t lg_block_size = host_sb_header(sb_id[i]).m_lg_block_size;

      if ( lg_block_size != 0 ) {
        printf( "SB_ID    BLOCK ID    USED_BLOCKS\n" );

        uint32_t block_size_id = lg_block_size - LG_MIN_BLOCK_SIZE;
        uint32_t pages_per_sb = m_blocksize_info[block_size_id].m_pages_per_sb;

        for ( uint32_t j = 0; j < pages_per_sb; ++j ) {
          unsigned start_pos = ( sb_id[i] << m_lg_max_sb_blocks ) + j * BLOCKS_PER_PAGE;
          unsigned end_pos = start_pos + BLOCKS_PER_PAGE;
          uint32_t num_allocated_blocks = 0;

          for ( unsigned k = start_pos; k < end_pos; ++k ) {
            num_allocated_blocks += m_sb_blocks.test( k );
          }

          printf( "%5u    %8u    %11u\n", sb_id[i], j, num_allocated_blocks );
        }

        printf( "\n" );
      }
    }
*/
#endif

    printf( "   Used blocks: %10u / %10u = %10.6lf\n", used_blocks, total_blocks,
            percent_used_blocks );
    printf( "    Used pages: %10u / %10u = %10.6lf\n", used_pages, total_pages,
            percent_used_pages );
    printf( "       Used SB: %10zu / %10zu = %10.6lf\n", m_num_sb - num_empty_sb, m_num_sb,
            percent_used_sb );
    printf( "     Active SB: %10u\n", num_active_sb );
    printf( "      Empty SB: %10u\n", num_empty_sb );
    printf( "   Partfull SB: %10u\n", num_partfull_sb );
    printf( "       Full SB: %10lu\n",
            m_num_sb - num_active_sb - num_empty_sb - num_partfull_sb );
    printf( "Ave. SB Full %%: %10.6lf\n", ave_sb_full );
    printf( "\n" );
    fflush( stdout );

#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    fflush( stdout );
#endif
  }

  KOKKOS_INLINE_FUNCTION
  size_t get_min_block_size() const { return MIN_BLOCK_SIZE; }

  KOKKOS_INLINE_FUNCTION
  size_t get_mem_size() const { return m_data_size; }

private:
  /// \brief Returns the index into the active array for the given size.
  ///
  /// Computes log2 of the largest power of two >= the given size
  /// ( ie ceil( log2(size) ) ) shifted by LG_MIN_BLOCK_SIZE.
  KOKKOS_FORCEINLINE_FUNCTION
  int get_block_size_index( const size_t size ) const
  {
    // We know the size fits in a 32 bit unsigned because the size of a
    // superblock is limited to 2^31, so casting to an unsigned is safe.

    // Find the most significant nonzero bit.
    uint32_t first_nonzero_bit =
      Kokkos::Impl::bit_scan_reverse( static_cast<unsigned>( size ) );

    // If size is an integral power of 2, ceil( log2(size) ) is equal to the
    // most significant nonzero bit.  Otherwise, you need to add 1.  Since the
    // minimum block size is MIN_BLOCK_SIZE, make sure ceil( log2(size) ) is at
    // least LG_MIN_BLOCK_SIZE.
    uint32_t lg2_size = first_nonzero_bit + !Kokkos::Impl::is_integral_power_of_two( size );
    lg2_size = lg2_size > LG_MIN_BLOCK_SIZE ? lg2_size : LG_MIN_BLOCK_SIZE;

    // Return ceil( log2(size) ) shifted so that the value for MIN_BLOCK_SIZE
    // is 0.
    return lg2_size - LG_MIN_BLOCK_SIZE;
  }

  /// \brief Finds a superblock with free space to become a new active superblock.
  ///
  /// If this function is called, the current active superblock needs to be replaced
  /// because it is full.  Initially, only the thread that sets the active superblock
  /// to full calls this function.  Other threads can still allocate from the "full"
  /// active superblock because a full superblock still has locations available.  If
  /// a thread tries to allocate from the active superblock when it has no free
  /// locations, then that thread will call this function, too, and spin on a lock
  /// waiting until the active superblock has been replaced.
  KOKKOS_FUNCTION
  uint32_t find_superblock( int block_size_id, uint32_t old_sb ) const
  {
    // Try to grab the lock on the head.
    uint32_t lock_sb =
      Kokkos::atomic_compare_exchange( &m_active(block_size_id), old_sb, SUPERBLOCK_LOCK );

    load_fence();

    // Initialize the new superblock to be the previous one so the previous
    // superblock is returned if a new superblock can't be found.
    uint32_t new_sb = lock_sb;

    if ( lock_sb == old_sb ) {
      // This thread has the lock.

      // 1. Look for a partially filled superblock that is of the right block
      //    size.

      size_t max_tries = m_ceil_num_sb >> LG_BLOCKS_PER_PAGE;
      size_t tries = 0;
      bool search_done = false;

      // Set the starting search position to the beginning of this block
      // size's bitset.
      unsigned pos = block_size_id * m_ceil_num_sb;

      while ( !search_done ) {
        bool success = false;
        unsigned prev_val = 0;

        Kokkos::tie( success, prev_val ) = m_partfull_sb.reset_any_in_word( pos );

        if ( !success ) {
          if ( ++tries >= max_tries ) {
            // Exceeded number of words for this block size's bitset.
            search_done = true;
          }
          else {
            pos += BLOCKS_PER_PAGE;
          }
        }
        else {
          // Found a superblock.

          // It is possible that the newly found superblock is the same as the
          // old superblock.  In this case putting the old value back in yields
          // correct behavior.  This could happen as follows.  This thread
          // grabs the lock and transitions the superblock to the full state.
          // Before it searches for a new superblock, other threads perform
          // enough deallocations to transition the superblock to the partially
          // full state.  This thread then searches for a partially full
          // superblock and finds the one it removed.  There's potential for
          // this to cause a performance issue if the same superblock keeps
          // being removed and added due to the right mix and ordering of
          // allocations and deallocations.
          search_done = true;
          new_sb = pos - block_size_id * m_ceil_num_sb;

          // Set the head status for the superblock.
          volatile_store( &m_sb_header(new_sb).m_is_active, uint32_t(true) );

          // If there was a previous active superblock, mark it as not active.
          // It is now in the full category and as such isn't tracked.
          if ( lock_sb != INVALID_SUPERBLOCK ) {
            volatile_store( &m_sb_header(lock_sb).m_is_active, uint32_t(false) );
          }

          store_fence();
        }
      }

      // 2. Look for an empty superblock.
      if ( new_sb == lock_sb ) {
        tries = 0;
        search_done = false;

        // Set the starting search position to the beginning of this block
        // size's bitset.
        pos = 0;

        while ( !search_done ) {
          bool success = false;
          unsigned prev_val = 0;

          Kokkos::tie( success, prev_val ) = m_empty_sb.reset_any_in_word( pos );

          if ( !success ) {
            if ( ++tries >= max_tries ) {
              // Exceeded number of words for this block size's bitset.
              search_done = true;
            }
            else {
              pos += BLOCKS_PER_PAGE;
            }
          }
          else {
            // Found a superblock.

            // It is possible that the newly found superblock is the same as
            // the old superblock.  In this case putting the old value back in
            // yields correct behavior.  This could happen as follows.  This
            // thread grabs the lock and transitions the superblock to the full
            // state.  Before it searches for a new superblock, other threads
            // perform enough deallocations to transition the superblock to the
            // partially full state and then the empty state.  This thread then
            // searches for a partially full superblock and none exist.  This
            // thread then searches for an empty superblock and finds the one
            // it removed.  The likelihood of this happening is so remote that
            // the potential for this to cause a performance issue is
            // infinitesimal.
            search_done = true;
            new_sb = pos;

            // Set the empty pages, block size, and head status for the
            // superblock.
            volatile_store( &m_sb_header(new_sb).m_empty_pages,
                            m_blocksize_info[block_size_id].m_pages_per_sb );
            volatile_store( &m_sb_header(new_sb).m_lg_block_size,
                            block_size_id + LG_MIN_BLOCK_SIZE );
            volatile_store( &m_sb_header(new_sb).m_is_active, uint32_t(true) );

            // If there was a previous active superblock, mark it as not active.
            // It is now in the full category and as such isn't tracked.
            if ( lock_sb != INVALID_SUPERBLOCK ) {
              volatile_store( &m_sb_header(lock_sb).m_is_active, uint32_t(false) );
            }

            store_fence();
          }
        }
      }

      // Write the new active superblock to release the lock.
      atomic_exchange( &m_active(block_size_id), new_sb );
    }
    else {
      // Either another thread has the lock and is switching the active
      // superblock for this block size or another thread has already changed
      // the active superblock since this thread read its value.  Keep
      // atomically reading the active superblock until it isn't locked to get
      // the new active superblock.
      do {
        new_sb = atomic_fetch_or( &m_active(block_size_id), uint32_t(0) );
      } while ( new_sb == SUPERBLOCK_LOCK );

      load_fence();

      // Assertions:
      //   1. An invalid superblock should never be found here.
      //   2. If the new superblock is the same as the previous superblock, the
      //      allocator is empty.
#ifdef KOKKOS_ENABLE_MEMPOOL_PRINTERR
      if ( new_sb == INVALID_SUPERBLOCK ) {
        printf( "\n** MemoryPool::find_superblock() FOUND_INACTIVE_SUPERBLOCK **\n" );
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
        fflush( stdout );
#endif
        Kokkos::abort( "" );
      }
#endif
    }

    return new_sb;
  }

  /// Returns 64 bits from a clock register.
  KOKKOS_FORCEINLINE_FUNCTION
  uint64_t get_clock_register(void) const
  {
#if defined( __CUDA_ARCH__ )
    // Return value of 64-bit hi-res clock register.
    return clock64();
#elif defined( __i386__ ) || defined( __x86_64 )
    // Return value of 64-bit hi-res clock register.
    unsigned a = 0, d = 0;

    __asm__ volatile( "rdtsc" : "=a" (a), "=d" (d) );

    return ( (uint64_t) a ) | ( ( (uint64_t) d ) << 32 );
#elif defined( __powerpc )   || defined( __powerpc__ ) || defined( __powerpc64__ ) || \
      defined( __POWERPC__ ) || defined( __ppc__ )     || defined( __ppc64__ )
  unsigned int cycles = 0;

  asm volatile( "mftb %0" : "=r" (cycles) );

  return (uint64_t) cycles;
#else
    const uint64_t ticks =
      std::chrono::high_resolution_clock::now().time_since_epoch().count();

    return ticks;
#endif
  }
};

} // namespace Experimental
} // namespace Kokkos

#ifdef KOKKOS_ENABLE_MEMPOOL_PRINTERR
#undef KOKKOS_ENABLE_MEMPOOL_PRINTERR
#endif

#ifdef KOKKOS_ENABLE_MEMPOOL_PRINT_INFO
#undef KOKKOS_ENABLE_MEMPOOL_PRINT_INFO
#endif

#ifdef KOKKOS_ENABLE_MEMPOOL_PRINT_BLOCKSIZE_INFO
#undef KOKKOS_ENABLE_MEMPOOL_PRINT_BLOCKSIZE_INFO
#endif

#ifdef KOKKOS_ENABLE_MEMPOOL_PRINT_SUPERBLOCK_INFO
#undef KOKKOS_ENABLE_MEMPOOL_PRINT_SUPERBLOCK_INFO
#endif

#ifdef KOKKOS_ENABLE_MEMPOOL_PRINT_PAGE_INFO
#undef KOKKOS_ENABLE_MEMPOOL_PRINT_PAGE_INFO
#endif

#ifdef KOKKOS_ENABLE_MEMPOOL_PRINT_INDIVIDUAL_PAGE_INFO
#undef KOKKOS_ENABLE_MEMPOOL_PRINT_INDIVIDUAL_PAGE_INFO
#endif

#endif // KOKKOS_MEMORYPOOL_HPP
