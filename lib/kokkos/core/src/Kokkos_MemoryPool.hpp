//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif
#ifndef KOKKOS_MEMORYPOOL_HPP
#define KOKKOS_MEMORYPOOL_HPP

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_Atomic.hpp>
#include <impl/Kokkos_ConcurrentBitset.hpp>
#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_SharedAlloc.hpp>

namespace Kokkos {
namespace Impl {
/* Report violation of size constraints:
 *   min_block_alloc_size <= max_block_alloc_size
 *   max_block_alloc_size <= min_superblock_size
 *   min_superblock_size  <= max_superblock_size
 *   min_superblock_size  <= min_total_alloc_size
 *   min_superblock_size  <= min_block_alloc_size *
 *                           max_block_per_superblock
 */
void memory_pool_bounds_verification(size_t min_block_alloc_size,
                                     size_t max_block_alloc_size,
                                     size_t min_superblock_size,
                                     size_t max_superblock_size,
                                     size_t max_block_per_superblock,
                                     size_t min_total_alloc_size);
}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos {

namespace Impl {

void _print_memory_pool_state(std::ostream &s, uint32_t const *sb_state_ptr,
                              int32_t sb_count, uint32_t sb_size_lg2,
                              uint32_t sb_state_size, uint32_t state_shift,
                              uint32_t state_used_mask);

}  // end namespace Impl

template <typename DeviceType>
class MemoryPool {
 private:
  using CB = Kokkos::Impl::concurrent_bitset;

  enum : uint32_t { bits_per_int_lg2 = CB::bits_per_int_lg2 };
  enum : uint32_t { state_shift = CB::state_shift };
  enum : uint32_t { state_used_mask = CB::state_used_mask };
  enum : uint32_t { state_header_mask = CB::state_header_mask };
  enum : uint32_t { max_bit_count_lg2 = CB::max_bit_count_lg2 };
  enum : uint32_t { max_bit_count = CB::max_bit_count };

  enum : uint32_t { HINT_PER_BLOCK_SIZE = 2 };

  /*  Each superblock has a concurrent bitset state
   *  which is an array of uint32_t integers.
   *    [ { block_count_lg2  : state_shift bits
   *      , used_block_count : ( 32 - state_shift ) bits
   *      }
   *    , { block allocation bit set }* ]
   *
   *  As superblocks are assigned (allocated) to a block size
   *  and released (deallocated) back to empty the superblock state
   *  is concurrently updated.
   */

  /*  Mapping between block_size <-> block_state
   *
   *  block_state = ( m_sb_size_lg2 - block_size_lg2 ) << state_shift
   *  block_size  = m_sb_size_lg2 - ( block_state >> state_shift )
   *
   *  Thus A_block_size < B_block_size  <=>  A_block_state > B_block_state
   */

  using base_memory_space = typename DeviceType::memory_space;

  enum {
    accessible = Kokkos::Impl::MemorySpaceAccess<Kokkos::HostSpace,
                                                 base_memory_space>::accessible
  };

  using Tracker = Kokkos::Impl::SharedAllocationTracker;
  using Record  = Kokkos::Impl::SharedAllocationRecord<base_memory_space>;

  Tracker m_tracker;
  uint32_t *m_sb_state_array;
  uint32_t m_sb_state_size;
  uint32_t m_sb_size_lg2;
  uint32_t m_max_block_size_lg2;
  uint32_t m_min_block_size_lg2;
  int32_t m_sb_count;
  int32_t m_hint_offset;  // Offset to K * #block_size array of hints
  int32_t m_data_offset;  // Offset to 0th superblock data
  int32_t m_unused_padding;

 public:
  using memory_space = typename DeviceType::memory_space;

  /**\brief  The maximum size of a superblock and block */
  enum : uint32_t { max_superblock_size = 1LU << 31 /* 2 gigabytes */ };
  enum : uint32_t { max_block_per_superblock = max_bit_count };

  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  bool operator==(MemoryPool const &other) const {
    return m_sb_state_array == other.m_sb_state_array;
  }

  KOKKOS_INLINE_FUNCTION
  size_t capacity() const noexcept {
    return size_t(m_sb_count) << m_sb_size_lg2;
  }

  KOKKOS_INLINE_FUNCTION
  size_t min_block_size() const noexcept {
    return (1LU << m_min_block_size_lg2);
  }

  KOKKOS_INLINE_FUNCTION
  size_t max_block_size() const noexcept {
    return (1LU << m_max_block_size_lg2);
  }

  struct usage_statistics {
    size_t capacity_bytes;        ///<  Capacity in bytes
    size_t superblock_bytes;      ///<  Superblock size in bytes
    size_t max_block_bytes;       ///<  Maximum block size in bytes
    size_t min_block_bytes;       ///<  Minimum block size in bytes
    size_t capacity_superblocks;  ///<  Number of superblocks
    size_t consumed_superblocks;  ///<  Superblocks assigned to allocations
    size_t consumed_blocks;       ///<  Number of allocations
    size_t consumed_bytes;        ///<  Bytes allocated
    size_t reserved_blocks;  ///<  Unallocated blocks in assigned superblocks
    size_t reserved_bytes;   ///<  Unallocated bytes in assigned superblocks
  };

  void get_usage_statistics(usage_statistics &stats) const {
    Kokkos::HostSpace host;

    const size_t alloc_size = m_hint_offset * sizeof(uint32_t);

    uint32_t *const sb_state_array =
        accessible ? m_sb_state_array : (uint32_t *)host.allocate(alloc_size);

    if (!accessible) {
      Kokkos::Impl::DeepCopy<Kokkos::HostSpace, base_memory_space>(
          sb_state_array, m_sb_state_array, alloc_size);
      Kokkos::fence(
          "MemoryPool::get_usage_statistics(): fence after copying state "
          "array to HostSpace");
    }

    stats.superblock_bytes     = (1LU << m_sb_size_lg2);
    stats.max_block_bytes      = (1LU << m_max_block_size_lg2);
    stats.min_block_bytes      = (1LU << m_min_block_size_lg2);
    stats.capacity_bytes       = stats.superblock_bytes * m_sb_count;
    stats.capacity_superblocks = m_sb_count;
    stats.consumed_superblocks = 0;
    stats.consumed_blocks      = 0;
    stats.consumed_bytes       = 0;
    stats.reserved_blocks      = 0;
    stats.reserved_bytes       = 0;

    const uint32_t *sb_state_ptr = sb_state_array;

    for (int32_t i = 0; i < m_sb_count; ++i, sb_state_ptr += m_sb_state_size) {
      const uint32_t block_count_lg2 = (*sb_state_ptr) >> state_shift;

      if (block_count_lg2) {
        const uint32_t block_count    = 1u << block_count_lg2;
        const uint32_t block_size_lg2 = m_sb_size_lg2 - block_count_lg2;
        const uint32_t block_size     = 1u << block_size_lg2;
        const uint32_t block_used     = (*sb_state_ptr) & state_used_mask;

        stats.consumed_superblocks++;
        stats.consumed_blocks += block_used;
        stats.consumed_bytes += block_used * block_size;
        stats.reserved_blocks += block_count - block_used;
        stats.reserved_bytes += (block_count - block_used) * block_size;
      }
    }

    if (!accessible) {
      host.deallocate(sb_state_array, alloc_size);
    }
  }

  void print_state(std::ostream &s) const {
    Kokkos::HostSpace host;

    const size_t alloc_size = m_hint_offset * sizeof(uint32_t);

    uint32_t *const sb_state_array =
        accessible ? m_sb_state_array : (uint32_t *)host.allocate(alloc_size);

    if (!accessible) {
      Kokkos::Impl::DeepCopy<Kokkos::HostSpace, base_memory_space>(
          sb_state_array, m_sb_state_array, alloc_size);
      Kokkos::fence(
          "MemoryPool::print_state(): fence after copying state array to "
          "HostSpace");
    }

    Impl::_print_memory_pool_state(s, sb_state_array, m_sb_count, m_sb_size_lg2,
                                   m_sb_state_size, state_shift,
                                   state_used_mask);

    if (!accessible) {
      host.deallocate(sb_state_array, alloc_size);
    }
  }

  //--------------------------------------------------------------------------

  KOKKOS_DEFAULTED_FUNCTION MemoryPool(MemoryPool &&)      = default;
  KOKKOS_DEFAULTED_FUNCTION MemoryPool(const MemoryPool &) = default;
  KOKKOS_DEFAULTED_FUNCTION MemoryPool &operator=(MemoryPool &&) = default;
  KOKKOS_DEFAULTED_FUNCTION MemoryPool &operator=(const MemoryPool &) = default;

  KOKKOS_INLINE_FUNCTION MemoryPool()
      : m_tracker(),
        m_sb_state_array(nullptr),
        m_sb_state_size(0),
        m_sb_size_lg2(0),
        m_max_block_size_lg2(0),
        m_min_block_size_lg2(0),
        m_sb_count(0),
        m_hint_offset(0),
        m_data_offset(0),
        m_unused_padding(0) {}

  /**\brief  Allocate a memory pool from 'memspace'.
   *
   *  The memory pool will have at least 'min_total_alloc_size' bytes
   *  of memory to allocate divided among superblocks of at least
   *  'min_superblock_size' bytes.  A single allocation must fit
   *  within a single superblock, so 'min_superblock_size' must be
   *  at least as large as the maximum single allocation.
   *  Both 'min_total_alloc_size' and 'min_superblock_size'
   *  are rounded up to the smallest power-of-two value that
   *  contains the corresponding sizes.
   *  Individual allocations will always consume a block of memory that
   *  is also a power-of-two.  These roundings are made to enable
   *  significant runtime performance improvements.
   */
  MemoryPool(const base_memory_space &memspace,
             const size_t min_total_alloc_size, size_t min_block_alloc_size = 0,
             size_t max_block_alloc_size = 0, size_t min_superblock_size = 0)
      : m_tracker(),
        m_sb_state_array(nullptr),
        m_sb_state_size(0),
        m_sb_size_lg2(0),
        m_max_block_size_lg2(0),
        m_min_block_size_lg2(0),
        m_sb_count(0),
        m_hint_offset(0),
        m_data_offset(0),
        m_unused_padding(0) {
    const uint32_t int_align_lg2               = 3; /* align as int[8] */
    const uint32_t int_align_mask              = (1u << int_align_lg2) - 1;
    const uint32_t default_min_block_size      = 1u << 6;  /* 64 bytes */
    const uint32_t default_max_block_size      = 1u << 12; /* 4k bytes */
    const uint32_t default_min_superblock_size = 1u << 20; /* 1M bytes */

    //--------------------------------------------------
    // Default block and superblock sizes:

    if (0 == min_block_alloc_size) {
      // Default all sizes:

      min_superblock_size =
          std::min(size_t(default_min_superblock_size), min_total_alloc_size);

      min_block_alloc_size =
          std::min(size_t(default_min_block_size), min_superblock_size);

      max_block_alloc_size =
          std::min(size_t(default_max_block_size), min_superblock_size);
    } else if (0 == min_superblock_size) {
      // Choose superblock size as minimum of:
      //   max_block_per_superblock * min_block_size
      //   max_superblock_size
      //   min_total_alloc_size

      const size_t max_superblock =
          min_block_alloc_size * max_block_per_superblock;

      min_superblock_size =
          std::min(max_superblock,
                   std::min(size_t(max_superblock_size), min_total_alloc_size));
    }

    if (0 == max_block_alloc_size) {
      max_block_alloc_size = min_superblock_size;
    }

    //--------------------------------------------------

    /* Enforce size constraints:
     *   min_block_alloc_size <= max_block_alloc_size
     *   max_block_alloc_size <= min_superblock_size
     *   min_superblock_size  <= max_superblock_size
     *   min_superblock_size  <= min_total_alloc_size
     *   min_superblock_size  <= min_block_alloc_size *
     *                           max_block_per_superblock
     */

    Kokkos::Impl::memory_pool_bounds_verification(
        min_block_alloc_size, max_block_alloc_size, min_superblock_size,
        max_superblock_size, max_block_per_superblock, min_total_alloc_size);

    //--------------------------------------------------
    // Block and superblock size is power of two:
    // Maximum value is 'max_superblock_size'

    m_min_block_size_lg2 =
        Kokkos::Impl::integral_power_of_two_that_contains(min_block_alloc_size);

    m_max_block_size_lg2 =
        Kokkos::Impl::integral_power_of_two_that_contains(max_block_alloc_size);

    m_sb_size_lg2 =
        Kokkos::Impl::integral_power_of_two_that_contains(min_superblock_size);

    {
      // number of superblocks is multiple of superblock size that
      // can hold min_total_alloc_size.

      const uint64_t sb_size_mask = (1LU << m_sb_size_lg2) - 1;

      m_sb_count = (min_total_alloc_size + sb_size_mask) >> m_sb_size_lg2;
    }

    {
      // Any superblock can be assigned to the smallest size block
      // Size the block bitset to maximum number of blocks

      const uint32_t max_block_count_lg2 = m_sb_size_lg2 - m_min_block_size_lg2;

      m_sb_state_size =
          (CB::buffer_bound_lg2(max_block_count_lg2) + int_align_mask) &
          ~int_align_mask;
    }

    // Array of all superblock states

    const size_t all_sb_state_size =
        (m_sb_count * m_sb_state_size + int_align_mask) & ~int_align_mask;

    // Number of block sizes

    const int32_t number_block_sizes =
        1 + m_max_block_size_lg2 - m_min_block_size_lg2;

    // Array length for possible block sizes
    // Hint array is one uint32_t per block size

    const int32_t block_size_array_size =
        (number_block_sizes + int_align_mask) & ~int_align_mask;

    m_hint_offset = all_sb_state_size;
    m_data_offset = m_hint_offset + block_size_array_size * HINT_PER_BLOCK_SIZE;

    // Allocation:

    const size_t header_size = m_data_offset * sizeof(uint32_t);
    const size_t alloc_size =
        header_size + (size_t(m_sb_count) << m_sb_size_lg2);

    Record *rec = Record::allocate(memspace, "Kokkos::MemoryPool", alloc_size);

    m_tracker.assign_allocated_record_to_uninitialized(rec);

    m_sb_state_array = (uint32_t *)rec->data();

    Kokkos::HostSpace host;

    uint32_t *const sb_state_array =
        accessible ? m_sb_state_array : (uint32_t *)host.allocate(header_size);

    for (int32_t i = 0; i < m_data_offset; ++i) sb_state_array[i] = 0;

    // Initial assignment of empty superblocks to block sizes:

    for (int32_t i = 0; i < number_block_sizes; ++i) {
      const uint32_t block_size_lg2  = i + m_min_block_size_lg2;
      const uint32_t block_count_lg2 = m_sb_size_lg2 - block_size_lg2;
      const uint32_t block_state     = block_count_lg2 << state_shift;
      const uint32_t hint_begin      = m_hint_offset + i * HINT_PER_BLOCK_SIZE;

      // for block size index 'i':
      //   sb_id_hint  = sb_state_array[ hint_begin ];
      //   sb_id_begin = sb_state_array[ hint_begin + 1 ];

      const int32_t jbeg = (i * m_sb_count) / number_block_sizes;
      const int32_t jend = ((i + 1) * m_sb_count) / number_block_sizes;

      sb_state_array[hint_begin]     = uint32_t(jbeg);
      sb_state_array[hint_begin + 1] = uint32_t(jbeg);

      for (int32_t j = jbeg; j < jend; ++j) {
        sb_state_array[j * m_sb_state_size] = block_state;
      }
    }

    // Write out initialized state:

    if (!accessible) {
      Kokkos::Impl::DeepCopy<base_memory_space, Kokkos::HostSpace>(
          m_sb_state_array, sb_state_array, header_size);
      Kokkos::fence(
          "MemoryPool::MemoryPool(): fence after copying state array from "
          "HostSpace");

      host.deallocate(sb_state_array, header_size);
    } else {
      Kokkos::memory_fence();
    }
  }

  //--------------------------------------------------------------------------

 private:
  /* Given a size 'n' get the block size in which it can be allocated.
   * Restrict lower bound to minimum block size.
   */
  KOKKOS_FORCEINLINE_FUNCTION
  uint32_t get_block_size_lg2(uint32_t n) const noexcept {
    const unsigned i = Kokkos::Impl::integral_power_of_two_that_contains(n);

    return i < m_min_block_size_lg2 ? m_min_block_size_lg2 : i;
  }

 public:
  /* Return 0 for invalid block size */
  KOKKOS_INLINE_FUNCTION
  uint32_t allocate_block_size(uint64_t alloc_size) const noexcept {
    return alloc_size <= (1UL << m_max_block_size_lg2)
               ? (1UL << get_block_size_lg2(uint32_t(alloc_size)))
               : 0;
  }

  //--------------------------------------------------------------------------
  /**\brief  Allocate a block of memory that is at least 'alloc_size'
   *
   *  The block of memory is aligned to the minimum block size,
   *  currently is 64 bytes, will never be less than 32 bytes.
   *
   *  If concurrent allocations and deallocations are taking place
   *  then a single allocation attempt may fail due to lack of available space.
   *  The allocation attempt will try up to 'attempt_limit' times.
   */
  KOKKOS_FUNCTION
  void *allocate(size_t alloc_size, int32_t attempt_limit = 1) const noexcept {
    if (size_t(1LU << m_max_block_size_lg2) < alloc_size) {
      Kokkos::abort(
          "Kokkos MemoryPool allocation request exceeded specified maximum "
          "allocation size");
    }

    if (0 == alloc_size) return nullptr;

    void *p = nullptr;

    const uint32_t block_size_lg2 = get_block_size_lg2(alloc_size);

    // Allocation will fit within a superblock
    // that has block sizes ( 1 << block_size_lg2 )

    const uint32_t block_count_lg2 = m_sb_size_lg2 - block_size_lg2;
    const uint32_t block_state     = block_count_lg2 << state_shift;
    const uint32_t block_count     = 1u << block_count_lg2;

    // Superblock hints for this block size:
    //   hint_sb_id_ptr[0] is the dynamically changing hint
    //   hint_sb_id_ptr[1] is the static start point

    volatile uint32_t *const hint_sb_id_ptr =
        m_sb_state_array      /* memory pool state array */
        + m_hint_offset       /* offset to hint portion of array */
        + HINT_PER_BLOCK_SIZE /* number of hints per block size */
              * (block_size_lg2 - m_min_block_size_lg2); /* block size id */

    const int32_t sb_id_begin = int32_t(hint_sb_id_ptr[1]);

    // Fast query clock register 'tic' to pseudo-randomize
    // the guess for which block within a superblock should
    // be claimed.  If not available then a search occurs.
#if defined(KOKKOS_ENABLE_SYCL) && !defined(KOKKOS_ARCH_INTEL_GPU)
    const uint32_t block_id_hint = alloc_size;
#else
    const uint32_t block_id_hint =
        (uint32_t)(Kokkos::Impl::clock_tic()
#ifdef __CUDA_ARCH__  // FIXME_CUDA
                   // Spread out potentially concurrent access
                   // by threads within a warp or thread block.
                   + (threadIdx.x + blockDim.x * threadIdx.y)
#endif
        );
#endif

    // expected state of superblock for allocation
    uint32_t sb_state = block_state;

    int32_t sb_id = -1;

    volatile uint32_t *sb_state_array = nullptr;

    while (attempt_limit) {
      int32_t hint_sb_id = -1;

      if (sb_id < 0) {
        // No superblock specified, try the hint for this block size

        sb_id = hint_sb_id = int32_t(*hint_sb_id_ptr);

        sb_state_array = m_sb_state_array + (sb_id * m_sb_state_size);
      }

      // Require:
      //   0 <= sb_id
      //   sb_state_array == m_sb_state_array + m_sb_state_size * sb_id

      if (sb_state == (state_header_mask & *sb_state_array)) {
        // This superblock state is as expected, for the moment.
        // Attempt to claim a bit.  The attempt updates the state
        // so have already made sure the state header is as expected.

        const uint32_t count_lg2 = sb_state >> state_shift;
        const uint32_t mask      = (1u << count_lg2) - 1;

        const Kokkos::pair<int, int> result = CB::acquire_bounded_lg2(
            sb_state_array, count_lg2, block_id_hint & mask, sb_state);

        // If result.first < 0 then failed to acquire
        // due to either full or buffer was wrong state.
        // Could be wrong state if a deallocation raced the
        // superblock to empty before the acquire could succeed.

        if (0 <= result.first) {  // acquired a bit

          const uint32_t size_lg2 = m_sb_size_lg2 - count_lg2;

          // Set the allocated block pointer

          p = ((char *)(m_sb_state_array + m_data_offset)) +
              (uint64_t(sb_id) << m_sb_size_lg2)       // superblock memory
              + (uint64_t(result.first) << size_lg2);  // block memory

          break;  // Success
        }
      }
      //------------------------------------------------------------------
      //  Arrive here if failed to acquire a block.
      //  Must find a new superblock.

      //  Start searching at designated index for this block size.
      //  Look for superblock that, in preferential order,
      //  1) part-full superblock of this block size
      //  2) empty superblock to claim for this block size
      //  3) part-full superblock of the next larger block size

      sb_state = block_state;  // Expect to find the desired state
      sb_id    = -1;

      bool update_hint        = false;
      int32_t sb_id_empty     = -1;
      int32_t sb_id_large     = -1;
      uint32_t sb_state_large = 0;

      sb_state_array = m_sb_state_array + sb_id_begin * m_sb_state_size;

      for (int32_t i = 0, id = sb_id_begin; i < m_sb_count; ++i) {
        //  Query state of the candidate superblock.
        //  Note that the state may change at any moment
        //  as concurrent allocations and deallocations occur.

        const uint32_t full_state = *sb_state_array;
        const uint32_t used       = full_state & state_used_mask;
        const uint32_t state      = full_state & state_header_mask;

        if (state == block_state) {
          //  Superblock is assigned to this block size

          if (used < block_count) {
            // There is room to allocate one block

            sb_id = id;

            // Is there room to allocate more than one block?

            update_hint = used + 1 < block_count;

            break;
          }
        } else if (0 == used) {
          // Superblock is empty

          if (-1 == sb_id_empty) {
            // Superblock is not assigned to this block size
            // and is the first empty superblock encountered.
            // Save this id to use if a partfull superblock is not found.

            sb_id_empty = id;
          }
        } else if ((-1 == sb_id_empty /* have not found an empty */) &&
                   (-1 == sb_id_large /* have not found a larger */) &&
                   (state < block_state /* a larger block */) &&
                   // is not full:
                   (used < (1u << (state >> state_shift)))) {
          //  First superblock encountered that is
          //  larger than this block size and
          //  has room for an allocation.
          //  Save this id to use of partfull or empty superblock not found
          sb_id_large    = id;
          sb_state_large = state;
        }

        // Iterate around the superblock array:

        if (++id < m_sb_count) {
          sb_state_array += m_sb_state_size;
        } else {
          id             = 0;
          sb_state_array = m_sb_state_array;
        }
      }

      // printf("  search m_sb_count(%d) sb_id(%d) sb_id_empty(%d)
      // sb_id_large(%d)\n" , m_sb_count , sb_id , sb_id_empty , sb_id_large);

      if (sb_id < 0) {
        //  Did not find a partfull superblock for this block size.

        if (0 <= sb_id_empty) {
          //  Found first empty superblock following designated superblock
          //  Attempt to claim it for this block size.
          //  If the claim fails assume that another thread claimed it
          //  for this block size and try to use it anyway,
          //  but do not update hint.

          sb_id = sb_id_empty;

          sb_state_array = m_sb_state_array + (sb_id * m_sb_state_size);

          //  If successfully changed assignment of empty superblock 'sb_id'
          //  to this block_size then update the hint.

          const uint32_t state_empty = state_header_mask & *sb_state_array;

          // If this thread claims the empty block then update the hint
          update_hint =
              state_empty == Kokkos::atomic_compare_exchange(
                                 sb_state_array, state_empty, block_state);
        } else if (0 <= sb_id_large) {
          // Found a larger superblock with space available

          sb_id    = sb_id_large;
          sb_state = sb_state_large;

          sb_state_array = m_sb_state_array + (sb_id * m_sb_state_size);
        } else {
          // Did not find a potentially usable superblock
          --attempt_limit;
        }
      }

      if (update_hint) {
        Kokkos::atomic_compare_exchange(hint_sb_id_ptr, uint32_t(hint_sb_id),
                                        uint32_t(sb_id));
      }
    }  // end allocation attempt loop
    //--------------------------------------------------------------------

    return p;
  }
  // end allocate
  //--------------------------------------------------------------------------

  /**\brief  Return an allocated block of memory to the pool.
   *
   *  Requires: p is return value from allocate( alloc_size );
   *
   *  For now the alloc_size is ignored.
   */
  KOKKOS_INLINE_FUNCTION
  void deallocate(void *p, size_t /* alloc_size */) const noexcept {
    if (nullptr == p) return;

    // Determine which superblock and block
    const ptrdiff_t d =
        static_cast<char *>(p) -
        reinterpret_cast<char *>(m_sb_state_array + m_data_offset);

    // Verify contained within the memory pool's superblocks:
    const int ok_contains =
        (0 <= d) && (size_t(d) < (size_t(m_sb_count) << m_sb_size_lg2));

    int ok_block_aligned = 0;
    int ok_dealloc_once  = 0;

    if (ok_contains) {
      const int sb_id = d >> m_sb_size_lg2;

      // State array for the superblock.
      volatile uint32_t *const sb_state_array =
          m_sb_state_array + (sb_id * m_sb_state_size);

      const uint32_t block_state = (*sb_state_array) & state_header_mask;
      const uint32_t block_size_lg2 =
          m_sb_size_lg2 - (block_state >> state_shift);

      ok_block_aligned = 0 == (d & ((1UL << block_size_lg2) - 1));

      if (ok_block_aligned) {
        // Map address to block's bit
        // mask into superblock and then shift down for block index

        const uint32_t bit =
            (d & (ptrdiff_t(1LU << m_sb_size_lg2) - 1)) >> block_size_lg2;

        const int result = CB::release(sb_state_array, bit, block_state);

        ok_dealloc_once = 0 <= result;
      }
    }

    if (!ok_contains || !ok_block_aligned || !ok_dealloc_once) {
      Kokkos::abort("Kokkos MemoryPool::deallocate given erroneous pointer");
    }
  }
  // end deallocate
  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  int number_of_superblocks() const noexcept { return m_sb_count; }

  KOKKOS_INLINE_FUNCTION
  void superblock_state(int sb_id, int &block_size, int &block_count_capacity,
                        int &block_count_used) const noexcept {
    block_size           = 0;
    block_count_capacity = 0;
    block_count_used     = 0;

    bool can_access_state_array = []() {
      KOKKOS_IF_ON_HOST(
          (return SpaceAccessibility<DefaultHostExecutionSpace,
                                     base_memory_space>::accessible;))
      KOKKOS_IF_ON_DEVICE(
          (return SpaceAccessibility<DefaultExecutionSpace,
                                     base_memory_space>::accessible;))
    }();

    if (can_access_state_array) {
      // Can access the state array

      const uint32_t state =
          ((uint32_t volatile *)m_sb_state_array)[sb_id * m_sb_state_size];

      const uint32_t block_count_lg2 = state >> state_shift;
      const uint32_t block_used      = state & state_used_mask;

      block_size           = 1LU << (m_sb_size_lg2 - block_count_lg2);
      block_count_capacity = 1LU << block_count_lg2;
      block_count_used     = block_used;
    }
  }
};

}  // namespace Kokkos

#endif /* #ifndef KOKKOS_MEMORYPOOL_HPP */
