/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//              Copyright (2019) Sandia Corporation
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
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

#ifndef KOKKOS_IMPL_KOKKOS_FIXEDBUFFERMEMORYPOOL_HPP
#define KOKKOS_IMPL_KOKKOS_FIXEDBUFFERMEMORYPOOL_HPP

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_Atomic.hpp>

#include <Kokkos_PointerOwnership.hpp>
#include <impl/Kokkos_SimpleTaskScheduler.hpp>

namespace Kokkos {
namespace Impl {

template <class DeviceType, size_t Size, size_t Align = 1,
          class SizeType = typename DeviceType::execution_space::size_type>
class FixedBlockSizeMemoryPool
    : private MemorySpaceInstanceStorage<typename DeviceType::memory_space> {
 public:
  using memory_space = typename DeviceType::memory_space;
  using size_type    = SizeType;

 private:
  using memory_space_storage_base =
      MemorySpaceInstanceStorage<typename DeviceType::memory_space>;
  using tracker_type = Kokkos::Impl::SharedAllocationTracker;
  using record_type  = Kokkos::Impl::SharedAllocationRecord<memory_space>;

  struct alignas(Align) Block {
    union {
      char ignore;
      char data[Size];
    };
  };

  static constexpr auto actual_size = sizeof(Block);

  // TODO shared allocation tracker
  // TODO @optimization put the index values on different cache lines (CPU) or
  // pages (GPU)?

  tracker_type m_tracker                         = {};
  size_type m_num_blocks                         = 0;
  size_type m_first_free_idx                     = 0;
  size_type m_last_free_idx                      = 0;
  Kokkos::OwningRawPtr<Block> m_first_block      = nullptr;
  Kokkos::OwningRawPtr<size_type> m_free_indices = nullptr;

  enum : size_type { IndexInUse = ~size_type(0) };

 public:
  FixedBlockSizeMemoryPool(memory_space const& mem_space, size_type num_blocks)
      : memory_space_storage_base(mem_space),
        m_tracker(),
        m_num_blocks(num_blocks),
        m_first_free_idx(0),
        m_last_free_idx(num_blocks) {
    // TODO alignment?
    auto block_record = record_type::allocate(
        mem_space, "FixedBlockSizeMemPool_blocks", num_blocks * sizeof(Block));
    KOKKOS_ASSERT(intptr_t(block_record->data()) % Align == 0);
    m_tracker.assign_allocated_record_to_uninitialized(block_record);
    m_first_block = (Block*)block_record->data();

    auto idx_record =
        record_type::allocate(mem_space, "FixedBlockSizeMemPool_blocks",
                              num_blocks * sizeof(size_type));
    KOKKOS_ASSERT(intptr_t(idx_record->data()) % alignof(size_type) == 0);
    m_tracker.assign_allocated_record_to_uninitialized(idx_record);
    m_free_indices = (size_type*)idx_record->data();

    for (size_type i = 0; i < num_blocks; ++i) {
      m_free_indices[i] = i;
    }

    Kokkos::memory_fence();
  }

  // For compatibility with MemoryPool<>
  FixedBlockSizeMemoryPool(memory_space const& mem_space,
                           size_t mempool_capacity, unsigned, unsigned,
                           unsigned)
      : FixedBlockSizeMemoryPool(
            mem_space, mempool_capacity /
                           actual_size) { /* forwarding ctor, must be empty */
  }

  KOKKOS_DEFAULTED_FUNCTION FixedBlockSizeMemoryPool() = default;
  KOKKOS_DEFAULTED_FUNCTION FixedBlockSizeMemoryPool(
      FixedBlockSizeMemoryPool&&) = default;
  KOKKOS_DEFAULTED_FUNCTION FixedBlockSizeMemoryPool(
      FixedBlockSizeMemoryPool const&)                        = default;
  KOKKOS_DEFAULTED_FUNCTION FixedBlockSizeMemoryPool& operator=(
      FixedBlockSizeMemoryPool&&) = default;
  KOKKOS_DEFAULTED_FUNCTION FixedBlockSizeMemoryPool& operator=(
      FixedBlockSizeMemoryPool const&) = default;

  KOKKOS_INLINE_FUNCTION
  void* allocate(size_type alloc_size) const noexcept {
    (void)alloc_size;
    KOKKOS_EXPECTS(alloc_size <= Size);
    auto free_idx_counter = Kokkos::atomic_fetch_add(
        (volatile size_type*)&m_first_free_idx, size_type(1));
    auto free_idx_idx = free_idx_counter % m_num_blocks;

    // We don't have exclusive access to m_free_indices[free_idx_idx] because
    // the allocate counter might have lapped us since we incremented it
    auto current_free_idx = m_free_indices[free_idx_idx];
    size_type free_idx    = IndexInUse;
    free_idx = Kokkos::atomic_compare_exchange(&m_free_indices[free_idx_idx],
                                               current_free_idx, free_idx);
    Kokkos::memory_fence();

    // TODO figure out how to decrement here?

    if (free_idx == IndexInUse) {
      return nullptr;
    } else {
      return (void*)&m_first_block[free_idx];
    }
  }

  KOKKOS_INLINE_FUNCTION
  void deallocate(void* ptr, size_type /*alloc_size*/) const noexcept {
    // figure out which block we are
    auto offset = intptr_t(ptr) - intptr_t(m_first_block);

    KOKKOS_EXPECTS(offset % actual_size == 0 &&
                   offset / actual_size < m_num_blocks);

    Kokkos::memory_fence();
    auto last_idx_idx = Kokkos::atomic_fetch_add(
        (volatile size_type*)&m_last_free_idx, size_type(1));
    last_idx_idx %= m_num_blocks;
    m_free_indices[last_idx_idx] = offset / actual_size;
  }
};

#if 0
template <
  class DeviceType,
  size_t Size,
  size_t Align=1,
  class SizeType = typename DeviceType::execution_space::size_type
>
class FixedBlockSizeChaseLevMemoryPool
  : private MemorySpaceInstanceStorage<typename DeviceType::memory_space>
{
public:

  using memory_space = typename DeviceType::memory_space;
  using size_type = SizeType;

private:

  using memory_space_storage_base = MemorySpaceInstanceStorage<typename DeviceType::memory_space>;
  using tracker_type = Kokkos::Impl::SharedAllocationTracker;
  using record_type = Kokkos::Impl::SharedAllocationRecord<memory_space>;

  struct alignas(Align) Block { union { char ignore; char data[Size]; }; };

  static constexpr auto actual_size = sizeof(Block);

  tracker_type m_tracker = { };
  size_type m_num_blocks = 0;
  size_type m_first_free_idx = 0;
  size_type m_last_free_idx = 0;


  enum : size_type { IndexInUse = ~size_type(0) };

public:

  FixedBlockSizeMemoryPool(
    memory_space const& mem_space,
    size_type num_blocks
  ) : memory_space_storage_base(mem_space),
    m_tracker(),
    m_num_blocks(num_blocks),
    m_first_free_idx(0),
    m_last_free_idx(num_blocks)
  {
    // TODO alignment?
    auto block_record = record_type::allocate(
      mem_space, "FixedBlockSizeMemPool_blocks", num_blocks * sizeof(Block)
    );
    KOKKOS_ASSERT(intptr_t(block_record->data()) % Align == 0);
    m_tracker.assign_allocated_record_to_uninitialized(block_record);
    m_first_block = (Block*)block_record->data();

    auto idx_record = record_type::allocate(
      mem_space, "FixedBlockSizeMemPool_blocks", num_blocks * sizeof(size_type)
    );
    KOKKOS_ASSERT(intptr_t(idx_record->data()) % alignof(size_type) == 0);
    m_tracker.assign_allocated_record_to_uninitialized(idx_record);
    m_free_indices = (size_type*)idx_record->data();

    for(size_type i = 0; i < num_blocks; ++i) {
      m_free_indices[i] = i;
    }

    Kokkos::memory_fence();
  }

  // For compatibility with MemoryPool<>
  FixedBlockSizeMemoryPool(
    memory_space const& mem_space,
    size_t mempool_capacity,
    unsigned, unsigned, unsigned
  ) : FixedBlockSizeMemoryPool(mem_space, mempool_capacity / actual_size)
  { /* forwarding ctor, must be empty */ }

  KOKKOS_DEFAULTED_FUNCTION FixedBlockSizeMemoryPool() = default;
  KOKKOS_DEFAULTED_FUNCTION FixedBlockSizeMemoryPool(FixedBlockSizeMemoryPool&&) = default;
  KOKKOS_DEFAULTED_FUNCTION FixedBlockSizeMemoryPool(FixedBlockSizeMemoryPool const&) = default;
  KOKKOS_DEFAULTED_FUNCTION FixedBlockSizeMemoryPool& operator=(FixedBlockSizeMemoryPool&&) = default;
  KOKKOS_DEFAULTED_FUNCTION FixedBlockSizeMemoryPool& operator=(FixedBlockSizeMemoryPool const&) = default;


  KOKKOS_INLINE_FUNCTION
  void* allocate(size_type alloc_size) const noexcept
  {
    KOKKOS_EXPECTS(alloc_size <= Size);
    auto free_idx_counter = Kokkos::atomic_fetch_add((volatile size_type*)&m_first_free_idx, size_type(1));
    auto free_idx_idx = free_idx_counter % m_num_blocks;

    // We don't have exclusive access to m_free_indices[free_idx_idx] because
    // the allocate counter might have lapped us since we incremented it
    auto current_free_idx = m_free_indices[free_idx_idx];
    size_type free_idx = IndexInUse;
    free_idx =
      Kokkos::atomic_compare_exchange(&m_free_indices[free_idx_idx], current_free_idx, free_idx);
    Kokkos::memory_fence();

    // TODO figure out how to decrement here?

    if(free_idx == IndexInUse) {
      return nullptr;
    }
    else {
      return (void*)&m_first_block[free_idx];
    }
  }

  KOKKOS_INLINE_FUNCTION
  void deallocate(void* ptr, size_type alloc_size) const noexcept
  {
    // figure out which block we are
    auto offset = intptr_t(ptr) - intptr_t(m_first_block);

    KOKKOS_EXPECTS(offset % actual_size == 0 && offset/actual_size < m_num_blocks);

    Kokkos::memory_fence();
    auto last_idx_idx = Kokkos::atomic_fetch_add((volatile size_type*)&m_last_free_idx, size_type(1));
    last_idx_idx %= m_num_blocks;
    m_free_indices[last_idx_idx] = offset / actual_size;
  }

};
#endif

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_IMPL_KOKKOS_FIXEDBUFFERMEMORYPOOL_HPP
