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
        record_type::allocate(mem_space, "Kokkos::FixedBlockSizeMemPool_blocks",
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
