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

// Experimental unified task-data parallel manycore LDRD

#ifndef KOKKOS_IMPL_MEMORYPOOLALLOCATOR_HPP
#define KOKKOS_IMPL_MEMORYPOOLALLOCATOR_HPP

#include <Kokkos_Macros.hpp>

#include <Kokkos_Core_fwd.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
namespace Kokkos {
namespace Impl {

template <class MemoryPool, class T>
class MemoryPoolAllocator {
 public:
  using memory_pool = MemoryPool;

 private:
  memory_pool m_pool;

 public:
  KOKKOS_DEFAULTED_FUNCTION
  MemoryPoolAllocator() = default;
  KOKKOS_DEFAULTED_FUNCTION
  MemoryPoolAllocator(MemoryPoolAllocator const&) = default;
  KOKKOS_DEFAULTED_FUNCTION
  MemoryPoolAllocator(MemoryPoolAllocator&&) = default;
  KOKKOS_DEFAULTED_FUNCTION
  MemoryPoolAllocator& operator=(MemoryPoolAllocator const&) = default;
  KOKKOS_DEFAULTED_FUNCTION
  MemoryPoolAllocator& operator=(MemoryPoolAllocator&&) = default;
  KOKKOS_DEFAULTED_FUNCTION
  ~MemoryPoolAllocator() = default;

  KOKKOS_INLINE_FUNCTION
  explicit MemoryPoolAllocator(memory_pool const& arg_pool)
      : m_pool(arg_pool) {}
  KOKKOS_INLINE_FUNCTION
  explicit MemoryPoolAllocator(memory_pool&& arg_pool)
      : m_pool(std::move(arg_pool)) {}

 public:
  using value_type      = T;
  using pointer         = T*;
  using size_type       = typename MemoryPool::memory_space::size_type;
  using difference_type = std::make_signed_t<size_type>;

  template <class U>
  struct rebind {
    using other = MemoryPoolAllocator<MemoryPool, U>;
  };

  KOKKOS_INLINE_FUNCTION
  pointer allocate(size_t n) {
    void* rv = m_pool.allocate(n * sizeof(T));
    if (rv == nullptr) {
      Kokkos::abort("Kokkos MemoryPool allocator failed to allocate memory");
    }
    return reinterpret_cast<T*>(rv);
  }

  KOKKOS_INLINE_FUNCTION
  void deallocate(T* ptr, size_t n) { m_pool.deallocate(ptr, n * sizeof(T)); }

  KOKKOS_INLINE_FUNCTION
  size_type max_size() const { return m_pool.max_block_size(); }

  KOKKOS_INLINE_FUNCTION
  bool operator==(MemoryPoolAllocator const& other) const {
    return m_pool == other.m_pool;
  }

  KOKKOS_INLINE_FUNCTION
  bool operator!=(MemoryPoolAllocator const& other) const {
    return !(*this == other);
  }
};

}  // end namespace Impl
}  // end namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_IMPL_MEMORYPOOLALLOCATOR_HPP */
