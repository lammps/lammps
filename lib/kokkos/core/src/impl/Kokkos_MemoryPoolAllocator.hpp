/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
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
  using difference_type = typename std::make_signed<size_type>::type;

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
