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

#ifndef KOKKOS_IMPL_TASKQUEUEMEMORYMANAGER_HPP
#define KOKKOS_IMPL_TASKQUEUEMEMORYMANAGER_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_TASKDAG)

#include <Kokkos_TaskScheduler_fwd.hpp>
#include <Kokkos_Core_fwd.hpp>

#include <Kokkos_MemoryPool.hpp>

#include <impl/Kokkos_TaskBase.hpp>
#include <impl/Kokkos_TaskResult.hpp>

#include <Kokkos_Atomic.hpp>
#include <impl/Kokkos_OptionalRef.hpp>
#include <impl/Kokkos_LIFO.hpp>

#include <string>
#include <typeinfo>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class ExecSpace, class MemorySpace,
          class MemoryPool =
              Kokkos::MemoryPool<Kokkos::Device<ExecSpace, MemorySpace>>>
class TaskQueueMemoryManager : public TaskQueueBase {
 public:
  using execution_space      = ExecSpace;
  using memory_space         = MemorySpace;
  using device_type          = Kokkos::Device<execution_space, memory_space>;
  using memory_pool          = MemoryPool;
  using allocation_size_type = size_t;

 private:
  memory_pool m_pool;
  // TODO @tasking @generalization DSH re-enable this with a flag in the type
  // long m_accum_alloc = 0;
  int m_count_alloc = 0;
  int m_max_alloc   = 0;

  struct _allocation_result {
    bool success;
    void* pointer;
  };

  KOKKOS_INLINE_FUNCTION
  _allocation_result _do_pool_allocate(allocation_size_type requested_size) {
    // KOKKOS_EXPECTS(requested_size >= 0); generates a warning when
    // allocation_size_type is unsigned
    if (requested_size == 0) {
      return {true, nullptr};
    } else {
      void* data = m_pool.allocate(static_cast<size_t>(requested_size));

      desul::atomic_inc(
          &m_count_alloc, desul::MemoryOrderSeqCst(),
          desul::MemoryScopeDevice());  // TODO? memory_order_relaxed
      // TODO @tasking @minor DSH make this thread safe? (otherwise, it's just
      // an approximation, which is probably fine...)
      if (m_max_alloc < m_count_alloc) m_max_alloc = m_count_alloc;

      return {data != nullptr, data};
    }
  }

  template <class T, class... Args>
  KOKKOS_INLINE_FUNCTION T* _do_contruct(void* allocated,
                                         allocation_size_type allocated_size,
                                         Args&&... args) {
    static_assert(std::is_base_of<PoolAllocatedObjectBase<int32_t>, T>::value,
                  "TaskQueueMemoryManager can only allocate objects with "
                  "PoolAllocatedObjectBase base class");

    // TODO @tasking DSH figure out why this isn't working
    // static_assert(
    //  std::is_constructible<T, Args..., int32_t>::value,
    //  "TaskQueueMemoryManager can't construct object of the requested type
    //  from the " " allocation size and the given arguments"
    //);

    auto rv = new (allocated) T(std::forward<Args>(args)..., allocated_size);

    // It feels like there should be a way to check this at compile-time
    KOKKOS_ASSERT(
        (intptr_t)(rv) ==
            (intptr_t)(static_cast<PoolAllocatedObjectBase<int32_t>*>(rv)) &&
        "PoolAllocatedObjectBase must be the first base class of the allocated "
        "type");

    return rv;
  }

 public:
  explicit TaskQueueMemoryManager(memory_pool const& pool) : m_pool(pool) {}

  template <class T, class... Args>
  KOKKOS_FUNCTION T* allocate_and_construct(Args&&... args)
  // requires
  //   std::is_base_of_v<PoolAllocatedObjectBase<typename
  //   memory_pool::size_type>, T>
  //     && std::is_constructible_v<T, Args&&..., allocation_size_type>
  {
    constexpr auto allocation_size = sizeof(T);

    auto result = _do_pool_allocate(allocation_size);

    KOKKOS_ASSERT(result.success && "Memory allocation failure");

    auto rv = _do_contruct<T>(result.pointer, allocation_size,
                              std::forward<Args>(args)...);

    KOKKOS_ENSURES(intptr_t(rv) % alignof(T) == 0 &&
                   "alignment not preserved!");

    return rv;
  }

  template <class T, class VLAValueType, class... Args>
  KOKKOS_INLINE_FUNCTION T* allocate_and_construct_with_vla_emulation(
      allocation_size_type n_vla_entries, Args&&... args)
  // requires
  //   std::is_base_of_v<PoolAllocatedObjectBase<typename
  //   memory_pool::size_type>, T>
  //     && std::is_base_of<ObjectWithVLAEmulation<T, VLAValueType>, T>::value
  //     && std::is_constructible_v<T, allocation_size_type, Args&&...>
  {
    static_assert(
        std::is_base_of<ObjectWithVLAEmulation<T, VLAValueType>, T>::value,
        "Can't append emulated variable length array of type with greater "
        "alignment than"
        "  the type to which the VLA is being appended");

    using vla_emulation_base = ObjectWithVLAEmulation<T, VLAValueType>;

    auto const allocation_size =
        vla_emulation_base::required_allocation_size(n_vla_entries);
    auto result = _do_pool_allocate(allocation_size);

    KOKKOS_ASSERT(result.success && "Memory allocation failure");

    auto rv = _do_contruct<T>(result.pointer, allocation_size,
                              std::forward<Args>(args)...);

    KOKKOS_ENSURES(intptr_t(rv) % alignof(T) == 0);

    return rv;
  }

  template <class CountType>
  KOKKOS_INLINE_FUNCTION void deallocate(
      PoolAllocatedObjectBase<CountType>&& obj) {
    m_pool.deallocate((void*)&obj, 1);
    desul::atomic_dec(
        &m_count_alloc, desul::MemoryOrderSeqCst(),
        desul::MemoryScopeDevice());  // TODO? memory_order_relaxed
  }

  KOKKOS_INLINE_FUNCTION
  memory_pool& get_memory_pool() { return m_pool; }
  KOKKOS_INLINE_FUNCTION
  memory_pool const& get_memory_pool() const { return m_pool; }

  KOKKOS_INLINE_FUNCTION
  int allocation_count() const noexcept { return m_count_alloc; }
};

} /* namespace Impl */
} /* namespace Kokkos */

////////////////////////////////////////////////////////////////////////////////
// END OLD CODE
////////////////////////////////////////////////////////////////////////////////

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_TASKDAG ) */
#endif /* #ifndef KOKKOS_IMPL_TASKQUEUEMEMORYMANAGER_HPP */
