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

#ifndef KOKKOS_IMPL_LOCKFREEDEQUE_HPP
#define KOKKOS_IMPL_LOCKFREEDEQUE_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_TASKDAG

#include <Kokkos_Core_fwd.hpp>

#include <Kokkos_PointerOwnership.hpp>
#include <impl/Kokkos_OptionalRef.hpp>
#include <impl/Kokkos_Error.hpp>           // KOKKOS_EXPECTS
#include <impl/Kokkos_LinkedListNode.hpp>  // KOKKOS_EXPECTS

#include <Kokkos_Atomic.hpp>  // atomic_compare_exchange, atomic_fence
#include <impl/Kokkos_LIFO.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template <class NodeType, size_t CircularBufferSize, class SizeType = size_t>
struct fixed_size_circular_buffer {
 public:
  using node_type = NodeType;
  using size_type = SizeType;

 private:
  node_type* m_buffer[CircularBufferSize] = {nullptr};

 public:
  fixed_size_circular_buffer()                                  = default;
  fixed_size_circular_buffer(fixed_size_circular_buffer const&) = delete;
  fixed_size_circular_buffer(fixed_size_circular_buffer&&)      = default;
  fixed_size_circular_buffer& operator=(fixed_size_circular_buffer const&) =
      delete;
  fixed_size_circular_buffer& operator=(fixed_size_circular_buffer&&) = default;
  ~fixed_size_circular_buffer()                                       = default;

  KOKKOS_FORCEINLINE_FUNCTION
  static constexpr size_type size() noexcept {
    return size_type(CircularBufferSize);
  }

  KOKKOS_FORCEINLINE_FUNCTION
  node_type* operator[](size_type idx) const noexcept {
    return m_buffer[idx % size()];
  }

  KOKKOS_FORCEINLINE_FUNCTION
  node_type*& operator[](size_type idx) noexcept {
    return m_buffer[idx % size()];
  }
};

template <class NodeType, class SizeType = size_t>
struct non_owning_variable_size_circular_buffer {
 public:
  using node_type = NodeType;
  using size_type = SizeType;

 private:
  ObservingRawPtr<node_type*> m_buffer = nullptr;
  size_type m_size                     = 0;

 public:
  KOKKOS_INLINE_FUNCTION
  non_owning_variable_size_circular_buffer(ObservingRawPtr<node_type*> buffer,
                                           size_type arg_size) noexcept
      : m_buffer(buffer), m_size(arg_size) {}

  non_owning_variable_size_circular_buffer() = default;
  non_owning_variable_size_circular_buffer(
      non_owning_variable_size_circular_buffer const&) = delete;
  non_owning_variable_size_circular_buffer(
      non_owning_variable_size_circular_buffer&&)      = default;
  non_owning_variable_size_circular_buffer& operator   =(
      non_owning_variable_size_circular_buffer const&) = delete;
  non_owning_variable_size_circular_buffer& operator   =(
      non_owning_variable_size_circular_buffer&&) = default;
  ~non_owning_variable_size_circular_buffer()          = default;

  KOKKOS_FORCEINLINE_FUNCTION
  constexpr size_type size() const noexcept { return m_size; }

  KOKKOS_FORCEINLINE_FUNCTION
  node_type* operator[](size_type idx) const noexcept {
    return m_buffer[idx % size()];
  }

  KOKKOS_FORCEINLINE_FUNCTION
  node_type*& operator[](size_type idx) noexcept {
    return m_buffer[idx % size()];
  }
};

/** Based on "Correct and Efficient Work-Stealing for Weak Memory Models,"
 * PPoPP '13, https://www.di.ens.fr/~zappa/readings/ppopp13.pdf
 *
 */
template <class T, class CircularBufferT, class SizeType = int32_t>
struct ChaseLevDeque {
 public:
  using size_type  = SizeType;
  using value_type = T;
  // Still using intrusive linked list for waiting queue
  using node_type = SimpleSinglyLinkedListNode<>;

 private:
  // TODO @tasking @new_feature DSH variable size circular buffer?

  CircularBufferT m_array;
  size_type m_top    = 0;
  size_type m_bottom = 0;

 public:
  template <class _ignore = void,
            class         = std::enable_if_t<
                std::is_default_constructible<CircularBufferT>::value>>
  ChaseLevDeque() : m_array() {}

  explicit ChaseLevDeque(CircularBufferT buffer) : m_array(std::move(buffer)) {}

  KOKKOS_INLINE_FUNCTION
  bool empty() const {
    // TODO @tasking @memory_order DSH memory order
    return m_top > m_bottom - 1;
  }

  KOKKOS_INLINE_FUNCTION
  OptionalRef<T> pop() {
    auto b   = m_bottom - 1;  // atomic load relaxed
    auto& a  = m_array;       // atomic load relaxed
    m_bottom = b;             // atomic store relaxed
    Kokkos::memory_fence();   // memory order seq_cst
    auto t = m_top;           // atomic load relaxed
    OptionalRef<T> return_value;
    if (t <= b) {
      /* non-empty queue */
      return_value = *static_cast<T*>(a[b]);  // relaxed load
      if (t == b) {
        /* single last element in the queue. */
#ifdef _WIN32
        Kokkos::memory_fence();
        bool const success =
            Kokkos::atomic_compare_exchange_strong(&m_top, t, t + 1);
        Kokkos::memory_fence();
        if (!success) {
          return_value = nullptr;
        }
#else
        if (!Impl::atomic_compare_exchange_strong(
                &m_top, t, t + 1, desul::MemoryOrderSeqCst(),
                desul::MemoryOrderRelaxed())) {
          /* failed race, someone else stole it */
          return_value = nullptr;
        }
#endif
        m_bottom = b + 1;  // memory order relaxed
      }
    } else {
      /* empty queue */
      m_bottom = b + 1;  // memory order relaxed
    }
    return return_value;
  }

  KOKKOS_INLINE_FUNCTION
  bool push(node_type&& node) {
    // Just forward to the lvalue version
    return push(node);
  }

  KOKKOS_INLINE_FUNCTION
  bool push(node_type& node) {
    auto b  = m_bottom;  // memory order relaxed
    auto t  = Impl::atomic_load(&m_top, desul::MemoryOrderAcquire());
    auto& a = m_array;
    if (b - t > a.size() - 1) {
      /* queue is full, resize */
      // m_array = a->grow();
      // a = m_array;
      return false;
    }
    a[b] = &node;  // relaxed
    Impl::atomic_store(&m_bottom, b + 1, desul::MemoryOrderRelease());
    return true;
  }

  KOKKOS_INLINE_FUNCTION
  OptionalRef<T> steal() {
    auto t = m_top;  // TODO @tasking @memory_order DSH: atomic load acquire
    Kokkos::memory_fence();  // seq_cst fence, so why does the above need to be
                             // acquire?
    auto b = Impl::atomic_load(&m_bottom, desul::MemoryOrderAcquire());
    OptionalRef<T> return_value;
    if (t < b) {
      /* Non-empty queue */
      auto& a = m_array;     // TODO @tasking @memory_order DSH: technically
                             // consume ordered, but acquire should be fine
      Kokkos::load_fence();  // TODO @tasking @memory_order DSH memory order
                             // instead of fence
      return_value = *static_cast<T*>(a[t]);  // relaxed
#ifdef _WIN32
      Kokkos::memory_fence();
      bool const success =
          Kokkos::atomic_compare_exchange_strong(&m_top, t, t + 1);
      Kokkos::memory_fence();
      if (!success) {
        return_value = nullptr;
      }
#else
      if (!Impl::atomic_compare_exchange_strong(&m_top, t, t + 1,
                                                desul::MemoryOrderSeqCst(),
                                                desul::MemoryOrderRelaxed())) {
        return_value = nullptr;
      }
#endif
    }
    return return_value;
  }
};

/*
      // The atomicity of this load was more important in the paper's version
      // because that version had a circular buffer that could grow.  We're
      // essentially using the memory order in this version as a fence, which
      // may be unnecessary
      auto buffer_ptr = (node_type***)&m_array.buffer;
      auto a = Impl::atomic_load(buffer_ptr, desul::MemoryOrderAcquire()); //
   technically consume ordered, but acquire should be fine return_value =
   *static_cast<T*>(a[t % m_array->size]); // relaxed; we'd have to replace the
   m_array->size if we ever allow growth
*/

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template <size_t CircularBufferSize>
struct TaskQueueTraitsChaseLev {
  template <class Task>
  using ready_queue_type =
      ChaseLevDeque<Task,
                    fixed_size_circular_buffer<SimpleSinglyLinkedListNode<>,
                                               CircularBufferSize, int32_t>,
                    int32_t>;

  template <class Task>
  using waiting_queue_type = SingleConsumeOperationLIFO<Task>;

  template <class Task>
  using intrusive_task_base_type = typename ready_queue_type<Task>::node_type;

  static constexpr auto ready_queue_insertion_may_fail = true;
};

}  // end namespace Impl
}  // end namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* defined KOKKOS_ENABLE_TASKDAG */
#endif /* #ifndef KOKKOS_IMPL_LOCKFREEDEQUE_HPP */
