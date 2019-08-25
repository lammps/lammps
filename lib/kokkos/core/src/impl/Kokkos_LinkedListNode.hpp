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
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

// Experimental unified task-data parallel manycore LDRD

#ifndef KOKKOS_IMPL_LINKEDLISTNODE_HPP
#define KOKKOS_IMPL_LINKEDLISTNODE_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_TASKDAG // Note: implies CUDA_VERSION >= 8000 if using CUDA

#include <Kokkos_Core_fwd.hpp>

#include <Kokkos_PointerOwnership.hpp>
#include <impl/Kokkos_OptionalRef.hpp>
#include <impl/Kokkos_Error.hpp> // KOKKOS_EXPECTS

#include <Kokkos_Atomic.hpp>  // atomic_compare_exchange, atomic_fence

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

struct LinkedListNodeAccess;

template <
  uintptr_t NotEnqueuedValue = 0,
  template <class> class PointerTemplate = std::add_pointer
>
struct SimpleSinglyLinkedListNode
{

private:

  using pointer_type = typename PointerTemplate<SimpleSinglyLinkedListNode>::type;

  pointer_type m_next = reinterpret_cast<pointer_type>(NotEnqueuedValue);

  // These are private because they are an implementation detail of the queue
  // and should not get added to the value type's interface via the intrusive
  // wrapper.

  KOKKOS_INLINE_FUNCTION
  void mark_as_not_enqueued() noexcept {
    // TODO @tasking @memory_order DSH make this an atomic store with memory order
    m_next = (pointer_type)NotEnqueuedValue;
  }

  KOKKOS_INLINE_FUNCTION
  void mark_as_not_enqueued() volatile noexcept {
    // TODO @tasking @memory_order DSH make this an atomic store with memory order
    m_next = (pointer_type)NotEnqueuedValue;
  }

  KOKKOS_INLINE_FUNCTION
  pointer_type& _next_ptr() noexcept {
    return m_next;
  }

  KOKKOS_INLINE_FUNCTION
  pointer_type volatile& _next_ptr() volatile noexcept {
    return m_next;
  }

  KOKKOS_INLINE_FUNCTION
  pointer_type const& _next_ptr() const noexcept {
    return m_next;
  }
  
  KOKKOS_INLINE_FUNCTION
  pointer_type const volatile& _next_ptr() const volatile noexcept {
    return m_next;
  }

  friend struct LinkedListNodeAccess;

public:

  // KOKKOS_CONSTEXPR_14
  KOKKOS_INLINE_FUNCTION
  bool is_enqueued() const noexcept {
    // TODO @tasking @memory_order DSH make this an atomic load with memory order
    return m_next != reinterpret_cast<pointer_type>(NotEnqueuedValue);
  }

  // KOKKOS_CONSTEXPR_14
  KOKKOS_INLINE_FUNCTION
  bool is_enqueued() const volatile noexcept {
    // TODO @tasking @memory_order DSH make this an atomic load with memory order
    return m_next != reinterpret_cast<pointer_type>(NotEnqueuedValue);
  }

};

/// Attorney for LinkedListNode, since user types inherit from it
struct LinkedListNodeAccess
{

  template <class Node>
  KOKKOS_INLINE_FUNCTION
  static void mark_as_not_enqueued(Node& node) noexcept {
    node.mark_as_not_enqueued();
  }

  template <class Node>
  KOKKOS_INLINE_FUNCTION
  static void mark_as_not_enqueued(Node volatile& node) noexcept {
    node.mark_as_not_enqueued();
  }

  template <class Node>
  KOKKOS_INLINE_FUNCTION
  static
  typename Node::pointer_type&
  next_ptr(Node& node) noexcept {
    return node._next_ptr();
  }

  template <class Node>
  KOKKOS_INLINE_FUNCTION
  static
  typename Node::pointer_type&
  next_ptr(Node volatile& node) noexcept {
    return node._next_ptr();
  }

  template <class Node>
  KOKKOS_INLINE_FUNCTION
  static
  typename Node::pointer_type&
  next_ptr(Node const& node) noexcept {
    return node._next_ptr();
  }

  template <class Node>
  KOKKOS_INLINE_FUNCTION
  static
  typename Node::pointer_type&
  prev_ptr(Node& node) noexcept {
    return node._prev_ptr();
  }

  template <class Node>
  KOKKOS_INLINE_FUNCTION
  static
  typename Node::pointer_type&
  prev_ptr(Node const& node) noexcept {
    return node._prev_ptr();
  }

};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // end namespace Impl
} // end namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* defined KOKKOS_ENABLE_TASKDAG */
#endif /* #ifndef KOKKOS_IMPL_LINKEDLISTNODE_HPP */

