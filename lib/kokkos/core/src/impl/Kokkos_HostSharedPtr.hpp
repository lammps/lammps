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

#ifndef KOKKOS_IMPL_HOST_SHARED_PTR_HPP
#define KOKKOS_IMPL_HOST_SHARED_PTR_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Atomic.hpp>
#include <impl/Kokkos_Error.hpp>

#include <functional>

namespace Kokkos {
namespace Impl {

template <typename T>
class HostSharedPtr {
 public:
  using element_type = T;

  KOKKOS_DEFAULTED_FUNCTION constexpr HostSharedPtr() = default;
  KOKKOS_FUNCTION constexpr HostSharedPtr(std::nullptr_t) {}

  explicit HostSharedPtr(T* element_ptr)
      : HostSharedPtr(element_ptr, [](T* const t) { delete t; }) {}

  template <class Deleter>
  HostSharedPtr(T* element_ptr, const Deleter& deleter)
      : m_element_ptr(element_ptr) {
    static_assert(std::is_invocable_v<Deleter, T*> &&
                  std::is_copy_constructible_v<Deleter>);
    if (element_ptr) {
      try {
        m_control = new Control{deleter, 1};
      } catch (...) {
        deleter(element_ptr);
        throw;
      }
    }
  }

  KOKKOS_FUNCTION HostSharedPtr(HostSharedPtr&& other) noexcept
      : m_element_ptr(other.m_element_ptr), m_control(other.m_control) {
    other.m_element_ptr = nullptr;
    other.m_control     = nullptr;
  }

  KOKKOS_FUNCTION HostSharedPtr(const HostSharedPtr& other) noexcept
      : m_element_ptr(other.m_element_ptr), m_control(other.m_control) {
    KOKKOS_IF_ON_HOST(
        (if (m_control) Kokkos::atomic_add(&(m_control->m_counter), 1);))
    KOKKOS_IF_ON_DEVICE(m_control = nullptr;)
  }

  KOKKOS_FUNCTION HostSharedPtr& operator=(HostSharedPtr&& other) noexcept {
    if (&other != this) {
      cleanup();
      m_element_ptr       = other.m_element_ptr;
      other.m_element_ptr = nullptr;
      m_control           = other.m_control;
      other.m_control     = nullptr;
    }
    return *this;
  }

  KOKKOS_FUNCTION HostSharedPtr& operator=(
      const HostSharedPtr& other) noexcept {
    if (&other != this) {
      cleanup();
      m_element_ptr = other.m_element_ptr;
      m_control     = other.m_control;
      KOKKOS_IF_ON_HOST(
          (if (m_control) Kokkos::atomic_add(&(m_control->m_counter), 1);))
      KOKKOS_IF_ON_DEVICE(m_control = nullptr;)
    }
    return *this;
  }

  KOKKOS_FUNCTION ~HostSharedPtr() { cleanup(); }

  // returns the stored pointer
  KOKKOS_FUNCTION T* get() const noexcept { return m_element_ptr; }
  // dereferences the stored pointer
  KOKKOS_FUNCTION T& operator*() const noexcept {
    KOKKOS_EXPECTS(bool(*this));
    return *get();
  }
  // dereferences the stored pointer
  KOKKOS_FUNCTION T* operator->() const noexcept {
    KOKKOS_EXPECTS(bool(*this));
    return get();
  }

  // checks if the stored pointer is not null
  KOKKOS_FUNCTION explicit operator bool() const noexcept {
    return get() != nullptr;
  }

  // returns the number of HostSharedPtr instances managing the current object
  // or 0 if there is no managed object.
  int use_count() const noexcept {
    return m_control ? m_control->m_counter : 0;
  }

 private:
  KOKKOS_FUNCTION void cleanup() noexcept {
    KOKKOS_IF_ON_HOST((
        // If m_counter is set, then this instance is responsible for managing
        // the object pointed to by m_counter and m_element_ptr.
        if (m_control) {
          int const count =
              Kokkos::atomic_fetch_sub(&(m_control->m_counter), 1);
          // atomic_fetch_sub might have memory order relaxed, so we need to
          // force synchronization to avoid multiple threads doing the cleanup.
          Kokkos::memory_fence();
          if (count == 1) {
            (m_control->m_deleter)(m_element_ptr);
            m_element_ptr = nullptr;
            delete m_control;
            m_control = nullptr;
          }
        }))
  }

  struct Control {
    std::function<void(T*)> m_deleter;
    int m_counter;
  };

  T* m_element_ptr   = nullptr;
  Control* m_control = nullptr;
};
}  // namespace Impl
}  // namespace Kokkos

#endif
