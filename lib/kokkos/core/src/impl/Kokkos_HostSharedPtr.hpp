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
#ifdef KOKKOS_ENABLE_CXX17
    static_assert(std::is_invocable_v<Deleter, T*> &&
                  std::is_copy_constructible_v<Deleter>);
#endif
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
