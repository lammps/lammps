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

#ifndef KOKKOS_IMPL_OPTIONALREF_HPP
#define KOKKOS_IMPL_OPTIONALREF_HPP

#include <Kokkos_Macros.hpp>

#include <Kokkos_Core_fwd.hpp>

#include <Kokkos_PointerOwnership.hpp>
#include <impl/Kokkos_Error.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
namespace Kokkos {
namespace Impl {

struct InPlaceTag {};

template <class T>
struct OptionalRef {
 private:
  ObservingRawPtr<T> m_value = nullptr;

 public:
  using value_type = T;

  KOKKOS_DEFAULTED_FUNCTION
  OptionalRef() = default;

  KOKKOS_DEFAULTED_FUNCTION
  OptionalRef(OptionalRef const&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  OptionalRef(OptionalRef&&) = default;

  KOKKOS_INLINE_FUNCTION
  // MSVC requires that this copy constructor is not defaulted
  // if there exists a (non-defaulted) volatile one.
  OptionalRef& operator=(OptionalRef const& other) noexcept {
    m_value = other.m_value;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  // Can't return a reference to volatile OptionalRef, since GCC issues a
  // warning about reference to volatile not accessing the underlying value
  void operator=(OptionalRef const volatile& other) volatile noexcept {
    m_value = other.m_value;
  }

  KOKKOS_DEFAULTED_FUNCTION
  OptionalRef& operator=(OptionalRef&&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  ~OptionalRef() = default;

  KOKKOS_INLINE_FUNCTION
  explicit OptionalRef(T& arg_value) : m_value(&arg_value) {}

  KOKKOS_INLINE_FUNCTION
  explicit OptionalRef(std::nullptr_t) : m_value(nullptr) {}

  KOKKOS_INLINE_FUNCTION
  OptionalRef& operator=(T& arg_value) {
    m_value = &arg_value;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  OptionalRef& operator=(std::nullptr_t) {
    m_value = nullptr;
    return *this;
  }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  OptionalRef<std::add_volatile_t<T>> as_volatile() volatile noexcept {
    return OptionalRef<std::add_volatile_t<T>>(*(*this));
  }

  KOKKOS_INLINE_FUNCTION
  OptionalRef<std::add_volatile_t<std::add_const_t<T>>> as_volatile() const
      volatile noexcept {
    return OptionalRef<std::add_volatile_t<std::add_const_t<T>>>(*(*this));
  }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  T& operator*() & {
    KOKKOS_EXPECTS(this->has_value());
    return *m_value;
  }

  KOKKOS_INLINE_FUNCTION
  T const& operator*() const& {
    KOKKOS_EXPECTS(this->has_value());
    return *m_value;
  }

  KOKKOS_INLINE_FUNCTION
  T volatile& operator*() volatile& {
    KOKKOS_EXPECTS(this->has_value());
    return *m_value;
  }

  KOKKOS_INLINE_FUNCTION
  T const volatile& operator*() const volatile& {
    KOKKOS_EXPECTS(this->has_value());
    return *m_value;
  }

  KOKKOS_INLINE_FUNCTION
  T&& operator*() && {
    KOKKOS_EXPECTS(this->has_value());
    return std::move(*m_value);
  }

  KOKKOS_INLINE_FUNCTION
  T* operator->() {
    KOKKOS_EXPECTS(this->has_value());
    return m_value;
  }

  KOKKOS_INLINE_FUNCTION
  T const* operator->() const {
    KOKKOS_EXPECTS(this->has_value());
    return m_value;
  }

  KOKKOS_INLINE_FUNCTION
  T volatile* operator->() volatile {
    KOKKOS_EXPECTS(this->has_value());
    return m_value;
  }

  KOKKOS_INLINE_FUNCTION
  T const volatile* operator->() const volatile {
    KOKKOS_EXPECTS(this->has_value());
    return m_value;
  }

  KOKKOS_INLINE_FUNCTION
  T* get() { return m_value; }

  KOKKOS_INLINE_FUNCTION
  T const* get() const { return m_value; }

  KOKKOS_INLINE_FUNCTION
  T volatile* get() volatile { return m_value; }

  KOKKOS_INLINE_FUNCTION
  T const volatile* get() const volatile { return m_value; }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  operator bool() { return m_value != nullptr; }

  KOKKOS_INLINE_FUNCTION
  operator bool() const { return m_value != nullptr; }

  KOKKOS_INLINE_FUNCTION
  operator bool() volatile { return m_value != nullptr; }

  KOKKOS_INLINE_FUNCTION
  operator bool() const volatile { return m_value != nullptr; }

  KOKKOS_INLINE_FUNCTION
  bool has_value() { return m_value != nullptr; }

  KOKKOS_INLINE_FUNCTION
  bool has_value() const { return m_value != nullptr; }

  KOKKOS_INLINE_FUNCTION
  bool has_value() volatile { return m_value != nullptr; }

  KOKKOS_INLINE_FUNCTION
  bool has_value() const volatile { return m_value != nullptr; }
};

}  // end namespace Impl
}  // end namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_IMPL_OPTIONALREF_HPP */
