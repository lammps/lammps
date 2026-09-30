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
#pragma once

#include "macros.hpp"
#include "trait_backports.hpp"

#if !defined(MDSPAN_IMPL_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)
#  include "no_unique_address.hpp"
#endif

namespace MDSPAN_IMPL_STANDARD_NAMESPACE {
namespace detail {

// For no unique address emulation, this is the case taken when neither are empty.
// For real `[[no_unique_address]]`, this case is always taken.
template <class T1, class T2, class Enable = void> struct impl_compressed_pair {
  MDSPAN_IMPL_NO_UNIQUE_ADDRESS T1 m_t1_val{};
  MDSPAN_IMPL_NO_UNIQUE_ADDRESS T2 m_t2_val{};
  MDSPAN_FORCE_INLINE_FUNCTION MDSPAN_IMPL_CONSTEXPR_14 T1 &first() noexcept { return m_t1_val; }
  MDSPAN_FORCE_INLINE_FUNCTION constexpr T1 const &first() const noexcept {
    return m_t1_val;
  }
  MDSPAN_FORCE_INLINE_FUNCTION MDSPAN_IMPL_CONSTEXPR_14 T2 &second() noexcept { return m_t2_val; }
  MDSPAN_FORCE_INLINE_FUNCTION constexpr T2 const &second() const noexcept {
    return m_t2_val;
  }

  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr impl_compressed_pair() = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr impl_compressed_pair(impl_compressed_pair const &) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr impl_compressed_pair(impl_compressed_pair &&) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  MDSPAN_IMPL_CONSTEXPR_14_DEFAULTED impl_compressed_pair &
  operator=(impl_compressed_pair const &) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  MDSPAN_IMPL_CONSTEXPR_14_DEFAULTED impl_compressed_pair &
  operator=(impl_compressed_pair &&) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  ~impl_compressed_pair() = default;
  template <class T1Like, class T2Like>
  MDSPAN_INLINE_FUNCTION constexpr impl_compressed_pair(T1Like &&t1, T2Like &&t2)
      : m_t1_val((T1Like &&) t1), m_t2_val((T2Like &&) t2) {}
};

#if !defined(MDSPAN_IMPL_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)

// First empty.
template <class T1, class T2>
struct impl_compressed_pair<
    T1, T2,
    std::enable_if_t<MDSPAN_IMPL_TRAIT(std::is_empty, T1) && !MDSPAN_IMPL_TRAIT(std::is_empty, T2)>>
    : private T1 {
  T2 m_t2_val{};
  MDSPAN_FORCE_INLINE_FUNCTION MDSPAN_IMPL_CONSTEXPR_14 T1 &first() noexcept {
    return *static_cast<T1 *>(this);
  }
  MDSPAN_FORCE_INLINE_FUNCTION constexpr T1 const &first() const noexcept {
    return *static_cast<T1 const *>(this);
  }
  MDSPAN_FORCE_INLINE_FUNCTION MDSPAN_IMPL_CONSTEXPR_14 T2 &second() noexcept { return m_t2_val; }
  MDSPAN_FORCE_INLINE_FUNCTION constexpr T2 const &second() const noexcept {
    return m_t2_val;
  }

  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr impl_compressed_pair() = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr impl_compressed_pair(impl_compressed_pair const &) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr impl_compressed_pair(impl_compressed_pair &&) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  MDSPAN_IMPL_CONSTEXPR_14_DEFAULTED impl_compressed_pair &
  operator=(impl_compressed_pair const &) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  MDSPAN_IMPL_CONSTEXPR_14_DEFAULTED impl_compressed_pair &
  operator=(impl_compressed_pair &&) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  ~impl_compressed_pair() = default;
  template <class T1Like, class T2Like>
  MDSPAN_INLINE_FUNCTION constexpr impl_compressed_pair(T1Like &&t1, T2Like &&t2)
      : T1((T1Like &&) t1), m_t2_val((T2Like &&) t2) {}
};

// Second empty.
template <class T1, class T2>
struct impl_compressed_pair<
    T1, T2,
    std::enable_if_t<!MDSPAN_IMPL_TRAIT(std::is_empty, T1) && MDSPAN_IMPL_TRAIT(std::is_empty, T2)>>
    : private T2 {
  T1 m_t1_val{};
  MDSPAN_FORCE_INLINE_FUNCTION MDSPAN_IMPL_CONSTEXPR_14 T1 &first() noexcept { return m_t1_val; }
  MDSPAN_FORCE_INLINE_FUNCTION constexpr T1 const &first() const noexcept {
    return m_t1_val;
  }
  MDSPAN_FORCE_INLINE_FUNCTION MDSPAN_IMPL_CONSTEXPR_14 T2 &second() noexcept {
    return *static_cast<T2 *>(this);
  }
  MDSPAN_FORCE_INLINE_FUNCTION constexpr T2 const &second() const noexcept {
    return *static_cast<T2 const *>(this);
  }

  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr impl_compressed_pair() = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr impl_compressed_pair(impl_compressed_pair const &) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr impl_compressed_pair(impl_compressed_pair &&) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  MDSPAN_IMPL_CONSTEXPR_14_DEFAULTED impl_compressed_pair &
  operator=(impl_compressed_pair const &) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  MDSPAN_IMPL_CONSTEXPR_14_DEFAULTED impl_compressed_pair &
  operator=(impl_compressed_pair &&) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  ~impl_compressed_pair() = default;

  template <class T1Like, class T2Like>
  MDSPAN_INLINE_FUNCTION constexpr impl_compressed_pair(T1Like &&t1, T2Like &&t2)
      : T2((T2Like &&) t2), m_t1_val((T1Like &&) t1) {}
};

// Both empty.
template <class T1, class T2>
struct impl_compressed_pair<
    T1, T2,
    std::enable_if_t<MDSPAN_IMPL_TRAIT(std::is_empty, T1) && MDSPAN_IMPL_TRAIT(std::is_empty, T2)>>
    // We need to use the no_unique_address_emulation wrapper here to avoid
    // base class ambiguities.
#ifdef MDSPAN_IMPL_COMPILER_MSVC
// MSVC doesn't allow you to access public static member functions of a type
// when you *happen* to privately inherit from that type.
    : protected no_unique_address_emulation<T1, 0>,
      protected no_unique_address_emulation<T2, 1>
#else
    : private no_unique_address_emulation<T1, 0>,
      private no_unique_address_emulation<T2, 1>
#endif
{
  using first_base_t = no_unique_address_emulation<T1, 0>;
  using second_base_t = no_unique_address_emulation<T2, 1>;

  MDSPAN_FORCE_INLINE_FUNCTION MDSPAN_IMPL_CONSTEXPR_14 T1 &first() noexcept {
    return this->first_base_t::ref();
  }
  MDSPAN_FORCE_INLINE_FUNCTION constexpr T1 const &first() const noexcept {
    return this->first_base_t::ref();
  }
  MDSPAN_FORCE_INLINE_FUNCTION MDSPAN_IMPL_CONSTEXPR_14 T2 &second() noexcept {
    return this->second_base_t::ref();
  }
  MDSPAN_FORCE_INLINE_FUNCTION constexpr T2 const &second() const noexcept {
    return this->second_base_t::ref();
  }

  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr impl_compressed_pair() = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr impl_compressed_pair(impl_compressed_pair const &) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr impl_compressed_pair(impl_compressed_pair &&) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  MDSPAN_IMPL_CONSTEXPR_14_DEFAULTED impl_compressed_pair &
  operator=(impl_compressed_pair const &) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  MDSPAN_IMPL_CONSTEXPR_14_DEFAULTED impl_compressed_pair &
  operator=(impl_compressed_pair &&) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  ~impl_compressed_pair() = default;
  template <class T1Like, class T2Like>
  MDSPAN_INLINE_FUNCTION constexpr impl_compressed_pair(T1Like &&t1, T2Like &&t2) noexcept
    : first_base_t(T1((T1Like &&) t1)),
      second_base_t(T2((T2Like &&) t2))
  { }
};

#endif // !defined(MDSPAN_IMPL_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)

} // end namespace detail
} // end namespace MDSPAN_IMPL_STANDARD_NAMESPACE
