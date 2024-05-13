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

#ifndef KOKKOS_ARRAY_HPP
#define KOKKOS_ARRAY_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_ARRAY
#endif

#include <Kokkos_Macros.hpp>
#include <Kokkos_Swap.hpp>
#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_StringManipulation.hpp>

#include <type_traits>
#include <algorithm>
#include <utility>
#include <limits>
#include <cstddef>

namespace Kokkos {

#ifdef KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK
namespace Impl {
template <typename Integral, bool Signed = std::is_signed<Integral>::value>
struct ArrayBoundsCheck;

template <typename Integral>
struct ArrayBoundsCheck<Integral, true> {
  KOKKOS_INLINE_FUNCTION
  constexpr ArrayBoundsCheck(Integral i, size_t N) {
    if (i < 0) {
      char err[128] = "Kokkos::Array: index ";
      to_chars_i(err + strlen(err), err + 128, i);
      strcat(err, " < 0");
      Kokkos::abort(err);
    }
    ArrayBoundsCheck<Integral, false>(i, N);
  }
};

template <typename Integral>
struct ArrayBoundsCheck<Integral, false> {
  KOKKOS_INLINE_FUNCTION
  constexpr ArrayBoundsCheck(Integral i, size_t N) {
    if (size_t(i) >= N) {
      char err[128] = "Kokkos::Array: index ";
      to_chars_i(err + strlen(err), err + 128, i);
      strcat(err, " >= ");
      to_chars_i(err + strlen(err), err + 128, N);
      Kokkos::abort(err);
    }
  }
};
}  // end namespace Impl

#define KOKKOS_ARRAY_BOUNDS_CHECK(i, N) \
  Kokkos::Impl::ArrayBoundsCheck<decltype(i)>(i, N)

#else  // !defined( KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK )

#define KOKKOS_ARRAY_BOUNDS_CHECK(i, N) (void)0

#endif  // !defined( KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK )

/**\brief  Derived from the C++17 'std::array'.
 *         Dropping the iterator interface.
 */
template <class T = void, size_t N = KOKKOS_INVALID_INDEX, class Proxy = void>
struct Array {
 public:
  /**
   * The elements of this C array shall not be accessed directly. The data
   * member has to be declared public to enable aggregate initialization as for
   * std::array. We mark it as private in the documentation.
   * @private
   */
  T m_internal_implementation_private_member_data[N];

 public:
  using reference       = T&;
  using const_reference = std::add_const_t<T>&;
  using size_type       = size_t;
  using difference_type = ptrdiff_t;
  using value_type      = T;
  using pointer         = T*;
  using const_pointer   = std::add_const_t<T>*;

  KOKKOS_INLINE_FUNCTION static constexpr size_type size() { return N; }
  KOKKOS_INLINE_FUNCTION static constexpr bool empty() { return false; }
  KOKKOS_INLINE_FUNCTION constexpr size_type max_size() const { return N; }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION constexpr reference operator[](const iType& i) {
    static_assert(
        (std::is_integral<iType>::value || std::is_enum<iType>::value),
        "Must be integral argument");
    KOKKOS_ARRAY_BOUNDS_CHECK(i, N);
    return m_internal_implementation_private_member_data[i];
  }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION constexpr const_reference operator[](
      const iType& i) const {
    static_assert(
        (std::is_integral<iType>::value || std::is_enum<iType>::value),
        "Must be integral argument");
    KOKKOS_ARRAY_BOUNDS_CHECK(i, N);
    return m_internal_implementation_private_member_data[i];
  }

  KOKKOS_INLINE_FUNCTION constexpr pointer data() {
    return &m_internal_implementation_private_member_data[0];
  }
  KOKKOS_INLINE_FUNCTION constexpr const_pointer data() const {
    return &m_internal_implementation_private_member_data[0];
  }
};

template <class T, class Proxy>
struct Array<T, 0, Proxy> {
 public:
  using reference       = T&;
  using const_reference = std::add_const_t<T>&;
  using size_type       = size_t;
  using difference_type = ptrdiff_t;
  using value_type      = T;
  using pointer         = T*;
  using const_pointer   = std::add_const_t<T>*;

  KOKKOS_INLINE_FUNCTION static constexpr size_type size() { return 0; }
  KOKKOS_INLINE_FUNCTION static constexpr bool empty() { return true; }
  KOKKOS_INLINE_FUNCTION constexpr size_type max_size() const { return 0; }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION reference operator[](const iType&) {
    static_assert(
        (std::is_integral<iType>::value || std::is_enum<iType>::value),
        "Must be integer argument");
    Kokkos::abort("Unreachable code");
    return *reinterpret_cast<pointer>(-1);
  }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION const_reference operator[](const iType&) const {
    static_assert(
        (std::is_integral<iType>::value || std::is_enum<iType>::value),
        "Must be integer argument");
    Kokkos::abort("Unreachable code");
    return *reinterpret_cast<const_pointer>(-1);
  }

  KOKKOS_INLINE_FUNCTION pointer data() { return nullptr; }
  KOKKOS_INLINE_FUNCTION const_pointer data() const { return nullptr; }

  KOKKOS_DEFAULTED_FUNCTION ~Array()            = default;
  KOKKOS_DEFAULTED_FUNCTION Array()             = default;
  KOKKOS_DEFAULTED_FUNCTION Array(const Array&) = default;
  KOKKOS_DEFAULTED_FUNCTION Array& operator=(const Array&) = default;

  // Some supported compilers are not sufficiently C++11 compliant
  // for default move constructor and move assignment operator.
  // Array( Array && ) = default ;
  // Array & operator = ( Array && ) = default ;
};

template <>
struct Array<void, KOKKOS_INVALID_INDEX, void> {
  struct contiguous {};
  struct strided {};
};

template <class T>
struct Array<T, KOKKOS_INVALID_INDEX, Array<>::contiguous> {
 private:
  T* m_elem;
  size_t m_size;

 public:
  using reference       = T&;
  using const_reference = std::add_const_t<T>&;
  using size_type       = size_t;
  using difference_type = ptrdiff_t;
  using value_type      = T;
  using pointer         = T*;
  using const_pointer   = std::add_const_t<T>*;

  KOKKOS_INLINE_FUNCTION constexpr size_type size() const { return m_size; }
  KOKKOS_INLINE_FUNCTION constexpr bool empty() const { return 0 == m_size; }
  KOKKOS_INLINE_FUNCTION constexpr size_type max_size() const { return m_size; }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION reference operator[](const iType& i) {
    static_assert(
        (std::is_integral<iType>::value || std::is_enum<iType>::value),
        "Must be integral argument");
    KOKKOS_ARRAY_BOUNDS_CHECK(i, m_size);
    return m_elem[i];
  }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION const_reference operator[](const iType& i) const {
    static_assert(
        (std::is_integral<iType>::value || std::is_enum<iType>::value),
        "Must be integral argument");
    KOKKOS_ARRAY_BOUNDS_CHECK(i, m_size);
    return m_elem[i];
  }

  KOKKOS_INLINE_FUNCTION pointer data() { return m_elem; }
  KOKKOS_INLINE_FUNCTION const_pointer data() const { return m_elem; }

  KOKKOS_DEFAULTED_FUNCTION ~Array()                     = default;
  KOKKOS_INLINE_FUNCTION_DELETED Array()                 = delete;
  KOKKOS_INLINE_FUNCTION_DELETED Array(const Array& rhs) = delete;

  // Some supported compilers are not sufficiently C++11 compliant
  // for default move constructor and move assignment operator.
  // Array( Array && rhs ) = default ;
  // Array & operator = ( Array && rhs ) = delete ;

  KOKKOS_INLINE_FUNCTION
  Array& operator=(const Array& rhs) {
    const size_t n = size() < rhs.size() ? size() : rhs.size();
    for (size_t i = 0; i < n; ++i) m_elem[i] = rhs[i];
    return *this;
  }

  template <size_t N, class P>
  KOKKOS_INLINE_FUNCTION Array& operator=(const Array<T, N, P>& rhs) {
    const size_t n = size() < rhs.size() ? size() : rhs.size();
    for (size_t i = 0; i < n; ++i) m_elem[i] = rhs[i];
    return *this;
  }

  KOKKOS_INLINE_FUNCTION constexpr Array(pointer arg_ptr, size_type arg_size,
                                         size_type = 0)
      : m_elem(arg_ptr), m_size(arg_size) {}
};

template <class T>
struct Array<T, KOKKOS_INVALID_INDEX, Array<>::strided> {
 private:
  T* m_elem;
  size_t m_size;
  size_t m_stride;

 public:
  using reference       = T&;
  using const_reference = std::add_const_t<T>&;
  using size_type       = size_t;
  using difference_type = ptrdiff_t;
  using value_type      = T;
  using pointer         = T*;
  using const_pointer   = std::add_const_t<T>*;

  KOKKOS_INLINE_FUNCTION constexpr size_type size() const { return m_size; }
  KOKKOS_INLINE_FUNCTION constexpr bool empty() const { return 0 == m_size; }
  KOKKOS_INLINE_FUNCTION constexpr size_type max_size() const { return m_size; }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION reference operator[](const iType& i) {
    static_assert(
        (std::is_integral<iType>::value || std::is_enum<iType>::value),
        "Must be integral argument");
    KOKKOS_ARRAY_BOUNDS_CHECK(i, m_size);
    return m_elem[i * m_stride];
  }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION const_reference operator[](const iType& i) const {
    static_assert(
        (std::is_integral<iType>::value || std::is_enum<iType>::value),
        "Must be integral argument");
    KOKKOS_ARRAY_BOUNDS_CHECK(i, m_size);
    return m_elem[i * m_stride];
  }

  KOKKOS_INLINE_FUNCTION pointer data() { return m_elem; }
  KOKKOS_INLINE_FUNCTION const_pointer data() const { return m_elem; }

  KOKKOS_DEFAULTED_FUNCTION ~Array()                 = default;
  KOKKOS_INLINE_FUNCTION_DELETED Array()             = delete;
  KOKKOS_INLINE_FUNCTION_DELETED Array(const Array&) = delete;

  // Some supported compilers are not sufficiently C++11 compliant
  // for default move constructor and move assignment operator.
  // Array( Array && rhs ) = default ;
  // Array & operator = ( Array && rhs ) = delete ;

  KOKKOS_INLINE_FUNCTION
  Array& operator=(const Array& rhs) {
    const size_t n = size() < rhs.size() ? size() : rhs.size();
    for (size_t i = 0; i < n; ++i) m_elem[i * m_stride] = rhs[i];
    return *this;
  }

  template <size_t N, class P>
  KOKKOS_INLINE_FUNCTION Array& operator=(const Array<T, N, P>& rhs) {
    const size_t n = size() < rhs.size() ? size() : rhs.size();
    for (size_t i = 0; i < n; ++i) m_elem[i * m_stride] = rhs[i];
    return *this;
  }

  KOKKOS_INLINE_FUNCTION constexpr Array(pointer arg_ptr, size_type arg_size,
                                         size_type arg_stride)
      : m_elem(arg_ptr), m_size(arg_size), m_stride(arg_stride) {}
};

template <typename T, typename... Us>
Array(T, Us...)->Array<T, 1 + sizeof...(Us)>;

}  // namespace Kokkos

//<editor-fold desc="Support for structured binding">
template <class T, std::size_t N>
struct std::tuple_size<Kokkos::Array<T, N>>
    : std::integral_constant<std::size_t, N> {};

template <std::size_t I, class T, std::size_t N>
struct std::tuple_element<I, Kokkos::Array<T, N>> {
  using type = T;
};

namespace Kokkos {

template <std::size_t I, class T, std::size_t N>
KOKKOS_FUNCTION constexpr T& get(Array<T, N>& a) noexcept {
  return a[I];
}

template <std::size_t I, class T, std::size_t N>
KOKKOS_FUNCTION constexpr T const& get(Array<T, N> const& a) noexcept {
  return a[I];
}

template <std::size_t I, class T, std::size_t N>
KOKKOS_FUNCTION constexpr T&& get(Array<T, N>&& a) noexcept {
  return std::move(a[I]);
}

template <std::size_t I, class T, std::size_t N>
KOKKOS_FUNCTION constexpr T const&& get(Array<T, N> const&& a) noexcept {
  return std::move(a[I]);
}

}  // namespace Kokkos
//</editor-fold>

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_ARRAY
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_ARRAY
#endif
#endif /* #ifndef KOKKOS_ARRAY_HPP */
