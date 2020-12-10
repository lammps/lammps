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

#ifndef KOKKOS_ARRAY_HPP
#define KOKKOS_ARRAY_HPP

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_Error.hpp>

#include <type_traits>
#include <algorithm>
#include <limits>
#include <cstddef>
#include <string>

namespace Kokkos {

#ifdef KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK
namespace Impl {
template <typename Integral, bool Signed = std::is_signed<Integral>::value>
struct ArrayBoundsCheck;

template <typename Integral>
struct ArrayBoundsCheck<Integral, true> {
  KOKKOS_INLINE_FUNCTION
  ArrayBoundsCheck(Integral i, size_t N) {
    if (i < 0) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
      std::string s = "Kokkos::Array: index ";
      s += std::to_string(i);
      s += " < 0";
      Kokkos::Impl::throw_runtime_exception(s);
#else
      Kokkos::abort("Kokkos::Array: negative index in device code");
#endif
    }
    ArrayBoundsCheck<Integral, false>(i, N);
  }
};

template <typename Integral>
struct ArrayBoundsCheck<Integral, false> {
  KOKKOS_INLINE_FUNCTION
  ArrayBoundsCheck(Integral i, size_t N) {
    if (size_t(i) >= N) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
      std::string s = "Kokkos::Array: index ";
      s += std::to_string(i);
      s += " >= ";
      s += std::to_string(N);
      Kokkos::Impl::throw_runtime_exception(s);
#else
      Kokkos::abort("Kokkos::Array: index >= size");
#endif
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
  using const_reference = typename std::add_const<T>::type&;
  using size_type       = size_t;
  using difference_type = ptrdiff_t;
  using value_type      = T;
  using pointer         = T*;
  using const_pointer   = typename std::add_const<T>::type*;

  KOKKOS_INLINE_FUNCTION static constexpr size_type size() { return N; }
  KOKKOS_INLINE_FUNCTION static constexpr bool empty() { return false; }
  KOKKOS_INLINE_FUNCTION constexpr size_type max_size() const { return N; }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION reference operator[](const iType& i) {
    static_assert(
        (std::is_integral<iType>::value || std::is_enum<iType>::value),
        "Must be integral argument");
    KOKKOS_ARRAY_BOUNDS_CHECK(i, N);
    return m_internal_implementation_private_member_data[i];
  }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION const_reference operator[](const iType& i) const {
    static_assert(
        (std::is_integral<iType>::value || std::is_enum<iType>::value),
        "Must be integral argument");
    KOKKOS_ARRAY_BOUNDS_CHECK(i, N);
    return m_internal_implementation_private_member_data[i];
  }

  KOKKOS_INLINE_FUNCTION pointer data() {
    return &m_internal_implementation_private_member_data[0];
  }
  KOKKOS_INLINE_FUNCTION const_pointer data() const {
    return &m_internal_implementation_private_member_data[0];
  }
};

template <class T, class Proxy>
struct Array<T, 0, Proxy> {
 public:
  using reference       = T&;
  using const_reference = typename std::add_const<T>::type&;
  using size_type       = size_t;
  using difference_type = ptrdiff_t;
  using value_type      = T;
  using pointer         = T*;
  using const_pointer   = typename std::add_const<T>::type*;

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

  KOKKOS_INLINE_FUNCTION pointer data() { return pointer(0); }
  KOKKOS_INLINE_FUNCTION const_pointer data() const { return const_pointer(0); }

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
  using const_reference = typename std::add_const<T>::type&;
  using size_type       = size_t;
  using difference_type = ptrdiff_t;
  using value_type      = T;
  using pointer         = T*;
  using const_pointer   = typename std::add_const<T>::type*;

  KOKKOS_INLINE_FUNCTION constexpr size_type size() const { return m_size; }
  KOKKOS_INLINE_FUNCTION constexpr bool empty() const { return 0 != m_size; }
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
    const size_t n = std::min(m_size, rhs.size());
    for (size_t i = 0; i < n; ++i) m_elem[i] = rhs[i];
    return *this;
  }

  template <size_t N, class P>
  KOKKOS_INLINE_FUNCTION Array& operator=(const Array<T, N, P>& rhs) {
    const size_t n = std::min(m_size, rhs.size());
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
  using const_reference = typename std::add_const<T>::type&;
  using size_type       = size_t;
  using difference_type = ptrdiff_t;
  using value_type      = T;
  using pointer         = T*;
  using const_pointer   = typename std::add_const<T>::type*;

  KOKKOS_INLINE_FUNCTION constexpr size_type size() const { return m_size; }
  KOKKOS_INLINE_FUNCTION constexpr bool empty() const { return 0 != m_size; }
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
    const size_t n = std::min(m_size, rhs.size());
    for (size_t i = 0; i < n; ++i) m_elem[i] = rhs[i];
    return *this;
  }

  template <size_t N, class P>
  KOKKOS_INLINE_FUNCTION Array& operator=(const Array<T, N, P>& rhs) {
    const size_t n = std::min(m_size, rhs.size());
    for (size_t i = 0; i < n; ++i) m_elem[i] = rhs[i];
    return *this;
  }

  KOKKOS_INLINE_FUNCTION constexpr Array(pointer arg_ptr, size_type arg_size,
                                         size_type arg_stride)
      : m_elem(arg_ptr), m_size(arg_size), m_stride(arg_stride) {}
};

}  // namespace Kokkos

#endif /* #ifndef KOKKOS_ARRAY_HPP */
