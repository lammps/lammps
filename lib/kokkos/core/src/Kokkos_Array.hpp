// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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
#include <cstddef>

namespace Kokkos {

#ifdef KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK
namespace Impl {

struct ArrayBoundsCheck {
  KOKKOS_INLINE_FUNCTION
  constexpr ArrayBoundsCheck(size_t i, size_t N) {
    if (i >= N) {
      char err[128] = "Kokkos::Array: index ";
      to_chars_i(err + strlen(err), err + 128, i);
      strcat(err, " >= ");
      to_chars_i(err + strlen(err), err + 128, N);
      Kokkos::abort(err);
    }
  }
};
}  // end namespace Impl

#define KOKKOS_ARRAY_BOUNDS_CHECK(i, N) Kokkos::Impl::ArrayBoundsCheck(i, N)

#else  // !defined( KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK )

#define KOKKOS_ARRAY_BOUNDS_CHECK(i, N) (void)0

#endif  // !defined( KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK )

/**\brief  Derived from the C++17 'std::array'.
 *         Dropping the iterator interface.
 */
template <class T, size_t N>
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

  KOKKOS_INLINE_FUNCTION static constexpr size_type size() noexcept {
    return N;
  }
  KOKKOS_INLINE_FUNCTION static constexpr bool empty() noexcept {
    return false;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type max_size() const noexcept {
    return N;
  }

  KOKKOS_INLINE_FUNCTION constexpr reference operator[](size_type i) {
    KOKKOS_ARRAY_BOUNDS_CHECK(i, N);
    return m_internal_implementation_private_member_data[i];
  }

  KOKKOS_INLINE_FUNCTION constexpr const_reference operator[](
      size_type i) const {
    KOKKOS_ARRAY_BOUNDS_CHECK(i, N);
    return m_internal_implementation_private_member_data[i];
  }

  KOKKOS_INLINE_FUNCTION constexpr pointer data() noexcept {
    return &m_internal_implementation_private_member_data[0];
  }
  KOKKOS_INLINE_FUNCTION constexpr const_pointer data() const noexcept {
    return &m_internal_implementation_private_member_data[0];
  }

  KOKKOS_INLINE_FUNCTION constexpr pointer begin() noexcept { return data(); }

  KOKKOS_INLINE_FUNCTION constexpr const_pointer begin() const noexcept {
    return data();
  }

  KOKKOS_INLINE_FUNCTION constexpr pointer end() noexcept { return data() + N; }

  KOKKOS_INLINE_FUNCTION constexpr const_pointer end() const noexcept {
    return data() + N;
  }

  KOKKOS_INLINE_FUNCTION constexpr const_pointer cbegin() const noexcept {
    return data();
  }

  KOKKOS_INLINE_FUNCTION constexpr const_pointer cend() const noexcept {
    return data() + N;
  }

  friend KOKKOS_FUNCTION constexpr bool operator==(Array const& lhs,
                                                   Array const& rhs) noexcept {
    for (size_t i = 0; i != N; ++i)
      if (lhs[i] != rhs[i]) return false;
    return true;
  }

  friend KOKKOS_FUNCTION constexpr bool operator!=(Array const& lhs,
                                                   Array const& rhs) noexcept {
    return !(lhs == rhs);
  }

 private:
  template <class U = T>
  friend KOKKOS_INLINE_FUNCTION constexpr std::enable_if_t<
      Impl::is_swappable<U>::value>
  kokkos_swap(Array<T, N>& a,
              Array<T, N>& b) noexcept(Impl::is_nothrow_swappable_v<U>) {
    for (std::size_t i = 0; i < N; ++i) {
      kokkos_swap(a[i], b[i]);
    }
  }
};

template <class T>
struct Array<T, 0> {
 public:
  using reference       = T&;
  using const_reference = std::add_const_t<T>&;
  using size_type       = size_t;
  using difference_type = ptrdiff_t;
  using value_type      = T;
  using pointer         = T*;
  using const_pointer   = std::add_const_t<T>*;

  KOKKOS_INLINE_FUNCTION static constexpr size_type size() noexcept {
    return 0;
  }
  KOKKOS_INLINE_FUNCTION static constexpr bool empty() noexcept { return true; }
  KOKKOS_INLINE_FUNCTION constexpr size_type max_size() const noexcept {
    return 0;
  }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION reference operator[](const iType&) {
    static_assert((std::is_integral_v<iType> || std::is_enum_v<iType>),
                  "Must be integer argument");
    Kokkos::abort("Unreachable code");
    return *reinterpret_cast<pointer>(-1);
  }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION const_reference operator[](const iType&) const {
    static_assert((std::is_integral_v<iType> || std::is_enum_v<iType>),
                  "Must be integer argument");
    Kokkos::abort("Unreachable code");
    return *reinterpret_cast<const_pointer>(-1);
  }

  KOKKOS_INLINE_FUNCTION constexpr pointer data() noexcept { return nullptr; }
  KOKKOS_INLINE_FUNCTION constexpr const_pointer data() const noexcept {
    return nullptr;
  }

  KOKKOS_INLINE_FUNCTION constexpr pointer begin() noexcept { return data(); }

  KOKKOS_INLINE_FUNCTION constexpr const_pointer begin() const noexcept {
    return data();
  }

  KOKKOS_INLINE_FUNCTION constexpr pointer end() noexcept { return data(); }

  KOKKOS_INLINE_FUNCTION constexpr const_pointer end() const noexcept {
    return data();
  }

  KOKKOS_INLINE_FUNCTION constexpr const_pointer cbegin() const noexcept {
    return data();
  }

  KOKKOS_INLINE_FUNCTION constexpr const_pointer cend() const noexcept {
    return data();
  }

  friend KOKKOS_FUNCTION constexpr bool operator==(Array const&,
                                                   Array const&) noexcept {
    return true;
  }
  friend KOKKOS_FUNCTION constexpr bool operator!=(Array const&,
                                                   Array const&) noexcept {
    return false;
  }

 private:
  friend KOKKOS_INLINE_FUNCTION constexpr void kokkos_swap(
      Array<T, 0>&, Array<T, 0>&) noexcept {}
};

template <typename T, typename... Us>
Array(T, Us...) -> Array<T, 1 + sizeof...(Us)>;

namespace Impl {

template <typename T, size_t N, size_t... I>
KOKKOS_FUNCTION constexpr Array<std::remove_cv_t<T>, N> to_array_impl(
    T (&a)[N], std::index_sequence<I...>) {
  return {{a[I]...}};
}

template <typename T, size_t N, size_t... I>
KOKKOS_FUNCTION constexpr Array<std::remove_cv_t<T>, N> to_array_impl(
    T (&&a)[N], std::index_sequence<I...>) {
  return {{std::move(a[I])...}};
}

}  // namespace Impl

template <typename T, size_t N>
KOKKOS_FUNCTION constexpr auto to_array(T (&a)[N]) {
  return Impl::to_array_impl(a, std::make_index_sequence<N>{});
}

template <typename T, size_t N>
KOKKOS_FUNCTION constexpr auto to_array(T (&&a)[N]) {
  return Impl::to_array_impl(std::move(a), std::make_index_sequence<N>{});
}

}  // namespace Kokkos

//<editor-fold desc="Support for structured binding">
template <class T, std::size_t N>
struct std::tuple_size<Kokkos::Array<T, N>>
    : std::integral_constant<std::size_t, N> {};

template <std::size_t I, class T, std::size_t N>
struct std::tuple_element<I, Kokkos::Array<T, N>> {
  static_assert(I < N);
  using type = T;
};

namespace Kokkos {

template <std::size_t I, class T, std::size_t N>
KOKKOS_FUNCTION constexpr T& get(Array<T, N>& a) noexcept {
  static_assert(I < N);
  return a[I];
}

template <std::size_t I, class T, std::size_t N>
KOKKOS_FUNCTION constexpr T const& get(Array<T, N> const& a) noexcept {
  static_assert(I < N);
  return a[I];
}

template <std::size_t I, class T, std::size_t N>
KOKKOS_FUNCTION constexpr T&& get(Array<T, N>&& a) noexcept {
  static_assert(I < N);
  return std::move(a[I]);
}

template <std::size_t I, class T, std::size_t N>
KOKKOS_FUNCTION constexpr T const&& get(Array<T, N> const&& a) noexcept {
  static_assert(I < N);
  return std::move(a[I]);
}

}  // namespace Kokkos
//</editor-fold>

//<editor-fold desc="Support for range-based for loop">
namespace Kokkos {

template <class T, std::size_t N>
KOKKOS_FUNCTION constexpr T const* begin(Array<T, N> const& a) noexcept {
  return a.data();
}

template <class T, std::size_t N>
KOKKOS_FUNCTION constexpr T* begin(Array<T, N>& a) noexcept {
  return a.data();
}

template <class T, std::size_t N>
KOKKOS_FUNCTION constexpr T const* end(Array<T, N> const& a) noexcept {
  return a.data() + a.size();
}

template <class T, std::size_t N>
KOKKOS_FUNCTION constexpr T* end(Array<T, N>& a) noexcept {
  return a.data() + a.size();
}

}  // namespace Kokkos
//</editor-fold>

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_ARRAY
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_ARRAY
#endif
#endif /* #ifndef KOKKOS_ARRAY_HPP */
