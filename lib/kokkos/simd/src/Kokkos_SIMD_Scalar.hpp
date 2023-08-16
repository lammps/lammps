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

#ifndef KOKKOS_SIMD_SCALAR_HPP
#define KOKKOS_SIMD_SCALAR_HPP

#include <type_traits>
#include <climits>
#include <cfloat>

#include <Kokkos_SIMD_Common.hpp>

namespace Kokkos {
namespace Experimental {

namespace simd_abi {

class scalar {};

}  // namespace simd_abi

template <class T>
class simd_mask<T, simd_abi::scalar> {
  bool m_value;

 public:
  using value_type                      = bool;
  using simd_type                       = simd<T, simd_abi::scalar>;
  using abi_type                        = simd_abi::scalar;
  using reference                       = value_type&;
  KOKKOS_DEFAULTED_FUNCTION simd_mask() = default;
  KOKKOS_FORCEINLINE_FUNCTION static constexpr std::size_t size() { return 1; }
  KOKKOS_FORCEINLINE_FUNCTION explicit simd_mask(value_type value)
      : m_value(value) {}
  template <class U>
  KOKKOS_FORCEINLINE_FUNCTION simd_mask(
      simd_mask<U, simd_abi::scalar> const& other)
      : m_value(static_cast<bool>(other)) {}
  KOKKOS_FORCEINLINE_FUNCTION constexpr explicit operator bool() const {
    return m_value;
  }
  KOKKOS_FORCEINLINE_FUNCTION reference operator[](std::size_t) {
    return m_value;
  }
  KOKKOS_FORCEINLINE_FUNCTION value_type operator[](std::size_t) const {
    return m_value;
  }
  KOKKOS_FORCEINLINE_FUNCTION simd_mask
  operator||(simd_mask const& other) const {
    return simd_mask(m_value || other.m_value);
  }
  KOKKOS_FORCEINLINE_FUNCTION simd_mask
  operator&&(simd_mask const& other) const {
    return simd_mask(m_value && other.m_value);
  }
  KOKKOS_FORCEINLINE_FUNCTION simd_mask operator!() const {
    return simd_mask(!m_value);
  }
  KOKKOS_FORCEINLINE_FUNCTION bool operator==(simd_mask const& other) const {
    return m_value == other.m_value;
  }
  KOKKOS_FORCEINLINE_FUNCTION bool operator!=(simd_mask const& other) const {
    return m_value != other.m_value;
  }
};

template <class T>
class simd<T, simd_abi::scalar> {
  T m_value;

 public:
  using value_type                            = T;
  using abi_type                              = simd_abi::scalar;
  using mask_type                             = simd_mask<T, abi_type>;
  using reference                             = value_type&;
  KOKKOS_DEFAULTED_FUNCTION simd()            = default;
  KOKKOS_DEFAULTED_FUNCTION simd(simd const&) = default;
  KOKKOS_DEFAULTED_FUNCTION simd(simd&&)      = default;
  KOKKOS_DEFAULTED_FUNCTION simd& operator=(simd const&) = default;
  KOKKOS_DEFAULTED_FUNCTION simd& operator=(simd&&) = default;
  KOKKOS_FORCEINLINE_FUNCTION static constexpr std::size_t size() { return 1; }
  template <class U, std::enable_if_t<std::is_convertible_v<U, value_type>,
                                      bool> = false>
  KOKKOS_FORCEINLINE_FUNCTION simd(U&& value) : m_value(value) {}
  template <class U, std::enable_if_t<std::is_convertible_v<U, value_type>,
                                      bool> = false>
  KOKKOS_FORCEINLINE_FUNCTION explicit simd(simd<U, abi_type> const& other)
      : m_value(static_cast<U>(other)) {}
  template <class G,
            std::enable_if_t<
                // basically, can you do { value_type r =
                // gen(std::integral_constant<std::size_t, i>()); }
                std::is_invocable_r_v<value_type, G,
                                      std::integral_constant<std::size_t, 0>>,
                bool> = false>
  KOKKOS_FORCEINLINE_FUNCTION simd(G&& gen)
      : m_value(gen(std::integral_constant<std::size_t, 0>())) {}
  KOKKOS_FORCEINLINE_FUNCTION simd operator-() const { return simd(-m_value); }
  KOKKOS_FORCEINLINE_FUNCTION simd operator>>(int rhs) const {
    return simd(m_value >> rhs);
  }
  KOKKOS_FORCEINLINE_FUNCTION simd
  operator>>(simd<int, abi_type> const& rhs) const {
    return simd(m_value >> static_cast<int>(rhs));
  }
  KOKKOS_FORCEINLINE_FUNCTION simd operator<<(int rhs) const {
    return simd(m_value << rhs);
  }
  KOKKOS_FORCEINLINE_FUNCTION simd
  operator<<(simd<int, abi_type> const& rhs) const {
    return simd(m_value << static_cast<int>(rhs));
  }
  KOKKOS_FORCEINLINE_FUNCTION simd operator&(simd const& other) const {
    return m_value & other.m_value;
  }
  KOKKOS_FORCEINLINE_FUNCTION simd operator|(simd const& other) const {
    return m_value | other.m_value;
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr explicit operator T() const {
    return m_value;
  }
  KOKKOS_FORCEINLINE_FUNCTION mask_type operator<(simd const& other) const {
    return mask_type(m_value < other.m_value);
  }
  KOKKOS_FORCEINLINE_FUNCTION mask_type operator>(simd const& other) const {
    return mask_type(m_value > other.m_value);
  }
  KOKKOS_FORCEINLINE_FUNCTION mask_type operator<=(simd const& other) const {
    return mask_type(m_value <= other.m_value);
  }
  KOKKOS_FORCEINLINE_FUNCTION mask_type operator>=(simd const& other) const {
    return mask_type(m_value >= other.m_value);
  }
  KOKKOS_FORCEINLINE_FUNCTION mask_type operator==(simd const& other) const {
    return mask_type(m_value == other.m_value);
  }
  KOKKOS_FORCEINLINE_FUNCTION mask_type operator!=(simd const& other) const {
    return mask_type(m_value != other.m_value);
  }
  KOKKOS_FORCEINLINE_FUNCTION void copy_from(T const* ptr,
                                             element_aligned_tag) {
    m_value = *ptr;
  }
  KOKKOS_FORCEINLINE_FUNCTION void copy_to(T* ptr, element_aligned_tag) const {
    *ptr = m_value;
  }
  KOKKOS_FORCEINLINE_FUNCTION reference operator[](std::size_t) {
    return m_value;
  }
  KOKKOS_FORCEINLINE_FUNCTION value_type operator[](std::size_t) const {
    return m_value;
  }
};

template <class T>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION simd<T, simd_abi::scalar> operator*(
    simd<T, simd_abi::scalar> const& lhs,
    simd<T, simd_abi::scalar> const& rhs) {
  return simd<T, simd_abi::scalar>(static_cast<T>(lhs) * static_cast<T>(rhs));
}

template <class T>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION simd<T, simd_abi::scalar> operator/(
    simd<T, simd_abi::scalar> const& lhs,
    simd<T, simd_abi::scalar> const& rhs) {
  return simd<T, simd_abi::scalar>(static_cast<T>(lhs) / static_cast<T>(rhs));
}

template <class T>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION simd<T, simd_abi::scalar> operator+(
    simd<T, simd_abi::scalar> const& lhs,
    simd<T, simd_abi::scalar> const& rhs) {
  return simd<T, simd_abi::scalar>(static_cast<T>(lhs) + static_cast<T>(rhs));
}

template <class T>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION simd<T, simd_abi::scalar> operator-(
    simd<T, simd_abi::scalar> const& lhs,
    simd<T, simd_abi::scalar> const& rhs) {
  return simd<T, simd_abi::scalar>(static_cast<T>(lhs) - static_cast<T>(rhs));
}

template <class T>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION simd<T, simd_abi::scalar> abs(
    simd<T, simd_abi::scalar> const& a) {
  return simd<T, simd_abi::scalar>(std::abs(static_cast<T>(a)));
}

template <class T>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION simd<T, simd_abi::scalar> sqrt(
    simd<T, simd_abi::scalar> const& a) {
  return simd<T, simd_abi::scalar>(std::sqrt(static_cast<T>(a)));
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION simd<T, simd_abi::scalar> fma(
    simd<T, simd_abi::scalar> const& x, simd<T, simd_abi::scalar> const& y,
    simd<T, simd_abi::scalar> const& z) {
  return simd<T, simd_abi::scalar>((static_cast<T>(x) * static_cast<T>(y)) +
                                   static_cast<T>(z));
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION simd<T, simd_abi::scalar> condition(
    desul::Impl::dont_deduce_this_parameter_t<
        simd_mask<T, simd_abi::scalar>> const& a,
    simd<T, simd_abi::scalar> const& b, simd<T, simd_abi::scalar> const& c) {
  return simd<T, simd_abi::scalar>(static_cast<bool>(a) ? static_cast<T>(b)
                                                        : static_cast<T>(c));
}

template <class T, class Abi>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION simd<T, Abi> copysign(
    simd<T, Abi> const& a, simd<T, Abi> const& b) {
  return std::copysign(static_cast<T>(a), static_cast<T>(b));
}

template <class T>
class const_where_expression<simd_mask<T, simd_abi::scalar>,
                             simd<T, simd_abi::scalar>> {
 public:
  using abi_type   = simd_abi::scalar;
  using value_type = simd<T, abi_type>;
  using mask_type  = simd_mask<T, abi_type>;

 protected:
  value_type& m_value;
  mask_type const& m_mask;

 public:
  KOKKOS_FORCEINLINE_FUNCTION
  const_where_expression(mask_type const& mask_arg, value_type const& value_arg)
      : m_value(const_cast<value_type&>(value_arg)), m_mask(mask_arg) {}

  KOKKOS_FORCEINLINE_FUNCTION
  void copy_to(T* mem, element_aligned_tag) const {
    if (static_cast<bool>(m_mask)) *mem = static_cast<T>(m_value);
  }
  template <class Integral>
  KOKKOS_FORCEINLINE_FUNCTION std::enable_if_t<std::is_integral_v<Integral>>
  scatter_to(T* mem, simd<Integral, simd_abi::scalar> const& index) const {
    if (static_cast<bool>(m_mask))
      mem[static_cast<Integral>(index)] = static_cast<T>(m_value);
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type const&
  impl_get_value() const {
    return m_value;
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type const&
  impl_get_mask() const {
    return m_mask;
  }
};

template <class T>
class where_expression<simd_mask<T, simd_abi::scalar>,
                       simd<T, simd_abi::scalar>>
    : public const_where_expression<simd_mask<T, simd_abi::scalar>,
                                    simd<T, simd_abi::scalar>> {
  using base_type = const_where_expression<simd_mask<T, simd_abi::scalar>,
                                           simd<T, simd_abi::scalar>>;

 public:
  using typename base_type::value_type;
  KOKKOS_FORCEINLINE_FUNCTION
  where_expression(simd_mask<T, simd_abi::scalar> const& mask_arg,
                   simd<T, simd_abi::scalar>& value_arg)
      : base_type(mask_arg, value_arg) {}
  KOKKOS_FORCEINLINE_FUNCTION
  void copy_from(T const* mem, element_aligned_tag) {
    if (static_cast<bool>(this->m_mask)) this->m_value = *mem;
  }
  template <class Integral>
  KOKKOS_FORCEINLINE_FUNCTION std::enable_if_t<std::is_integral_v<Integral>>
  gather_from(T const* mem, simd<Integral, simd_abi::scalar> const& index) {
    if (static_cast<bool>(this->m_mask))
      this->m_value = mem[static_cast<Integral>(index)];
  }
  template <class U, std::enable_if_t<
                         std::is_convertible_v<U, simd<T, simd_abi::scalar>>,
                         bool> = false>
  KOKKOS_FORCEINLINE_FUNCTION void operator=(U&& x) {
    if (static_cast<bool>(this->m_mask))
      this->m_value =
          static_cast<simd<T, simd_abi::scalar>>(std::forward<U>(x));
  }
};

template <class T>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION
    where_expression<simd_mask<T, Kokkos::Experimental::simd_abi::scalar>,
                     simd<T, Kokkos::Experimental::simd_abi::scalar>>
    where(typename simd<
              T, Kokkos::Experimental::simd_abi::scalar>::mask_type const& mask,
          simd<T, Kokkos::Experimental::simd_abi::scalar>& value) {
  return where_expression(mask, value);
}

template <class T>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION
    const_where_expression<simd_mask<T, Kokkos::Experimental::simd_abi::scalar>,
                           simd<T, Kokkos::Experimental::simd_abi::scalar>>
    where(typename simd<
              T, Kokkos::Experimental::simd_abi::scalar>::mask_type const& mask,
          simd<T, Kokkos::Experimental::simd_abi::scalar> const& value) {
  return const_where_expression(mask, value);
}

template <class T>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION bool all_of(
    simd_mask<T, Kokkos::Experimental::simd_abi::scalar> const& a) {
  return a == simd_mask<T, Kokkos::Experimental::simd_abi::scalar>(true);
}

template <class T>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION bool any_of(
    simd_mask<T, Kokkos::Experimental::simd_abi::scalar> const& a) {
  return a != simd_mask<T, Kokkos::Experimental::simd_abi::scalar>(false);
}

template <class T>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION bool none_of(
    simd_mask<T, Kokkos::Experimental::simd_abi::scalar> const& a) {
  return a == simd_mask<T, Kokkos::Experimental::simd_abi::scalar>(false);
}

template <class T>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION T
reduce(const_where_expression<simd_mask<T, simd_abi::scalar>,
                              simd<T, simd_abi::scalar>> const& x,
       T identity_element, std::plus<>) {
  return static_cast<bool>(x.impl_get_mask())
             ? static_cast<T>(x.impl_get_value())
             : identity_element;
}

template <class T>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION T
hmax(const_where_expression<simd_mask<T, simd_abi::scalar>,
                            simd<T, simd_abi::scalar>> const& x) {
  return static_cast<bool>(x.impl_get_mask())
             ? static_cast<T>(x.impl_get_value())
             : Kokkos::reduction_identity<T>::max();
}

template <class T>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION T
hmin(const_where_expression<simd_mask<T, simd_abi::scalar>,
                            simd<T, simd_abi::scalar>> const& x) {
  return static_cast<bool>(x.impl_get_mask())
             ? static_cast<T>(x.impl_get_value())
             : Kokkos::reduction_identity<T>::min();
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
