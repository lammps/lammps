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

#ifndef KOKKOSTRAITS_HPP
#define KOKKOSTRAITS_HPP

#include <cstddef>
#include <cstdint>
#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_BitOps.hpp>
#include <string>
#include <type_traits>

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
// Help with C++11 variadic argument packs

template <unsigned I, typename... Pack>
struct get_type {
  using type = void;
};

template <typename T, typename... Pack>
struct get_type<0, T, Pack...> {
  using type = T;
};

template <unsigned I, typename T, typename... Pack>
struct get_type<I, T, Pack...> {
  using type = typename get_type<I - 1, Pack...>::type;
};

template <typename T, typename... Pack>
struct has_type {
  enum : bool { value = false };
};

template <typename T, typename S, typename... Pack>
struct has_type<T, S, Pack...> {
 private:
  enum { self_value = std::is_same<T, S>::value };

  using next = has_type<T, Pack...>;

  static_assert(
      !(self_value && next::value),
      "Error: more than one member of the argument pack matches the type");

 public:
  enum : bool { value = self_value || next::value };
};

template <typename DefaultType, template <typename> class Condition,
          typename... Pack>
struct has_condition {
  enum : bool { value = false };
  using type = DefaultType;
};

template <typename DefaultType, template <typename> class Condition, typename S,
          typename... Pack>
struct has_condition<DefaultType, Condition, S, Pack...> {
 private:
  enum { self_value = Condition<S>::value };

  using next = has_condition<DefaultType, Condition, Pack...>;

  static_assert(
      !(self_value && next::value),
      "Error: more than one member of the argument pack satisfies condition");

 public:
  enum : bool { value = self_value || next::value };

  using type = std::conditional_t<self_value, S, typename next::type>;
};

template <class... Args>
struct are_integral {
  enum : bool { value = true };
};

template <typename T, class... Args>
struct are_integral<T, Args...> {
  enum {
    value =
        // Accept std::is_integral OR std::is_enum as an integral value
        // since a simple enum value is automically convertible to an
        // integral value.
    (std::is_integral<T>::value || std::is_enum<T>::value) &&
    are_integral<Args...>::value
  };
};

//----------------------------------------------------------------------------
}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Other traits

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
// if_

template <bool Cond, typename TrueType, typename FalseType>
struct if_c {
  enum : bool { value = Cond };

  using type = FalseType;

  using value_type = std::remove_const_t<std::remove_reference_t<type>>;

  using const_value_type = std::add_const_t<value_type>;

  static KOKKOS_INLINE_FUNCTION const_value_type& select(const_value_type& v) {
    return v;
  }

  static KOKKOS_INLINE_FUNCTION value_type& select(value_type& v) { return v; }

  template <class T>
  static KOKKOS_INLINE_FUNCTION value_type& select(const T&) {
    value_type* ptr(0);
    return *ptr;
  }

  template <class T>
  static KOKKOS_INLINE_FUNCTION const_value_type& select(const T&,
                                                         const_value_type& v) {
    return v;
  }

  template <class T>
  static KOKKOS_INLINE_FUNCTION value_type& select(const T&, value_type& v) {
    return v;
  }
};

template <typename TrueType, typename FalseType>
struct if_c<true, TrueType, FalseType> {
  enum : bool { value = true };

  using type = TrueType;

  using value_type = std::remove_const_t<std::remove_reference_t<type>>;

  using const_value_type = std::add_const_t<value_type>;

  static KOKKOS_INLINE_FUNCTION const_value_type& select(const_value_type& v) {
    return v;
  }

  static KOKKOS_INLINE_FUNCTION value_type& select(value_type& v) { return v; }

  template <class T>
  static KOKKOS_INLINE_FUNCTION value_type& select(const T&) {
    value_type* ptr(0);
    return *ptr;
  }

  template <class F>
  static KOKKOS_INLINE_FUNCTION const_value_type& select(const_value_type& v,
                                                         const F&) {
    return v;
  }

  template <class F>
  static KOKKOS_INLINE_FUNCTION value_type& select(value_type& v, const F&) {
    return v;
  }
};

template <typename TrueType>
struct if_c<false, TrueType, void> {
  enum : bool { value = false };

  using type       = void;
  using value_type = void;
};

template <typename FalseType>
struct if_c<true, void, FalseType> {
  enum : bool { value = true };

  using type       = void;
  using value_type = void;
};

//----------------------------------------------------------------------------
// These 'constexpr'functions can be used as
// both regular functions and meta-function.

/**\brief  There exists integral 'k' such that N = 2^k */
KOKKOS_INLINE_FUNCTION
constexpr bool is_integral_power_of_two(const size_t N) {
  return (0 < N) && (0 == (N & (N - 1)));
}

/**\brief  Return integral 'k' such that N = 2^k, assuming valid.  */
KOKKOS_INLINE_FUNCTION
constexpr unsigned integral_power_of_two_assume_valid(const size_t N) {
  return N == 1 ? 0 : 1 + integral_power_of_two_assume_valid(N >> 1);
}

/**\brief  Return integral 'k' such that N = 2^k, if exists.
 *         If does not exist return ~0u.
 */
KOKKOS_INLINE_FUNCTION
constexpr unsigned integral_power_of_two(const size_t N) {
  return is_integral_power_of_two(N) ? integral_power_of_two_assume_valid(N)
                                     : ~0u;
}

/** \brief  If power of two then return power,
 *          otherwise return ~0u.
 */
KOKKOS_FORCEINLINE_FUNCTION
unsigned power_of_two_if_valid(const unsigned N) {
  unsigned p = ~0u;
  if (is_integral_power_of_two(N)) {
    p = bit_scan_forward(N);
  }
  return p;
}

//----------------------------------------------------------------------------

template <typename T, T v, bool NonZero = (v != T(0))>
struct integral_nonzero_constant {
  // Declaration of 'static const' causes an unresolved linker symbol in debug
  // static const T value = v ;
  enum { value = T(v) };
  using value_type = T;
  using type       = integral_nonzero_constant<T, v>;
  KOKKOS_INLINE_FUNCTION integral_nonzero_constant(const T&) {}
};

template <typename T, T zero>
struct integral_nonzero_constant<T, zero, false> {
  const T value;
  using value_type = T;
  using type       = integral_nonzero_constant<T, 0>;
  KOKKOS_INLINE_FUNCTION integral_nonzero_constant(const T& v) : value(v) {}
};

//----------------------------------------------------------------------------

template <class T>
struct make_all_extents_into_pointers {
  using type = T;
};

template <class T, unsigned N>
struct make_all_extents_into_pointers<T[N]> {
  using type = typename make_all_extents_into_pointers<T>::type*;
};

template <class T>
struct make_all_extents_into_pointers<T*> {
  using type = typename make_all_extents_into_pointers<T>::type*;
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSTRAITS_HPP */
