// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSTRAITS_HPP
#define KOKKOSTRAITS_HPP

#include <cstddef>
#include <cstdint>
#include <Kokkos_Macros.hpp>
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
  enum { self_value = std::is_same_v<T, S> };

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
    (std::is_integral_v<T> || std::is_enum_v<T>)&&are_integral<Args...>::value
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
