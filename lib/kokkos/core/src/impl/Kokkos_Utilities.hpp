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

#ifndef KOKKOS_CORE_IMPL_UTILITIES_HPP
#define KOKKOS_CORE_IMPL_UTILITIES_HPP

#include <Kokkos_Macros.hpp>
#include <cstdint>
#include <type_traits>
#include <initializer_list>  // in-order comma operator fold emulation
#include <utility>           // integer_sequence and friends

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

// same as std::integral_constant but with __host__ __device__ annotations on
// the implicit conversion function and the call operator
template <class T, T v>
struct integral_constant {
  using value_type         = T;
  using type               = integral_constant<T, v>;
  static constexpr T value = v;
  KOKKOS_FUNCTION constexpr operator value_type() const noexcept {
    return value;
  }
  KOKKOS_FUNCTION constexpr value_type operator()() const noexcept {
    return value;
  }
};

//==============================================================================

template <typename... Is>
struct always_true : std::true_type {};

//==============================================================================

#if defined(__cpp_lib_type_identity)
// since C++20
using std::type_identity;
using std::type_identity_t;
#else
template <typename T>
struct type_identity {
  using type = T;
};

template <typename T>
using type_identity_t = typename type_identity<T>::type;
#endif

#if defined(__cpp_lib_remove_cvref)
// since C++20
using std::remove_cvref;
using std::remove_cvref_t;
#else
template <class T>
struct remove_cvref {
  using type = std::remove_cv_t<std::remove_reference_t<T>>;
};

template <class T>
using remove_cvref_t = typename remove_cvref<T>::type;
#endif

//==============================================================================
// <editor-fold desc="is_specialization_of"> {{{1

template <class Type, template <class...> class Template, class Enable = void>
struct is_specialization_of : std::false_type {};

template <template <class...> class Template, class... Args>
struct is_specialization_of<Template<Args...>, Template> : std::true_type {};

// </editor-fold> end is_specialization_of }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="type_list"> {{{1

// An intentionally uninstantiateable type_list for metaprogramming purposes
template <class...>
struct type_list;

//------------------------------------------------------------------------------
// <editor-fold desc="type_list_remove_first"> {{{2

// Currently linear complexity; if we use this a lot, maybe make it better?

template <class Entry, class InList, class OutList>
struct _type_list_remove_first_impl;

template <class Entry, class T, class... Ts, class... OutTs>
struct _type_list_remove_first_impl<Entry, type_list<T, Ts...>,
                                    type_list<OutTs...>>
    : _type_list_remove_first_impl<Entry, type_list<Ts...>,
                                   type_list<OutTs..., T>> {};

template <class Entry, class... Ts, class... OutTs>
struct _type_list_remove_first_impl<Entry, type_list<Entry, Ts...>,
                                    type_list<OutTs...>>
    : _type_list_remove_first_impl<Entry, type_list<>,
                                   type_list<OutTs..., Ts...>> {};

template <class Entry, class... OutTs>
struct _type_list_remove_first_impl<Entry, type_list<>, type_list<OutTs...>>
    : type_identity<type_list<OutTs...>> {};

template <class Entry, class List>
struct type_list_remove_first
    : _type_list_remove_first_impl<Entry, List, type_list<>> {};

// </editor-fold> end type_list_remove_first }}}2
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// <editor-fold desc="type_list_any"> {{{2

template <template <class> class UnaryPred, class List>
struct type_list_any;

template <template <class> class UnaryPred, class... Ts>
struct type_list_any<UnaryPred, type_list<Ts...>>
    : std::bool_constant<(UnaryPred<Ts>::value || ...)> {};

// </editor-fold> end type_list_any }}}2
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// <editor-fold desc="concat_type_list"> {{{2
//  concat_type_list combines types in multiple type_lists

// forward declaration
template <typename... T>
struct concat_type_list;

// alias
template <typename... T>
using concat_type_list_t = typename concat_type_list<T...>::type;

// final instantiation
template <typename... T>
struct concat_type_list<type_list<T...>> {
  using type = type_list<T...>;
};

// combine consecutive type_lists
template <typename... T, typename... U, typename... Tail>
struct concat_type_list<type_list<T...>, type_list<U...>, Tail...>
    : concat_type_list<type_list<T..., U...>, Tail...> {};
// </editor-fold> end concat_type_list }}}2
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// <editor-fold desc="filter_type_list"> {{{2
//  filter_type_list generates type-list of types which satisfy
//  PredicateT<T>::value == ValueT

template <template <typename> class PredicateT, typename TypeListT,
          bool ValueT = true>
struct filter_type_list;

template <template <typename> class PredicateT, typename... T, bool ValueT>
struct filter_type_list<PredicateT, type_list<T...>, ValueT> {
  using type =
      concat_type_list_t<std::conditional_t<PredicateT<T>::value == ValueT,
                                            type_list<T>, type_list<>>...>;
};

template <template <typename> class PredicateT, typename T, bool ValueT = true>
using filter_type_list_t =
    typename filter_type_list<PredicateT, T, ValueT>::type;

// </editor-fold> end filter_type_list }}}2
//------------------------------------------------------------------------------

// </editor-fold> end type_list }}}1
//==============================================================================

//==============================================================================
// The weird !sizeof(F*) to express false is to make the
// expression dependent on the type of F, and thus only applicable
// at instantiation and not first-pass semantic analysis of the
// template definition.
template <typename T>
constexpr bool dependent_false_v = !sizeof(T*);
//==============================================================================

}  // namespace Impl
}  // namespace Kokkos

#endif  // KOKKOS_CORE_IMPL_UTILITIES_HPP
