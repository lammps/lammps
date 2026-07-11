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
#include "macros.hpp"

#include "trait_backports.hpp" // make_index_sequence

namespace MDSPAN_IMPL_STANDARD_NAMESPACE {

//==============================================================================

namespace detail {

template <class... Ts> struct type_list { static constexpr auto size = sizeof...(Ts); };

// Implementation of type_list at() that's heavily optimized for small typelists
template <size_t, class> struct type_at;
template <size_t, class Seq, class=std::make_index_sequence<Seq::size>> struct type_at_large_impl;

template <size_t I, size_t Idx, class T>
struct type_at_entry { };

template <class Result>
struct type_at_assign_op_ignore_rest {
  template <class T>
  type_at_assign_op_ignore_rest<Result> operator=(T&&);
  using type = Result;
};

struct type_at_assign_op_impl {
  template <size_t I, size_t Idx, class T>
  type_at_assign_op_impl operator=(type_at_entry<I, Idx, T>&&);
  template <size_t I, class T>
  type_at_assign_op_ignore_rest<T> operator=(type_at_entry<I, I, T>&&);
};

template <size_t I, class... Ts, size_t... Idxs>
struct type_at_large_impl<I, type_list<Ts...>, std::integer_sequence<size_t, Idxs...>>
  : decltype(
      MDSPAN_IMPL_FOLD_ASSIGN_LEFT(type_at_assign_op_impl{}, /* = ... = */ type_at_entry<I, Idxs, Ts>{})
    )
{ };

template <size_t I, class... Ts>
struct type_at<I, type_list<Ts...>>
    : type_at_large_impl<I, type_list<Ts...>>
{ };

template <class T0, class... Ts>
struct type_at<0, type_list<T0, Ts...>> {
  using type = T0;
};

template <class T0, class T1, class... Ts>
struct type_at<1, type_list<T0, T1, Ts...>> {
  using type = T1;
};

template <class T0, class T1, class T2, class... Ts>
struct type_at<2, type_list<T0, T1, T2, Ts...>> {
  using type = T2;
};

template <class T0, class T1, class T2, class T3, class... Ts>
struct type_at<3, type_list<T0, T1, T2, T3, Ts...>> {
  using type = T3;
};


} // namespace detail

//==============================================================================

} // end namespace MDSPAN_IMPL_STANDARD_NAMESPACE
