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

#ifndef KOKKOS_KOKKOS_GRAPHIMPL_UTILITIES_HPP
#define KOKKOS_KOKKOS_GRAPHIMPL_UTILITIES_HPP

#include <Kokkos_Macros.hpp>

#include <Kokkos_Graph_fwd.hpp>

#include <type_traits>

namespace Kokkos {
namespace Impl {

//==============================================================================
// <editor-fold desc="is_compatible_type_erasure"> {{{1

template <class Src, class Dst, class Enable = void>
struct is_compatible_type_erasure : std::false_type {};

template <class T>
struct is_compatible_type_erasure<T, Kokkos::Experimental::TypeErasedTag>
    : std::true_type {};

template <>
struct is_compatible_type_erasure<Kokkos::Experimental::TypeErasedTag,
                                  Kokkos::Experimental::TypeErasedTag>
    : std::true_type {};

template <class T>
struct is_compatible_type_erasure<T, T> : std::true_type {};

// So there are a couple of ways we could do this, but I didn't want to set up
// all of the machinery to do a lazy instantiation of the convertibility
// condition in the converting constructor of GraphNodeRef, so I'm going with
// this for now:
// TODO @desul-integration make this variadic once we have a meta-conjunction
template <template <class, class, class> class Template, class TSrc, class USrc,
          class VSrc, class TDst, class UDst, class VDst>
struct is_compatible_type_erasure<
    Template<TSrc, USrc, VSrc>, Template<TDst, UDst, VDst>,
    // Because gcc thinks this is ambiguous, we need to add this:
    std::enable_if_t<!std::is_same<TSrc, TDst>::value ||
                     !std::is_same<USrc, UDst>::value ||
                     !std::is_same<VSrc, VDst>::value>>
    : std::bool_constant<is_compatible_type_erasure<TSrc, TDst>::value &&
                         is_compatible_type_erasure<USrc, UDst>::value &&
                         is_compatible_type_erasure<VSrc, VDst>::value> {};

// </editor-fold> end is_compatible_type_erasure }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="is_more_type_erased"> {{{1

template <class T, class U>
struct is_more_type_erased : std::false_type {};

template <class T>
struct is_more_type_erased<Kokkos::Experimental::TypeErasedTag, T>
    : std::true_type {};

template <>
struct is_more_type_erased<Kokkos::Experimental::TypeErasedTag,
                           Kokkos::Experimental::TypeErasedTag>
    : std::false_type {};

// TODO @desul-integration variadic version of this, like the above

// </editor-fold> end is_more_type_erased }}}1
//==============================================================================

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_GRAPHIMPL_UTILITIES_HPP
