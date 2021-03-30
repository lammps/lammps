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
    : std::integral_constant<
          bool, is_compatible_type_erasure<TSrc, TDst>::value &&
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
