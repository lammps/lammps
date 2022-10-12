/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//              Copyright (2019) Sandia Corporation
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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
#ifndef KOKKOS_ENABLE_DEPRECATED_CODE_3
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#else
KOKKOS_IMPL_WARNING("Including non-public Kokkos header files is not allowed.")
#endif
#endif
#ifndef KOKKOS_KOKKOS_EXTENTS_HPP
#define KOKKOS_KOKKOS_EXTENTS_HPP

#include <cstddef>
#include <type_traits>
#include <Kokkos_Macros.hpp>

namespace Kokkos {
namespace Experimental {

constexpr ptrdiff_t dynamic_extent = -1;

template <ptrdiff_t... ExtentSpecs>
struct Extents {
  /* TODO @enhancement flesh this out more */
};

template <class Exts, ptrdiff_t NewExtent>
struct PrependExtent;

template <ptrdiff_t... Exts, ptrdiff_t NewExtent>
struct PrependExtent<Extents<Exts...>, NewExtent> {
  using type = Extents<NewExtent, Exts...>;
};

template <class Exts, ptrdiff_t NewExtent>
struct AppendExtent;

template <ptrdiff_t... Exts, ptrdiff_t NewExtent>
struct AppendExtent<Extents<Exts...>, NewExtent> {
  using type = Extents<Exts..., NewExtent>;
};

}  // end namespace Experimental

namespace Impl {

namespace _parse_view_extents_impl {

template <class T>
struct _all_remaining_extents_dynamic : std::true_type {};

template <class T>
struct _all_remaining_extents_dynamic<T*> : _all_remaining_extents_dynamic<T> {
};

template <class T, unsigned N>
struct _all_remaining_extents_dynamic<T[N]> : std::false_type {};

template <class T, class Result, class = void>
struct _parse_impl {
  using type = Result;
};

// We have to treat the case of int**[x] specially, since it *doesn't* go
// backwards
template <class T, ptrdiff_t... ExtentSpec>
struct _parse_impl<T*, Kokkos::Experimental::Extents<ExtentSpec...>,
                   std::enable_if_t<_all_remaining_extents_dynamic<T>::value>>
    : _parse_impl<T, Kokkos::Experimental::Extents<
                         Kokkos::Experimental::dynamic_extent, ExtentSpec...>> {
};

// int*(*[x])[y] should still work also (meaning int[][x][][y])
template <class T, ptrdiff_t... ExtentSpec>
struct _parse_impl<
    T*, Kokkos::Experimental::Extents<ExtentSpec...>,
    std::enable_if_t<!_all_remaining_extents_dynamic<T>::value>> {
  using _next = Kokkos::Experimental::AppendExtent<
      typename _parse_impl<T, Kokkos::Experimental::Extents<ExtentSpec...>,
                           void>::type,
      Kokkos::Experimental::dynamic_extent>;
  using type = typename _next::type;
};

template <class T, ptrdiff_t... ExtentSpec, unsigned N>
struct _parse_impl<T[N], Kokkos::Experimental::Extents<ExtentSpec...>, void>
    : _parse_impl<
          T, Kokkos::Experimental::Extents<ExtentSpec...,
                                           ptrdiff_t(N)>  // TODO @pedantic this
                                                          // could be a
                                                          // narrowing cast
          > {};

}  // end namespace _parse_view_extents_impl

template <class DataType>
struct ParseViewExtents {
  using type = typename _parse_view_extents_impl ::_parse_impl<
      DataType, Kokkos::Experimental::Extents<>>::type;
};

template <class ValueType, ptrdiff_t Ext>
struct ApplyExtent {
  using type = ValueType[Ext];
};

template <class ValueType>
struct ApplyExtent<ValueType, Kokkos::Experimental::dynamic_extent> {
  using type = ValueType*;
};

template <class ValueType, unsigned N, ptrdiff_t Ext>
struct ApplyExtent<ValueType[N], Ext> {
  using type = typename ApplyExtent<ValueType, Ext>::type[N];
};

template <class ValueType, ptrdiff_t Ext>
struct ApplyExtent<ValueType*, Ext> {
  using type = ValueType * [Ext];
};

template <class ValueType>
struct ApplyExtent<ValueType*, Kokkos::Experimental::dynamic_extent> {
  using type =
      typename ApplyExtent<ValueType,
                           Kokkos::Experimental::dynamic_extent>::type*;
};

template <class ValueType, unsigned N>
struct ApplyExtent<ValueType[N], Kokkos::Experimental::dynamic_extent> {
  using type =
      typename ApplyExtent<ValueType,
                           Kokkos::Experimental::dynamic_extent>::type[N];
};

}  // end namespace Impl

}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_EXTENTS_HPP
