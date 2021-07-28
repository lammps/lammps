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

template <typename T>
struct identity {
  using type = T;
};

template <typename T>
using identity_t = typename identity<T>::type;

struct not_a_type {
  not_a_type()                  = delete;
  ~not_a_type()                 = delete;
  not_a_type(not_a_type const&) = delete;
  void operator=(not_a_type const&) = delete;
};

#if defined(__cpp_lib_void_t)
// since C++17
using std::void_t;
#else
template <class...>
using void_t = void;
#endif

//==============================================================================
// <editor-fold desc="remove_cvref_t"> {{{1

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

// </editor-fold> end remove_cvref_t }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="is_specialization_of"> {{{1

template <class Type, template <class...> class Template, class Enable = void>
struct is_specialization_of : std::false_type {};

template <template <class...> class Template, class... Args>
struct is_specialization_of<Template<Args...>, Template> : std::true_type {};

// </editor-fold> end is_specialization_of }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="Folding emulation"> {{{1

// acts like void for comma fold emulation
struct _fold_comma_emulation_return {};

template <class... Ts>
constexpr KOKKOS_INLINE_FUNCTION _fold_comma_emulation_return
emulate_fold_comma_operator(Ts&&...) noexcept {
  return _fold_comma_emulation_return{};
}

#define KOKKOS_IMPL_FOLD_COMMA_OPERATOR(expr)                                \
  ::Kokkos::Impl::emulate_fold_comma_operator(                               \
      ::std::initializer_list<::Kokkos::Impl::_fold_comma_emulation_return>{ \
          ((expr), ::Kokkos::Impl::_fold_comma_emulation_return{})...})

// </editor-fold> end Folding emulation }}}1
//==============================================================================

//==============================================================================
// destruct_delete is a unique_ptr deleter for objects
// created by placement new into already allocated memory
// by only calling the destructor on the object.
//
// Because unique_ptr never calls its deleter with a nullptr value,
// no need to check if p == nullptr.
//
// Note:  This differs in interface from std::default_delete in that the
// function call operator is templated instead of the class, to make
// it easier to use and disallow specialization.
struct destruct_delete {
  template <typename T>
  KOKKOS_INLINE_FUNCTION constexpr void operator()(T* p) const noexcept {
    p->~T();
  }
};
//==============================================================================

//==============================================================================
// <editor-fold desc="type_list"> {{{1

// An intentionally uninstantiateable type_list for metaprogramming purposes
template <class...>
struct type_list;

// </editor-fold> end type_list }}}1
//==============================================================================

}  // namespace Impl
}  // namespace Kokkos

#endif  // KOKKOS_CORE_IMPL_UTILITIES_HPP
