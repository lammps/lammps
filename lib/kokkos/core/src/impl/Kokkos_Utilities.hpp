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

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------
// C++14 integer sequence
template <typename T, T... Ints>
struct integer_sequence {
  using value_type = T;
  static constexpr std::size_t size() noexcept { return sizeof...(Ints); }
};

template <typename T, std::size_t N>
struct make_integer_sequence_helper;

template <typename T, T N>
using make_integer_sequence = typename make_integer_sequence_helper<T, N>::type;

template <typename T>
struct make_integer_sequence_helper<T, 0> {
  using type = integer_sequence<T>;
};

template <typename T>
struct make_integer_sequence_helper<T, 1> {
  using type = integer_sequence<T, 0>;
};

template <typename T>
struct make_integer_sequence_helper<T, 2> {
  using type = integer_sequence<T, 0, 1>;
};

template <typename T>
struct make_integer_sequence_helper<T, 3> {
  using type = integer_sequence<T, 0, 1, 2>;
};

template <typename T>
struct make_integer_sequence_helper<T, 4> {
  using type = integer_sequence<T, 0, 1, 2, 3>;
};

template <typename T>
struct make_integer_sequence_helper<T, 5> {
  using type = integer_sequence<T, 0, 1, 2, 3, 4>;
};

template <typename T>
struct make_integer_sequence_helper<T, 6> {
  using type = integer_sequence<T, 0, 1, 2, 3, 4, 5>;
};

template <typename T>
struct make_integer_sequence_helper<T, 7> {
  using type = integer_sequence<T, 0, 1, 2, 3, 4, 5, 6>;
};

template <typename T>
struct make_integer_sequence_helper<T, 8> {
  using type = integer_sequence<T, 0, 1, 2, 3, 4, 5, 6, 7>;
};

template <typename X, typename Y>
struct make_integer_sequence_concat;

template <typename T, T... x, T... y>
struct make_integer_sequence_concat<integer_sequence<T, x...>,
                                    integer_sequence<T, y...>> {
  using type = integer_sequence<T, x..., (sizeof...(x) + y)...>;
};

template <typename T, std::size_t N>
struct make_integer_sequence_helper {
  using type = typename make_integer_sequence_concat<
      typename make_integer_sequence_helper<T, N / 2>::type,
      typename make_integer_sequence_helper<T, N - N / 2>::type>::type;
};

//----------------------------------------

template <std::size_t... Indices>
using index_sequence = integer_sequence<std::size_t, Indices...>;

template <std::size_t N>
using make_index_sequence = make_integer_sequence<std::size_t, N>;

//----------------------------------------

template <unsigned I, typename IntegerSequence>
struct integer_sequence_at;

template <unsigned I, typename T, T h0, T... tail>
struct integer_sequence_at<I, integer_sequence<T, h0, tail...>>
    : public integer_sequence_at<I - 1u, integer_sequence<T, tail...>> {
  static_assert(8 <= I, "Reasoning Error");
  static_assert(I < integer_sequence<T, h0, tail...>::size(),
                "Error: Index out of bounds");
};

template <typename T, T h0, T... tail>
struct integer_sequence_at<0u, integer_sequence<T, h0, tail...>> {
  using type               = T;
  static constexpr T value = h0;
};

template <typename T, T h0, T h1, T... tail>
struct integer_sequence_at<1u, integer_sequence<T, h0, h1, tail...>> {
  using type               = T;
  static constexpr T value = h1;
};

template <typename T, T h0, T h1, T h2, T... tail>
struct integer_sequence_at<2u, integer_sequence<T, h0, h1, h2, tail...>> {
  using type               = T;
  static constexpr T value = h2;
};

template <typename T, T h0, T h1, T h2, T h3, T... tail>
struct integer_sequence_at<3u, integer_sequence<T, h0, h1, h2, h3, tail...>> {
  using type               = T;
  static constexpr T value = h3;
};

template <typename T, T h0, T h1, T h2, T h3, T h4, T... tail>
struct integer_sequence_at<4u,
                           integer_sequence<T, h0, h1, h2, h3, h4, tail...>> {
  using type               = T;
  static constexpr T value = h4;
};

template <typename T, T h0, T h1, T h2, T h3, T h4, T h5, T... tail>
struct integer_sequence_at<
    5u, integer_sequence<T, h0, h1, h2, h3, h4, h5, tail...>> {
  using type               = T;
  static constexpr T value = h5;
};

template <typename T, T h0, T h1, T h2, T h3, T h4, T h5, T h6, T... tail>
struct integer_sequence_at<
    6u, integer_sequence<T, h0, h1, h2, h3, h4, h5, h6, tail...>> {
  using type               = T;
  static constexpr T value = h6;
};

template <typename T, T h0, T h1, T h2, T h3, T h4, T h5, T h6, T h7, T... tail>
struct integer_sequence_at<
    7u, integer_sequence<T, h0, h1, h2, h3, h4, h5, h6, h7, tail...>> {
  using type               = T;
  static constexpr T value = h7;
};

//----------------------------------------

template <typename T>
constexpr T at(const unsigned, integer_sequence<T>) noexcept {
  return ~static_cast<T>(0);
}

template <typename T, T h0, T... tail>
constexpr T at(const unsigned i, integer_sequence<T, h0>) noexcept {
  return i == 0u ? h0 : ~static_cast<T>(0);
}

template <typename T, T h0, T h1>
constexpr T at(const unsigned i, integer_sequence<T, h0, h1>) noexcept {
  return i == 0u ? h0 : i == 1u ? h1 : ~static_cast<T>(0);
}

template <typename T, T h0, T h1, T h2>
constexpr T at(const unsigned i, integer_sequence<T, h0, h1, h2>) noexcept {
  return i == 0u ? h0 : i == 1u ? h1 : i == 2u ? h2 : ~static_cast<T>(0);
}

template <typename T, T h0, T h1, T h2, T h3>
constexpr T at(const unsigned i, integer_sequence<T, h0, h1, h2, h3>) noexcept {
  return i == 0u
             ? h0
             : i == 1u ? h1 : i == 2u ? h2 : i == 3u ? h3 : ~static_cast<T>(0);
}

template <typename T, T h0, T h1, T h2, T h3, T h4>
constexpr T at(const unsigned i,
               integer_sequence<T, h0, h1, h2, h3, h4>) noexcept {
  return i == 0u
             ? h0
             : i == 1u
                   ? h1
                   : i == 2u ? h2
                             : i == 3u ? h3 : i == 4u ? h4 : ~static_cast<T>(0);
}

template <typename T, T h0, T h1, T h2, T h3, T h4, T h5>
constexpr T at(const unsigned i,
               integer_sequence<T, h0, h1, h2, h3, h4, h5>) noexcept {
  return i == 0u
             ? h0
             : i == 1u
                   ? h1
                   : i == 2u ? h2
                             : i == 3u ? h3
                                       : i == 4u ? h4
                                                 : i == 5u ? h5
                                                           : ~static_cast<T>(0);
}

template <typename T, T h0, T h1, T h2, T h3, T h4, T h5, T h6>
constexpr T at(const unsigned i,
               integer_sequence<T, h0, h1, h2, h3, h4, h5, h6>) noexcept {
  return i == 0u
             ? h0
             : i == 1u
                   ? h1
                   : i == 2u
                         ? h2
                         : i == 3u
                               ? h3
                               : i == 4u
                                     ? h4
                                     : i == 5u
                                           ? h5
                                           : i == 6u ? h6 : ~static_cast<T>(0);
}

template <typename T, T h0, T h1, T h2, T h3, T h4, T h5, T h6, T h7, T... tail>
constexpr T at(
    const unsigned i,
    integer_sequence<T, h0, h1, h2, h3, h4, h5, h6, h7, tail...>) noexcept {
  return i == 0u
             ? h0
             : i == 1u
                   ? h1
                   : i == 2u
                         ? h2
                         : i == 3u
                               ? h3
                               : i == 4u
                                     ? h4
                                     : i == 5u
                                           ? h5
                                           : i == 6u
                                                 ? h6
                                                 : i == 7u
                                                       ? h7
                                                       : at(i - 8u,
                                                            integer_sequence<
                                                                T, tail...>{});
}

//----------------------------------------

template <typename IntegerSequence,
          typename ResultSequence =
              integer_sequence<typename IntegerSequence::value_type>>
struct reverse_integer_sequence_helper;

template <typename T, T h0, T... tail, T... results>
struct reverse_integer_sequence_helper<integer_sequence<T, h0, tail...>,
                                       integer_sequence<T, results...>>
    : public reverse_integer_sequence_helper<
          integer_sequence<T, tail...>, integer_sequence<T, h0, results...>> {};

template <typename T, T... results>
struct reverse_integer_sequence_helper<integer_sequence<T>,
                                       integer_sequence<T, results...>> {
  using type = integer_sequence<T, results...>;
};

template <typename IntegerSequence>
using reverse_integer_sequence =
    typename reverse_integer_sequence_helper<IntegerSequence>::type;

//----------------------------------------

template <typename IntegerSequence, typename Result,
          typename ResultSequence =
              integer_sequence<typename IntegerSequence::value_type>>
struct exclusive_scan_integer_sequence_helper;

template <typename T, T h0, T... tail, typename Result, T... results>
struct exclusive_scan_integer_sequence_helper<
    integer_sequence<T, h0, tail...>, Result, integer_sequence<T, results...>>
    : public exclusive_scan_integer_sequence_helper<
          integer_sequence<T, tail...>,
          std::integral_constant<T, Result::value + h0>,
          integer_sequence<T, 0, (results + h0)...>> {};

template <typename T, typename Result, T... results>
struct exclusive_scan_integer_sequence_helper<integer_sequence<T>, Result,
                                              integer_sequence<T, results...>> {
  using type               = integer_sequence<T, results...>;
  static constexpr T value = Result::value;
};

template <typename IntegerSequence>
struct exclusive_scan_integer_sequence {
  using value_type = typename IntegerSequence::value_type;
  using helper     = exclusive_scan_integer_sequence_helper<
      reverse_integer_sequence<IntegerSequence>,
      std::integral_constant<value_type, 0>>;
  using type                        = typename helper::type;
  static constexpr value_type value = helper::value;
};

//----------------------------------------

template <typename IntegerSequence, typename Result,
          typename ResultSequence =
              integer_sequence<typename IntegerSequence::value_type>>
struct inclusive_scan_integer_sequence_helper;

template <typename T, T h0, T... tail, typename Result, T... results>
struct inclusive_scan_integer_sequence_helper<
    integer_sequence<T, h0, tail...>, Result, integer_sequence<T, results...>>
    : public inclusive_scan_integer_sequence_helper<
          integer_sequence<T, tail...>,
          std::integral_constant<T, Result::value + h0>,
          integer_sequence<T, h0, (results + h0)...>> {};

template <typename T, typename Result, T... results>
struct inclusive_scan_integer_sequence_helper<integer_sequence<T>, Result,
                                              integer_sequence<T, results...>> {
  using type               = integer_sequence<T, results...>;
  static constexpr T value = Result::value;
};

template <typename IntegerSequence>
struct inclusive_scan_integer_sequence {
  using value_type = typename IntegerSequence::value_type;
  using helper     = inclusive_scan_integer_sequence_helper<
      reverse_integer_sequence<IntegerSequence>,
      std::integral_constant<value_type, 0>>;
  using type                        = typename helper::type;
  static constexpr value_type value = helper::value;
};

template <typename T>
struct identity {
  using type = T;
};

template <typename T>
using identity_t = typename identity<T>::type;

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

}  // namespace Impl
}  // namespace Kokkos

#endif  // KOKKOS_CORE_IMPL_UTILITIES_HPP
