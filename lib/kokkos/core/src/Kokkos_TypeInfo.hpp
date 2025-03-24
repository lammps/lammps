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

#ifndef KOKKOS_TYPE_INFO_HPP
#define KOKKOS_TYPE_INFO_HPP

#include <array>
#include <string_view>
#include <utility>

#include <Kokkos_Macros.hpp>

// Intel C++ Compiler Classic version 2021.2.0 works but 2021.1.2 doesn't
// Both have __INTEL_COMPILER defined to 2021 so using
// __INTEL_COMPILER_BUILD_DATE to discriminate.
// Experimenting on the compiler explorer gave
//     icc version | __INTEL_COMPILER | __INTEL_COMPILER_BUILD_DATE
//     2021.1.2    | 2021             | 20201208
//     2021.2.0    | 2021             | 20210228
// NVCC versions less than 11.3.0 segfault when that header is included
// NVCC+MSVC doesn't work at all - it simply reports "T" inside type_name
#if (!defined(KOKKOS_COMPILER_INTEL) ||                                   \
     (__INTEL_COMPILER_BUILD_DATE >= 20210228)) &&                        \
    (!defined(KOKKOS_COMPILER_NVCC) || (KOKKOS_COMPILER_NVCC >= 1130)) && \
    (!(defined(KOKKOS_COMPILER_NVCC) && defined(KOKKOS_COMPILER_MSVC)))

#define KOKKOS_ENABLE_IMPL_TYPEINFO

namespace Kokkos::Impl {

template <size_t N>
constexpr std::array<char, N> to_array(std::string_view src) {
  std::array<char, N> dst{};
  for (size_t i = 0; i < N; ++i) {
    dst[i] = src[i];
  }
  return dst;
}

template <class T>
constexpr auto type_name() {
#if defined(__clang__)
  constexpr std::string_view func = __PRETTY_FUNCTION__;
  constexpr std::string_view prefix{"[T = "};
  constexpr std::string_view suffix{"]"};
#elif defined(__GNUC__)
  constexpr std::string_view func = __PRETTY_FUNCTION__;
  constexpr std::string_view prefix{"[with T = "};
  constexpr std::string_view suffix{"]"};
#elif defined(_MSC_VER)
  constexpr std::string_view func = __FUNCSIG__;
  constexpr std::string_view prefix{"type_name<"};
  constexpr std::string_view suffix{">(void)"};
#else
#error bug
#endif
  constexpr auto beg = func.find(prefix) + prefix.size();
  constexpr auto end = func.rfind(suffix);
  static_assert(beg != std::string_view::npos);
  static_assert(end != std::string_view::npos);
  return to_array<end - beg>(func.substr(beg, end));
}

template <class T>
class TypeInfo {
  static constexpr auto value_ = type_name<T>();

 public:
  static constexpr std::string_view name() noexcept {
    return {value_.data(), value_.size()};
  }
};

}  // namespace Kokkos::Impl

#else  // out of luck, using Intel C++ Compiler Classic

namespace Kokkos::Impl {

template <class T>
class TypeInfo {
 public:
  static constexpr std::string_view name() noexcept { return "not supported"; }
};

}  // namespace Kokkos::Impl

#endif

#endif
