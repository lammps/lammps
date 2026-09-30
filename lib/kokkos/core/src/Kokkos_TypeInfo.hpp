// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TYPE_INFO_HPP
#define KOKKOS_TYPE_INFO_HPP

#include <array>
#include <string_view>
#include <utility>

#include <Kokkos_Macros.hpp>

// NVCC+MSVC doesn't work at all - it simply reports "T" inside type_name
#if !(defined(KOKKOS_COMPILER_NVCC) && defined(KOKKOS_COMPILER_MSVC))

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
