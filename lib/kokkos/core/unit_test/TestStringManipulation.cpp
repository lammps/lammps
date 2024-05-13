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

#include <impl/Kokkos_StringManipulation.hpp>
#include <climits>

namespace {

KOKKOS_FUNCTION constexpr bool test_strlen() {
  using Kokkos::Impl::strlen;
  constexpr char str[] = "How many characters does this string contain?";
  static_assert(strlen(str) == 45);  // without null character
  static_assert(sizeof str == 46);   // with null character
  static_assert(strlen("") == 0);
  return true;
}
static_assert(test_strlen());

KOKKOS_FUNCTION constexpr bool test_strcmp() {
  using Kokkos::Impl::strcmp;
  constexpr char cat1[] = "Heathcliff";
  constexpr char cat2[] = "Snagglepuss";
  constexpr char cat3[] = "Hobbes";
  constexpr char cat4[] = "Garfield";
  static_assert(strcmp(cat1, cat1) == 0);
#if (!defined(KOKKOS_COMPILER_NVCC) ||                                 \
     ((__CUDACC_VER_MAJOR__ >= 11) && (__CUDACC_VER_MINOR__ >= 3))) && \
    (!defined(__INTEL_COMPILER_BUILD_DATE) ||                          \
     (__INTEL_COMPILER_BUILD_DATE >= 20210228))
  static_assert(strcmp(cat1, cat2) < 0);
  static_assert(strcmp(cat1, cat3) < 0);
#endif
  static_assert(strcmp(cat1, cat4) > 0);
  static_assert(strcmp(cat2, cat2) == 0);
  static_assert(strcmp(cat2, cat3) > 0);
  static_assert(strcmp(cat2, cat4) > 0);
  static_assert(strcmp(cat3, cat3) == 0);
  static_assert(strcmp(cat3, cat4) > 0);
  static_assert(strcmp(cat4, cat4) == 0);
  return true;
}
static_assert(test_strcmp());

KOKKOS_FUNCTION constexpr bool test_strncmp() {
  using Kokkos::Impl::strncmp;
  constexpr char greet1[] = "Hello, world!";
  constexpr char greet2[] = "Hello, everybody!";
  constexpr char greet3[] = "Hello, somebody!";
  static_assert(strncmp(greet1, greet2, 13) > 0);
  static_assert(strncmp(greet2, greet1, 13) < 0);
  static_assert(strncmp(greet2, greet1, 7) == 0);
  static_assert(strncmp(greet2 + 12, greet3 + 11, 5) == 0);
  static_assert(strncmp(greet1, greet2, 0) == 0);
  return true;
}
static_assert(test_strncmp());

KOKKOS_FUNCTION constexpr bool strcpy_helper(const char* dest, const char* src,
                                             const char* ref) {
  using Kokkos::Impl::strcmp;
  using Kokkos::Impl::strcpy;
  char buffer[50] = {};
  strcpy(buffer, dest);
  strcpy(buffer, src);
  return strcmp(buffer, ref) == 0;
}

KOKKOS_FUNCTION constexpr bool test_strcpy() {
  static_assert(strcpy_helper("abcdef", "hi", "hi\0\0\0f"));
  return true;
}
static_assert(test_strcpy());

KOKKOS_FUNCTION constexpr bool strncpy_helper(const char* dest, const char* src,
                                              std::size_t count,
                                              const char* ref) {
  using Kokkos::Impl::strcmp;
  using Kokkos::Impl::strlen;
  using Kokkos::Impl::strncpy;
  char buffer[50] = {};
  strncpy(buffer, dest, strlen(dest));
  strncpy(buffer, src, count);
  return strcmp(buffer, ref) == 0;
}

KOKKOS_FUNCTION constexpr bool test_strncpy() {
  static_assert(strncpy_helper("abcdef", "hi", 5, "hi\0\0\0f"));
  static_assert(strncpy_helper("abcdef", "hi", 0, "abcdef"));
  return true;
}
static_assert(test_strncpy());

KOKKOS_FUNCTION constexpr bool strcat_helper(const char* dest, const char* src,
                                             const char* ref) {
  using Kokkos::Impl::strcat;
  using Kokkos::Impl::strcmp;
  char buffer[50] = {};
  strcat(buffer, dest);
  strcat(buffer, src);
  return strcmp(buffer, ref) == 0;
}

KOKKOS_FUNCTION constexpr bool test_strcat() {
  static_assert(strcat_helper("Hello ", "World!", "Hello World!"));
  static_assert(strcat_helper("Hello World!", " Goodbye World!",
                              "Hello World! Goodbye World!"));
  return true;
}
static_assert(test_strcat());

KOKKOS_FUNCTION constexpr bool strncat_helper(const char* dest, const char* src,
                                              std::size_t count,
                                              const char* ref) {
  using Kokkos::Impl::strcmp;
  using Kokkos::Impl::strlen;
  using Kokkos::Impl::strncat;
  char buffer[50] = {};
  strncat(buffer, dest, strlen(dest));
  strncat(buffer, src, count);
  return strcmp(buffer, ref) == 0;
}

KOKKOS_FUNCTION constexpr bool test_strncat() {
  static_assert(
      strncat_helper("Hello World!", " Goodbye World!", 3, "Hello World! Go"));
  static_assert(
      strncat_helper("Hello World!", " Goodbye World!", 0, "Hello World!"));
  return true;
}
static_assert(test_strncat());

template <class Integral>
KOKKOS_FUNCTION constexpr bool to_chars_helper(Integral val, char const* ref) {
  using Kokkos::Impl::strcmp;
  using Kokkos::Impl::strlen;
  using Kokkos::Impl::to_chars_i;
  constexpr int BUFFER_SIZE = 21;
  char buffer[BUFFER_SIZE]  = {};
  return (buffer + strlen(ref) ==
          to_chars_i(buffer, buffer + BUFFER_SIZE, val).ptr) &&
         (strcmp(buffer, ref) == 0);
}

KOKKOS_FUNCTION constexpr bool test_to_chars() {
  static_assert(to_chars_helper(0, "0"));
  static_assert(to_chars_helper(123, "123"));
  static_assert(to_chars_helper(-456, "-456"));
  static_assert(to_chars_helper(INT_MAX, "2147483647"));
  static_assert(to_chars_helper(INT_MIN, "-2147483648"));

  static_assert(to_chars_helper(0u, "0"));
  static_assert(to_chars_helper(78u, "78"));
  static_assert(to_chars_helper(UINT_MAX, "4294967295"));

  static_assert(to_chars_helper(0ll, "0"));
  static_assert(to_chars_helper(LLONG_MAX, "9223372036854775807"));
  static_assert(to_chars_helper(LLONG_MIN, "-9223372036854775808"));

  static_assert(to_chars_helper(0ull, "0"));
  static_assert(to_chars_helper(ULLONG_MAX, "18446744073709551615"));

  return true;
}
static_assert(test_to_chars());

}  // namespace
