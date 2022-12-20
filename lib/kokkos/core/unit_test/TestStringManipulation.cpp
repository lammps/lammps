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

#include <impl/Kokkos_StringManipulation.hpp>
#include <climits>

namespace {

#define STATIC_ASSERT(cond) static_assert(cond, "")

KOKKOS_FUNCTION constexpr bool test_strlen() {
  using Kokkos::Impl::strlen;
  constexpr char str[] = "How many characters does this string contain?";
  STATIC_ASSERT(strlen(str) == 45);  // without null character
  STATIC_ASSERT(sizeof str == 46);   // with null character
  STATIC_ASSERT(strlen("") == 0);
  return true;
}
STATIC_ASSERT(test_strlen());

KOKKOS_FUNCTION constexpr bool test_strcmp() {
  using Kokkos::Impl::strcmp;
  constexpr char cat1[] = "Heathcliff";
  constexpr char cat2[] = "Snagglepuss";
  constexpr char cat3[] = "Hobbes";
  constexpr char cat4[] = "Garfield";
  STATIC_ASSERT(strcmp(cat1, cat1) == 0);
#if (!defined(KOKKOS_COMPILER_NVCC) ||                                 \
     ((__CUDACC_VER_MAJOR__ >= 11) && (__CUDACC_VER_MINOR__ >= 3))) && \
    (!defined(__INTEL_COMPILER_BUILD_DATE) ||                          \
     (__INTEL_COMPILER_BUILD_DATE >= 20210228))
  STATIC_ASSERT(strcmp(cat1, cat2) < 0);
  STATIC_ASSERT(strcmp(cat1, cat3) < 0);
#endif
  STATIC_ASSERT(strcmp(cat1, cat4) > 0);
  STATIC_ASSERT(strcmp(cat2, cat2) == 0);
  STATIC_ASSERT(strcmp(cat2, cat3) > 0);
  STATIC_ASSERT(strcmp(cat2, cat4) > 0);
  STATIC_ASSERT(strcmp(cat3, cat3) == 0);
  STATIC_ASSERT(strcmp(cat3, cat4) > 0);
  STATIC_ASSERT(strcmp(cat4, cat4) == 0);
  return true;
}
STATIC_ASSERT(test_strcmp());

KOKKOS_FUNCTION constexpr bool test_strncmp() {
  using Kokkos::Impl::strncmp;
  constexpr char greet1[] = "Hello, world!";
  constexpr char greet2[] = "Hello, everybody!";
  constexpr char greet3[] = "Hello, somebody!";
  STATIC_ASSERT(strncmp(greet1, greet2, 13) > 0);
  STATIC_ASSERT(strncmp(greet2, greet1, 13) < 0);
  STATIC_ASSERT(strncmp(greet2, greet1, 7) == 0);
#if defined(KOKKOS_COMPILER_GNU) && (KOKKOS_COMPILER_GNU < 610)
  (void)greet3;
#elif defined(KOKKOS_COMPILER_GNU) && (KOKKOS_COMPILER_GNU < 710)
  STATIC_ASSERT(strncmp(&greet2[12], &greet3[11], 5) == 0);
#else
  STATIC_ASSERT(strncmp(greet2 + 12, greet3 + 11, 5) == 0);
#endif
  STATIC_ASSERT(strncmp(greet1, greet2, 0) == 0);
  return true;
}
STATIC_ASSERT(test_strncmp());

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
  STATIC_ASSERT(strcpy_helper("abcdef", "hi", "hi\0\0\0f"));
  return true;
}
STATIC_ASSERT(test_strcpy());

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
  STATIC_ASSERT(strncpy_helper("abcdef", "hi", 5, "hi\0\0\0f"));
  STATIC_ASSERT(strncpy_helper("abcdef", "hi", 0, "abcdef"));
  return true;
}
STATIC_ASSERT(test_strncpy());

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
  STATIC_ASSERT(strcat_helper("Hello ", "World!", "Hello World!"));
  STATIC_ASSERT(strcat_helper("Hello World!", " Goodbye World!",
                              "Hello World! Goodbye World!"));
  return true;
}
STATIC_ASSERT(test_strcat());

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
  STATIC_ASSERT(
      strncat_helper("Hello World!", " Goodbye World!", 3, "Hello World! Go"));
  STATIC_ASSERT(
      strncat_helper("Hello World!", " Goodbye World!", 0, "Hello World!"));
  return true;
}
STATIC_ASSERT(test_strncat());

#if !defined(KOKKOS_COMPILER_GNU) && (KOKKOS_COMPILER_GNU >= 540)
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
  STATIC_ASSERT(to_chars_helper(0, "0"));
  STATIC_ASSERT(to_chars_helper(123, "123"));
  STATIC_ASSERT(to_chars_helper(-456, "-456"));
  STATIC_ASSERT(to_chars_helper(INT_MAX, "2147483647"));
  STATIC_ASSERT(to_chars_helper(INT_MIN, "-2147483648"));

  STATIC_ASSERT(to_chars_helper(0u, "0"));
  STATIC_ASSERT(to_chars_helper(78u, "78"));
  STATIC_ASSERT(to_chars_helper(UINT_MAX, "4294967295"));

  STATIC_ASSERT(to_chars_helper(0ll, "0"));
  STATIC_ASSERT(to_chars_helper(LLONG_MAX, "9223372036854775807"));
  STATIC_ASSERT(to_chars_helper(LLONG_MIN, "-9223372036854775808"));

  STATIC_ASSERT(to_chars_helper(0ull, "0"));
  STATIC_ASSERT(to_chars_helper(ULLONG_MAX, "18446744073709551615"));

  return true;
}
STATIC_ASSERT(test_to_chars());
#endif

}  // namespace
