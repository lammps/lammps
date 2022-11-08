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

#ifndef KOKKOS_STRING_MANIPULATION_HPP
#define KOKKOS_STRING_MANIPULATION_HPP

#include <Kokkos_Macros.hpp>
#include <cstddef>
#include <type_traits>

namespace Kokkos {
namespace Impl {

// This header provides a subset of the functionality from <cstring>.  In
// contrast to the standard library header, functions are usable on the device
// and in constant expressions.  It also includes functionality from <charconv>
// to convert an integer value to a character sequence.

//<editor-fold desc="String examination">
// returns the length of a given string
KOKKOS_INLINE_FUNCTION constexpr std::size_t strlen(const char *str) {
  std::size_t i = 0;
  while (str[i] != '\0') {
    ++i;
  }
  return i;
}

// compares two strings
KOKKOS_INLINE_FUNCTION constexpr int strcmp(const char *lhs, const char *rhs) {
  while (*lhs == *rhs++) {
    if (*lhs++ == '\0') {
      return 0;
    }
  }
  return static_cast<unsigned int>(*lhs) -
         static_cast<unsigned int>(*(rhs - 1));
}

// compares a certain number of characters from two strings
KOKKOS_INLINE_FUNCTION constexpr int strncmp(const char *lhs, const char *rhs,
                                             std::size_t count) {
  for (std::size_t i = 0; i < count; ++i) {
    if (lhs[i] != rhs[i]) {
      return lhs[i] < rhs[i] ? -1 : 1;
    } else if (lhs[i] == '\0') {
      return 0;
    }
  }
  return 0;
}
//</editor-fold>

//<editor-fold desc="String manipulation">
// copies one string to another
KOKKOS_INLINE_FUNCTION constexpr char *strcpy(char *dest, const char *src) {
  char *d = dest;
  for (; (*d = *src) != '\0'; ++d, ++src) {
  }
  return dest;
}

// copies a certain amount of characters from one string to another
KOKKOS_INLINE_FUNCTION constexpr char *strncpy(char *dest, const char *src,
                                               std::size_t count) {
  if (count != 0) {
    char *d = dest;
    do {
      if ((*d++ = *src++) == '\0') {
        while (--count != 0) {
          *d++ = '\0';
        }
        break;
      }
    } while (--count != 0);
  }
  return dest;
}

// concatenates two strings
KOKKOS_INLINE_FUNCTION constexpr char *strcat(char *dest, const char *src) {
  char *d = dest;
  for (; *d != '\0'; ++d) {
  }
  while ((*d++ = *src++) != '\0') {
  }
  return dest;
}

// concatenates a certain amount of characters of two strings
KOKKOS_INLINE_FUNCTION constexpr char *strncat(char *dest, const char *src,
                                               std::size_t count) {
  if (count != 0) {
    char *d = dest;
    for (; *d != '\0'; ++d) {
    }
    do {
      if ((*d = *src++) == '\0') {
        break;
      }
      d++;
    } while (--count != 0);
    *d = '\0';
  }
  return dest;
}
//</editor-fold>

//<editor-fold desc="Character conversions">
template <class Unsigned>
KOKKOS_FUNCTION constexpr unsigned int to_chars_len(Unsigned val) {
  unsigned int const base = 10;
  static_assert(std::is_integral<Unsigned>::value, "implementation bug");
  static_assert(std::is_unsigned<Unsigned>::value, "implementation bug");
  unsigned int n = 1;
  while (val >= base) {
    val /= base;
    ++n;
  }
  return n;
}
template <class Unsigned>
KOKKOS_FUNCTION constexpr void to_chars_impl(char *first, unsigned int len,
                                             Unsigned val) {
  unsigned int const base = 10;
  static_assert(std::is_integral<Unsigned>::value, "implementation bug");
  static_assert(std::is_unsigned<Unsigned>::value, "implementation bug");
  unsigned int pos = len - 1;
  while (val > 0) {
    auto const num = val % base;
    val /= base;
    first[pos] = '0' + num;
    --pos;
  }
}

// define values of portable error conditions that correspond to the POSIX error
// codes
enum class errc {
  value_too_large = 75  // equivalent POSIX error is EOVERFLOW
};
struct to_chars_result {
  char *ptr;
  errc ec;
};

// converts an integer value to a character sequence
template <class Integral>
KOKKOS_FUNCTION constexpr to_chars_result to_chars_i(char *first, char *last,
                                                     Integral value) {
  using Unsigned = std::conditional_t<sizeof(Integral) <= sizeof(unsigned int),
                                      unsigned int, unsigned long long>;
  Unsigned unsigned_val = value;
  if (value == 0) {
    *first = '0';
    return {first + 1, {}};
  } else if
#ifdef KOKKOS_ENABLE_CXX17
      constexpr
#endif
      (std::is_signed<Integral>::value) {
    if (value < 0) {
      *first++     = '-';
      unsigned_val = Unsigned(~value) + Unsigned(1);
    }
  }
  unsigned int const len = to_chars_len(unsigned_val);
  if (last - first < len) {
    return {last, errc::value_too_large};
  }
  to_chars_impl(first, len, unsigned_val);
  return {first + len, {}};
}
//</editor-fold>

}  // namespace Impl
}  // namespace Kokkos

#endif
