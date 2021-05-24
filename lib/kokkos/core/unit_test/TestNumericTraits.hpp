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

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <type_traits>
#include "Kokkos_NumericTraits.hpp"
#include "Kokkos_ExecPolicy.hpp"

struct extrema {
#define DEFINE_EXTREMA(T, m, M)                 \
  KOKKOS_FUNCTION static T min(T) { return m; } \
  KOKKOS_FUNCTION static T max(T) { return M; }

  DEFINE_EXTREMA(char, CHAR_MIN, CHAR_MAX);
  DEFINE_EXTREMA(signed char, SCHAR_MIN, SCHAR_MAX);
  DEFINE_EXTREMA(unsigned char, 0, UCHAR_MAX);
  DEFINE_EXTREMA(short, SHRT_MIN, SHRT_MAX);
  DEFINE_EXTREMA(unsigned short, 0, USHRT_MAX);
  DEFINE_EXTREMA(int, INT_MIN, INT_MAX);
  DEFINE_EXTREMA(unsigned, 0U, UINT_MAX);
  DEFINE_EXTREMA(long, LONG_MIN, LONG_MAX);
  DEFINE_EXTREMA(unsigned long, 0UL, ULONG_MAX);
  DEFINE_EXTREMA(long long, LLONG_MIN, LLONG_MAX);
  DEFINE_EXTREMA(unsigned long long, 0ULL, ULLONG_MAX);

  DEFINE_EXTREMA(float, -FLT_MAX, FLT_MAX);
  DEFINE_EXTREMA(double, -DBL_MAX, DBL_MAX);
  DEFINE_EXTREMA(long double, -LDBL_MAX, LDBL_MAX);

#undef DEFINE_EXTREMA
};

// clang-format off
struct Infinity { template <class T> using trait = Kokkos::Experimental::infinity<T>; };
struct Epsilon { template <class T> using trait = Kokkos::Experimental::epsilon<T>; };
struct FiniteMin { template <class T> using trait = Kokkos::Experimental::finite_min<T>; };
struct FiniteMax { template <class T> using trait = Kokkos::Experimental::finite_max<T>; };
struct RoundError { template <class T> using trait = Kokkos::Experimental::round_error<T>; };
struct NormMin { template <class T> using trait = Kokkos::Experimental::norm_min<T>; };
struct Digits { template <class T> using trait = Kokkos::Experimental::digits<T>; };
struct Digits10 { template <class T> using trait = Kokkos::Experimental::digits10<T>; };
struct MaxDigits10 { template <class T> using trait = Kokkos::Experimental::max_digits10<T>; };
struct Radix { template <class T> using trait = Kokkos::Experimental::radix<T>; };
struct MinExponent { template <class T> using trait = Kokkos::Experimental::min_exponent<T>; };
struct MaxExponent { template <class T> using trait = Kokkos::Experimental::max_exponent<T>; };
struct MinExponent10 { template <class T> using trait = Kokkos::Experimental::min_exponent10<T>; };
struct MaxExponent10 { template <class T> using trait = Kokkos::Experimental::max_exponent10<T>; };
// clang-format on

template <class T>
KOKKOS_FUNCTION T* take_address_of(T& arg) {
  return &arg;
}

template <class T>
KOKKOS_FUNCTION void take_by_value(T) {}

template <class Space, class T, class Tag>
struct TestNumericTraits {
  template <class U>
  using trait = typename Tag::template trait<U>;

  Kokkos::View<T, Space> compare;
  TestNumericTraits() {
    compare = Kokkos::View<T, Space>("C");
    run();
  }

  void run() const {
    int errors = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Space, Tag>(0, 1), *this,
                            errors);
    ASSERT_EQ(errors, 0);
    (void)take_address_of(trait<T>::value);  // use on host
  }

  KOKKOS_FUNCTION void operator()(Infinity, int, int& e) const {
    using Kokkos::Experimental::infinity;
    auto const inf  = infinity<T>::value;
    auto const zero = T(0);
    e += (int)!(inf + inf == inf);
    e += (int)!(inf != zero);
    use_on_device();
  }

  KOKKOS_FUNCTION void operator()(Epsilon, int, int& e) const {
    using Kokkos::Experimental::epsilon;
    auto const eps = epsilon<T>::value;
    auto const one = T(1);
    // Avoid higher precision intermediate representation
    compare() = one + eps;
    e += (int)!(compare() != one);
    compare() = one + eps / 2;
    e += (int)!(compare() == one);
    use_on_device();
  }

  KOKKOS_FUNCTION void operator()(FiniteMin, int, int& e) const {
    using Kokkos::Experimental::finite_max;
    using Kokkos::Experimental::finite_min;
    auto const min = finite_min<T>::value;
    auto const max = finite_max<T>::value;
    e += (int)!(min == extrema::min(T{}));
    e += (int)!(max == extrema::max(T{}));
    use_on_device();
  }

  // clang-format off
  KOKKOS_FUNCTION void operator()(FiniteMax, int, int&) const { use_on_device(); }
  KOKKOS_FUNCTION void operator()(RoundError, int, int&) const { use_on_device(); }
  KOKKOS_FUNCTION void operator()(NormMin, int, int&) const { use_on_device(); }
  KOKKOS_FUNCTION void operator()(Digits, int, int&) const { use_on_device(); }
  KOKKOS_FUNCTION void operator()(Digits10, int, int&) const { use_on_device(); }
  KOKKOS_FUNCTION void operator()(MaxDigits10, int, int&) const { use_on_device(); }
  KOKKOS_FUNCTION void operator()(Radix, int, int&) const { use_on_device(); }
  KOKKOS_FUNCTION void operator()(MinExponent, int, int&) const { use_on_device(); }
  KOKKOS_FUNCTION void operator()(MaxExponent, int, int&) const { use_on_device(); }
  KOKKOS_FUNCTION void operator()(MinExponent10, int, int&) const { use_on_device(); }
  KOKKOS_FUNCTION void operator()(MaxExponent10, int, int&) const { use_on_device(); }
  // clang-format on

  KOKKOS_FUNCTION void use_on_device() const {
#if defined(KOKKOS_COMPILER_NVCC) || defined(KOKKOS_ENABLE_OPENMPTARGET)
    take_by_value(trait<T>::value);
#else
    (void)take_address_of(trait<T>::value);
#endif
  }
};

#if defined(KOKKOS_COMPILER_NVCC) || defined(KOKKOS_ENABLE_SYCL) || \
    defined(KOKKOS_ENABLE_OPENMPTARGET)
template <class Tag>
struct TestNumericTraits<
#if defined(KOKKOS_ENABLE_CUDA)
    Kokkos::Cuda,
#elif defined(KOKKOS_ENABLE_SYCL)
    Kokkos::Experimental::SYCL,
#else
    Kokkos::Experimental::OpenMPTarget,
#endif
    long double, Tag> {
  template <class T>
  using trait = typename Tag::template trait<T>;
  TestNumericTraits() {
    (void)take_address_of(trait<long double>::value);
    // Do nothing on the device.
    // According to the doc
    // https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#constexpr-variables
    // the traits member constant value cannot be directly used in device code.
  }
};
#endif

TEST(TEST_CATEGORY, numeric_traits_infinity) {
  TestNumericTraits<TEST_EXECSPACE, float, Infinity>();
  TestNumericTraits<TEST_EXECSPACE, double, Infinity>();
  TestNumericTraits<TEST_EXECSPACE, long double, Infinity>();
}

TEST(TEST_CATEGORY, numeric_traits_epsilon) {
  TestNumericTraits<TEST_EXECSPACE, float, Epsilon>();
  TestNumericTraits<TEST_EXECSPACE, double, Epsilon>();
#ifndef KOKKOS_COMPILER_IBM  // fails with XL 16.1.1
  TestNumericTraits<TEST_EXECSPACE, long double, Epsilon>();
#endif
}

TEST(TEST_CATEGORY, numeric_traits_round_error) {
  TestNumericTraits<TEST_EXECSPACE, float, RoundError>();
  TestNumericTraits<TEST_EXECSPACE, double, RoundError>();
  TestNumericTraits<TEST_EXECSPACE, long double, RoundError>();
}

TEST(TEST_CATEGORY, numeric_traits_norm_min) {
  TestNumericTraits<TEST_EXECSPACE, float, NormMin>();
  TestNumericTraits<TEST_EXECSPACE, double, NormMin>();
  TestNumericTraits<TEST_EXECSPACE, long double, NormMin>();
}

TEST(TEST_CATEGORY, numeric_traits_finite_min_max) {
  TestNumericTraits<TEST_EXECSPACE, char, FiniteMin>();
  TestNumericTraits<TEST_EXECSPACE, char, FiniteMax>();
  TestNumericTraits<TEST_EXECSPACE, signed char, FiniteMin>();
  TestNumericTraits<TEST_EXECSPACE, signed char, FiniteMax>();
  TestNumericTraits<TEST_EXECSPACE, unsigned char, FiniteMin>();
  TestNumericTraits<TEST_EXECSPACE, unsigned char, FiniteMax>();

  TestNumericTraits<TEST_EXECSPACE, short, FiniteMin>();
  TestNumericTraits<TEST_EXECSPACE, short, FiniteMax>();
  TestNumericTraits<TEST_EXECSPACE, unsigned short, FiniteMin>();
  TestNumericTraits<TEST_EXECSPACE, unsigned short, FiniteMax>();

  TestNumericTraits<TEST_EXECSPACE, int, FiniteMin>();
  TestNumericTraits<TEST_EXECSPACE, int, FiniteMax>();
  TestNumericTraits<TEST_EXECSPACE, unsigned int, FiniteMin>();
  TestNumericTraits<TEST_EXECSPACE, unsigned int, FiniteMax>();

  TestNumericTraits<TEST_EXECSPACE, long, FiniteMin>();
  TestNumericTraits<TEST_EXECSPACE, long, FiniteMax>();
  TestNumericTraits<TEST_EXECSPACE, unsigned long, FiniteMin>();
  TestNumericTraits<TEST_EXECSPACE, unsigned long, FiniteMax>();

  TestNumericTraits<TEST_EXECSPACE, long long, FiniteMin>();
  TestNumericTraits<TEST_EXECSPACE, long long, FiniteMax>();
  TestNumericTraits<TEST_EXECSPACE, unsigned long long, FiniteMin>();
  TestNumericTraits<TEST_EXECSPACE, unsigned long long, FiniteMax>();

  TestNumericTraits<TEST_EXECSPACE, float, FiniteMin>();
  TestNumericTraits<TEST_EXECSPACE, float, FiniteMax>();
  TestNumericTraits<TEST_EXECSPACE, double, FiniteMin>();
  TestNumericTraits<TEST_EXECSPACE, double, FiniteMax>();
  TestNumericTraits<TEST_EXECSPACE, long double, FiniteMin>();
  TestNumericTraits<TEST_EXECSPACE, long double, FiniteMax>();
}

TEST(TEST_CATEGORY, numeric_traits_digits) {
  TestNumericTraits<TEST_EXECSPACE, bool, Digits>();
  TestNumericTraits<TEST_EXECSPACE, char, Digits>();
  TestNumericTraits<TEST_EXECSPACE, signed char, Digits>();
  TestNumericTraits<TEST_EXECSPACE, unsigned char, Digits>();
  TestNumericTraits<TEST_EXECSPACE, short, Digits>();
  TestNumericTraits<TEST_EXECSPACE, unsigned short, Digits>();
  TestNumericTraits<TEST_EXECSPACE, int, Digits>();
  TestNumericTraits<TEST_EXECSPACE, unsigned int, Digits>();
  TestNumericTraits<TEST_EXECSPACE, long int, Digits>();
  TestNumericTraits<TEST_EXECSPACE, unsigned long int, Digits>();
  TestNumericTraits<TEST_EXECSPACE, long long int, Digits>();
  TestNumericTraits<TEST_EXECSPACE, unsigned long long int, Digits>();
  TestNumericTraits<TEST_EXECSPACE, float, Digits>();
  TestNumericTraits<TEST_EXECSPACE, double, Digits>();
  TestNumericTraits<TEST_EXECSPACE, long double, Digits>();
}

TEST(TEST_CATEGORY, numeric_traits_digits10) {
  TestNumericTraits<TEST_EXECSPACE, bool, Digits10>();
  TestNumericTraits<TEST_EXECSPACE, char, Digits10>();
  TestNumericTraits<TEST_EXECSPACE, signed char, Digits10>();
  TestNumericTraits<TEST_EXECSPACE, unsigned char, Digits10>();
  TestNumericTraits<TEST_EXECSPACE, short, Digits10>();
  TestNumericTraits<TEST_EXECSPACE, unsigned short, Digits10>();
  TestNumericTraits<TEST_EXECSPACE, int, Digits10>();
  TestNumericTraits<TEST_EXECSPACE, unsigned int, Digits10>();
  TestNumericTraits<TEST_EXECSPACE, long int, Digits10>();
  TestNumericTraits<TEST_EXECSPACE, unsigned long int, Digits10>();
  TestNumericTraits<TEST_EXECSPACE, long long int, Digits10>();
  TestNumericTraits<TEST_EXECSPACE, unsigned long long int, Digits10>();
  TestNumericTraits<TEST_EXECSPACE, float, Digits10>();
  TestNumericTraits<TEST_EXECSPACE, double, Digits10>();
  TestNumericTraits<TEST_EXECSPACE, long double, Digits10>();
}

TEST(TEST_CATEGORY, numeric_traits_max_digits10) {
  TestNumericTraits<TEST_EXECSPACE, float, MaxDigits10>();
  TestNumericTraits<TEST_EXECSPACE, double, MaxDigits10>();
  TestNumericTraits<TEST_EXECSPACE, long double, MaxDigits10>();
}

TEST(TEST_CATEGORY, numeric_traits_radix) {
  TestNumericTraits<TEST_EXECSPACE, bool, Radix>();
  TestNumericTraits<TEST_EXECSPACE, char, Radix>();
  TestNumericTraits<TEST_EXECSPACE, signed char, Radix>();
  TestNumericTraits<TEST_EXECSPACE, unsigned char, Radix>();
  TestNumericTraits<TEST_EXECSPACE, short, Radix>();
  TestNumericTraits<TEST_EXECSPACE, unsigned short, Radix>();
  TestNumericTraits<TEST_EXECSPACE, int, Radix>();
  TestNumericTraits<TEST_EXECSPACE, unsigned int, Radix>();
  TestNumericTraits<TEST_EXECSPACE, long int, Radix>();
  TestNumericTraits<TEST_EXECSPACE, unsigned long int, Radix>();
  TestNumericTraits<TEST_EXECSPACE, long long int, Radix>();
  TestNumericTraits<TEST_EXECSPACE, unsigned long long int, Radix>();
  TestNumericTraits<TEST_EXECSPACE, float, Radix>();
  TestNumericTraits<TEST_EXECSPACE, double, Radix>();
  TestNumericTraits<TEST_EXECSPACE, long double, Radix>();
}

TEST(TEST_CATEGORY, numeric_traits_min_max_exponent) {
  TestNumericTraits<TEST_EXECSPACE, float, MinExponent>();
  TestNumericTraits<TEST_EXECSPACE, float, MaxExponent>();
  TestNumericTraits<TEST_EXECSPACE, double, MinExponent>();
  TestNumericTraits<TEST_EXECSPACE, double, MaxExponent>();
  TestNumericTraits<TEST_EXECSPACE, long double, MinExponent>();
  TestNumericTraits<TEST_EXECSPACE, long double, MaxExponent>();
}

TEST(TEST_CATEGORY, numeric_traits_min_max_exponent10) {
  TestNumericTraits<TEST_EXECSPACE, float, MinExponent10>();
  TestNumericTraits<TEST_EXECSPACE, float, MaxExponent10>();
  TestNumericTraits<TEST_EXECSPACE, double, MinExponent10>();
  TestNumericTraits<TEST_EXECSPACE, double, MaxExponent10>();
  TestNumericTraits<TEST_EXECSPACE, long double, MinExponent10>();
  TestNumericTraits<TEST_EXECSPACE, long double, MaxExponent10>();
}
