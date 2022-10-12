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
#include <limits>
#include "Kokkos_NumericTraits.hpp"

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
struct DenormMin { template <class T> using trait = Kokkos::Experimental::denorm_min<T>; };
struct ReciprocalOverflowThreshold { template <class T> using trait = Kokkos::Experimental::reciprocal_overflow_threshold<T>; };
struct Digits { template <class T> using trait = Kokkos::Experimental::digits<T>; };
struct Digits10 { template <class T> using trait = Kokkos::Experimental::digits10<T>; };
struct MaxDigits10 { template <class T> using trait = Kokkos::Experimental::max_digits10<T>; };
struct Radix { template <class T> using trait = Kokkos::Experimental::radix<T>; };
struct MinExponent { template <class T> using trait = Kokkos::Experimental::min_exponent<T>; };
struct MaxExponent { template <class T> using trait = Kokkos::Experimental::max_exponent<T>; };
struct MinExponent10 { template <class T> using trait = Kokkos::Experimental::min_exponent10<T>; };
struct MaxExponent10 { template <class T> using trait = Kokkos::Experimental::max_exponent10<T>; };
struct QuietNaN { template <class T> using trait = Kokkos::Experimental::quiet_NaN<T>; };
struct SignalingNaN { template <class T> using trait = Kokkos::Experimental::signaling_NaN<T>; };
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

  KOKKOS_FUNCTION void operator()(ReciprocalOverflowThreshold, int,
                                  int& e) const {
    using Kokkos::Experimental::reciprocal_overflow_threshold;
    auto const inv = 1 / reciprocal_overflow_threshold<T>::value;
    if (inv + inv == inv && inv != 0) {
      KOKKOS_IMPL_DO_NOT_USE_PRINTF(
          "inverse of reciprocal overflow threshold is inf\n");
      ++e;
    }
    use_on_device();
  }

  // clang-format off
  KOKKOS_FUNCTION void operator()(FiniteMax, int, int&) const { use_on_device(); }
  KOKKOS_FUNCTION void operator()(RoundError, int, int&) const { use_on_device(); }
  KOKKOS_FUNCTION void operator()(NormMin, int, int&) const { use_on_device(); }
  KOKKOS_FUNCTION void operator()(DenormMin, int, int&) const { use_on_device(); }
  KOKKOS_FUNCTION void operator()(Digits, int, int&) const { use_on_device(); }
  KOKKOS_FUNCTION void operator()(Digits10, int, int&) const { use_on_device(); }
  KOKKOS_FUNCTION void operator()(MaxDigits10, int, int&) const { use_on_device(); }
  KOKKOS_FUNCTION void operator()(Radix, int, int&) const { use_on_device(); }
  KOKKOS_FUNCTION void operator()(MinExponent, int, int&) const { use_on_device(); }
  KOKKOS_FUNCTION void operator()(MaxExponent, int, int&) const { use_on_device(); }
  KOKKOS_FUNCTION void operator()(MinExponent10, int, int&) const { use_on_device(); }
  KOKKOS_FUNCTION void operator()(MaxExponent10, int, int&) const { use_on_device(); }
  // clang-format on
  KOKKOS_FUNCTION void operator()(QuietNaN, int, int& e) const {
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC
    using Kokkos::Experimental::quiet_NaN;
    constexpr auto nan  = quiet_NaN<T>::value;
    constexpr auto zero = T(0);
    e += (int)!(nan != nan);
    e += (int)!(nan != zero);
#else
    (void)e;
#endif
    use_on_device();
  }
  KOKKOS_FUNCTION void operator()(SignalingNaN, int, int& e) const {
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC
    using Kokkos::Experimental::signaling_NaN;
    constexpr auto nan  = signaling_NaN<T>::value;
    constexpr auto zero = T(0);
    e += (int)!(nan != nan);
    e += (int)!(nan != zero);
#else
    (void)e;
#endif
    use_on_device();
  }

  KOKKOS_FUNCTION void use_on_device() const {
#if defined(KOKKOS_COMPILER_NVCC) || defined(KOKKOS_ENABLE_OPENMPTARGET)
    take_by_value(trait<T>::value);
#else
    (void)take_address_of(trait<T>::value);
#endif
  }
};

#if (defined(KOKKOS_COMPILER_NVCC) && defined(KOKKOS_ENABLE_CUDA)) || \
    defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ENABLE_OPENMPTARGET)
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
  // fails with XL 16.1.1 see issue #4100
  // FIXME_NVHPC long double not supported
#if !defined(KOKKOS_COMPILER_IBM) && !defined(KOKKOS_COMPILER_NVHPC)
  TestNumericTraits<TEST_EXECSPACE, long double, Infinity>();
#endif
}

TEST(TEST_CATEGORY, numeric_traits_epsilon) {
  TestNumericTraits<TEST_EXECSPACE, float, Epsilon>();
  TestNumericTraits<TEST_EXECSPACE, double, Epsilon>();
  // fails with XL 16.1.1 see issue #4100
  // FIXME_NVHPC long double not supported
#if !defined(KOKKOS_COMPILER_IBM) && !defined(KOKKOS_COMPILER_NVHPC)
  TestNumericTraits<TEST_EXECSPACE, long double, Epsilon>();
#endif
}

TEST(TEST_CATEGORY, numeric_traits_round_error) {
  TestNumericTraits<TEST_EXECSPACE, float, RoundError>();
  TestNumericTraits<TEST_EXECSPACE, double, RoundError>();
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC:
  // nvc++-Fatal-/home/projects/x86-64/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/bin/tools/cpp2
  // TERMINATED by signal 11
  TestNumericTraits<TEST_EXECSPACE, long double, RoundError>();
#endif
}

TEST(TEST_CATEGORY, numeric_traits_norm_min) {
  TestNumericTraits<TEST_EXECSPACE, float, NormMin>();
  TestNumericTraits<TEST_EXECSPACE, double, NormMin>();
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC:
  // nvc++-Fatal-/home/projects/x86-64/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/bin/tools/cpp2
  // TERMINATED by signal 11
  TestNumericTraits<TEST_EXECSPACE, long double, NormMin>();
#endif
}

TEST(TEST_CATEGORY, numeric_traits_denorm_min) {
  TestNumericTraits<TEST_EXECSPACE, float, DenormMin>();
  TestNumericTraits<TEST_EXECSPACE, double, DenormMin>();
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC:
  // nvc++-Fatal-/home/projects/x86-64/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/bin/tools/cpp2
  // TERMINATED by signal 11
  TestNumericTraits<TEST_EXECSPACE, long double, DenormMin>();
#endif
}

TEST(TEST_CATEGORY, numeric_traits_reciprocal_overflow_threshold) {
  TestNumericTraits<TEST_EXECSPACE, float, ReciprocalOverflowThreshold>();
  TestNumericTraits<TEST_EXECSPACE, double, ReciprocalOverflowThreshold>();
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC:
  // nvc++-Fatal-/home/projects/x86-64/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/bin/tools/cpp2
  // TERMINATED by signal 11
  TestNumericTraits<TEST_EXECSPACE, long double, ReciprocalOverflowThreshold>();
#endif
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
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC:
  // nvc++-Fatal-/home/projects/x86-64/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/bin/tools/cpp2
  // TERMINATED by signal 11
  TestNumericTraits<TEST_EXECSPACE, long double, FiniteMin>();
  TestNumericTraits<TEST_EXECSPACE, long double, FiniteMax>();
#endif
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
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC:
  // nvc++-Fatal-/home/projects/x86-64/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/bin/tools/cpp2
  // TERMINATED by signal 11
  TestNumericTraits<TEST_EXECSPACE, long double, Digits>();
#endif
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
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC:
  // nvc++-Fatal-/home/projects/x86-64/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/bin/tools/cpp2
  // TERMINATED by signal 11
  TestNumericTraits<TEST_EXECSPACE, long double, Digits10>();
#endif
}

TEST(TEST_CATEGORY, numeric_traits_max_digits10) {
  TestNumericTraits<TEST_EXECSPACE, float, MaxDigits10>();
  TestNumericTraits<TEST_EXECSPACE, double, MaxDigits10>();
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC:
  // nvc++-Fatal-/home/projects/x86-64/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/bin/tools/cpp2
  // TERMINATED by signal 11
  TestNumericTraits<TEST_EXECSPACE, long double, MaxDigits10>();
#endif
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
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC:
  // nvc++-Fatal-/home/projects/x86-64/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/bin/tools/cpp2
  // TERMINATED by signal 11
  TestNumericTraits<TEST_EXECSPACE, long double, Radix>();
#endif
}

TEST(TEST_CATEGORY, numeric_traits_min_max_exponent) {
  TestNumericTraits<TEST_EXECSPACE, float, MinExponent>();
  TestNumericTraits<TEST_EXECSPACE, float, MaxExponent>();
  TestNumericTraits<TEST_EXECSPACE, double, MinExponent>();
  TestNumericTraits<TEST_EXECSPACE, double, MaxExponent>();
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC:
  // nvc++-Fatal-/home/projects/x86-64/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/bin/tools/cpp2
  // TERMINATED by signal 11
  TestNumericTraits<TEST_EXECSPACE, long double, MinExponent>();
  TestNumericTraits<TEST_EXECSPACE, long double, MaxExponent>();
#endif
}

TEST(TEST_CATEGORY, numeric_traits_min_max_exponent10) {
  TestNumericTraits<TEST_EXECSPACE, float, MinExponent10>();
  TestNumericTraits<TEST_EXECSPACE, float, MaxExponent10>();
  TestNumericTraits<TEST_EXECSPACE, double, MinExponent10>();
  TestNumericTraits<TEST_EXECSPACE, double, MaxExponent10>();
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC:
  // nvc++-Fatal-/home/projects/x86-64/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/bin/tools/cpp2
  // TERMINATED by signal 11
  TestNumericTraits<TEST_EXECSPACE, long double, MinExponent10>();
  TestNumericTraits<TEST_EXECSPACE, long double, MaxExponent10>();
#endif
}
TEST(TEST_CATEGORY, numeric_traits_quiet_and_signaling_nan) {
  TestNumericTraits<TEST_EXECSPACE, float, QuietNaN>();
  TestNumericTraits<TEST_EXECSPACE, float, SignalingNaN>();
  TestNumericTraits<TEST_EXECSPACE, double, QuietNaN>();
  TestNumericTraits<TEST_EXECSPACE, double, SignalingNaN>();
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC:
  // Unsupported unknown data type 38.
  // Unsupported unknown data type 38.
  // Unsupported unknown data type 38.
  // nvc++-Fatal-/home/projects/x86-64/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/bin/tools/cpp2
  // TERMINATED by signal 11
  TestNumericTraits<TEST_EXECSPACE, long double, QuietNaN>();
  TestNumericTraits<TEST_EXECSPACE, long double, SignalingNaN>();
#endif
}

namespace NumericTraitsSFINAE {

struct HasNoSpecialization {};

#define CHECK_TRAIT_IS_SFINAE_FRIENDLY(TRAIT)                              \
  template <class T>                                                       \
  using TRAIT##_value_t = decltype(Kokkos::Experimental::TRAIT<T>::value); \
  template <class T>                                                       \
  using has_##TRAIT = Kokkos::is_detected<TRAIT##_value_t, T>;             \
  static_assert(!has_##TRAIT<HasNoSpecialization>::value, "");

CHECK_TRAIT_IS_SFINAE_FRIENDLY(infinity)
CHECK_TRAIT_IS_SFINAE_FRIENDLY(finite_min)
CHECK_TRAIT_IS_SFINAE_FRIENDLY(finite_max)
CHECK_TRAIT_IS_SFINAE_FRIENDLY(epsilon)
CHECK_TRAIT_IS_SFINAE_FRIENDLY(round_error)
CHECK_TRAIT_IS_SFINAE_FRIENDLY(norm_min)
CHECK_TRAIT_IS_SFINAE_FRIENDLY(denorm_min)
CHECK_TRAIT_IS_SFINAE_FRIENDLY(quiet_NaN)
CHECK_TRAIT_IS_SFINAE_FRIENDLY(signaling_NaN)

CHECK_TRAIT_IS_SFINAE_FRIENDLY(digits)
CHECK_TRAIT_IS_SFINAE_FRIENDLY(digits10)
CHECK_TRAIT_IS_SFINAE_FRIENDLY(max_digits10)
CHECK_TRAIT_IS_SFINAE_FRIENDLY(radix)
CHECK_TRAIT_IS_SFINAE_FRIENDLY(min_exponent)
CHECK_TRAIT_IS_SFINAE_FRIENDLY(min_exponent10)
CHECK_TRAIT_IS_SFINAE_FRIENDLY(max_exponent)
CHECK_TRAIT_IS_SFINAE_FRIENDLY(max_exponent10)

}  // namespace NumericTraitsSFINAE

// Example detecting presence or absence of values
template <class T>
using infinity_value_t = decltype(Kokkos::Experimental::infinity<T>::value);

template <class T>
using has_infinity = Kokkos::is_detected<infinity_value_t, T>;

template <class T, std::enable_if_t<has_infinity<T>::value>* = nullptr>
constexpr T legacy_std_numeric_limits_infinity() {
  return Kokkos::Experimental::infinity<T>::value;
}

template <class T, std::enable_if_t<!has_infinity<T>::value>* = nullptr>
constexpr T legacy_std_numeric_limits_infinity() {
  return T();
}

TEST(TEST_CATEGORY, numeric_traits_sfinae_friendly) {
  ASSERT_EQ(legacy_std_numeric_limits_infinity<int>(), 0);
}

// Compare to std::numeric_limits
template <int V1, int V2>
struct AssertIntEquality {
  static constexpr bool value = false;
};
template <int V>
struct AssertIntEquality<V, V> {
  static constexpr bool value = true;
};
#define CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(T, TRAIT)           \
  static_assert(AssertIntEquality<Kokkos::Experimental::TRAIT<T>::value, \
                                  std::numeric_limits<T>::TRAIT>::value, \
                "")
#define CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_FUNCTION(T, TRAIT) \
  static_assert(Kokkos::Experimental::TRAIT<T>::value ==       \
                    std::numeric_limits<T>::TRAIT(),           \
                "")

CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_FUNCTION(float, infinity);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_FUNCTION(double, infinity);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_FUNCTION(long double, infinity);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_FUNCTION(float, epsilon);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_FUNCTION(double, epsilon);
#ifndef KOKKOS_COMPILER_IBM  // fails with XL 16.1.1
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_FUNCTION(long double, epsilon);
#endif
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_FUNCTION(float, round_error);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_FUNCTION(double, round_error);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_FUNCTION(long double, round_error);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_FUNCTION(float, denorm_min);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_FUNCTION(double, denorm_min);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_FUNCTION(long double, denorm_min);
// NOTE reciprocal_overflow_threshold purposefully omitted since it does not
// exist in std::numeric_limits
// clang-format off
static_assert(Kokkos::Experimental::norm_min<float      >::value == std::numeric_limits<      float>::min(), "");
static_assert(Kokkos::Experimental::norm_min<double     >::value == std::numeric_limits<     double>::min(), "");
static_assert(Kokkos::Experimental::norm_min<long double>::value == std::numeric_limits<long double>::min(), "");
// integer types
static_assert(Kokkos::Experimental::finite_min<char                  >::value == std::numeric_limits<                  char>::min(), "");
static_assert(Kokkos::Experimental::finite_min<signed char           >::value == std::numeric_limits<           signed char>::min(), "");
static_assert(Kokkos::Experimental::finite_min<unsigned char         >::value == std::numeric_limits<         unsigned char>::min(), "");
static_assert(Kokkos::Experimental::finite_min<short                 >::value == std::numeric_limits<                 short>::min(), "");
static_assert(Kokkos::Experimental::finite_min<unsigned short        >::value == std::numeric_limits<        unsigned short>::min(), "");
static_assert(Kokkos::Experimental::finite_min<int                   >::value == std::numeric_limits<                   int>::min(), "");
static_assert(Kokkos::Experimental::finite_min<unsigned int          >::value == std::numeric_limits<          unsigned int>::min(), "");
static_assert(Kokkos::Experimental::finite_min<long int              >::value == std::numeric_limits<              long int>::min(), "");
static_assert(Kokkos::Experimental::finite_min<unsigned long int     >::value == std::numeric_limits<     unsigned long int>::min(), "");
static_assert(Kokkos::Experimental::finite_min<long long int         >::value == std::numeric_limits<         long long int>::min(), "");
static_assert(Kokkos::Experimental::finite_min<unsigned long long int>::value == std::numeric_limits<unsigned long long int>::min(), "");
static_assert(Kokkos::Experimental::finite_max<char                  >::value == std::numeric_limits<                  char>::max(), "");
static_assert(Kokkos::Experimental::finite_max<signed char           >::value == std::numeric_limits<           signed char>::max(), "");
static_assert(Kokkos::Experimental::finite_max<unsigned char         >::value == std::numeric_limits<         unsigned char>::max(), "");
static_assert(Kokkos::Experimental::finite_max<short                 >::value == std::numeric_limits<                 short>::max(), "");
static_assert(Kokkos::Experimental::finite_max<unsigned short        >::value == std::numeric_limits<        unsigned short>::max(), "");
static_assert(Kokkos::Experimental::finite_max<int                   >::value == std::numeric_limits<                   int>::max(), "");
static_assert(Kokkos::Experimental::finite_max<unsigned int          >::value == std::numeric_limits<          unsigned int>::max(), "");
static_assert(Kokkos::Experimental::finite_max<long int              >::value == std::numeric_limits<              long int>::max(), "");
static_assert(Kokkos::Experimental::finite_max<unsigned long int     >::value == std::numeric_limits<     unsigned long int>::max(), "");
static_assert(Kokkos::Experimental::finite_max<long long int         >::value == std::numeric_limits<         long long int>::max(), "");
static_assert(Kokkos::Experimental::finite_max<unsigned long long int>::value == std::numeric_limits<unsigned long long int>::max(), "");
// floating point types
static_assert(Kokkos::Experimental::finite_min<float      >::value == -std::numeric_limits<      float>::max(), "");
static_assert(Kokkos::Experimental::finite_min<double     >::value == -std::numeric_limits<     double>::max(), "");
static_assert(Kokkos::Experimental::finite_min<long double>::value == -std::numeric_limits<long double>::max(), "");
static_assert(Kokkos::Experimental::finite_max<float      >::value ==  std::numeric_limits<      float>::max(), "");
static_assert(Kokkos::Experimental::finite_max<double     >::value ==  std::numeric_limits<     double>::max(), "");
static_assert(Kokkos::Experimental::finite_max<long double>::value ==  std::numeric_limits<long double>::max(), "");
// clang-format on

CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(bool, digits);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(char, digits);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(signed char, digits);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(unsigned char, digits);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(short, digits);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(unsigned short, digits);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(int, digits);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(unsigned int, digits);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(long int, digits);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(unsigned long int, digits);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(long long int, digits);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(unsigned long long int, digits);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(float, digits);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(double, digits);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(long double, digits);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(bool, digits10);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(char, digits10);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(signed char, digits10);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(unsigned char, digits10);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(short, digits10);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(unsigned short, digits10);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(int, digits10);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(unsigned int, digits10);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(long int, digits10);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(unsigned long int, digits10);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(long long int, digits10);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(unsigned long long int, digits10);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(float, digits10);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(double, digits10);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(long double, digits10);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(float, max_digits10);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(double, max_digits10);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(long double, max_digits10);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(bool, radix);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(char, radix);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(signed char, radix);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(unsigned char, radix);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(short, radix);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(unsigned short, radix);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(int, radix);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(unsigned int, radix);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(long int, radix);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(unsigned long int, radix);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(long long int, radix);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(unsigned long long int, radix);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(float, radix);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(double, radix);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(long double, radix);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(float, min_exponent);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(float, max_exponent);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(double, min_exponent);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(double, max_exponent);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(long double, min_exponent);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(long double, max_exponent);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(float, min_exponent10);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(float, max_exponent10);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(double, min_exponent10);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(double, max_exponent10);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(long double, min_exponent10);
CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT(long double, max_exponent10);

#undef CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_FUNCTION
#undef CHECK_SAME_AS_NUMERIC_LIMITS_MEMBER_CONSTANT

#define CHECK_NAN_SAME_AS_NUMERIC_LIMITS_MEMBER_FUNCTION(T, TRAIT)             \
  static_assert(Kokkos::Experimental::TRAIT<T>::value !=                       \
                    Kokkos::Experimental::TRAIT<T>::value,                     \
                "");                                                           \
  static_assert(                                                               \
      std::numeric_limits<T>::TRAIT() != std::numeric_limits<T>::TRAIT(), ""); \
  static_assert(Kokkos::Experimental::TRAIT<T>::value !=                       \
                    std::numeric_limits<T>::TRAIT(),                           \
                "")

// Workaround compiler issue error: expression must have a constant value
// See kokkos/kokkos#4574
// There is the same bug with CUDA 11.6
// FIXME_NVHPC FIXME_CUDA FIXME_NVCC
#if !defined(KOKKOS_COMPILER_NVHPC) && (CUDA_VERSION < 11060) && \
    !(defined(KOKKOS_COMPILER_NVCC) && !defined(KOKKOS_ENABLE_CUDA))
CHECK_NAN_SAME_AS_NUMERIC_LIMITS_MEMBER_FUNCTION(float, quiet_NaN);
CHECK_NAN_SAME_AS_NUMERIC_LIMITS_MEMBER_FUNCTION(double, quiet_NaN);
CHECK_NAN_SAME_AS_NUMERIC_LIMITS_MEMBER_FUNCTION(long double, quiet_NaN);
CHECK_NAN_SAME_AS_NUMERIC_LIMITS_MEMBER_FUNCTION(float, signaling_NaN);
CHECK_NAN_SAME_AS_NUMERIC_LIMITS_MEMBER_FUNCTION(double, signaling_NaN);
CHECK_NAN_SAME_AS_NUMERIC_LIMITS_MEMBER_FUNCTION(long double, signaling_NaN);
#endif

#undef CHECK_NAN_SAME_AS_NUMERIC_LIMITS_MEMBER_FUNCTION

#define CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES(T, TRAIT)              \
  static_assert(Kokkos::Experimental::TRAIT<T const>::value ==          \
                    Kokkos::Experimental::TRAIT<T>::value,              \
                "");                                                    \
  static_assert(Kokkos::Experimental::TRAIT<T volatile>::value ==       \
                    Kokkos::Experimental::TRAIT<T>::value,              \
                "");                                                    \
  static_assert(Kokkos::Experimental::TRAIT<T const volatile>::value == \
                    Kokkos::Experimental::TRAIT<T>::value,              \
                "")

#define CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_FLOATING_POINT(TRAIT) \
  CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES(float, TRAIT);              \
  CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES(double, TRAIT);             \
  CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES(long double, TRAIT)

#define CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_INTEGRAL(TRAIT)      \
  CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES(bool, TRAIT);              \
  CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES(char, TRAIT);              \
  CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES(signed char, TRAIT);       \
  CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES(unsigned char, TRAIT);     \
  CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES(short, TRAIT);             \
  CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES(unsigned short, TRAIT);    \
  CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES(int, TRAIT);               \
  CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES(unsigned int, TRAIT);      \
  CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES(long int, TRAIT);          \
  CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES(unsigned long int, TRAIT); \
  CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES(long long int, TRAIT);     \
  CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES(unsigned long long int, TRAIT)

CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_FLOATING_POINT(infinity);
CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_FLOATING_POINT(finite_min);
CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_INTEGRAL(finite_min);
CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_FLOATING_POINT(finite_max);
CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_INTEGRAL(finite_max);
CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_FLOATING_POINT(epsilon);
CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_FLOATING_POINT(round_error);
CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_FLOATING_POINT(norm_min);
CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_FLOATING_POINT(
    reciprocal_overflow_threshold);

CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_FLOATING_POINT(digits);
CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_INTEGRAL(digits);
CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_FLOATING_POINT(digits10);
CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_INTEGRAL(digits10);
CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_FLOATING_POINT(max_digits10);
CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_FLOATING_POINT(radix);
CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_INTEGRAL(radix);
CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_FLOATING_POINT(min_exponent);
CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_FLOATING_POINT(min_exponent10);
CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_FLOATING_POINT(max_exponent);
CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_FLOATING_POINT(max_exponent10);

#undef CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_INTEGRAL
#undef CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES_FLOATING_POINT
#undef CHECK_INSTANTIATED_ON_CV_QUALIFIED_TYPES

#define CHECK_NAN_INSTANTIATED_ON_CV_QUALIFIED_TYPES(T, TRAIT)          \
  static_assert(Kokkos::Experimental::TRAIT<T>::value !=                \
                    Kokkos::Experimental::TRAIT<T>::value,              \
                "");                                                    \
  static_assert(Kokkos::Experimental::TRAIT<T const>::value !=          \
                    Kokkos::Experimental::TRAIT<T>::value,              \
                "");                                                    \
  static_assert(Kokkos::Experimental::TRAIT<T volatile>::value !=       \
                    Kokkos::Experimental::TRAIT<T>::value,              \
                "");                                                    \
  static_assert(Kokkos::Experimental::TRAIT<T const volatile>::value != \
                    Kokkos::Experimental::TRAIT<T>::value,              \
                "")

#define CHECK_NAN_INSTANTIATED_ON_CV_QUALIFIED_TYPES_FLOATING_POINT(TRAIT) \
  CHECK_NAN_INSTANTIATED_ON_CV_QUALIFIED_TYPES(float, TRAIT);              \
  CHECK_NAN_INSTANTIATED_ON_CV_QUALIFIED_TYPES(double, TRAIT);             \
  CHECK_NAN_INSTANTIATED_ON_CV_QUALIFIED_TYPES(long double, TRAIT)

CHECK_NAN_INSTANTIATED_ON_CV_QUALIFIED_TYPES_FLOATING_POINT(quiet_NaN);
CHECK_NAN_INSTANTIATED_ON_CV_QUALIFIED_TYPES_FLOATING_POINT(signaling_NaN);

#undef CHECK_NAN_INSTANTIATED_ON_CV_QUALIFIED_TYPES_FLOATING_POINT
#undef CHECK_NAN_INSTANTIATED_ON_CV_QUALIFIED_TYPES
