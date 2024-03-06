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

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <algorithm>
#include <initializer_list>
#include <type_traits>

#include <cfloat>

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) ||          \
    defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ENABLE_OPENMPTARGET) || \
    defined(KOKKOS_ENABLE_OPENACC)
#else
#define MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
#endif

#if defined KOKKOS_COMPILER_INTEL ||                                  \
    (defined(KOKKOS_COMPILER_NVCC) && KOKKOS_COMPILER_NVCC >= 1130 && \
     !defined(KOKKOS_COMPILER_MSVC))
#define MATHEMATICAL_FUNCTIONS_TEST_UNREACHABLE __builtin_unreachable();
#else
#define MATHEMATICAL_FUNCTIONS_TEST_UNREACHABLE
#endif

namespace KE = Kokkos::Experimental;

// clang-format off
template <class>
struct math_unary_function_return_type;
// Floating-point types
#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
template <> struct math_unary_function_return_type<KE::half_t> { using type = KE::half_t; };
#endif // defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
#if defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT
template <> struct math_unary_function_return_type<KE::bhalf_t> { using type = KE::bhalf_t; };
#endif // defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT
template <> struct math_unary_function_return_type<      float> { using type =       float; };
template <> struct math_unary_function_return_type<     double> { using type =      double; };
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
template <> struct math_unary_function_return_type<long double> { using type = long double; };
#endif
// Integral types
template <> struct math_unary_function_return_type<              bool> { using type = double; };
template <> struct math_unary_function_return_type<             short> { using type = double; };
template <> struct math_unary_function_return_type<    unsigned short> { using type = double; };
template <> struct math_unary_function_return_type<               int> { using type = double; };
template <> struct math_unary_function_return_type<      unsigned int> { using type = double; };
template <> struct math_unary_function_return_type<              long> { using type = double; };
template <> struct math_unary_function_return_type<     unsigned long> { using type = double; };
template <> struct math_unary_function_return_type<         long long> { using type = double; };
template <> struct math_unary_function_return_type<unsigned long long> { using type = double; };
template <class T>
using math_unary_function_return_type_t = typename math_unary_function_return_type<T>::type;
template <class, class>
struct math_binary_function_return_type;
#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
template <> struct math_binary_function_return_type<KE::half_t, KE::half_t> { using type = KE::half_t; };
template <> struct math_binary_function_return_type<short, KE::half_t> { using type = double; };
template <> struct math_binary_function_return_type<unsigned short, KE::half_t> { using type = double; };
template <> struct math_binary_function_return_type<int, KE::half_t> { using type = double; };
template <> struct math_binary_function_return_type<unsigned int, KE::half_t> { using type = double; };
template <> struct math_binary_function_return_type<long, KE::half_t> { using type = double; };
template <> struct math_binary_function_return_type<unsigned long, KE::half_t> { using type = double; };
template <> struct math_binary_function_return_type<long long, KE::half_t> { using type = double; };
template <> struct math_binary_function_return_type<unsigned long long, KE::half_t> { using type = double; };
#endif // defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
#if defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT
template <> struct math_binary_function_return_type<KE::bhalf_t, KE::bhalf_t> { using type = KE::bhalf_t; };
template <> struct math_binary_function_return_type<KE::half_t, KE::bhalf_t> { using type = KE::half_t; };
template <> struct math_binary_function_return_type<short, KE::bhalf_t> { using type = double; };
template <> struct math_binary_function_return_type<unsigned short, KE::bhalf_t> { using type = double; };
template <> struct math_binary_function_return_type<int, KE::bhalf_t> { using type = double; };
template <> struct math_binary_function_return_type<unsigned int, KE::bhalf_t> { using type = double; };
template <> struct math_binary_function_return_type<long, KE::bhalf_t> { using type = double; };
template <> struct math_binary_function_return_type<unsigned long, KE::bhalf_t> { using type = double; };
template <> struct math_binary_function_return_type<long long, KE::bhalf_t> { using type = double; };
template <> struct math_binary_function_return_type<unsigned long long, KE::bhalf_t> { using type = double; };
#endif // defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT
template <> struct math_binary_function_return_type<             float,              float> { using type =       float; };
template <> struct math_binary_function_return_type<             float,             double> { using type =      double; };
template <> struct math_binary_function_return_type<             float,               bool> { using type =      double; };
template <> struct math_binary_function_return_type<             float,              short> { using type =      double; };
template <> struct math_binary_function_return_type<             float,                int> { using type =      double; };
template <> struct math_binary_function_return_type<             float,               long> { using type =      double; };
template <> struct math_binary_function_return_type<             float,          long long> { using type =      double; };
template <> struct math_binary_function_return_type<             float,     unsigned short> { using type =      double; };
template <> struct math_binary_function_return_type<             float,       unsigned int> { using type =      double; };
template <> struct math_binary_function_return_type<             float,      unsigned long> { using type =      double; };
template <> struct math_binary_function_return_type<             float, unsigned long long> { using type =      double; };
template <> struct math_binary_function_return_type<            double,              float> { using type =      double; };
template <> struct math_binary_function_return_type<            double,             double> { using type =      double; };
template <> struct math_binary_function_return_type<            double,               bool> { using type =      double; };
template <> struct math_binary_function_return_type<            double,              short> { using type =      double; };
template <> struct math_binary_function_return_type<            double,                int> { using type =      double; };
template <> struct math_binary_function_return_type<            double,               long> { using type =      double; };
template <> struct math_binary_function_return_type<            double,          long long> { using type =      double; };
template <> struct math_binary_function_return_type<            double,     unsigned short> { using type =      double; };
template <> struct math_binary_function_return_type<            double,       unsigned int> { using type =      double; };
template <> struct math_binary_function_return_type<            double,      unsigned long> { using type =      double; };
template <> struct math_binary_function_return_type<            double, unsigned long long> { using type =      double; };
template <> struct math_binary_function_return_type<             short,              float> { using type =      double; };
template <> struct math_binary_function_return_type<             short,             double> { using type =      double; };
template <> struct math_binary_function_return_type<             short,               bool> { using type =      double; };
template <> struct math_binary_function_return_type<             short,              short> { using type =      double; };
template <> struct math_binary_function_return_type<             short,                int> { using type =      double; };
template <> struct math_binary_function_return_type<             short,               long> { using type =      double; };
template <> struct math_binary_function_return_type<             short,          long long> { using type =      double; };
template <> struct math_binary_function_return_type<             short,     unsigned short> { using type =      double; };
template <> struct math_binary_function_return_type<             short,       unsigned int> { using type =      double; };
template <> struct math_binary_function_return_type<             short,      unsigned long> { using type =      double; };
template <> struct math_binary_function_return_type<             short, unsigned long long> { using type =      double; };
template <> struct math_binary_function_return_type<               int,              float> { using type =      double; };
template <> struct math_binary_function_return_type<               int,             double> { using type =      double; };
template <> struct math_binary_function_return_type<               int,               bool> { using type =      double; };
template <> struct math_binary_function_return_type<               int,              short> { using type =      double; };
template <> struct math_binary_function_return_type<               int,                int> { using type =      double; };
template <> struct math_binary_function_return_type<               int,               long> { using type =      double; };
template <> struct math_binary_function_return_type<               int,          long long> { using type =      double; };
template <> struct math_binary_function_return_type<               int,     unsigned short> { using type =      double; };
template <> struct math_binary_function_return_type<               int,       unsigned int> { using type =      double; };
template <> struct math_binary_function_return_type<               int,      unsigned long> { using type =      double; };
template <> struct math_binary_function_return_type<               int, unsigned long long> { using type =      double; };
template <> struct math_binary_function_return_type<              long,              float> { using type =      double; };
template <> struct math_binary_function_return_type<              long,             double> { using type =      double; };
template <> struct math_binary_function_return_type<              long,               bool> { using type =      double; };
template <> struct math_binary_function_return_type<              long,              short> { using type =      double; };
template <> struct math_binary_function_return_type<              long,                int> { using type =      double; };
template <> struct math_binary_function_return_type<              long,               long> { using type =      double; };
template <> struct math_binary_function_return_type<              long,          long long> { using type =      double; };
template <> struct math_binary_function_return_type<              long,     unsigned short> { using type =      double; };
template <> struct math_binary_function_return_type<              long,       unsigned int> { using type =      double; };
template <> struct math_binary_function_return_type<              long,      unsigned long> { using type =      double; };
template <> struct math_binary_function_return_type<              long, unsigned long long> { using type =      double; };
template <> struct math_binary_function_return_type<         long long,              float> { using type =      double; };
template <> struct math_binary_function_return_type<         long long,             double> { using type =      double; };
template <> struct math_binary_function_return_type<         long long,               bool> { using type =      double; };
template <> struct math_binary_function_return_type<         long long,              short> { using type =      double; };
template <> struct math_binary_function_return_type<         long long,                int> { using type =      double; };
template <> struct math_binary_function_return_type<         long long,               long> { using type =      double; };
template <> struct math_binary_function_return_type<         long long,          long long> { using type =      double; };
template <> struct math_binary_function_return_type<         long long,     unsigned short> { using type =      double; };
template <> struct math_binary_function_return_type<         long long,       unsigned int> { using type =      double; };
template <> struct math_binary_function_return_type<         long long,      unsigned long> { using type =      double; };
template <> struct math_binary_function_return_type<         long long, unsigned long long> { using type =      double; };
template <> struct math_binary_function_return_type<    unsigned short,              float> { using type =      double; };
template <> struct math_binary_function_return_type<    unsigned short,             double> { using type =      double; };
template <> struct math_binary_function_return_type<    unsigned short,               bool> { using type =      double; };
template <> struct math_binary_function_return_type<    unsigned short,              short> { using type =      double; };
template <> struct math_binary_function_return_type<    unsigned short,                int> { using type =      double; };
template <> struct math_binary_function_return_type<    unsigned short,               long> { using type =      double; };
template <> struct math_binary_function_return_type<    unsigned short,          long long> { using type =      double; };
template <> struct math_binary_function_return_type<    unsigned short,     unsigned short> { using type =      double; };
template <> struct math_binary_function_return_type<    unsigned short,       unsigned int> { using type =      double; };
template <> struct math_binary_function_return_type<    unsigned short,      unsigned long> { using type =      double; };
template <> struct math_binary_function_return_type<    unsigned short, unsigned long long> { using type =      double; };
template <> struct math_binary_function_return_type<      unsigned int,              float> { using type =      double; };
template <> struct math_binary_function_return_type<      unsigned int,             double> { using type =      double; };
template <> struct math_binary_function_return_type<      unsigned int,               bool> { using type =      double; };
template <> struct math_binary_function_return_type<      unsigned int,              short> { using type =      double; };
template <> struct math_binary_function_return_type<      unsigned int,                int> { using type =      double; };
template <> struct math_binary_function_return_type<      unsigned int,               long> { using type =      double; };
template <> struct math_binary_function_return_type<      unsigned int,          long long> { using type =      double; };
template <> struct math_binary_function_return_type<      unsigned int,     unsigned short> { using type =      double; };
template <> struct math_binary_function_return_type<      unsigned int,       unsigned int> { using type =      double; };
template <> struct math_binary_function_return_type<      unsigned int,      unsigned long> { using type =      double; };
template <> struct math_binary_function_return_type<      unsigned int, unsigned long long> { using type =      double; };
template <> struct math_binary_function_return_type<     unsigned long,              float> { using type =      double; };
template <> struct math_binary_function_return_type<     unsigned long,             double> { using type =      double; };
template <> struct math_binary_function_return_type<     unsigned long,               bool> { using type =      double; };
template <> struct math_binary_function_return_type<     unsigned long,              short> { using type =      double; };
template <> struct math_binary_function_return_type<     unsigned long,                int> { using type =      double; };
template <> struct math_binary_function_return_type<     unsigned long,               long> { using type =      double; };
template <> struct math_binary_function_return_type<     unsigned long,          long long> { using type =      double; };
template <> struct math_binary_function_return_type<     unsigned long,     unsigned short> { using type =      double; };
template <> struct math_binary_function_return_type<     unsigned long,       unsigned int> { using type =      double; };
template <> struct math_binary_function_return_type<     unsigned long,      unsigned long> { using type =      double; };
template <> struct math_binary_function_return_type<     unsigned long, unsigned long long> { using type =      double; };
template <> struct math_binary_function_return_type<unsigned long long,              float> { using type =      double; };
template <> struct math_binary_function_return_type<unsigned long long,             double> { using type =      double; };
template <> struct math_binary_function_return_type<unsigned long long,               bool> { using type =      double; };
template <> struct math_binary_function_return_type<unsigned long long,              short> { using type =      double; };
template <> struct math_binary_function_return_type<unsigned long long,                int> { using type =      double; };
template <> struct math_binary_function_return_type<unsigned long long,               long> { using type =      double; };
template <> struct math_binary_function_return_type<unsigned long long,          long long> { using type =      double; };
template <> struct math_binary_function_return_type<unsigned long long,     unsigned short> { using type =      double; };
template <> struct math_binary_function_return_type<unsigned long long,       unsigned int> { using type =      double; };
template <> struct math_binary_function_return_type<unsigned long long,      unsigned long> { using type =      double; };
template <> struct math_binary_function_return_type<unsigned long long, unsigned long long> { using type =      double; };
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
template <> struct math_binary_function_return_type<             float,        long double> { using type = long double; };
template <> struct math_binary_function_return_type<            double,        long double> { using type = long double; };
template <> struct math_binary_function_return_type<       long double,              float> { using type = long double; };
template <> struct math_binary_function_return_type<       long double,             double> { using type = long double; };
template <> struct math_binary_function_return_type<       long double,        long double> { using type = long double; };
template <> struct math_binary_function_return_type<       long double,               bool> { using type = long double; };
template <> struct math_binary_function_return_type<       long double,              short> { using type = long double; };
template <> struct math_binary_function_return_type<       long double,                int> { using type = long double; };
template <> struct math_binary_function_return_type<       long double,               long> { using type = long double; };
template <> struct math_binary_function_return_type<       long double,          long long> { using type = long double; };
template <> struct math_binary_function_return_type<       long double,     unsigned short> { using type = long double; };
template <> struct math_binary_function_return_type<       long double,       unsigned int> { using type = long double; };
template <> struct math_binary_function_return_type<       long double,      unsigned long> { using type = long double; };
template <> struct math_binary_function_return_type<       long double, unsigned long long> { using type = long double; };
template <> struct math_binary_function_return_type<             short,        long double> { using type = long double; };
template <> struct math_binary_function_return_type<               int,        long double> { using type = long double; };
template <> struct math_binary_function_return_type<              long,        long double> { using type = long double; };
template <> struct math_binary_function_return_type<         long long,        long double> { using type = long double; };
template <> struct math_binary_function_return_type<    unsigned short,        long double> { using type = long double; };
template <> struct math_binary_function_return_type<      unsigned int,        long double> { using type = long double; };
template <> struct math_binary_function_return_type<     unsigned long,        long double> { using type = long double; };
template <> struct math_binary_function_return_type<unsigned long long,        long double> { using type = long double; };
#endif
template <class T, class U>
using math_binary_function_return_type_t = typename math_binary_function_return_type<T, U>::type;
// clang-format on
template <class T, class U, class V>
using math_ternary_function_return_type_t = math_binary_function_return_type_t<
    T, math_binary_function_return_type_t<U, V>>;

struct FloatingPointComparison {
 private:
  template <class T>
  KOKKOS_FUNCTION double eps(T) const {
    return DBL_EPSILON;
  }
#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
  KOKKOS_FUNCTION
  KE::half_t eps(KE::half_t) const {
// FIXME_NVHPC compile-time error
#ifdef KOKKOS_COMPILER_NVHPC
    return 0.0009765625F;
#else
    return KE::epsilon<KE::half_t>::value;
#endif
  }
#endif
#if defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT
  KOKKOS_FUNCTION
  KE::bhalf_t eps(KE::bhalf_t) const {
// FIXME_NVHPC compile-time error
#ifdef KOKKOS_COMPILER_NVHPC
    return 0.0078125;
#else
    return KE::epsilon<KE::bhalf_t>::value;
#endif
  }
#endif
  KOKKOS_FUNCTION
  double eps(float) const { return FLT_EPSILON; }
// POWER9 gives unexpected values with LDBL_EPSILON issues
// https://stackoverflow.com/questions/68960416/ppc64-long-doubles-machine-epsilon-calculation
#if defined(KOKKOS_ARCH_POWER9) || defined(KOKKOS_ARCH_POWER8)
  KOKKOS_FUNCTION
  double eps(long double) const { return DBL_EPSILON; }
#else
  KOKKOS_FUNCTION
  double eps(long double) const { return LDBL_EPSILON; }
#endif
  // Using absolute here instead of abs, since we actually test abs ...
  template <class T>
  KOKKOS_FUNCTION std::enable_if_t<std::is_signed<T>::value, T> absolute(
      T val) const {
    return val < T(0) ? -val : val;
  }

  template <class T>
  KOKKOS_FUNCTION std::enable_if_t<!std::is_signed<T>::value, T> absolute(
      T val) const {
    return val;
  }

 public:
  template <class FPT>
  KOKKOS_FUNCTION bool compare_near_zero(FPT const& fpv, double ulp) const {
    auto abs_tol = eps(fpv) * ulp;

    bool ar = absolute(fpv) < abs_tol;
    if (!ar) {
      Kokkos::printf("absolute value exceeds tolerance [|%e| > %e]\n",
                     (double)fpv, abs_tol);
    }

    return ar;
  }

  template <class Lhs, class Rhs>
  KOKKOS_FUNCTION bool compare(Lhs const& lhs, Rhs const& rhs,
                               double ulp) const {
    if (lhs == 0) {
      return compare_near_zero(rhs, ulp);
    } else if (rhs == 0) {
      return compare_near_zero(lhs, ulp);
    } else {
      auto rel_tol     = (eps(lhs) < eps(rhs) ? eps(lhs) : eps(rhs)) * ulp;
      double abs_diff  = static_cast<double>(rhs > lhs ? rhs - lhs : lhs - rhs);
      double min_denom = static_cast<double>(
          absolute(rhs) < absolute(lhs) ? absolute(rhs) : absolute(lhs));
      double rel_diff = abs_diff / min_denom;
      bool ar         = abs_diff == 0 || rel_diff < rel_tol;
      if (!ar) {
        Kokkos::printf("relative difference exceeds tolerance [%e > %e]\n",
                       (double)rel_diff, rel_tol);
      }

      return ar;
    }
  }
};

template <class>
struct math_function_name;

#define DEFINE_UNARY_FUNCTION_EVAL(FUNC, ULP_FACTOR)                    \
  struct MathUnaryFunction_##FUNC {                                     \
    template <typename T>                                               \
    static KOKKOS_FUNCTION auto eval(T x) {                             \
      static_assert(                                                    \
          std::is_same<decltype(Kokkos::FUNC((T)0)),                    \
                       math_unary_function_return_type_t<T>>::value);   \
      return Kokkos::FUNC(x);                                           \
    }                                                                   \
    template <typename T>                                               \
    static auto eval_std(T x) {                                         \
      if constexpr (std::is_same<T, KE::half_t>::value ||               \
                    std::is_same<T, KE::bhalf_t>::value) {              \
        return std::FUNC(static_cast<float>(x));                        \
      } else {                                                          \
        static_assert(                                                  \
            std::is_same<decltype(std::FUNC((T)0)),                     \
                         math_unary_function_return_type_t<T>>::value); \
        return std::FUNC(x);                                            \
      }                                                                 \
      MATHEMATICAL_FUNCTIONS_TEST_UNREACHABLE                           \
    }                                                                   \
    static KOKKOS_FUNCTION double ulp_factor() { return ULP_FACTOR; }   \
  };                                                                    \
  using kk_##FUNC = MathUnaryFunction_##FUNC;                           \
  template <>                                                           \
  struct math_function_name<MathUnaryFunction_##FUNC> {                 \
    static constexpr char name[] = #FUNC;                               \
  };                                                                    \
  constexpr char math_function_name<MathUnaryFunction_##FUNC>::name[]

#define DEFINE_UNARY_FUNCTION_EVAL_CUSTOM(FUNC, ULP_FACTOR, REF_FUNC) \
  struct MathUnaryFunction_##FUNC {                                   \
    template <typename T>                                             \
    static KOKKOS_FUNCTION auto eval(T x) {                           \
      static_assert(                                                  \
          std::is_same<decltype(Kokkos::FUNC((T)0)),                  \
                       math_unary_function_return_type_t<T>>::value); \
      return Kokkos::FUNC(x);                                         \
    }                                                                 \
    template <typename T>                                             \
    static auto eval_std(T x) {                                       \
      static_assert(                                                  \
          std::is_same<decltype(REF_FUNC),                            \
                       math_unary_function_return_type_t<T>>::value); \
      return REF_FUNC;                                                \
    }                                                                 \
    static KOKKOS_FUNCTION double ulp_factor() { return ULP_FACTOR; } \
  };                                                                  \
  using kk_##FUNC = MathUnaryFunction_##FUNC;                         \
  template <>                                                         \
  struct math_function_name<MathUnaryFunction_##FUNC> {               \
    static constexpr char name[] = #FUNC;                             \
  };                                                                  \
  constexpr char math_function_name<MathUnaryFunction_##FUNC>::name[]

#ifndef KOKKOS_MATHEMATICAL_FUNCTIONS_SKIP_3
// Generally the expected ULP error should come from here:
// https://www.gnu.org/software/libc/manual/html_node/Errors-in-Math-Functions.html
// For now 1s largely seem to work ...
DEFINE_UNARY_FUNCTION_EVAL(exp, 2);
DEFINE_UNARY_FUNCTION_EVAL(exp2, 2);
DEFINE_UNARY_FUNCTION_EVAL(expm1, 2);
DEFINE_UNARY_FUNCTION_EVAL(log, 2);
DEFINE_UNARY_FUNCTION_EVAL(log10, 2);
DEFINE_UNARY_FUNCTION_EVAL(log2, 2);
DEFINE_UNARY_FUNCTION_EVAL(log1p, 2);
#endif

#ifndef KOKKOS_MATHEMATICAL_FUNCTIONS_SKIP_2
DEFINE_UNARY_FUNCTION_EVAL(sqrt, 2);
DEFINE_UNARY_FUNCTION_EVAL(cbrt, 2);
#endif

#ifndef KOKKOS_MATHEMATICAL_FUNCTIONS_SKIP_1
DEFINE_UNARY_FUNCTION_EVAL(sin, 2);
DEFINE_UNARY_FUNCTION_EVAL(cos, 2);
DEFINE_UNARY_FUNCTION_EVAL(tan, 2);
DEFINE_UNARY_FUNCTION_EVAL(asin, 2);
DEFINE_UNARY_FUNCTION_EVAL(acos, 2);
DEFINE_UNARY_FUNCTION_EVAL(atan, 2);

DEFINE_UNARY_FUNCTION_EVAL(sinh, 2);
DEFINE_UNARY_FUNCTION_EVAL(cosh, 2);
DEFINE_UNARY_FUNCTION_EVAL(tanh, 2);
DEFINE_UNARY_FUNCTION_EVAL(asinh, 4);
DEFINE_UNARY_FUNCTION_EVAL(acosh, 2);
DEFINE_UNARY_FUNCTION_EVAL(atanh, 2);

// non-standard math functions
DEFINE_UNARY_FUNCTION_EVAL_CUSTOM(rsqrt, 2,
                                  decltype(std::sqrt(x))(1) / std::sqrt(x));
#endif

#ifndef KOKKOS_MATHEMATICAL_FUNCTIONS_SKIP_2
#if defined(__APPLE__)
// Apple's standard library implementation seems to have a poor implementation
DEFINE_UNARY_FUNCTION_EVAL(erf, 5);
#else
DEFINE_UNARY_FUNCTION_EVAL(erf, 2);
#endif

DEFINE_UNARY_FUNCTION_EVAL(erfc, 5);
// has a larger error due to some impls doing integer exact.
// We cast always to double leading to larger difference when comparing our
// tgamma to std::tgamma on the host.
DEFINE_UNARY_FUNCTION_EVAL(tgamma, 200);
DEFINE_UNARY_FUNCTION_EVAL(lgamma, 2);

DEFINE_UNARY_FUNCTION_EVAL(ceil, 2);
DEFINE_UNARY_FUNCTION_EVAL(floor, 2);
DEFINE_UNARY_FUNCTION_EVAL(trunc, 2);
DEFINE_UNARY_FUNCTION_EVAL(round, 1);
#ifndef KOKKOS_ENABLE_SYCL
DEFINE_UNARY_FUNCTION_EVAL(nearbyint, 2);
#endif

DEFINE_UNARY_FUNCTION_EVAL(logb, 2);
#endif

#undef DEFINE_UNARY_FUNCTION_EVAL

#define DEFINE_BINARY_FUNCTION_EVAL(FUNC, ULP_FACTOR)                          \
  struct MathBinaryFunction_##FUNC {                                           \
    template <typename T, typename U>                                          \
    static KOKKOS_FUNCTION auto eval(T x, U y) {                               \
      static_assert(                                                           \
          std::is_same<decltype(Kokkos::FUNC((T)0, (U)0)),                     \
                       math_binary_function_return_type_t<T, U>>::value);      \
      return Kokkos::FUNC(x, y);                                               \
    }                                                                          \
    template <typename T, typename U>                                          \
    static auto eval_std(T x, U y) {                                           \
      constexpr bool const x_is_half =                                         \
          (KE::Impl::is_float16<T>::value || KE::Impl::is_bfloat16<T>::value); \
      constexpr bool const y_is_half =                                         \
          (KE::Impl::is_float16<U>::value || KE::Impl::is_bfloat16<U>::value); \
      if constexpr (x_is_half && y_is_half)                                    \
        return std::FUNC(static_cast<float>(x), static_cast<float>(y));        \
      else if constexpr (x_is_half)                                            \
        return std::FUNC(static_cast<float>(x), y);                            \
      else if constexpr (y_is_half)                                            \
        return std::FUNC(x, static_cast<float>(y));                            \
      else {                                                                   \
        static_assert(                                                         \
            std::is_same<decltype(std::FUNC((T)0, (U)0)),                      \
                         math_binary_function_return_type_t<T, U>>::value);    \
        return std::FUNC(x, y);                                                \
      }                                                                        \
      MATHEMATICAL_FUNCTIONS_TEST_UNREACHABLE                                  \
    }                                                                          \
    static KOKKOS_FUNCTION double ulp_factor() { return ULP_FACTOR; }          \
  };                                                                           \
  using kk_##FUNC = MathBinaryFunction_##FUNC;                                 \
  template <>                                                                  \
  struct math_function_name<MathBinaryFunction_##FUNC> {                       \
    static constexpr char name[] = #FUNC;                                      \
  };                                                                           \
  constexpr char math_function_name<MathBinaryFunction_##FUNC>::name[]

#ifndef KOKKOS_MATHEMATICAL_FUNCTIONS_SKIP_2
DEFINE_BINARY_FUNCTION_EVAL(pow, 2);
DEFINE_BINARY_FUNCTION_EVAL(hypot, 2);
DEFINE_BINARY_FUNCTION_EVAL(nextafter, 1);
DEFINE_BINARY_FUNCTION_EVAL(copysign, 1);
#endif

#undef DEFINE_BINARY_FUNCTION_EVAL

#define DEFINE_TERNARY_FUNCTION_EVAL(FUNC, ULP_FACTOR)                        \
  struct MathTernaryFunction_##FUNC {                                         \
    template <typename T, typename U, typename V>                             \
    static KOKKOS_FUNCTION auto eval(T x, U y, V z) {                         \
      static_assert(                                                          \
          std::is_same<decltype(Kokkos::FUNC((T)0, (U)0, (V)0)),              \
                       math_ternary_function_return_type_t<T, U, V>>::value); \
      return Kokkos::FUNC(x, y, z);                                           \
    }                                                                         \
    template <typename T, typename U, typename V>                             \
    static auto eval_std(T x, U y, V z) {                                     \
      static_assert(                                                          \
          std::is_same<decltype(std::FUNC((T)0, (U)0, (V)0)),                 \
                       math_ternary_function_return_type_t<T, U, V>>::value); \
      return std::FUNC(x, y, z);                                              \
    }                                                                         \
    static KOKKOS_FUNCTION double ulp_factor() { return ULP_FACTOR; }         \
  };                                                                          \
  using kk3_##FUNC = MathTernaryFunction_##FUNC;                              \
  template <>                                                                 \
  struct math_function_name<MathTernaryFunction_##FUNC> {                     \
    static constexpr char name[] = #FUNC;                                     \
  };                                                                          \
  constexpr char math_function_name<MathTernaryFunction_##FUNC>::name[]

#ifndef KOKKOS_MATHEMATICAL_FUNCTIONS_SKIP_2
DEFINE_TERNARY_FUNCTION_EVAL(hypot, 2);
DEFINE_TERNARY_FUNCTION_EVAL(fma, 2);
#endif

#undef DEFINE_TERNARY_FUNCTION_EVAL

// clang-format off
template <class>
struct type_helper;
#define DEFINE_TYPE_NAME(T) \
template <> struct type_helper<T> { static char const * name() { return #T; } };
DEFINE_TYPE_NAME(bool)
DEFINE_TYPE_NAME(int)
DEFINE_TYPE_NAME(long)
DEFINE_TYPE_NAME(long long)
DEFINE_TYPE_NAME(unsigned int)
DEFINE_TYPE_NAME(unsigned long)
DEFINE_TYPE_NAME(unsigned long long)
#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
DEFINE_TYPE_NAME(KE::half_t)
#endif
#if defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT
DEFINE_TYPE_NAME(KE::bhalf_t)
#endif
DEFINE_TYPE_NAME(float)
DEFINE_TYPE_NAME(double)
DEFINE_TYPE_NAME(long double)
#undef DEFINE_TYPE_NAME
// clang-format on

template <class Space, class Func, class Arg, std::size_t N,
          class Ret = math_unary_function_return_type_t<Arg>>
struct TestMathUnaryFunction : FloatingPointComparison {
  Arg val_[N];
  Ret res_[N];
  TestMathUnaryFunction(const Arg (&val)[N]) {
    std::copy(val, val + N, val_);
    std::transform(val, val + N, res_,
                   [](auto x) { return Func::eval_std(x); });
    run();
  }
  void run() {
    int errors = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Space>(0, N), *this, errors);
    ASSERT_EQ(errors, 0) << "Failed check no error for "
                         << math_function_name<Func>::name << "("
                         << type_helper<Arg>::name() << ")";
  }
  KOKKOS_FUNCTION void operator()(int i, int& e) const {
    bool ar = compare(Func::eval(val_[i]), res_[i], Func::ulp_factor());
    if (!ar) {
      ++e;
      Kokkos::printf("value at %f which is %f was expected to be %f\n",
                     (double)val_[i], (double)Func::eval(val_[i]),
                     (double)res_[i]);
    }
  }
};

template <class Space, class... Func, class Arg, std::size_t N>
void do_test_math_unary_function(const Arg (&x)[N]) {
  (void)std::initializer_list<int>{
      (TestMathUnaryFunction<Space, Func, Arg, N>(x), 0)...};

  // test if potentially device specific math functions also work on host
  if constexpr (!std::is_same_v<Space, Kokkos::DefaultHostExecutionSpace>)
    (void)std::initializer_list<int>{
        (TestMathUnaryFunction<Kokkos::DefaultHostExecutionSpace, Func, Arg, N>(
             x),
         0)...};
}

#define TEST_MATH_FUNCTION(FUNC) \
  do_test_math_unary_function<TEST_EXECSPACE, MathUnaryFunction_##FUNC>

template <class Half, class Space, class... Func, class Arg, std::size_t N>
void do_test_half_math_unary_function(const Arg (&x)[N]) {
  Half y[N];
  std::copy(x, x + N, y);  // cast to array of half type
  (void)std::initializer_list<int>{
      (TestMathUnaryFunction<Space, Func, Half, N>(y), 0)...};

  // test if potentially device specific math functions also work on host
  if constexpr (!std::is_same_v<Space, Kokkos::DefaultHostExecutionSpace>)
    (void)std::initializer_list<int>{(
        TestMathUnaryFunction<Kokkos::DefaultHostExecutionSpace, Func, Half, N>(
            y),
        0)...};
}

#define TEST_HALF_MATH_FUNCTION(FUNC, T) \
  do_test_half_math_unary_function<T, TEST_EXECSPACE, MathUnaryFunction_##FUNC>

template <class Space, class Func, class Arg1, class Arg2,
          class Ret = math_binary_function_return_type_t<Arg1, Arg2>>
struct TestMathBinaryFunction : FloatingPointComparison {
  Arg1 val1_;
  Arg2 val2_;
  Ret res_;
  TestMathBinaryFunction(Arg1 val1, Arg2 val2)
      : val1_(val1), val2_(val2), res_(Func::eval_std(val1, val2)) {
    run();
  }
  void run() {
    int errors = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Space>(0, 1), *this, errors);
    ASSERT_EQ(errors, 0) << "Failed check no error for "
                         << math_function_name<Func>::name << "("
                         << type_helper<Arg1>::name() << ", "
                         << type_helper<Arg2>::name() << ")";
  }
  KOKKOS_FUNCTION void operator()(int, int& e) const {
    bool ar = compare(Func::eval(val1_, val2_), res_, Func::ulp_factor());
    if (!ar) {
      ++e;
      Kokkos::printf("value at %f, %f which is %f was expected to be %f\n",
                     (double)val1_, (double)val2_,
                     (double)Func::eval(val1_, val2_), (double)res_);
    }
  }
};

template <class Space, class... Func, class Arg1, class Arg2>
void do_test_math_binary_function(Arg1 arg1, Arg2 arg2) {
  (void)std::initializer_list<int>{
      (TestMathBinaryFunction<Space, Func, Arg1, Arg2>(arg1, arg2), 0)...};
}

template <class Space, class Func, class Arg1, class Arg2, class Arg3,
          class Ret = math_ternary_function_return_type_t<Arg1, Arg2, Arg3>>
struct TestMathTernaryFunction : FloatingPointComparison {
  Arg1 val1_;
  Arg2 val2_;
  Arg3 val3_;
  Ret res_;
  TestMathTernaryFunction(Arg1 val1, Arg2 val2, Arg3 val3)
      : val1_(val1),
        val2_(val2),
        val3_(val3),
        res_(Func::eval_std(val1, val2, val3)) {
    run();
  }
  void run() {
    int errors = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Space>(0, 1), *this, errors);
    ASSERT_EQ(errors, 0) << "Failed check no error for "
                         << math_function_name<Func>::name << "("
                         << type_helper<Arg1>::name() << ", "
                         << type_helper<Arg1>::name() << ", "
                         << type_helper<Arg3>::name() << ")";
  }
  KOKKOS_FUNCTION void operator()(int, int& e) const {
    bool ar =
        compare(Func::eval(val1_, val2_, val3_), res_, Func::ulp_factor());
    if (!ar) {
      ++e;
      Kokkos::printf("value at %f, %f, %f which is %f was expected to be %f\n",
                     (double)val1_, (double)val2_, (double)val3_,
                     (double)Func::eval(val1_, val2_, val3_), (double)res_);
    }
  }
};

template <class Space, class... Func, class Arg1, class Arg2, class Arg3>
void do_test_math_ternary_function(Arg1 arg1, Arg2 arg2, Arg3 arg3) {
  (void)std::initializer_list<int>{
      (TestMathTernaryFunction<Space, Func, Arg1, Arg2, Arg3>(arg1, arg2, arg3),
       0)...};
}

#ifndef KOKKOS_MATHEMATICAL_FUNCTIONS_SKIP_1

TEST(TEST_CATEGORY, mathematical_functions_trigonometric_functions) {
  TEST_MATH_FUNCTION(sin)({true, false});
  TEST_MATH_FUNCTION(sin)({-3, -2, -1, 0, 1});
  TEST_MATH_FUNCTION(sin)({-3l, -2l, -1l, 0l, 1l});
  TEST_MATH_FUNCTION(sin)({-3ll, -2ll, -1ll, 0ll, 1ll});
  TEST_MATH_FUNCTION(sin)({2u, 3u, 4u, 5u, 6u});
  TEST_MATH_FUNCTION(sin)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_MATH_FUNCTION(sin)({2ull, 3ull, 4ull, 5ull, 6ull});
  TEST_HALF_MATH_FUNCTION(sin, KE::half_t)({.1f, .2f, .3f});
  TEST_HALF_MATH_FUNCTION(sin, KE::bhalf_t)({.1f, .2f, .3f});
  TEST_MATH_FUNCTION(sin)({.1f, .2f, .3f});
  TEST_MATH_FUNCTION(sin)({.4, .5, .6});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(sin)({.7l, .8l, .9l});
#endif

  TEST_MATH_FUNCTION(cos)({true, false});
  TEST_MATH_FUNCTION(cos)({-3, -2, -1, 0, 1});
  TEST_MATH_FUNCTION(cos)({-3l, -2l, -1l, 0l, 1l});
  TEST_MATH_FUNCTION(cos)({-3ll, -2ll, -1ll, 0ll, 1ll});
  TEST_MATH_FUNCTION(cos)({2u, 3u, 4u, 5u, 6u});
  TEST_MATH_FUNCTION(cos)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_MATH_FUNCTION(cos)({2ull, 3ull, 4ull, 5ull, 6ull});
  TEST_HALF_MATH_FUNCTION(cos, KE::half_t)({.1f, .2f, .3f});
  TEST_HALF_MATH_FUNCTION(cos, KE::bhalf_t)({.1f, .2f, .3f});
  TEST_MATH_FUNCTION(cos)({.1f, .2f, .3f});
  TEST_MATH_FUNCTION(cos)({.4, .5, .6});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(cos)({.7l, .8l, .9l});
#endif

  TEST_MATH_FUNCTION(tan)({true, false});
  TEST_MATH_FUNCTION(tan)({-3, -2, -1, 0, 1});
  TEST_MATH_FUNCTION(tan)({-3l, -2l, -1l, 0l, 1l});
  TEST_MATH_FUNCTION(tan)({-3ll, -2ll, -1ll, 0ll, 1ll});
  TEST_MATH_FUNCTION(tan)({2u, 3u, 4u, 5u, 6u});
  TEST_MATH_FUNCTION(tan)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_MATH_FUNCTION(tan)({2ull, 3ull, 4ull, 5ull, 6ull});
  TEST_HALF_MATH_FUNCTION(tan, KE::half_t)({.1f, .2f, .3f});
  TEST_HALF_MATH_FUNCTION(tan, KE::bhalf_t)({.1f, .2f, .3f});
  TEST_MATH_FUNCTION(tan)({.1f, .2f, .3f});
  TEST_MATH_FUNCTION(tan)({.4, .5, .6});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(tan)({.7l, .8l, .9l});
#endif

  TEST_MATH_FUNCTION(asin)({true, false});
  TEST_MATH_FUNCTION(asin)({-1, 0, 1});
  TEST_MATH_FUNCTION(asin)({-1l, 0l, 1l});
  TEST_MATH_FUNCTION(asin)({-1ll, 0ll, 1ll});
  TEST_MATH_FUNCTION(asin)({0u, 1u});
  TEST_MATH_FUNCTION(asin)({0ul, 1ul});
  TEST_MATH_FUNCTION(asin)({0ull, 1ull});
  TEST_HALF_MATH_FUNCTION(asin, KE::half_t)({-1.f, .9f, -.8f, .7f, -.6f});
  TEST_HALF_MATH_FUNCTION(asin, KE::bhalf_t)({-1.f, .9f, -.8f, .7f, -.6f});
  TEST_MATH_FUNCTION(asin)({-1.f, .9f, -.8f, .7f, -.6f});
  TEST_MATH_FUNCTION(asin)({-.5, .4, -.3, .2, -.1, 0.});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(asin)({-.5l, .3l, 0.l, .2l, .4l, .6l});
#endif

  TEST_MATH_FUNCTION(acos)({true, false});
  TEST_MATH_FUNCTION(acos)({-1, 0, 1});
  TEST_MATH_FUNCTION(acos)({-1l, 0l, 1l});
  TEST_MATH_FUNCTION(acos)({-1ll, 0ll, 1ll});
  TEST_MATH_FUNCTION(acos)({0u, 1u});
  TEST_MATH_FUNCTION(acos)({0ul, 1ul});
  TEST_MATH_FUNCTION(acos)({0ull, 1ull});
  TEST_HALF_MATH_FUNCTION(acos, KE::half_t)({-1.f, .9f, -.8f, .7f, -.6f});
  TEST_HALF_MATH_FUNCTION(acos, KE::bhalf_t)({-1.f, .9f, -.8f, .7f, -.6f});
  TEST_MATH_FUNCTION(acos)({-1.f, .9f, -.8f, .7f, -.6f});
  TEST_MATH_FUNCTION(acos)({-.5, .4, -.3, .2, -.1, 0.});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(acos)({-.5l, .3l, 0.l, .2l, .4l, .6l});
#endif

  TEST_MATH_FUNCTION(atan)({true, false});
  TEST_MATH_FUNCTION(atan)({-1, 0, 1});
  TEST_MATH_FUNCTION(atan)({-1l, 0l, 1l});
  TEST_MATH_FUNCTION(atan)({-1ll, 0ll, 1ll});
  TEST_MATH_FUNCTION(atan)({0u, 1u});
  TEST_MATH_FUNCTION(atan)({0ul, 1ul});
  TEST_MATH_FUNCTION(atan)({0ull, 1ull});
  TEST_HALF_MATH_FUNCTION(atan, KE::half_t)
  ({-1.5f, 1.3f, -1.1f, .9f, -.7f, .5f});
  TEST_HALF_MATH_FUNCTION(atan, KE::bhalf_t)
  ({-1.5f, 1.3f, -1.1f, .9f, -.7f, .5f});
  TEST_MATH_FUNCTION(atan)({-1.5f, 1.3f, -1.1f, .9f, -.7f, .5f});
  TEST_MATH_FUNCTION(atan)({1.4, -1.2, 1., -.8, .6, -.4, .2, -0.});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(atan)({-.98l, .67l, -54.l, .34l, -.21l});
#endif

  // TODO atan2
}
#endif

#ifndef KOKKOS_MATHEMATICAL_FUNCTIONS_SKIP_2
TEST(TEST_CATEGORY, mathematical_functions_power_functions) {
  TEST_MATH_FUNCTION(sqrt)({0, 1, 2, 3, 5, 7, 11});
  TEST_MATH_FUNCTION(sqrt)({0l, 1l, 2l, 3l, 5l, 7l, 11l});
  TEST_MATH_FUNCTION(sqrt)({0ll, 1ll, 2ll, 3ll, 5ll, 7ll, 11ll});
  TEST_MATH_FUNCTION(sqrt)({0u, 1u, 2u, 3u, 5u, 7u});
  TEST_MATH_FUNCTION(sqrt)({0ul, 1ul, 2ul, 3ul, 5ul, 7ul});
  TEST_MATH_FUNCTION(sqrt)({0ull, 1ull, 2ull, 3ull, 5ull, 7ull});
  TEST_HALF_MATH_FUNCTION(sqrt, KE::half_t)({10.f, 20.f, 30.f, 40.f});
  TEST_HALF_MATH_FUNCTION(sqrt, KE::bhalf_t)({10.f, 20.f, 30.f, 40.f});
  TEST_MATH_FUNCTION(sqrt)({10.f, 20.f, 30.f, 40.f});
  TEST_MATH_FUNCTION(sqrt)({11.1, 22.2, 33.3, 44.4});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(sqrt)({10.l, 20.l, 30.l, 40.l});
#endif

  TEST_MATH_FUNCTION(cbrt)({-5, -3, -1, 2, 4, 6});
  TEST_MATH_FUNCTION(cbrt)({-5l, -3l, -1l, 2l, 4l, 6l});
  TEST_MATH_FUNCTION(cbrt)({-5ll, -3ll, -1ll, 2ll, 4ll, 6ll});
  TEST_MATH_FUNCTION(cbrt)({0u, 1u, 2u, 3u, 4u, 5u});
  TEST_MATH_FUNCTION(cbrt)({0ul, 1ul, 2ul, 3ul, 4ul, 5ul});
  TEST_MATH_FUNCTION(cbrt)({0ull, 1ull, 2ull, 3ull, 4ull, 5ull});
  TEST_HALF_MATH_FUNCTION(cbrt, KE::half_t)({-1.f, .2f, -3.f, .4f, -5.f});
  TEST_HALF_MATH_FUNCTION(cbrt, KE::bhalf_t)({-1.f, .2f, -3.f, .4f, -5.f});
  TEST_MATH_FUNCTION(cbrt)({-1.f, .2f, -3.f, .4f, -5.f});
  TEST_MATH_FUNCTION(cbrt)({11.1, -2.2, 33.3, -4.4, 55.5});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(cbrt)({-10.l, 20.l, -30.l, 40.l, -50.l});
#endif

  do_test_math_binary_function<TEST_EXECSPACE, kk_pow>(
      static_cast<KE::half_t>(2.f), static_cast<KE::half_t>(3.f));
  do_test_math_binary_function<TEST_EXECSPACE, kk_pow>(
      static_cast<KE::bhalf_t>(2.f), static_cast<KE::bhalf_t>(3.f));
  do_test_math_binary_function<TEST_EXECSPACE, kk_pow>(2.f, 3.f);
  do_test_math_binary_function<TEST_EXECSPACE, kk_pow>(2.f, 3.f);
  do_test_math_binary_function<TEST_EXECSPACE, kk_pow>(2., 3.);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  do_test_math_binary_function<TEST_EXECSPACE, kk_pow>(2.l, 3.l);
#endif

  do_test_math_binary_function<TEST_EXECSPACE, kk_hypot>(
      static_cast<KE::half_t>(2.f), static_cast<KE::half_t>(3.f));
  do_test_math_binary_function<TEST_EXECSPACE, kk_hypot>(
      static_cast<KE::bhalf_t>(2.f), static_cast<KE::bhalf_t>(3.f));
  do_test_math_binary_function<TEST_EXECSPACE, kk_hypot>(2.f, 3.f);
  do_test_math_binary_function<TEST_EXECSPACE, kk_hypot>(2., 3.);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
// FIXME: fails with gcc on Power platforms
#if !(defined(KOKKOS_ARCH_POWER8) || defined(KOKKOS_ARCH_POWER9))
  do_test_math_binary_function<TEST_EXECSPACE, kk_hypot>(2.l, 3.l);
#endif
#endif

  do_test_math_ternary_function<TEST_EXECSPACE, kk3_hypot>(2.f, 3.f, 4.f);
  do_test_math_ternary_function<TEST_EXECSPACE, kk3_hypot>(2., 3., 4.);
  do_test_math_ternary_function<TEST_EXECSPACE, kk3_hypot>(2, 3.f, 4.);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
#if !(defined(KOKKOS_ARCH_POWER8) || defined(KOKKOS_ARCH_POWER9))
  do_test_math_ternary_function<TEST_EXECSPACE, kk3_hypot>(2.l, 3.l, 4.l);
#endif
#endif
}

TEST(TEST_CATEGORY, mathematical_functions_fma) {
  do_test_math_ternary_function<TEST_EXECSPACE, kk3_fma>(2.f, 3.f, 4.f);
  do_test_math_ternary_function<TEST_EXECSPACE, kk3_fma>(2., 3., 4.);
  do_test_math_ternary_function<TEST_EXECSPACE, kk3_fma>(2, 3.f, 4.);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  do_test_math_ternary_function<TEST_EXECSPACE, kk3_fma>(2.l, 3.l, 4.l);
#endif
}
#endif

#ifndef KOKKOS_MATHEMATICAL_FUNCTIONS_SKIP_3
TEST(TEST_CATEGORY, mathematical_functions_exponential_functions) {
  TEST_MATH_FUNCTION(exp)({-9, -8, -7, -6, -5, 4, 3, 2, 1, 0});
  TEST_MATH_FUNCTION(exp)({-9l, -8l, -7l, -6l, -5l, 4l, 3l, 2l, 1l, 0l});
  TEST_MATH_FUNCTION(exp)({-9ll, -8ll, -7ll, -6ll, -5ll, 4ll, 3ll, 2ll, 1ll});
  TEST_MATH_FUNCTION(exp)({0u, 1u, 2u, 3u, 4u, 5u});
  TEST_MATH_FUNCTION(exp)({0ul, 1ul, 2ul, 3ul, 4ul, 5ul});
  TEST_MATH_FUNCTION(exp)({0ull, 1ull, 2ull, 3ull, 4ull, 5ull});
  TEST_HALF_MATH_FUNCTION(exp, KE::half_t)
  ({-98.f, -7.6f, -.54f, 3.2f, 1.f, -0.f});
  TEST_HALF_MATH_FUNCTION(exp, KE::bhalf_t)
  ({-98.f, -7.6f, -.54f, 3.2f, 1.f, -0.f});
  TEST_MATH_FUNCTION(exp)({-98.f, -7.6f, -.54f, 3.2f, 1.f, -0.f});
  TEST_MATH_FUNCTION(exp)({-98., -7.6, -.54, 3.2, 1., -0.});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(exp)({-98.l, -7.6l, -.54l, 3.2l, 1.l, -0.l});
#endif

  TEST_MATH_FUNCTION(exp2)({-9, -8, -7, -6, -5, 4, 3, 2, 1, 0});
  TEST_MATH_FUNCTION(exp2)({-9l, -8l, -7l, -6l, -5l, 4l, 3l, 2l, 1l, 0l});
  TEST_MATH_FUNCTION(exp2)({-9ll, -8ll, -7ll, -6ll, -5ll, 4ll, 3ll, 2ll, 1ll});
  TEST_MATH_FUNCTION(exp2)({0u, 1u, 2u, 3u, 4u, 5u});
  TEST_MATH_FUNCTION(exp2)({0ul, 1ul, 2ul, 3ul, 4ul, 5ul});
  TEST_MATH_FUNCTION(exp2)({0ull, 1ull, 2ull, 3ull, 4ull, 5ull});
  TEST_HALF_MATH_FUNCTION(exp2, KE::half_t)
  ({-98.f, -7.6f, -.54f, 3.2f, 1.f, -0.f});
  TEST_HALF_MATH_FUNCTION(exp2, KE::bhalf_t)
  ({-98.f, -7.6f, -.54f, 3.2f, 1.f, -0.f});
  TEST_MATH_FUNCTION(exp2)({-98.f, -7.6f, -.54f, 3.2f, 1.f, -0.f});
  TEST_MATH_FUNCTION(exp2)({-98., -7.6, -.54, 3.2, 1., -0.});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(exp2)({-98.l, -7.6l, -.54l, 3.2l, 1.l, -0.l});
#endif

  TEST_MATH_FUNCTION(expm1)({-9, -8, -7, -6, -5, 4, 3, 2, 1, 0});
  TEST_MATH_FUNCTION(expm1)({-9l, -8l, -7l, -6l, -5l, 4l, 3l, 2l, 1l, 0l});
  TEST_MATH_FUNCTION(expm1)({-9ll, -8ll, -7ll, -6ll, -5ll, 4ll, 3ll, 2ll, 1ll});
  TEST_MATH_FUNCTION(expm1)({0u, 1u, 2u, 3u, 4u, 5u});
  TEST_MATH_FUNCTION(expm1)({0ul, 1ul, 2ul, 3ul, 4ul, 5ul});
  TEST_MATH_FUNCTION(expm1)({0ull, 1ull, 2ull, 3ull, 4ull, 5ull});
  TEST_HALF_MATH_FUNCTION(expm1, KE::half_t)
  ({-98.f, -7.6f, -.54f, 3.2f, 1.f, -0.f});
  TEST_HALF_MATH_FUNCTION(expm1, KE::bhalf_t)
  ({-98.f, -7.6f, -.54f, 3.2f, 1.f, -0.f});
  TEST_MATH_FUNCTION(expm1)({-98.f, -7.6f, -.54f, 3.2f, 1.f, -0.f});
  TEST_MATH_FUNCTION(expm1)({-98., -7.6, -.54, 3.2, 1., -0.});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(expm1)({-98.l, -7.6l, -.54l, 3.2l, 1.l, -0.l});
#endif

  TEST_MATH_FUNCTION(log)({1, 23, 456, 7890});
  TEST_MATH_FUNCTION(log)({1l, 23l, 456l, 7890l});
  TEST_MATH_FUNCTION(log)({1ll, 23ll, 456ll, 7890ll});
  TEST_MATH_FUNCTION(log)({1u, 23u, 456u, 7890u});
  TEST_MATH_FUNCTION(log)({1ul, 23ul, 456ul, 7890ul});
  TEST_MATH_FUNCTION(log)({1ull, 23ull, 456ull, 7890ull});
  TEST_HALF_MATH_FUNCTION(log, KE::half_t)({1234.f, 567.f, 89.f, .1f});
  TEST_HALF_MATH_FUNCTION(log, KE::bhalf_t)({1234.f, 567.f, 89.f, .1f});
  TEST_MATH_FUNCTION(log)({1234.f, 567.f, 89.f, .1f});
  TEST_MATH_FUNCTION(log)({1234., 567., 89., .02});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(log)({1234.l, 567.l, 89.l, .003l});
#endif

  TEST_MATH_FUNCTION(log10)({1, 23, 456, 7890});
  TEST_MATH_FUNCTION(log10)({1l, 23l, 456l, 7890l});
  TEST_MATH_FUNCTION(log10)({1ll, 23ll, 456ll, 7890ll});
  TEST_MATH_FUNCTION(log10)({1u, 23u, 456u, 7890u});
  TEST_MATH_FUNCTION(log10)({1ul, 23ul, 456ul, 7890ul});
  TEST_MATH_FUNCTION(log10)({1ull, 23ull, 456ull, 7890ull});
  TEST_HALF_MATH_FUNCTION(log10, KE::half_t)({1234.f, 567.f, 89.f, .1f});
  TEST_HALF_MATH_FUNCTION(log10, KE::bhalf_t)({1234.f, 567.f, 89.f, .1f});
  TEST_MATH_FUNCTION(log10)({1234.f, 567.f, 89.f, .1f});
  TEST_MATH_FUNCTION(log10)({1234., 567., 89., .02});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(log10)({1234.l, 567.l, 89.l, .003l});
#endif

// FIXME_OPENMPTARGET FIXME_AMD
#if defined(KOKKOS_ENABLE_OPENMPTARGET) &&                                 \
    (defined(KOKKOS_ARCH_AMD_GFX906) || defined(KOKKOS_ARCH_AMD_GFX908) || \
     defined(KOKKOS_ARCH_AMD_GFX90A) || defined(KOKKOS_ARCH_AMD_GFX942))

  TEST_MATH_FUNCTION(log2)({1, 23, 456, 7890});
#endif
  TEST_MATH_FUNCTION(log2)({1l, 23l, 456l, 7890l});
  TEST_MATH_FUNCTION(log2)({1ll, 23ll, 456ll, 7890ll});
  TEST_MATH_FUNCTION(log2)({1u, 23u, 456u, 7890u});
  TEST_MATH_FUNCTION(log2)({1ul, 23ul, 456ul, 7890ul});
  TEST_MATH_FUNCTION(log2)({1ull, 23ull, 456ull, 7890ull});
  TEST_HALF_MATH_FUNCTION(log2, KE::half_t)({1234.f, 567.f, 89.f, .1f});
  TEST_HALF_MATH_FUNCTION(log2, KE::bhalf_t)({1234.f, 567.f, 89.f, .1f});
  TEST_MATH_FUNCTION(log2)({1234.f, 567.f, 89.f, .1f});
  TEST_MATH_FUNCTION(log2)({1234., 567., 89., .02});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(log2)({1234.l, 567.l, 89.l, .003l});
#endif

  TEST_MATH_FUNCTION(log1p)({1, 23, 456, 7890, 0});
  TEST_MATH_FUNCTION(log1p)({1l, 23l, 456l, 7890l, 0l});
  TEST_MATH_FUNCTION(log1p)({1ll, 23ll, 456ll, 7890ll, 0ll});
  TEST_MATH_FUNCTION(log1p)({1u, 23u, 456u, 7890u, 0u});
  TEST_MATH_FUNCTION(log1p)({1ul, 23ul, 456ul, 7890ul, 0ul});
  TEST_MATH_FUNCTION(log1p)({1ull, 23ull, 456ull, 7890ull, 0ull});
  TEST_HALF_MATH_FUNCTION(log1p, KE::half_t)({1234.f, 567.f, 89.f, -.9f});
  TEST_HALF_MATH_FUNCTION(log1p, KE::bhalf_t)({1234.f, 567.f, 89.f, -.9f});
  TEST_MATH_FUNCTION(log1p)({1234.f, 567.f, 89.f, -.9f});
  TEST_MATH_FUNCTION(log1p)({1234., 567., 89., -.08});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(log1p)({1234.l, 567.l, 89.l, -.007l});
#endif
}
#endif

#ifndef KOKKOS_MATHEMATICAL_FUNCTIONS_SKIP_1
TEST(TEST_CATEGORY, mathematical_functions_hyperbolic_functions) {
  TEST_MATH_FUNCTION(sinh)({-3, -2, -1, 0, 1});
  TEST_MATH_FUNCTION(sinh)({-3l, -2l, -1l, 0l, 1l});
  TEST_MATH_FUNCTION(sinh)({-3ll, -2ll, -1ll, 0ll, 1ll});
  TEST_MATH_FUNCTION(sinh)({2u, 3u, 4u, 5u, 6u});
  TEST_MATH_FUNCTION(sinh)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_MATH_FUNCTION(sinh)({2ull, 3ull, 4ull, 5ull, 6ull});
  TEST_HALF_MATH_FUNCTION(sinh, KE::half_t)({.1f, -2.f, 3.f});
  TEST_HALF_MATH_FUNCTION(sinh, KE::bhalf_t)({.1f, -2.f, 3.f});
  TEST_MATH_FUNCTION(sinh)({.1f, -2.f, 3.f});
  TEST_MATH_FUNCTION(sinh)({-4., .5, -.6});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(sinh)({.7l, .8l, .9l});
#endif

  TEST_MATH_FUNCTION(cosh)({-3, -2, -1, 0, 1});
  TEST_MATH_FUNCTION(cosh)({-3l, -2l, -1l, 0l, 1l});
  TEST_MATH_FUNCTION(cosh)({-3ll, -2ll, -1ll, 0ll, 1ll});
  TEST_MATH_FUNCTION(cosh)({2u, 3u, 4u, 5u, 6u});
  TEST_MATH_FUNCTION(cosh)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_MATH_FUNCTION(cosh)({2ull, 3ull, 4ull, 5ull, 6ull});
  TEST_HALF_MATH_FUNCTION(cosh, KE::half_t)({.1f, -2.f, 3.f});
  TEST_HALF_MATH_FUNCTION(cosh, KE::bhalf_t)({.1f, -2.f, 3.f});
  TEST_MATH_FUNCTION(cosh)({.1f, -2.f, 3.f});
  TEST_MATH_FUNCTION(cosh)({-4., .5, -.6});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(cosh)({.7l, .8l, .9l});
#endif

  TEST_MATH_FUNCTION(tanh)({-3, -2, -1, 0, 1});
  TEST_MATH_FUNCTION(tanh)({-3l, -2l, -1l, 0l, 1l});
  TEST_MATH_FUNCTION(tanh)({-3ll, -2ll, -1ll, 0ll, 1ll});
  TEST_MATH_FUNCTION(tanh)({2u, 3u, 4u, 5u, 6u});
  TEST_MATH_FUNCTION(tanh)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_MATH_FUNCTION(tanh)({2ull, 3ull, 4ull, 5ull, 6ull});
  TEST_HALF_MATH_FUNCTION(tanh, KE::half_t)({.1f, -2.f, 3.f});
  TEST_HALF_MATH_FUNCTION(tanh, KE::bhalf_t)({.1f, -2.f, 3.f});
  TEST_MATH_FUNCTION(tanh)({.1f, -2.f, 3.f});
  TEST_MATH_FUNCTION(tanh)({-4., .5, -.6});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(tanh)({.7l, .8l, .9l});
#endif

  TEST_MATH_FUNCTION(asinh)({-3, -2, -1, 0, 1});
  TEST_MATH_FUNCTION(asinh)({-3l, -2l, -1l, 0l, 1l});
  TEST_MATH_FUNCTION(asinh)({-3ll, -2ll, -1ll, 0ll, 1ll});
  TEST_MATH_FUNCTION(asinh)({2u, 3u, 4u, 5u, 6u});
  TEST_MATH_FUNCTION(asinh)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_MATH_FUNCTION(asinh)({2ull, 3ull, 4ull, 5ull, 6ull});
  TEST_HALF_MATH_FUNCTION(asinh, KE::half_t)({.1f, -2.f, 3.f});
  TEST_HALF_MATH_FUNCTION(asinh, KE::bhalf_t)({.1f, -2.f, 3.f});
  TEST_MATH_FUNCTION(asinh)({.1f, -2.f, 3.f});
  TEST_MATH_FUNCTION(asinh)({-4., .5, -.6});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(asinh)({.7l, .8l, .9l});
#endif

  TEST_MATH_FUNCTION(acosh)({1, 2, 3, 4, 5, 6});
  TEST_MATH_FUNCTION(acosh)({1l, 2l, 3l, 4l, 5l, 6l});
  TEST_MATH_FUNCTION(acosh)({1ll, 2ll, 3ll, 4ll, 5ll, 6ll});
  TEST_MATH_FUNCTION(acosh)({1u, 2u, 3u, 4u, 5u, 6u});
  TEST_MATH_FUNCTION(acosh)({1ul, 2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_MATH_FUNCTION(acosh)({1ull, 2ull, 3ull, 4ull, 5ull, 6ull});
  TEST_HALF_MATH_FUNCTION(acosh, KE::half_t)({1.2f, 34.f, 56.f, 789.f});
  TEST_HALF_MATH_FUNCTION(acosh, KE::bhalf_t)({1.2f, 34.f, 56.f, 789.f});
  TEST_MATH_FUNCTION(acosh)({1.2f, 34.f, 56.f, 789.f});
  TEST_MATH_FUNCTION(acosh)({1.2, 34., 56., 789.});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(acosh)({1.2l, 34.l, 56.l, 789.l});
#endif

  TEST_MATH_FUNCTION(atanh)({0});
  TEST_MATH_FUNCTION(atanh)({0l});
  TEST_MATH_FUNCTION(atanh)({0ll});
  TEST_MATH_FUNCTION(atanh)({0u});
  TEST_MATH_FUNCTION(atanh)({0ul});
  TEST_MATH_FUNCTION(atanh)({0ull});
  TEST_HALF_MATH_FUNCTION(atanh, KE::half_t)
  ({-.97f, .86f, -.53f, .42f, -.1f, 0.f});
  TEST_HALF_MATH_FUNCTION(atanh, KE::bhalf_t)
  ({-.97f, .86f, -.53f, .42f, -.1f, 0.f});
  TEST_MATH_FUNCTION(atanh)({-.97f, .86f, -.53f, .42f, -.1f, 0.f});
  TEST_MATH_FUNCTION(atanh)({-.97, .86, -.53, .42, -.1, 0.});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(atanh)({-.97l, .86l, -.53l, .42l, -.1l, 0.l});
#endif
}

TEST(TEST_CATEGORY, mathematical_functions_non_standard) {
  TEST_MATH_FUNCTION(rsqrt)({1, 2, 3, 5, 7, 11});
  TEST_MATH_FUNCTION(rsqrt)({1l, 2l, 3l, 5l, 7l, 11l});
  TEST_MATH_FUNCTION(rsqrt)({1ll, 2ll, 3ll, 5ll, 7ll, 11ll});
  TEST_MATH_FUNCTION(rsqrt)({1u, 2u, 3u, 5u, 7u});
  TEST_MATH_FUNCTION(rsqrt)({1ul, 2ul, 3ul, 5ul, 7ul});
  TEST_MATH_FUNCTION(rsqrt)({1ull, 2ull, 3ull, 5ull, 7ull});
  TEST_MATH_FUNCTION(rsqrt)({10.f, 20.f, 30.f, 40.f});
  TEST_MATH_FUNCTION(rsqrt)({11.1, 22.2, 33.3, 44.4});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(rsqrt)({10.l, 20.l, 30.l, 40.l});
#endif
}
#endif

#ifndef KOKKOS_MATHEMATICAL_FUNCTIONS_SKIP_2

TEST(TEST_CATEGORY, mathematical_functions_error_and_gamma_functions) {
  TEST_MATH_FUNCTION(erf)({-3, -2, -1, 0, 1});
  TEST_MATH_FUNCTION(erf)({-3l, -2l, -1l, 0l, 1l});
  TEST_MATH_FUNCTION(erf)({-3ll, -2ll, -1ll, 0ll, 1ll});
  TEST_MATH_FUNCTION(erf)({2u, 3u, 4u, 5u, 6u});
  TEST_MATH_FUNCTION(erf)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_MATH_FUNCTION(erf)({2ull, 3ull, 4ull, 5ull, 6ull});
  TEST_HALF_MATH_FUNCTION(erf, KE::half_t)({.1f, -2.f, 3.f});
  TEST_HALF_MATH_FUNCTION(erf, KE::bhalf_t)({.1f, -2.f, 3.f});
  TEST_MATH_FUNCTION(erf)({.1f, -2.f, 3.f});
  TEST_MATH_FUNCTION(erf)({-4., .5, -.6});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(erf)({.7l, .8l, .9l});
#endif

  TEST_MATH_FUNCTION(erfc)({-3, -2, -1, 0, 1});
  TEST_MATH_FUNCTION(erfc)({-3l, -2l, -1l, 0l, 1l});
  TEST_MATH_FUNCTION(erfc)({-3ll, -2ll, -1ll, 0ll, 1ll});
  TEST_MATH_FUNCTION(erfc)({2u, 3u, 4u, 5u, 6u});
  TEST_MATH_FUNCTION(erfc)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_MATH_FUNCTION(erfc)({2ull, 3ull, 4ull, 5ull, 6ull});
  TEST_HALF_MATH_FUNCTION(erfc, KE::half_t)({.1f, -2.f, 3.f});
  TEST_HALF_MATH_FUNCTION(erfc, KE::bhalf_t)({.1f, -2.f, 3.f});
  TEST_MATH_FUNCTION(erfc)({.1f, -2.f, 3.f});
  TEST_MATH_FUNCTION(erfc)({-4., .5, -.6});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(erfc)({.7l, .8l, .9l});
#endif

  TEST_MATH_FUNCTION(tgamma)({1, 2, 3, 4, 56, 78});
  TEST_MATH_FUNCTION(tgamma)({1l, 2l, 3l, 4l, 56l, 78l});
  TEST_MATH_FUNCTION(tgamma)({1ll, 2ll, 3ll, 4ll, 56ll, 78ll});
  TEST_MATH_FUNCTION(tgamma)({1u, 2u, 3u, 4u, 56u, 78u});
  TEST_MATH_FUNCTION(tgamma)({1ul, 2ul, 3ul, 4ul, 56ul, 78ul});
  TEST_MATH_FUNCTION(tgamma)({1ull, 2ull, 3ull, 4ull, 56ull, 78ull});
  TEST_HALF_MATH_FUNCTION(tgamma, KE::half_t)({.1f, -2.2f, 3.f});
  TEST_HALF_MATH_FUNCTION(tgamma, KE::bhalf_t)({.1f, -2.2f, 3.f});
  TEST_MATH_FUNCTION(tgamma)({.1f, -2.2f, 3.f});
  TEST_MATH_FUNCTION(tgamma)({-4.4, .5, -.6});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(tgamma)({.7l, .8l, .9l});
#endif

  TEST_MATH_FUNCTION(lgamma)({1, 2, 3, 4, 56, 78});
  TEST_MATH_FUNCTION(lgamma)({1l, 2l, 3l, 4l, 56l, 78l});
  TEST_MATH_FUNCTION(lgamma)({1ll, 2ll, 3ll, 4ll, 56ll, 78ll});
  TEST_MATH_FUNCTION(lgamma)({1u, 2u, 3u, 4u, 56u, 78u});
  TEST_MATH_FUNCTION(lgamma)({1ul, 2ul, 3ul, 4ul, 56ul, 78ul});
  TEST_MATH_FUNCTION(lgamma)({1ull, 2ull, 3ull, 4ull, 56ull, 78ull});
  TEST_HALF_MATH_FUNCTION(lgamma, KE::half_t)({.1f, -2.2f, 3.f});
  TEST_HALF_MATH_FUNCTION(lgamma, KE::bhalf_t)({.1f, -2.2f, 3.f});
  TEST_MATH_FUNCTION(lgamma)({.1f, -2.2f, 3.f});
  TEST_MATH_FUNCTION(lgamma)({-4.4, .5, -.6});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(lgamma)({.7l, .8l, .9l});
#endif
}

TEST(TEST_CATEGORY,
     mathematical_functions_nearest_interger_floating_point_operations) {
  TEST_MATH_FUNCTION(ceil)({-3, -2, -1, 0, 1});
  TEST_MATH_FUNCTION(ceil)({-3l, -2l, -1l, 0l, 1l});
  TEST_MATH_FUNCTION(ceil)({-3ll, -2ll, -1ll, 0ll, 1ll});
  TEST_MATH_FUNCTION(ceil)({2u, 3u, 4u, 5u, 6u});
  TEST_MATH_FUNCTION(ceil)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_MATH_FUNCTION(ceil)({2ull, 3ull, 4ull, 5ull, 6ull});
  TEST_HALF_MATH_FUNCTION(ceil, KE::half_t)({-1.1f, 2.2f, -3.3f, 4.4f, -5.5f});
  TEST_HALF_MATH_FUNCTION(ceil, KE::bhalf_t)({-1.1f, 2.2f, -3.3f, 4.4f, -5.5f});
  TEST_MATH_FUNCTION(ceil)({-1.1f, 2.2f, -3.3f, 4.4f, -5.5f});
  TEST_MATH_FUNCTION(ceil)({-6.6, 7.7, -8.8, 9.9});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(ceil)({12.3l, 4.56l, 789.l});
#endif

  TEST_MATH_FUNCTION(floor)({-3, -2, -1, 0, 1});
  TEST_MATH_FUNCTION(floor)({-3l, -2l, -1l, 0l, 1l});
  TEST_MATH_FUNCTION(floor)({-3ll, -2ll, -1ll, 0ll, 1ll});
  TEST_MATH_FUNCTION(floor)({2u, 3u, 4u, 5u, 6u});
  TEST_MATH_FUNCTION(floor)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_MATH_FUNCTION(floor)({2ull, 3ull, 4ull, 5ull, 6ull});
  TEST_HALF_MATH_FUNCTION(floor, KE::half_t)({-1.1f, 2.2f, -3.3f, 4.4f, -5.5f});
  TEST_HALF_MATH_FUNCTION(floor, KE::bhalf_t)
  ({-1.1f, 2.2f, -3.3f, 4.4f, -5.5f});
  TEST_MATH_FUNCTION(floor)({-1.1f, 2.2f, -3.3f, 4.4f, -5.5f});
  TEST_MATH_FUNCTION(floor)({-6.6, 7.7, -8.8, 9.9});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(floor)({12.3l, 4.56l, 789.l});
#endif

  TEST_MATH_FUNCTION(trunc)({-3, -2, -1, 0, 1});
  TEST_MATH_FUNCTION(trunc)({-3l, -2l, -1l, 0l, 1l});
  TEST_MATH_FUNCTION(trunc)({-3ll, -2ll, -1ll, 0ll, 1ll});
  TEST_MATH_FUNCTION(trunc)({2u, 3u, 4u, 5u, 6u});
  TEST_MATH_FUNCTION(trunc)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_MATH_FUNCTION(trunc)({2ull, 3ull, 4ull, 5ull, 6ull});
  TEST_HALF_MATH_FUNCTION(trunc, KE::half_t)({-1.1f, 2.2f, -3.3f, 4.4f, -5.5f});
  TEST_HALF_MATH_FUNCTION(trunc, KE::bhalf_t)
  ({-1.1f, 2.2f, -3.3f, 4.4f, -5.5f});
  TEST_MATH_FUNCTION(trunc)({-1.1f, 2.2f, -3.3f, 4.4f, -5.5f});
  TEST_MATH_FUNCTION(trunc)({-6.6, 7.7, -8.8, 9.9});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(trunc)({12.3l, 4.56l, 789.l});
#endif

  TEST_MATH_FUNCTION(round)({-3, -2, -1, 0, 1});
  TEST_MATH_FUNCTION(round)({-3l, -2l, -1l, 0l, 1l});
  TEST_MATH_FUNCTION(round)({-3ll, -2ll, -1ll, 0ll, 1ll});
  TEST_MATH_FUNCTION(round)({2u, 3u, 4u, 5u, 6u});
  TEST_MATH_FUNCTION(round)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_MATH_FUNCTION(round)({2ull, 3ull, 4ull, 5ull, 6ull});
  TEST_HALF_MATH_FUNCTION(round, KE::half_t)
  ({2.3f, 2.5f, 2.7f, -2.3f, -2.5f, -2.7f, -0.0f});
  TEST_HALF_MATH_FUNCTION(round, KE::bhalf_t)
  ({2.3f, 2.5f, 2.7f, -2.3f, -2.5f, -2.7f, -0.0f});
  TEST_MATH_FUNCTION(round)({2.3f, 2.5f, 2.7f, -2.3f, -2.5f, -2.7f, -0.0f});
  TEST_MATH_FUNCTION(round)({2.3, 2.5, 2.7, -2.3, -2.5, -2.7, -0.0});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(round)({2.3l, 2.5l, 2.7l, -2.3l, -2.5l, -2.7l, -0.0l});
#endif

#ifndef KOKKOS_ENABLE_SYCL
  TEST_MATH_FUNCTION(nearbyint)({-3, -2, -1, 0, 1});
  TEST_MATH_FUNCTION(nearbyint)({-3l, -2l, -1l, 0l, 1l});
  TEST_MATH_FUNCTION(nearbyint)({-3ll, -2ll, -1ll, 0ll, 1ll});
  TEST_MATH_FUNCTION(nearbyint)({2u, 3u, 4u, 5u, 6u});
  TEST_MATH_FUNCTION(nearbyint)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_MATH_FUNCTION(nearbyint)({2ull, 3ull, 4ull, 5ull, 6ull});
  TEST_HALF_MATH_FUNCTION(nearbyint, KE::half_t)
  ({-1.1f, 2.2f, -3.3f, 4.4f, -5.5f});
  TEST_HALF_MATH_FUNCTION(nearbyint, KE::bhalf_t)
  ({-1.1f, 2.2f, -3.3f, 4.4f, -5.5f});
  TEST_MATH_FUNCTION(nearbyint)({-1.1f, 2.2f, -3.3f, 4.4f, -5.5f});
  TEST_MATH_FUNCTION(nearbyint)({-6.6, 7.7, -8.8, 9.9});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(nearbyint)({12.3l, 4.56l, 789.l});
#endif
#endif
}

TEST(TEST_CATEGORY,
     mathematical_functions_floating_point_manipulation_functions) {
  TEST_MATH_FUNCTION(logb)({2, 3, 4, 56, 789});
  TEST_MATH_FUNCTION(logb)({2l, 3l, 4l, 56l, 789l});
  TEST_MATH_FUNCTION(logb)({2ll, 3ll, 4ll, 56ll, 789ll});
  TEST_MATH_FUNCTION(logb)({2u, 3u, 4u, 5u, 6u});
  TEST_MATH_FUNCTION(logb)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_MATH_FUNCTION(logb)({2ull, 3ull, 4ull, 5ull, 6ull});
  TEST_HALF_MATH_FUNCTION(logb, KE::half_t)({123.45f, 6789.0f});
  TEST_HALF_MATH_FUNCTION(logb, KE::bhalf_t)({123.45f, 6789.0f});
  TEST_MATH_FUNCTION(logb)({123.45f, 6789.0f});
  TEST_MATH_FUNCTION(logb)({123.45, 6789.0});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(logb)({123.45l, 6789.0l});
#endif

  do_test_math_binary_function<TEST_EXECSPACE, kk_nextafter>(
      0, static_cast<KE::half_t>(1.f));
  do_test_math_binary_function<TEST_EXECSPACE, kk_nextafter>(
      1, static_cast<KE::half_t>(2.f));
  do_test_math_binary_function<TEST_EXECSPACE, kk_nextafter>(
      0, static_cast<KE::bhalf_t>(1.f));
  do_test_math_binary_function<TEST_EXECSPACE, kk_nextafter>(
      1, static_cast<KE::bhalf_t>(2.f));
  do_test_math_binary_function<TEST_EXECSPACE, kk_nextafter>(0, 1.f);
  do_test_math_binary_function<TEST_EXECSPACE, kk_nextafter>(1, 2.f);
  do_test_math_binary_function<TEST_EXECSPACE, kk_nextafter>(0.1, 0);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  do_test_math_binary_function<TEST_EXECSPACE, kk_nextafter>(1, 2.l);
  do_test_math_binary_function<TEST_EXECSPACE, kk_nextafter>(1.l, 2.l);
#endif

  do_test_math_binary_function<TEST_EXECSPACE, kk_copysign>(
      0, static_cast<KE::half_t>(1.f));
  do_test_math_binary_function<TEST_EXECSPACE, kk_copysign>(
      1, static_cast<KE::half_t>(2.f));
  do_test_math_binary_function<TEST_EXECSPACE, kk_copysign>(
      0, static_cast<KE::bhalf_t>(1.f));
  do_test_math_binary_function<TEST_EXECSPACE, kk_copysign>(
      1, static_cast<KE::bhalf_t>(2.f));
  do_test_math_binary_function<TEST_EXECSPACE, kk_copysign>(0, 1.f);
  do_test_math_binary_function<TEST_EXECSPACE, kk_copysign>(1, 2.f);
  do_test_math_binary_function<TEST_EXECSPACE, kk_copysign>(0.1, 0);
  do_test_math_binary_function<TEST_EXECSPACE, kk_copysign>(1.f, +2.f);
  do_test_math_binary_function<TEST_EXECSPACE, kk_copysign>(1.f, -2.f);
  do_test_math_binary_function<TEST_EXECSPACE, kk_copysign>(1., +2.);
  do_test_math_binary_function<TEST_EXECSPACE, kk_copysign>(1., -2.);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  do_test_math_binary_function<TEST_EXECSPACE, kk_copysign>(1, +2.l);
  do_test_math_binary_function<TEST_EXECSPACE, kk_copysign>(1.l, +2);
  do_test_math_binary_function<TEST_EXECSPACE, kk_copysign>(1.l, +2.l);
  do_test_math_binary_function<TEST_EXECSPACE, kk_copysign>(1.l, -2.l);
#endif
}

template <class Space>
struct TestAbsoluteValueFunction {
  TestAbsoluteValueFunction() { run(); }
  void run() const {
    int errors = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Space>(0, 1), *this, errors);
    ASSERT_EQ(errors, 0);
  }
  KOKKOS_FUNCTION void operator()(int, int& e) const {
    using Kokkos::abs;
    if (abs(1) != 1 || abs(-1) != 1) {
      ++e;
      Kokkos::printf("failed abs(int)\n");
    }
    if (abs(2l) != 2l || abs(-2l) != 2l) {
      ++e;
      Kokkos::printf("failed abs(long int)\n");
    }
    if (abs(3ll) != 3ll || abs(-3ll) != 3ll) {
      ++e;
      Kokkos::printf("failed abs(long long int)\n");
    }
    if (abs(4.f) != 4.f || abs(-4.f) != 4.f) {
      ++e;
      Kokkos::printf("failed abs(float)\n");
    }
    if (abs(static_cast<KE::half_t>(4.f)) != static_cast<KE::half_t>(4.f) ||
        abs(static_cast<KE::half_t>(-4.f)) != static_cast<KE::half_t>(4.f)) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed abs(KE::half_t)\n");
    }
    if (abs(static_cast<KE::bhalf_t>(4.f)) != static_cast<KE::bhalf_t>(4.f) ||
        abs(static_cast<KE::bhalf_t>(-4.f)) != static_cast<KE::bhalf_t>(4.f)) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed abs(KE::bhalf_t)\n");
    }
    if (abs(5.) != 5. || abs(-5.) != 5.) {
      ++e;
      Kokkos::printf("failed abs(double)\n");
    }
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    if (abs(6.l) != 6.l || abs(-6.l) != 6.l) {
      ++e;
      Kokkos::printf("failed abs(long double)\n");
    }
#endif
    // special values
    using Kokkos::isinf;
    using Kokkos::isnan;
    if (abs(-0.) != 0. || !isinf(abs(-INFINITY)) || !isnan(abs(-NAN))) {
      ++e;
      Kokkos::printf("failed abs(floating_point) special values\n");
    }

    static_assert(std::is_same<decltype(abs(1)), int>::value, "");
    static_assert(std::is_same<decltype(abs(2l)), long>::value, "");
    static_assert(std::is_same<decltype(abs(3ll)), long long>::value, "");
    static_assert(std::is_same<decltype(abs(static_cast<KE::half_t>(4.f))),
                               KE::half_t>::value,
                  "");
    static_assert(std::is_same<decltype(abs(static_cast<KE::bhalf_t>(4.f))),
                               KE::bhalf_t>::value,
                  "");
    static_assert(std::is_same<decltype(abs(4.f)), float>::value, "");
    static_assert(std::is_same<decltype(abs(5.)), double>::value, "");
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    static_assert(std::is_same<decltype(abs(6.l)), long double>::value, "");
#endif
  }
};

TEST(TEST_CATEGORY, mathematical_functions_absolute_value) {
  TestAbsoluteValueFunction<TEST_EXECSPACE>();
}

template <class Space>
struct TestFloatingPointAbsoluteValueFunction {
  TestFloatingPointAbsoluteValueFunction() { run(); }
  void run() const {
    int errors = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Space>(0, 1), *this, errors);
    ASSERT_EQ(errors, 0);
  }
  KOKKOS_FUNCTION void operator()(int, int& e) const {
    using Kokkos::fabs;
    if (fabs(4.f) != 4.f || fabs(-4.f) != 4.f) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed fabs(float)\n");
    }
    if (fabs(static_cast<KE::half_t>(4.f)) != static_cast<KE::half_t>(4.f) ||
        fabs(static_cast<KE::half_t>(-4.f)) != static_cast<KE::half_t>(4.f)) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed fabs(KE::half_t)\n");
    }
    if (fabs(static_cast<KE::bhalf_t>(4.f)) != static_cast<KE::bhalf_t>(4.f) ||
        fabs(static_cast<KE::bhalf_t>(-4.f)) != static_cast<KE::bhalf_t>(4.f)) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed fabs(KE::bhalf_t)\n");
    }
    if (fabs(5.) != 5. || fabs(-5.) != 5.) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed fabs(double)\n");
    }
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    if (fabs(6.l) != 6.l || fabs(-6.l) != 6.l) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed fabs(long double)\n");
    }
#endif
    // special values
    using Kokkos::isinf;
    using Kokkos::isnan;
    if (fabs(-0.) != 0. || !isinf(fabs(-INFINITY)) || !isnan(fabs(-NAN))) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF(
          "failed fabs(floating_point) special values\n");
    }

    static_assert(std::is_same<decltype(fabs(static_cast<KE::half_t>(4.f))),
                               KE::half_t>::value);
    static_assert(std::is_same<decltype(fabs(static_cast<KE::bhalf_t>(4.f))),
                               KE::bhalf_t>::value);
    static_assert(std::is_same<decltype(fabs(4.f)), float>::value);
    static_assert(std::is_same<decltype(fabs(5.)), double>::value);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    static_assert(std::is_same<decltype(fabs(6.l)), long double>::value);
#endif
  }
};

TEST(TEST_CATEGORY, mathematical_functions_floating_point_absolute_value) {
  TestFloatingPointAbsoluteValueFunction<TEST_EXECSPACE>();
}

template <class Space>
struct TestFloatingPointRemainderFunction : FloatingPointComparison {
  TestFloatingPointRemainderFunction() { run(); }
  void run() const {
    int errors = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Space>(0, 1), *this, errors);
    ASSERT_EQ(errors, 0);
  }
  KOKKOS_FUNCTION void operator()(int, int& e) const {
    using Kokkos::fmod;
    if (!compare(fmod(6.2f, 4.f), 2.2f, 1) &&
        !compare(fmod(-6.2f, 4.f), -2.2f, 1)) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed fmod(float)\n");
    }
    if (!compare(
            fmod(static_cast<KE::half_t>(6.2f), static_cast<KE::half_t>(4.f)),
            static_cast<KE::half_t>(2.2f), 1) &&
        !compare(
            fmod(static_cast<KE::half_t>(-6.2f), static_cast<KE::half_t>(4.f)),
            -static_cast<KE::half_t>(2.2f), 1)) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed fmod(KE::half_t)\n");
    }
    if (!compare(
            fmod(static_cast<KE::bhalf_t>(6.2f), static_cast<KE::bhalf_t>(4.f)),
            static_cast<KE::bhalf_t>(2.2f), 1) &&
        !compare(fmod(static_cast<KE::bhalf_t>(-6.2f),
                      static_cast<KE::bhalf_t>(4.f)),
                 -static_cast<KE::bhalf_t>(2.2f), 1)) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed fmod(KE::bhalf_t)\n");
    }
    if (!compare(fmod(6.2, 4.), 2.2, 1) && !compare(fmod(-6.2, 4.), -2.2, 1)) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed fmod(double)\n");
    }
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    if (!compare(fmod(6.2l, 4.l), 2.2l, 1) &&
        !compare(fmod(-6.2l, 4.l), -2.2l, 1)) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed fmod(long double)\n");
    }
#endif
    // special values
    using Kokkos::isinf;
    using Kokkos::isnan;
    if (!isinf(fmod(-KE::infinity<float>::value, 1.f)) &&
        !isnan(fmod(-KE::quiet_NaN<float>::value, 1.f))) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF(
          "failed fmod(floating_point) special values\n");
    }

    static_assert(std::is_same<decltype(fmod(static_cast<KE::half_t>(4.f),
                                             static_cast<KE::half_t>(4.f))),
                               KE::half_t>::value,
                  "");
    static_assert(std::is_same<decltype(fmod(static_cast<KE::bhalf_t>(4.f),
                                             static_cast<KE::bhalf_t>(4.f))),
                               KE::bhalf_t>::value,
                  "");
    static_assert(std::is_same<decltype(fmod(4.f, 4.f)), float>::value, "");
    static_assert(std::is_same<decltype(fmod(5., 5.)), double>::value, "");
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    static_assert(std::is_same<decltype(fmod(6.l, 6.l)), long double>::value,
                  "");
#endif
  }
};

TEST(TEST_CATEGORY, mathematical_functions_remainder_function) {
  TestFloatingPointRemainderFunction<TEST_EXECSPACE>();
}

#if 0
// TODO: Adjust expected values, see https://github.com/kokkos/kokkos/issues/6275
template <class Space>
struct TestIEEEFloatingPointRemainderFunction : FloatingPointComparison {
  TestIEEEFloatingPointRemainderFunction() { run(); }
  void run() const {
    int errors = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Space>(0, 1), *this, errors);
    ASSERT_EQ(errors, 0);
  }
  KOKKOS_FUNCTION void operator()(int, int& e) const {
    using Kokkos::remainder;
    if (!compare(remainder(6.2f, 4.f), 2.2f, 2) &&
        !compare(remainder(-6.2f, 4.f), 2.2f, 1)) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed remainder(float)\n");
    }
    if (!compare(remainder(static_cast<KE::half_t>(6.2f),
                           static_cast<KE::half_t>(4.f)),
                 static_cast<KE::half_t>(2.2f), 1) &&
        !compare(remainder(static_cast<KE::half_t>(-6.2f),
                           static_cast<KE::half_t>(4.f)),
                 -static_cast<KE::half_t>(2.2f), 1)) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed remainder(KE::half_t)\n");
    }
    if (!compare(remainder(static_cast<KE::bhalf_t>(6.2f),
                           static_cast<KE::bhalf_t>(4.f)),
                 static_cast<KE::bhalf_t>(2.2f), 1) &&
        !compare(remainder(static_cast<KE::bhalf_t>(-6.2f),
                           static_cast<KE::bhalf_t>(4.f)),
                 -static_cast<KE::bhalf_t>(2.2f), 1)) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed remainder(KE::bhalf_t)\n");
    }
    if (!compare(remainder(6.2, 4.), 2.2, 2) &&
        !compare(remainder(-6.2, 4.), 2.2, 1)) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed remainder(double)\n");
    }
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    if (!compare(remainder(6.2l, 4.l), 2.2l, 1) &&
        !compare(remainder(-6.2l, 4.l), -2.2l, 1)) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed remainder(long double)\n");
    }
#endif
    // special values
    using Kokkos::isinf;
    using Kokkos::isnan;
    if (!isinf(remainder(-KE::infinity<float>::value, 1.f)) &&
        !isnan(remainder(-KE::quiet_NaN<float>::value, 1.f))) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF(
          "failed remainder(floating_point) special values\n");
    }

    static_assert(
        std::is_same<decltype(remainder(static_cast<KE::half_t>(4.f),
                                        static_cast<KE::half_t>(4.f))),
                     KE::half_t>::value,
        "");
    static_assert(
        std::is_same<decltype(remainder(static_cast<KE::bhalf_t>(4.f),
                                        static_cast<KE::bhalf_t>(4.f))),
                     KE::bhalf_t>::value,
        "");
    static_assert(std::is_same<decltype(remainder(4.f, 4.f)), float>::value,
                  "");
    static_assert(std::is_same<decltype(remainder(5., 5.)), double>::value, "");
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    static_assert(
        std::is_same<decltype(remainder(6.l, 6.l)), long double>::value, "");
#endif
  }
};

TEST(TEST_CATEGORY, mathematical_functions_ieee_remainder_function) {
  TestIEEEFloatingPointRemainderFunction<TEST_EXECSPACE>();
}
#endif

// TODO: TestFpClassify, see https://github.com/kokkos/kokkos/issues/6279

#ifndef KOKKOS_MATHEMATICAL_FUNCTIONS_SKIP_2
template <class Space>
struct TestIsFinite {
  TestIsFinite() { run(); }
  void run() const {
    int errors = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Space>(0, 1), *this, errors);
    ASSERT_EQ(errors, 0);
  }
  KOKKOS_FUNCTION void operator()(int, int& e) const {
    using KE::infinity;
    using KE::quiet_NaN;
    using KE::signaling_NaN;
    using Kokkos::isfinite;
    if (!isfinite(1) || !isfinite(INT_MAX)) {
      ++e;
      Kokkos::printf("failed isfinite(integral)\n");
    }
    if (!isfinite(2.f) || isfinite(quiet_NaN<float>::value) ||
        isfinite(signaling_NaN<float>::value) ||
        isfinite(infinity<float>::value)) {
      ++e;
      Kokkos::printf("failed isfinite(float)\n");
    }
#if !(defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_COMPILER_MSVC))
    if (!isfinite(static_cast<KE::half_t>(2.f))
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC 23.7
        || isfinite(quiet_NaN<KE::half_t>::value) ||
        isfinite(signaling_NaN<KE::half_t>::value) ||
        isfinite(infinity<KE::half_t>::value)
#endif
    ) {
      ++e;
      Kokkos::printf("failed isfinite(KE::half_t)\n");
    }
    if (!isfinite(static_cast<KE::bhalf_t>(2.f))
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC 23.7
        || isfinite(quiet_NaN<KE::bhalf_t>::value) ||
        isfinite(signaling_NaN<KE::bhalf_t>::value) ||
        isfinite(infinity<KE::bhalf_t>::value)
#endif
    ) {
      ++e;
      Kokkos::printf("failed isfinite(KE::bhalf_t)\n");
    }
#endif
    if (!isfinite(3.)
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC 23.7
        || isfinite(quiet_NaN<double>::value) ||
        isfinite(signaling_NaN<double>::value) ||
        isfinite(infinity<double>::value)
#endif
    ) {
      ++e;
      Kokkos::printf("failed isfinite(double)\n");
    }
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    if (!isfinite(4.l) || isfinite(quiet_NaN<long double>::value) ||
        isfinite(signaling_NaN<long double>::value) ||
        isfinite(infinity<long double>::value)) {
      ++e;
      Kokkos::printf("failed isfinite(long double)\n");
    }
#endif
    // special values
    if (isfinite(INFINITY) || isfinite(NAN)) {
      ++e;
      Kokkos::printf("failed isfinite(floating_point) special values\n");
    }

    static_assert(std::is_same<decltype(isfinite(1)), bool>::value);
    static_assert(std::is_same<decltype(isfinite(2.f)), bool>::value);
    static_assert(std::is_same<decltype(isfinite(3.)), bool>::value);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    static_assert(std::is_same<decltype(isfinite(4.l)), bool>::value);
#endif
  }
};

TEST(TEST_CATEGORY, mathematical_functions_isfinite) {
  TestIsFinite<TEST_EXECSPACE>();
}

template <class Space>
struct TestIsInf {
  TestIsInf() { run(); }
  void run() const {
    int errors = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Space>(0, 1), *this, errors);
    ASSERT_EQ(errors, 0);
  }
  KOKKOS_FUNCTION void operator()(int, int& e) const {
    using KE::infinity;
    using KE::quiet_NaN;
    using KE::signaling_NaN;
    using Kokkos::isinf;
    if (isinf(1) || isinf(INT_MAX)) {
      ++e;
      Kokkos::printf("failed isinf(integral)\n");
    }
    if (isinf(2.f) || isinf(quiet_NaN<float>::value) ||
        isinf(signaling_NaN<float>::value) || !isinf(infinity<float>::value)) {
      ++e;
      Kokkos::printf("failed isinf(float)\n");
    }
#if !(defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_COMPILER_MSVC))
    if (isinf(static_cast<KE::half_t>(2.f))
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC 23.7
        || isinf(quiet_NaN<KE::half_t>::value) ||
        isinf(signaling_NaN<KE::half_t>::value) ||
        !isinf(infinity<KE::half_t>::value)
#endif
    ) {
      ++e;
      Kokkos::printf("failed isinf(KE::half_t)\n");
    }
    if (isinf(static_cast<KE::bhalf_t>(2.f))
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC 23.7
        || isinf(quiet_NaN<KE::bhalf_t>::value) ||
        isinf(signaling_NaN<KE::bhalf_t>::value) ||
        !isinf(infinity<KE::bhalf_t>::value)
#endif
    ) {
      ++e;
      Kokkos::printf("failed isinf(KE::bhalf_t)\n");
    }
#endif
    if (isinf(3.)
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC 23.7
        || isinf(quiet_NaN<double>::value) ||
        isinf(signaling_NaN<double>::value) || !isinf(infinity<double>::value)
#endif
    ) {
      ++e;
      Kokkos::printf("failed isinf(double)\n");
    }
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    if (isinf(4.l) || isinf(quiet_NaN<long double>::value) ||
        isinf(signaling_NaN<long double>::value) ||
        !isinf(infinity<long double>::value)) {
      ++e;
      Kokkos::printf("failed isinf(long double)\n");
    }
#endif
    // special values
    if (!isinf(INFINITY) || isinf(NAN)) {
      ++e;
      Kokkos::printf("failed isinf(floating_point) special values\n");
    }

    static_assert(std::is_same<decltype(isinf(1)), bool>::value);
    static_assert(std::is_same<decltype(isinf(2.f)), bool>::value);
    static_assert(std::is_same<decltype(isinf(3.)), bool>::value);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    static_assert(std::is_same<decltype(isinf(4.l)), bool>::value);
#endif
  }
};

TEST(TEST_CATEGORY, mathematical_functions_isinf) {
  TestIsInf<TEST_EXECSPACE>();
}

template <class Space>
struct TestIsNaN {
  TestIsNaN() { run(); }
  void run() const {
    int errors = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Space>(0, 1), *this, errors);
    ASSERT_EQ(errors, 0);
  }
  KOKKOS_FUNCTION void operator()(int, int& e) const {
    using KE::infinity;
    using KE::quiet_NaN;
    using KE::signaling_NaN;
    using Kokkos::isnan;
    if (isnan(1) || isnan(INT_MAX)) {
      ++e;
      Kokkos::printf("failed isnan(integral)\n");
    }
    if (isnan(2.f) || !isnan(quiet_NaN<float>::value) ||
        !isnan(signaling_NaN<float>::value) || isnan(infinity<float>::value)) {
      ++e;
      Kokkos::printf("failed isnan(float)\n");
    }
#if !(defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_COMPILER_MSVC))
    if (isnan(static_cast<KE::half_t>(2.f))
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC 23.7
        || !isnan(quiet_NaN<KE::half_t>::value) ||
        !isnan(signaling_NaN<KE::half_t>::value) ||
        isnan(infinity<KE::half_t>::value)
#endif
    ) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed isnan(KE::half_t)\n");
    }
    if (isnan(static_cast<KE::bhalf_t>(2.f))
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC 23.7
        || !isnan(quiet_NaN<KE::bhalf_t>::value) ||
        !isnan(signaling_NaN<KE::bhalf_t>::value) ||
        isnan(infinity<KE::bhalf_t>::value)
#endif
    ) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed isnan(KE::bhalf_t)\n");
    }
    if (isnan(3.)
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC 23.7
        || !isnan(quiet_NaN<double>::value) ||
        !isnan(signaling_NaN<double>::value) || isnan(infinity<double>::value)
#endif
    ) {
      ++e;
      Kokkos::printf("failed isnan(double)\n");
    }
#endif
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    if (isnan(4.l) || !isnan(quiet_NaN<long double>::value) ||
        !isnan(signaling_NaN<long double>::value) ||
        isnan(infinity<long double>::value)) {
      ++e;
      Kokkos::printf("failed isnan(long double)\n");
    }
#endif
    // special values
    if (isnan(INFINITY) || !isnan(NAN)) {
      ++e;
      Kokkos::printf("failed isnan(floating_point) special values\n");
    }

    static_assert(std::is_same<decltype(isnan(1)), bool>::value, "");
    static_assert(std::is_same<decltype(isnan(2.f)), bool>::value, "");
    static_assert(std::is_same<decltype(isnan(3.)), bool>::value, "");
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    static_assert(std::is_same<decltype(isnan(4.l)), bool>::value, "");
#endif
  }
};

TEST(TEST_CATEGORY, mathematical_functions_isnan) {
  TestIsNaN<TEST_EXECSPACE>();
}
#endif

// TODO: TestSignBit, see https://github.com/kokkos/kokkos/issues/6279
#endif
