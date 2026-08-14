// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include <impl/Kokkos_Half_FloatingPointWrapper.hpp>
#include <Kokkos_TypeInfo.hpp>

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <type_traits>
#include <cstdint>
#include <cfloat>

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || \
    defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ENABLE_OPENACC)
#else
#define MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
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
  KOKKOS_FUNCTION std::enable_if_t<std::is_signed_v<T>, T> absolute(
      T val) const {
    return val < T(0) ? -val : val;
  }

  template <class T>
  KOKKOS_FUNCTION std::enable_if_t<!std::is_signed_v<T>, T> absolute(
      T val) const {
    return val;
  }

 public:
  template <class FPT>
  KOKKOS_FUNCTION bool compare_near_zero(FPT const& fpv, int ulp) const {
    auto abs_tol = eps(fpv) * ulp;

    bool ar = absolute(fpv) <= abs_tol;
    if (!ar) {
      Kokkos::printf("absolute value exceeds tolerance [|%e| > %e]\n",
                     (double)fpv, (double)abs_tol);
    }

    return ar;
  }

  template <class Lhs, class Rhs>
  KOKKOS_FUNCTION bool compare(Lhs const& lhs, Rhs const& rhs, int ulp) const {
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
      bool ar         = rel_diff <= rel_tol;
      if (!ar) {
        Kokkos::printf("relative difference exceeds tolerance [%e > %e]\n",
                       (double)rel_diff, (double)rel_tol);
      }

      return ar;
    }
  }
};

struct IntegerComparison {
  template <class Lhs, class Rhs>
  KOKKOS_FUNCTION bool compare(Lhs const& lhs, Rhs const& rhs) const {
    static_assert(std::is_integral_v<Lhs>);
    static_assert(std::is_integral_v<Rhs>);
    return lhs == rhs;
  }
};

template <class Floating>
struct ConvertibleTo {
  operator Floating() const;
};

#define KOKKOS_TEST_STATIC_ASSERT_UNARY_PREDICATE(FUNC, FP_TYPE, RET_TYPE)   \
  static_assert(                                                             \
      std::is_same_v<decltype(FUNC(std::declval<ConvertibleTo<FP_TYPE>>())), \
                     RET_TYPE>)

template <class>
struct math_function_name;

#define DEFINE_UNARY_FUNCTION_EVAL(FUNC, ULP_FACTOR)                         \
  struct MathUnaryFunction_##FUNC {                                          \
    template <typename T>                                                    \
    static KOKKOS_FUNCTION auto eval(T x) {                                  \
      static_assert(std::is_same_v<decltype(Kokkos::FUNC((T)0)),             \
                                   math_unary_function_return_type_t<T>>);   \
      if constexpr (std::is_floating_point_v<T>) {                           \
        KOKKOS_TEST_STATIC_ASSERT_UNARY_PREDICATE(                           \
            Kokkos::FUNC, T, math_unary_function_return_type_t<T>);          \
      }                                                                      \
      return Kokkos::FUNC(x);                                                \
    }                                                                        \
    template <typename T>                                                    \
    static auto eval_std(T x) {                                              \
      if constexpr (std::is_same_v<T, KE::half_t> ||                         \
                    std::is_same_v<T, KE::bhalf_t>) {                        \
        return std::FUNC(static_cast<float>(x));                             \
      } else {                                                               \
        static_assert(std::is_same_v<decltype(std::FUNC((T)0)),              \
                                     math_unary_function_return_type_t<T>>); \
        if constexpr (std::is_floating_point_v<T>) {                         \
          KOKKOS_TEST_STATIC_ASSERT_UNARY_PREDICATE(                         \
              std::FUNC, T, math_unary_function_return_type_t<T>);           \
        }                                                                    \
        return std::FUNC(x);                                                 \
      }                                                                      \
    }                                                                        \
    static KOKKOS_FUNCTION int ulp_factor() { return ULP_FACTOR; }           \
  };                                                                         \
  using kk_##FUNC = MathUnaryFunction_##FUNC;                                \
  template <>                                                                \
  struct math_function_name<MathUnaryFunction_##FUNC> {                      \
    static constexpr char name[] = #FUNC;                                    \
  };                                                                         \
  constexpr char math_function_name<MathUnaryFunction_##FUNC>::name[]

#define DEFINE_UNARY_FUNCTION_EVAL_CUSTOM(FUNC, ULP_FACTOR, REF_FUNC)        \
  struct MathUnaryFunction_##FUNC {                                          \
    template <typename T>                                                    \
    static KOKKOS_FUNCTION auto eval(T x) {                                  \
      static_assert(std::is_same_v<decltype(Kokkos::FUNC((T)0)),             \
                                   math_unary_function_return_type_t<T>>);   \
      if constexpr (std::is_floating_point_v<T>) {                           \
        KOKKOS_TEST_STATIC_ASSERT_UNARY_PREDICATE(                           \
            Kokkos::FUNC, T, math_unary_function_return_type_t<T>);          \
      }                                                                      \
      return Kokkos::FUNC(x);                                                \
    }                                                                        \
    template <typename T>                                                    \
    static auto eval_std(T y) {                                              \
      if constexpr (std::is_same_v<T, KE::half_t> ||                         \
                    std::is_same_v<T, KE::bhalf_t>) {                        \
        auto x = static_cast<float>(y);                                      \
        return static_cast<T>(REF_FUNC);                                     \
      } else {                                                               \
        const T x = y;                                                       \
        static_assert(std::is_same_v<decltype(REF_FUNC),                     \
                                     math_unary_function_return_type_t<T>>); \
        return REF_FUNC;                                                     \
      }                                                                      \
    }                                                                        \
    static KOKKOS_FUNCTION int ulp_factor() { return ULP_FACTOR; }           \
  };                                                                         \
  using kk_##FUNC = MathUnaryFunction_##FUNC;                                \
  template <>                                                                \
  struct math_function_name<MathUnaryFunction_##FUNC> {                      \
    static constexpr char name[] = #FUNC;                                    \
  };                                                                         \
  constexpr char math_function_name<MathUnaryFunction_##FUNC>::name[]

#define DEFINE_UNARY_FUNCTION_EVAL_INT(FUNC)                           \
  struct MathUnaryFunction_##FUNC {                                    \
    template <typename T>                                              \
    static KOKKOS_FUNCTION auto eval(T x) {                            \
      static_assert(std::is_integral_v<decltype(Kokkos::FUNC((T)0))>); \
      if constexpr (std::is_floating_point_v<T>) {                     \
        static_assert(std::is_integral_v<decltype(Kokkos::FUNC(        \
                          std::declval<ConvertibleTo<T>>()))>);        \
      }                                                                \
      return Kokkos::FUNC(x);                                          \
    }                                                                  \
    template <typename T>                                              \
    static auto eval_std(T x) {                                        \
      if constexpr (std::is_same_v<T, KE::half_t> ||                   \
                    std::is_same_v<T, KE::bhalf_t>) {                  \
        return std::FUNC(static_cast<float>(x));                       \
      } else {                                                         \
        static_assert(std::is_integral_v<decltype(std::FUNC((T)0))>);  \
        if constexpr (std::is_floating_point_v<T>) {                   \
          static_assert(std::is_integral_v<decltype(std::FUNC(         \
                            std::declval<ConvertibleTo<T>>()))>);      \
        }                                                              \
        return std::FUNC(x);                                           \
      }                                                                \
    }                                                                  \
  };                                                                   \
  using kk_##FUNC = MathUnaryFunction_##FUNC;                          \
  template <>                                                          \
  struct math_function_name<MathUnaryFunction_##FUNC> {                \
    static constexpr char name[] = #FUNC;                              \
  };                                                                   \
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
DEFINE_UNARY_FUNCTION_EVAL_CUSTOM(rcp, 2,
                                  math_unary_function_return_type_t<T>(1) / x);
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
DEFINE_UNARY_FUNCTION_EVAL_INT(lround);
DEFINE_UNARY_FUNCTION_EVAL_INT(llround);
DEFINE_UNARY_FUNCTION_EVAL(nearbyint, 2);
#endif
DEFINE_UNARY_FUNCTION_EVAL(rint, 0);
#ifndef KOKKOS_ENABLE_SYCL
DEFINE_UNARY_FUNCTION_EVAL_INT(lrint);
DEFINE_UNARY_FUNCTION_EVAL_INT(llrint);
#endif

DEFINE_UNARY_FUNCTION_EVAL_INT(ilogb);
DEFINE_UNARY_FUNCTION_EVAL(logb, 2);
#endif

#undef DEFINE_UNARY_FUNCTION_EVAL

#define DEFINE_BINARY_FUNCTION_EVAL(FUNC, ULP_FACTOR)                          \
  struct MathBinaryFunction_##FUNC {                                           \
    template <typename T, typename U>                                          \
    static KOKKOS_FUNCTION auto eval(T x, U y) {                               \
      static_assert(std::is_same_v<decltype(Kokkos::FUNC((T)0, (U)0)),         \
                                   math_binary_function_return_type_t<T, U>>); \
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
            std::is_same_v<decltype(std::FUNC((T)0, (U)0)),                    \
                           math_binary_function_return_type_t<T, U>>);         \
        return std::FUNC(x, y);                                                \
      }                                                                        \
    }                                                                          \
    static KOKKOS_FUNCTION int ulp_factor() { return ULP_FACTOR; }             \
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
DEFINE_BINARY_FUNCTION_EVAL(fmax, 0);
DEFINE_BINARY_FUNCTION_EVAL(fmin, 0);
#endif

#undef DEFINE_BINARY_FUNCTION_EVAL

#define DEFINE_BINARY_PREDICATE_EVAL(FUNC)                                     \
  struct MathBinaryPredicate_##FUNC {                                          \
    template <typename T, typename U>                                          \
    static KOKKOS_FUNCTION auto eval(T x, U y) {                               \
      static_assert(std::is_same_v<decltype(Kokkos::FUNC((T)0, (U)0)), bool>); \
      return Kokkos::FUNC(x, y);                                               \
    }                                                                          \
    template <typename T, typename U>                                          \
    static auto eval_std(T x, U y) {                                           \
      using std::FUNC;                                                         \
      return FUNC(x, y);                                                       \
    }                                                                          \
  };                                                                           \
  using kk_##FUNC = MathBinaryPredicate_##FUNC;                                \
  template <>                                                                  \
  struct math_function_name<MathBinaryPredicate_##FUNC> {                      \
    static constexpr char name[] = #FUNC;                                      \
  };                                                                           \
  constexpr char math_function_name<MathBinaryPredicate_##FUNC>::name[]

#ifndef KOKKOS_MATHEMATICAL_FUNCTIONS_SKIP_2
DEFINE_BINARY_PREDICATE_EVAL(isgreater);
DEFINE_BINARY_PREDICATE_EVAL(isgreaterequal);
DEFINE_BINARY_PREDICATE_EVAL(isless);
DEFINE_BINARY_PREDICATE_EVAL(islessequal);
DEFINE_BINARY_PREDICATE_EVAL(islessgreater);
DEFINE_BINARY_PREDICATE_EVAL(isunordered);
#endif

#undef DEFINE_BINARY_PREDICATE_EVAL

#define DEFINE_BINARY_INT_FUNCTION_EVAL(FUNC, ULP_FACTOR)                  \
  struct MathBinaryIntFunction_##FUNC {                                    \
    template <typename T, typename U>                                      \
    static KOKKOS_FUNCTION auto eval(T x, U y) {                           \
      static_assert(std::is_same_v<decltype(Kokkos::FUNC((T)0, (U)0)),     \
                                   math_unary_function_return_type_t<T>>); \
      static_assert(std::is_integral_v<U>);                                \
      return Kokkos::FUNC(x, y);                                           \
    }                                                                      \
    template <typename T, typename U>                                      \
    static auto eval_std(T x, U y) {                                       \
      static_assert(std::is_same_v<decltype(std::FUNC((T)0, (U)0)),        \
                                   math_unary_function_return_type_t<T>>); \
      static_assert(std::is_integral_v<U>);                                \
      return std::FUNC(x, y);                                              \
    }                                                                      \
    static KOKKOS_FUNCTION int ulp_factor() { return ULP_FACTOR; }         \
  };                                                                       \
  using kk_##FUNC = MathBinaryIntFunction_##FUNC;                          \
  template <>                                                              \
  struct math_function_name<MathBinaryIntFunction_##FUNC> {                \
    static constexpr char name[] = #FUNC;                                  \
  };                                                                       \
  constexpr char math_function_name<MathBinaryIntFunction_##FUNC>::name[]

#ifndef KOKKOS_MATHEMATICAL_FUNCTIONS_SKIP_2
DEFINE_BINARY_INT_FUNCTION_EVAL(ldexp, 0);
#if !defined(KOKKOS_ENABLE_SYCL) || (defined(FLT_RADIX) && FLT_RADIX == 2)
DEFINE_BINARY_INT_FUNCTION_EVAL(scalbn, 0);
#endif
#ifndef KOKKOS_ENABLE_SYCL
DEFINE_BINARY_INT_FUNCTION_EVAL(scalbln, 0);
#endif
#endif

#undef DEFINE_BINARY_INT_FUNCTION_EVAL

#define DEFINE_BINARY_PTR_FUNCTION_EVAL(FUNC, ULP_FACTOR)          \
  struct MathBinaryPtrFunction_##FUNC {                            \
    template <typename T, typename U>                              \
    static KOKKOS_FUNCTION auto eval(T x, U* y) {                  \
      return Kokkos::FUNC(x, y);                                   \
    }                                                              \
    template <typename T, typename U>                              \
    static auto eval_std(T x, U* y) {                              \
      return std::FUNC(x, y);                                      \
    }                                                              \
    static KOKKOS_FUNCTION int ulp_factor() { return ULP_FACTOR; } \
  };                                                               \
  template <>                                                      \
  struct math_function_name<MathBinaryPtrFunction_##FUNC> {        \
    static constexpr char name[] = #FUNC;                          \
  };                                                               \
  constexpr char math_function_name<MathBinaryPtrFunction_##FUNC>::name[];

#ifndef KOKKOS_MATHEMATICAL_FUNCTIONS_SKIP_2
DEFINE_BINARY_PTR_FUNCTION_EVAL(modf, 0)
#endif

#undef DEFINE_BINARY_PTR_FUNCTION_EVAL

#define DEFINE_TERNARY_INT_PTR_FUNCTION_EVAL(FUNC, ULP_FACTOR)                 \
  struct MathTernaryIntPtrFunction_##FUNC {                                    \
    template <typename T, typename U>                                          \
    static KOKKOS_FUNCTION auto eval(T x, U y, int* z) {                       \
      static_assert(                                                           \
          std::is_same_v<decltype(Kokkos::FUNC((T)0, (U)0, (int*)nullptr)),    \
                         math_binary_function_return_type_t<T, U>>);           \
      return Kokkos::FUNC(x, y, z);                                            \
    }                                                                          \
    template <typename T, typename U>                                          \
    static auto eval_std(T x, U y, int* z) {                                   \
      constexpr bool const x_is_half =                                         \
          (KE::Impl::is_float16<T>::value || KE::Impl::is_bfloat16<T>::value); \
      constexpr bool const y_is_half =                                         \
          (KE::Impl::is_float16<U>::value || KE::Impl::is_bfloat16<U>::value); \
      if constexpr (x_is_half && y_is_half)                                    \
        return std::FUNC(static_cast<float>(x), static_cast<float>(y), z);     \
      else if constexpr (x_is_half)                                            \
        return std::FUNC(static_cast<float>(x), y, z);                         \
      else if constexpr (y_is_half)                                            \
        return std::FUNC(x, static_cast<float>(y), z);                         \
      else {                                                                   \
        static_assert(                                                         \
            std::is_same_v<decltype(std::FUNC((T)0, (U)0, (int*)nullptr)),     \
                           math_binary_function_return_type_t<T, U>>);         \
        return std::FUNC(x, y, z);                                             \
      }                                                                        \
    }                                                                          \
    static KOKKOS_FUNCTION int ulp_factor() { return ULP_FACTOR; }             \
  };                                                                           \
  using kk3_##FUNC = MathTernaryIntPtrFunction_##FUNC;                         \
  template <>                                                                  \
  struct math_function_name<MathTernaryIntPtrFunction_##FUNC> {                \
    static constexpr char name[] = #FUNC;                                      \
  };                                                                           \
  constexpr char math_function_name<MathTernaryIntPtrFunction_##FUNC>::name[]

#ifndef KOKKOS_MATHEMATICAL_FUNCTIONS_SKIP_2
DEFINE_TERNARY_INT_PTR_FUNCTION_EVAL(remquo, 0);
#endif

#undef DEFINE_TERNARY_INT_PTR_FUNCTION_EVAL

#define DEFINE_BINARY_INT_PTR_FUNCTION_EVAL(FUNC, ULP_FACTOR)              \
  struct MathBinaryIntPtrFunction_##FUNC {                                 \
    template <typename T, typename U>                                      \
    static KOKKOS_FUNCTION auto eval(T x, U* y) {                          \
      static_assert(std::is_same_v<decltype(Kokkos::FUNC((T)0, (U*)0)),    \
                                   math_unary_function_return_type_t<T>>); \
      static_assert(std::is_integral_v<U>);                                \
      return Kokkos::FUNC(x, y);                                           \
    }                                                                      \
    template <typename T, typename U>                                      \
    static auto eval_std(T x, U* y) {                                      \
      static_assert(std::is_same_v<decltype(std::FUNC((T)0, (U*)0)),       \
                                   math_unary_function_return_type_t<T>>); \
      static_assert(std::is_integral_v<U>);                                \
      return std::FUNC(x, y);                                              \
    }                                                                      \
    static KOKKOS_FUNCTION int ulp_factor() { return ULP_FACTOR; }         \
  };                                                                       \
  using kk_##FUNC = MathBinaryIntPtrFunction_##FUNC;                       \
  template <>                                                              \
  struct math_function_name<MathBinaryIntPtrFunction_##FUNC> {             \
    static constexpr char name[] = #FUNC;                                  \
  };                                                                       \
  constexpr char math_function_name<MathBinaryIntPtrFunction_##FUNC>::name[]

#ifndef KOKKOS_MATHEMATICAL_FUNCTIONS_SKIP_2
DEFINE_BINARY_INT_PTR_FUNCTION_EVAL(frexp, 0);
#endif

#undef DEFINE_BINARY_INT_PTR_FUNCTION_EVAL

#define DEFINE_TERNARY_FUNCTION_EVAL(FUNC, ULP_FACTOR)                   \
  struct MathTernaryFunction_##FUNC {                                    \
    template <typename T, typename U, typename V>                        \
    static KOKKOS_FUNCTION auto eval(T x, U y, V z) {                    \
      static_assert(                                                     \
          std::is_same_v<decltype(Kokkos::FUNC((T)0, (U)0, (V)0)),       \
                         math_ternary_function_return_type_t<T, U, V>>); \
      return Kokkos::FUNC(x, y, z);                                      \
    }                                                                    \
    template <typename T, typename U, typename V>                        \
    static auto eval_std(T x, U y, V z) {                                \
      static_assert(                                                     \
          std::is_same_v<decltype(std::FUNC((T)0, (U)0, (V)0)),          \
                         math_ternary_function_return_type_t<T, U, V>>); \
      return std::FUNC(x, y, z);                                         \
    }                                                                    \
    static KOKKOS_FUNCTION int ulp_factor() { return ULP_FACTOR; }       \
  };                                                                     \
  using kk3_##FUNC = MathTernaryFunction_##FUNC;                         \
  template <>                                                            \
  struct math_function_name<MathTernaryFunction_##FUNC> {                \
    static constexpr char name[] = #FUNC;                                \
  };                                                                     \
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

template <class Space, class Func, class Arg, std::size_t N>
struct TestIntMathUnaryFunction : IntegerComparison {
  Arg val_[N];
  int res_[N];
  TestIntMathUnaryFunction(const Arg (&val)[N]) {
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
    bool ar = compare(Func::eval(val_[i]), res_[i]);
    if (!ar) {
      ++e;
      Kokkos::printf("value at %f which is %f was expected to be %f\n",
                     (double)val_[i], (double)Func::eval(val_[i]),
                     (double)res_[i]);
    }
  }
};

template <class Space, class... Func, class Arg, std::size_t N>
void do_test_int_math_unary_function(const Arg (&x)[N]) {
  (void)std::initializer_list<int>{
      (TestIntMathUnaryFunction<Space, Func, Arg, N>(x), 0)...};

  // test if potentially device specific math functions also work on host
  if constexpr (!std::is_same_v<Space, Kokkos::DefaultHostExecutionSpace>)
    (void)std::initializer_list<int>{
        (TestIntMathUnaryFunction<Kokkos::DefaultHostExecutionSpace, Func, Arg,
                                  N>(x),
         0)...};
}

#define TEST_INT_MATH_FUNCTION(FUNC) \
  do_test_int_math_unary_function<TEST_EXECSPACE, MathUnaryFunction_##FUNC>

template <class Half, class Space, class... Func, class Arg, std::size_t N>
void do_test_int_half_math_unary_function(const Arg (&x)[N]) {
  Half y[N];
  std::copy(x, x + N, y);  // cast to array of half type
  (void)std::initializer_list<int>{
      (TestIntMathUnaryFunction<Space, Func, Half, N>(y), 0)...};

  // test if potentially device specific math functions also work on host
  if constexpr (!std::is_same_v<Space, Kokkos::DefaultHostExecutionSpace>)
    (void)std::initializer_list<int>{
        (TestIntMathUnaryFunction<Kokkos::DefaultHostExecutionSpace, Func, Half,
                                  N>(y),
         0)...};
}

#define TEST_INT_HALF_MATH_FUNCTION(FUNC, T)              \
  do_test_int_half_math_unary_function<T, TEST_EXECSPACE, \
                                       MathUnaryFunction_##FUNC>

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

template <class Space, class Func, class Arg1, class Arg2,
          class Ret = math_unary_function_return_type_t<Arg1>>
struct TestMathBinaryIntFunction : FloatingPointComparison {
  Arg1 val1_;
  Arg2 val2_;
  Ret res_;
  TestMathBinaryIntFunction(Arg1 val1, Arg2 val2)
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
void do_test_math_binary_int_function(Arg1 arg1, Arg2 arg2) {
  (void)std::initializer_list<int>{
      (TestMathBinaryIntFunction<Space, Func, Arg1, Arg2>(arg1, arg2), 0)...};
}

template <class Space, class Func, class Arg,
          class Ret = math_unary_function_return_type_t<Arg>>
struct TestMathBinaryPtrFunction : FloatingPointComparison {
  Arg val_;
  Ret res_frac_;
  Ret res_int_;
  const char* m_name;
  TestMathBinaryPtrFunction(Arg val)
      : val_(val), m_name(math_function_name<Func>::name) {
    res_frac_ = Func::eval_std(val_, &res_int_);
    run();
  }
  void run() {
    int errors = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Space>(0, 1), *this, errors);
    ASSERT_EQ(errors, 0) << "Failed " << m_name << " check for "
                         << type_helper<Arg>::name();
  }
  KOKKOS_FUNCTION void operator()(int, int& e) const {
    Ret iptr;
    Ret frac     = Func::eval(val_, &iptr);
    bool ar_frac = compare(frac, res_frac_, Func::ulp_factor());
    bool ar_int  = compare(iptr, res_int_, Func::ulp_factor());
    if (!ar_frac || !ar_int) {
      ++e;
      Kokkos::printf("%s failed: Val %f -> Frac %f (exp %f), Int %f (exp %f)\n",
                     m_name, (double)val_, (double)frac, (double)res_frac_,
                     (double)iptr, (double)res_int_);
    }
  }
};

template <class Space, class... Func, class Arg>
void do_test_math_binary_ptr_function(Arg x) {
  (void)std::initializer_list<int>{
      (TestMathBinaryPtrFunction<Space, Func, Arg>(x), 0)...};
  if constexpr (!std::is_same_v<Space, Kokkos::DefaultHostExecutionSpace>) {
    (void)std::initializer_list<int>{(
        TestMathBinaryPtrFunction<Kokkos::DefaultHostExecutionSpace, Func, Arg>(
            x),
        0)...};
  }
}

template <class Space, class Func, class Arg1, class Arg2>
struct TestMathBinaryPredicate {
  Arg1 val1_;
  Arg2 val2_;
  bool res_;
  TestMathBinaryPredicate(Arg1 val1, Arg2 val2)
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
    bool ar = Func::eval(val1_, val2_) == res_;
    if (!ar) {
      ++e;
      Kokkos::printf("value at %f, %f which is %s was expected to be %s\n",
                     (double)val1_, (double)val2_,
                     Func::eval(val1_, val2_) ? "true" : "false",
                     res_ ? "true" : "false");
    }
  }
};

template <class Space, class... Func, class Arg1, class Arg2>
void do_test_math_binary_predicate(Arg1 arg1, Arg2 arg2) {
  (void)std::initializer_list<int>{
      (TestMathBinaryPredicate<Space, Func, Arg1, Arg2>(arg1, arg2), 0)...};
}

template <class Space, class Func, class Arg,
          class Ret = math_unary_function_return_type_t<Arg>>
struct TestMathBinaryIntPtrFunction : FloatingPointComparison {
  Arg val_;
  int res1_;
  Ret res2_;
  TestMathBinaryIntPtrFunction(Arg val)
      : val_(val), res2_(Func::eval_std(val, &res1_)) {
    run();
  }
  void run() {
    int errors = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Space>(0, 1), *this, errors);
    ASSERT_EQ(errors, 0) << "Failed check no error for "
                         << math_function_name<Func>::name << "("
                         << type_helper<Arg>::name() << ")";
  }
  KOKKOS_FUNCTION void operator()(int, int& e) const {
    int res1;
    Ret res2 = Func::eval(val_, &res1);
    bool ar1 = (res1 == res1_);
    bool ar2 = compare(res2, res2_, Func::ulp_factor());
    if (!(ar1 && ar2)) {
      ++e;
      Kokkos::printf(
          "values at %f which are %f, %f were expected to be %f, %f\n",
          (double)val_, (double)res1, (double)res2, (double)res1_,
          (double)res2_);
    }
  }
};

template <class Space, class... Func, class Arg>
void do_test_math_binary_int_ptr_function(Arg x) {
  (void)std::initializer_list<int>{
      (TestMathBinaryIntPtrFunction<Space, Func, Arg>(x), 0)...};
  if constexpr (!std::is_same_v<Space, Kokkos::DefaultHostExecutionSpace>) {
    (void)std::initializer_list<int>{
        (TestMathBinaryIntPtrFunction<Kokkos::DefaultHostExecutionSpace, Func,
                                      Arg>(x),
         0)...};
  }
}

template <class Space, class Func, class Arg1, class Arg2,
          class Ret = math_binary_function_return_type_t<Arg1, Arg2>>
struct TestMathTernaryIntPtrFunction : FloatingPointComparison {
  Arg1 val1_;
  Arg2 val2_;
  int val_;
  Ret res_;
  TestMathTernaryIntPtrFunction(Arg1 val1, Arg2 val2)
      : val1_(val1), val2_(val2), res_(Func::eval_std(val1, val2, &val_)) {
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
    int val;
    auto res  = Func::eval(val1_, val2_, &val);
    bool ar_1 = compare(res, res_, Func::ulp_factor());
    bool ar_2 = (val_ == val);
    if (!(ar_1 && ar_2)) {
      ++e;
      Kokkos::printf(
          "value at %f, %f which is %f and %i was expected to be %f and %i\n",
          (double)val1_, (double)val2_, (double)res, val, (double)res_, val_);
    }
  }
};

template <class Space, class... Func, class Arg1, class Arg2>
void do_test_math_ternary_int_ptr_function(Arg1 arg1, Arg2 arg2) {
  (void)std::initializer_list<int>{
      (TestMathTernaryIntPtrFunction<Space, Func, Arg1, Arg2>(arg1, arg2),
       0)...};
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
                         << type_helper<Arg2>::name() << ", "
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

  do_test_math_binary_int_function<TEST_EXECSPACE, kk_ldexp>(42.765f, 3);
  do_test_math_binary_int_function<TEST_EXECSPACE, kk_ldexp>(-15.123, 3);
  do_test_math_binary_int_function<TEST_EXECSPACE, kk_ldexp>(15, 3);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  do_test_math_binary_int_function<TEST_EXECSPACE, kk_ldexp>(1234.5678l, 3);
#endif

#if !defined(KOKKOS_ENABLE_SYCL) || (defined(FLT_RADIX) && FLT_RADIX == 2)
  do_test_math_binary_int_function<TEST_EXECSPACE, kk_scalbn>(42.765f, -4);
  do_test_math_binary_int_function<TEST_EXECSPACE, kk_scalbn>(-15.123, -4);
  do_test_math_binary_int_function<TEST_EXECSPACE, kk_scalbn>(15, -4);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  do_test_math_binary_int_function<TEST_EXECSPACE, kk_scalbn>(1234.5678l, -4);
#endif
#endif

#ifndef KOKKOS_ENABLE_SYCL
  do_test_math_binary_int_function<TEST_EXECSPACE, kk_scalbln>(42.765f, -4l);
  do_test_math_binary_int_function<TEST_EXECSPACE, kk_scalbln>(-15.123, -4l);
  do_test_math_binary_int_function<TEST_EXECSPACE, kk_scalbln>(15, -4l);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  do_test_math_binary_int_function<TEST_EXECSPACE, kk_scalbln>(1234.5678l, -4l);
#endif
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

TEST(TEST_CATEGORY, mathematical_functions_modf) {
  using Func = MathBinaryPtrFunction_modf;

  do_test_math_binary_ptr_function<TEST_EXECSPACE, Func>(42.765f);
  do_test_math_binary_ptr_function<TEST_EXECSPACE, Func>(-15.123);
  do_test_math_binary_ptr_function<TEST_EXECSPACE, Func>(15);

#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  do_test_math_binary_ptr_function<TEST_EXECSPACE, Func>(1234.5678l);
#endif
}

TEST(TEST_CATEGORY, mathematical_functions_frexp) {
  using Func = MathBinaryIntPtrFunction_frexp;

  do_test_math_binary_int_ptr_function<TEST_EXECSPACE, Func>(42.765f);
  do_test_math_binary_int_ptr_function<TEST_EXECSPACE, Func>(-15.123);
  do_test_math_binary_int_ptr_function<TEST_EXECSPACE, Func>(15);

#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  do_test_math_binary_int_ptr_function<TEST_EXECSPACE, Func>(1234.5678l);
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

TEST(TEST_CATEGORY, mathematical_functions_remquo) {
  do_test_math_ternary_int_ptr_function<TEST_EXECSPACE, kk3_remquo>(
      static_cast<KE::half_t>(2.f), static_cast<KE::half_t>(3.f));
  do_test_math_ternary_int_ptr_function<TEST_EXECSPACE, kk3_remquo>(
      static_cast<KE::bhalf_t>(2.f), static_cast<KE::bhalf_t>(3.f));
  do_test_math_ternary_int_ptr_function<TEST_EXECSPACE, kk3_remquo>(2.f, 3.f);
  do_test_math_ternary_int_ptr_function<TEST_EXECSPACE, kk3_remquo>(2., 3.);
  do_test_math_ternary_int_ptr_function<TEST_EXECSPACE, kk3_remquo>(2, 3);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  do_test_math_ternary_int_ptr_function<TEST_EXECSPACE, kk3_remquo>(2.l, 3.l);
#endif
}

TEST(TEST_CATEGORY, mathematical_functions_fmax_fmin) {
  do_test_math_binary_function<TEST_EXECSPACE, kk_fmax>(
      static_cast<KE::half_t>(2.f), static_cast<KE::half_t>(3.f));
  do_test_math_binary_function<TEST_EXECSPACE, kk_fmin>(
      static_cast<KE::half_t>(2.f), static_cast<KE::half_t>(3.f));
  do_test_math_binary_function<TEST_EXECSPACE, kk_fmax>(
      static_cast<KE::bhalf_t>(2.f), static_cast<KE::bhalf_t>(3.f));
  do_test_math_binary_function<TEST_EXECSPACE, kk_fmin>(
      static_cast<KE::bhalf_t>(2.f), static_cast<KE::bhalf_t>(3.f));
  do_test_math_binary_function<TEST_EXECSPACE, kk_fmax>(2.f, 3.f);
  do_test_math_binary_function<TEST_EXECSPACE, kk_fmin>(2.f, 3.f);
  do_test_math_binary_function<TEST_EXECSPACE, kk_fmax>(2., 3.);
  do_test_math_binary_function<TEST_EXECSPACE, kk_fmin>(2., 3.);
  do_test_math_binary_function<TEST_EXECSPACE, kk_fmax>(2, 3.f);
  do_test_math_binary_function<TEST_EXECSPACE, kk_fmin>(2, 3.f);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  do_test_math_binary_function<TEST_EXECSPACE, kk_fmax>(2.l, 3.l);
  do_test_math_binary_function<TEST_EXECSPACE, kk_fmin>(2.l, 3.l);
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

  TEST_MATH_FUNCTION(log2)({1, 23, 456, 7890});
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
  TEST_HALF_MATH_FUNCTION(rsqrt, KE::half_t)({10.f, 20.f, 30.f, 40.f});
  TEST_HALF_MATH_FUNCTION(rsqrt, KE::bhalf_t)({10.f, 20.f, 30.f, 40.f});
  TEST_MATH_FUNCTION(rsqrt)({10.f, 20.f, 30.f, 40.f});
  TEST_MATH_FUNCTION(rsqrt)({11.1, 22.2, 33.3, 44.4});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(rsqrt)({10.l, 20.l, 30.l, 40.l});
#endif

  TEST_MATH_FUNCTION(rcp)({-13, -9, 1, 7, 11});
  TEST_MATH_FUNCTION(rcp)({-13l, -9l, 1l, 7l, 11l});
  TEST_MATH_FUNCTION(rcp)({-13ll, -9ll, 1ll, 7ll, 11ll});
  TEST_MATH_FUNCTION(rcp)({2u, 3u, 9u, 13u, 17u});
  TEST_MATH_FUNCTION(rcp)({2ul, 3ul, 9ul, 13ul, 17ul});
  TEST_MATH_FUNCTION(rcp)({2ull, 3ull, 9ull, 13ull, 17ull});
  TEST_HALF_MATH_FUNCTION(rcp, KE::half_t)({-13.f, -9.f, 1.f, 7.f, 11.f});
  TEST_HALF_MATH_FUNCTION(rcp, KE::bhalf_t)({-13.f, -9.f, 1.f, 7.f, 11.f});
  TEST_MATH_FUNCTION(rcp)({-13.1f, -9.2f, 1.3f, 7.4f, 11.5f});
  TEST_MATH_FUNCTION(rcp)({-13.1, -9.2, 1.3, 7.4, 11.5});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(rcp)({-13.1l, -9.2l, 1.3l, 7.4l, 11.5l});
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
  // NOTE there can be no domain error, since int has enough range to represent
  // any possible rounded half_t. Thus lround, llround can be implemented by
  // just casting from round. Thus they are implemented and tested.
  TEST_INT_HALF_MATH_FUNCTION(lround, KE::half_t)
  ({2.3f, 2.5f, 2.7f, -2.3f, -2.5f, -2.7f, -0.0f});
  TEST_INT_MATH_FUNCTION(lround)({-3, -2, -1, 0, 1});
  TEST_INT_MATH_FUNCTION(lround)({-3l, -2l, -1l, 0l, 1l});
  TEST_INT_MATH_FUNCTION(lround)({-3ll, -2ll, -1ll, 0ll, 1ll});
  TEST_INT_MATH_FUNCTION(lround)({2u, 3u, 4u, 5u, 6u});
  TEST_INT_MATH_FUNCTION(lround)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_INT_MATH_FUNCTION(lround)({2ull, 3ull, 4ull, 5ull, 6ull});
  TEST_INT_MATH_FUNCTION(lround)
  ({2.3f, 2.5f, 2.7f, -2.3f, -2.5f, -2.7f, -0.0f});
  TEST_INT_MATH_FUNCTION(lround)({2.3, 2.5, 2.7, -2.3, -2.5, -2.7, -0.0});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_INT_MATH_FUNCTION(lround)
  ({2.3l, 2.5l, 2.7l, -2.3l, -2.5l, -2.7l, -0.0l});
#endif

  TEST_INT_HALF_MATH_FUNCTION(llround, KE::half_t)
  ({2.3f, 2.5f, 2.7f, -2.3f, -2.5f, -2.7f, -0.0f});
  TEST_INT_MATH_FUNCTION(llround)({-3, -2, -1, 0, 1});
  TEST_INT_MATH_FUNCTION(llround)({-3l, -2l, -1l, 0l, 1l});
  TEST_INT_MATH_FUNCTION(llround)({-3ll, -2ll, -1ll, 0ll, 1ll});
  TEST_INT_MATH_FUNCTION(llround)({2u, 3u, 4u, 5u, 6u});
  TEST_INT_MATH_FUNCTION(llround)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_INT_MATH_FUNCTION(llround)({2ull, 3ull, 4ull, 5ull, 6ull});
  TEST_INT_MATH_FUNCTION(llround)
  ({2.3f, 2.5f, 2.7f, -2.3f, -2.5f, -2.7f, -0.0f});
  TEST_INT_MATH_FUNCTION(llround)({2.3, 2.5, 2.7, -2.3, -2.5, -2.7, -0.0});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_INT_MATH_FUNCTION(llround)
  ({2.3l, 2.5l, 2.7l, -2.3l, -2.5l, -2.7l, -0.0l});
#endif
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

  TEST_MATH_FUNCTION(rint)({-3, -2, -1, 0, 1});
  TEST_MATH_FUNCTION(rint)({-3l, -2l, -1l, 0l, 1l});
  TEST_MATH_FUNCTION(rint)({-3ll, -2ll, -1ll, 0ll, 1ll});
  TEST_MATH_FUNCTION(rint)({2u, 3u, 4u, 5u, 6u});
  TEST_MATH_FUNCTION(rint)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_MATH_FUNCTION(rint)({2ull, 3ull, 4ull, 5ull, 6ull});
  TEST_HALF_MATH_FUNCTION(rint, KE::half_t)
  ({2.3f, 2.5f, 2.7f, -2.3f, -2.5f, -2.7f, -0.0f});
  TEST_HALF_MATH_FUNCTION(rint, KE::bhalf_t)
  ({2.3f, 2.5f, 2.7f, -2.3f, -2.5f, -2.7f, -0.0f});
  TEST_MATH_FUNCTION(rint)({2.3f, 2.5f, 2.7f, -2.3f, -2.5f, -2.7f, -0.0f});
  TEST_MATH_FUNCTION(rint)({2.3, 2.5, 2.7, -2.3, -2.5, -2.7, -0.0});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(rint)({2.3l, 2.5l, 2.7l, -2.3l, -2.5l, -2.7l, -0.0l});
#endif

#ifndef KOKKOS_ENABLE_SYCL
  // NOTE there can be no domain error, since int has enough range to represent
  // any possible rounded half_t. Thus lrint, llrint can be implemented by just
  // casting from rint. Thus they are implemented and tested.
  TEST_INT_HALF_MATH_FUNCTION(lrint, KE::half_t)
  ({2.3f, 2.5f, 2.7f, -2.3f, -2.5f, -2.7f, -0.0f});
  TEST_INT_MATH_FUNCTION(lrint)({-3, -2, -1, 0, 1});
  TEST_INT_MATH_FUNCTION(lrint)({-3l, -2l, -1l, 0l, 1l});
  TEST_INT_MATH_FUNCTION(lrint)({-3ll, -2ll, -1ll, 0ll, 1ll});
  TEST_INT_MATH_FUNCTION(lrint)({2u, 3u, 4u, 5u, 6u});
  TEST_INT_MATH_FUNCTION(lrint)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_INT_MATH_FUNCTION(lrint)({2ull, 3ull, 4ull, 5ull, 6ull});
  TEST_INT_MATH_FUNCTION(lrint)({2.3f, 2.5f, 2.7f, -2.3f, -2.5f, -2.7f, -0.0f});
  TEST_INT_MATH_FUNCTION(lrint)({2.3, 2.5, 2.7, -2.3, -2.5, -2.7, -0.0});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_INT_MATH_FUNCTION(lrint)({2.3l, 2.5l, 2.7l, -2.3l, -2.5l, -2.7l, -0.0l});
#endif

  TEST_INT_HALF_MATH_FUNCTION(llrint, KE::half_t)
  ({2.3f, 2.5f, 2.7f, -2.3f, -2.5f, -2.7f, -0.0f});
  TEST_INT_MATH_FUNCTION(llrint)({-3, -2, -1, 0, 1});
  TEST_INT_MATH_FUNCTION(llrint)({-3l, -2l, -1l, 0l, 1l});
  TEST_INT_MATH_FUNCTION(llrint)({-3ll, -2ll, -1ll, 0ll, 1ll});
  TEST_INT_MATH_FUNCTION(llrint)({2u, 3u, 4u, 5u, 6u});
  TEST_INT_MATH_FUNCTION(llrint)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_INT_MATH_FUNCTION(llrint)({2ull, 3ull, 4ull, 5ull, 6ull});
  TEST_INT_MATH_FUNCTION(llrint)
  ({2.3f, 2.5f, 2.7f, -2.3f, -2.5f, -2.7f, -0.0f});
  TEST_INT_MATH_FUNCTION(llrint)({2.3, 2.5, 2.7, -2.3, -2.5, -2.7, -0.0});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_INT_MATH_FUNCTION(llrint)
  ({2.3l, 2.5l, 2.7l, -2.3l, -2.5l, -2.7l, -0.0l});
#endif
#endif
}

TEST(TEST_CATEGORY,
     mathematical_functions_floating_point_manipulation_functions) {
  TEST_INT_MATH_FUNCTION(ilogb)({1, 13, 132, 1282, 7839});
  TEST_INT_MATH_FUNCTION(ilogb)({1l, 13l, 132l, 1282l, 7839l});
  TEST_INT_MATH_FUNCTION(ilogb)({1ll, 13ll, 132ll, 1282ll, 7839ll});
  TEST_INT_MATH_FUNCTION(ilogb)({1u, 13u, 132u, 1282u, 7839u});
  TEST_INT_MATH_FUNCTION(ilogb)({1ul, 13ul, 132ul, 1282ul, 7839ul});
  TEST_INT_MATH_FUNCTION(ilogb)({1ull, 13ull, 132ull, 1282ull, 7839ull});
  TEST_INT_HALF_MATH_FUNCTION(ilogb, KE::half_t)
  ({0.3f, 13.7f, 132.7f, 1282.4f, 7839.9f});
  TEST_INT_HALF_MATH_FUNCTION(ilogb, KE::bhalf_t)
  ({0.3f, 13.7f, 132.7f, 1282.4f, 7839.9f});
  TEST_INT_MATH_FUNCTION(ilogb)({0.3f, 13.7f, 132.7f, 1282.4f, 7839.9f});
  TEST_INT_MATH_FUNCTION(ilogb)({0.3, 13.7, 132.7, 1282.4, 7839.9});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_INT_MATH_FUNCTION(ilogb)({0.3l, 13.7l, 132.7l, 1282.4l, 7839.9l});
#endif

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

#if defined(KOKKOS_HALF_T_IS_FLOAT) && KOKKOS_HALF_T_IS_FLOAT
  do_test_math_binary_function<TEST_EXECSPACE, kk_nextafter>(
      0, static_cast<KE::half_t>(1.f));
  do_test_math_binary_function<TEST_EXECSPACE, kk_nextafter>(
      1, static_cast<KE::half_t>(2.f));
#endif
#if defined(KOKKOS_BHALF_T_IS_FLOAT) && KOKKOS_BHALF_T_IS_FLOAT
  do_test_math_binary_function<TEST_EXECSPACE, kk_nextafter>(
      0, static_cast<KE::bhalf_t>(1.f));
  do_test_math_binary_function<TEST_EXECSPACE, kk_nextafter>(
      1, static_cast<KE::bhalf_t>(2.f));
#endif
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
      Kokkos::printf("failed abs(KE::half_t)\n");
    }
    if (abs(static_cast<KE::bhalf_t>(4.f)) != static_cast<KE::bhalf_t>(4.f) ||
        abs(static_cast<KE::bhalf_t>(-4.f)) != static_cast<KE::bhalf_t>(4.f)) {
      ++e;
      Kokkos::printf("failed abs(KE::bhalf_t)\n");
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
#if !__FINITE_MATH_ONLY__
    // special values
    using Kokkos::isinf;
    using Kokkos::isnan;
    if (abs(-0.) != 0. || !isinf(abs(-INFINITY)) || !isnan(abs(-NAN))) {
      ++e;
      Kokkos::printf("failed abs(floating_point) special values\n");
    }
#endif

    static_assert(std::is_same_v<decltype(abs(1)), int>);
    static_assert(std::is_same_v<decltype(abs(2l)), long>);
    static_assert(std::is_same_v<decltype(abs(3ll)), long long>);
    static_assert(std::is_same_v<decltype(abs(static_cast<KE::half_t>(4.f))),
                                 KE::half_t>);
    static_assert(std::is_same_v<decltype(abs(static_cast<KE::bhalf_t>(4.f))),
                                 KE::bhalf_t>);
    static_assert(std::is_same_v<decltype(abs(4.f)), float>);
    static_assert(std::is_same_v<decltype(abs(5.)), double>);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    static_assert(std::is_same_v<decltype(abs(6.l)), long double>);
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
      Kokkos::printf("failed fabs(float)\n");
    }
    if (fabs(static_cast<KE::half_t>(4.f)) != static_cast<KE::half_t>(4.f) ||
        fabs(static_cast<KE::half_t>(-4.f)) != static_cast<KE::half_t>(4.f)) {
      ++e;
      Kokkos::printf("failed fabs(KE::half_t)\n");
    }
    if (fabs(static_cast<KE::bhalf_t>(4.f)) != static_cast<KE::bhalf_t>(4.f) ||
        fabs(static_cast<KE::bhalf_t>(-4.f)) != static_cast<KE::bhalf_t>(4.f)) {
      ++e;
      Kokkos::printf("failed fabs(KE::bhalf_t)\n");
    }
    if (fabs(5.) != 5. || fabs(-5.) != 5.) {
      ++e;
      Kokkos::printf("failed fabs(double)\n");
    }
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    if (fabs(6.l) != 6.l || fabs(-6.l) != 6.l) {
      ++e;
      Kokkos::printf("failed fabs(long double)\n");
    }
#endif
#if !__FINITE_MATH_ONLY__
    // special values
    using Kokkos::isinf;
    using Kokkos::isnan;
    if (fabs(-0.) != 0. || !isinf(fabs(-INFINITY)) || !isnan(fabs(-NAN))) {
      ++e;
      Kokkos::printf("failed fabs(floating_point) special values\n");
    }
#endif

    static_assert(std::is_same_v<decltype(fabs(static_cast<KE::half_t>(4.f))),
                                 KE::half_t>);
    static_assert(std::is_same_v<decltype(fabs(static_cast<KE::bhalf_t>(4.f))),
                                 KE::bhalf_t>);
    static_assert(std::is_same_v<decltype(fabs(4.f)), float>);
    static_assert(std::is_same_v<decltype(fabs(5.)), double>);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    static_assert(std::is_same_v<decltype(fabs(6.l)), long double>);
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
    if (!compare(fmod(6.2f, 4.f), 2.2f, 1) ||
        !compare(fmod(-6.2f, 4.f), -2.2f, 1)) {
      ++e;
      Kokkos::printf("failed fmod(float)\n");
    }
    if (!compare(
            fmod(static_cast<KE::half_t>(6.2f), static_cast<KE::half_t>(4.f)),
            static_cast<KE::half_t>(2.2f), 1) ||
        !compare(
            fmod(static_cast<KE::half_t>(-6.2f), static_cast<KE::half_t>(4.f)),
            -static_cast<KE::half_t>(2.2f), 1)) {
      ++e;
      Kokkos::printf("failed fmod(KE::half_t)\n");
    }
    if (!compare(
            fmod(static_cast<KE::bhalf_t>(6.2f), static_cast<KE::bhalf_t>(4.f)),
            static_cast<KE::bhalf_t>(2.2f), 1) ||
        !compare(fmod(static_cast<KE::bhalf_t>(-6.2f),
                      static_cast<KE::bhalf_t>(4.f)),
                 -static_cast<KE::bhalf_t>(2.2f), 1)) {
      ++e;
      Kokkos::printf("failed fmod(KE::bhalf_t)\n");
    }
    if (!compare(fmod(6.2, 4.), 2.2, 1) || !compare(fmod(-6.2, 4.), -2.2, 1)) {
      ++e;
      Kokkos::printf("failed fmod(double)\n");
    }
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    if (!compare(fmod(6.2l, 4.l), 2.2l, 1) ||
        !compare(fmod(-6.2l, 4.l), -2.2l, 1)) {
      ++e;
      Kokkos::printf("failed fmod(long double)\n");
    }
#endif
#if !__FINITE_MATH_ONLY__
    // special values
    using Kokkos::isinf;
    using Kokkos::isnan;
    if (!isnan(fmod(-KE::infinity<float>::value, 1.f)) ||
        !(fmod(5.f, -KE::infinity<float>::value) == 5.f) ||
        !isnan(fmod(5.f, 0.f)) ||
        !isnan(fmod(-KE::quiet_NaN<float>::value, 1.f)) ||
        !isnan(fmod(1.f, -KE::quiet_NaN<float>::value))) {
      ++e;
      Kokkos::printf("failed fmod(floating_point) special values\n");
    }
#endif

    static_assert(std::is_same_v<decltype(fmod(static_cast<KE::half_t>(4.f),
                                               static_cast<KE::half_t>(4.f))),
                                 KE::half_t>);
    static_assert(std::is_same_v<decltype(fmod(static_cast<KE::bhalf_t>(4.f),
                                               static_cast<KE::bhalf_t>(4.f))),
                                 KE::bhalf_t>);
    static_assert(std::is_same_v<decltype(fmod(4.f, 4.f)), float>);
    static_assert(std::is_same_v<decltype(fmod(5., 5.)), double>);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    static_assert(std::is_same_v<decltype(fmod(6.l, 6.l)), long double>);
#endif
  }
};

TEST(TEST_CATEGORY, mathematical_functions_remainder_function) {
  TestFloatingPointRemainderFunction<TEST_EXECSPACE>();
}

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
    if (!compare(remainder(6.2f, 4.f), -1.8f, 2) ||
        !compare(remainder(-6.2f, 4.f), 1.8f, 2)) {
      ++e;
      Kokkos::printf("failed remainder(float)\n");
    }
    if (!compare(remainder(static_cast<KE::half_t>(6.2f),
                           static_cast<KE::half_t>(4.f)),
                 static_cast<KE::half_t>(-1.8f), 2) ||
        !compare(remainder(static_cast<KE::half_t>(-6.2f),
                           static_cast<KE::half_t>(4.f)),
                 static_cast<KE::half_t>(1.8f), 2)) {
      ++e;
      Kokkos::printf("failed remainder(KE::half_t)\n");
    }
    if (!compare(remainder(static_cast<KE::bhalf_t>(6.2f),
                           static_cast<KE::bhalf_t>(4.f)),
                 static_cast<KE::bhalf_t>(-1.8f), 2) ||
        !compare(remainder(static_cast<KE::bhalf_t>(-6.2f),
                           static_cast<KE::bhalf_t>(4.f)),
                 static_cast<KE::bhalf_t>(1.8f), 2)) {
      ++e;
      Kokkos::printf("failed remainder(KE::bhalf_t)\n");
    }
    if (!compare(remainder(6.2, 4.), -1.8, 2) ||
        !compare(remainder(-6.2, 4.), 1.8, 2)) {
      ++e;
      Kokkos::printf("failed remainder(double)\n");
    }
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    if (!compare(remainder(6.2l, 4.l), -1.8l, 2) ||
        !compare(remainder(-6.2l, 4.l), 1.8l, 2)) {
      ++e;
      Kokkos::printf("failed remainder(long double)\n");
    }
#endif
#if !__FINITE_MATH_ONLY__
    // special values
    using Kokkos::isinf;
    using Kokkos::isnan;
    if (!isnan(remainder(-KE::infinity<float>::value, 2.f)) ||
        !isnan(remainder(-KE::quiet_NaN<float>::value, 2.f))) {
      ++e;
      Kokkos::printf("failed remainder(floating_point) special values\n");
    }
#endif

    static_assert(
        std::is_same_v<decltype(remainder(static_cast<KE::half_t>(4.f),
                                          static_cast<KE::half_t>(4.f))),
                       KE::half_t>);
    static_assert(
        std::is_same_v<decltype(remainder(static_cast<KE::bhalf_t>(4.f),
                                          static_cast<KE::bhalf_t>(4.f))),
                       KE::bhalf_t>);
    static_assert(std::is_same_v<decltype(remainder(4.f, 4.f)), float>);
    static_assert(std::is_same_v<decltype(remainder(5., 5.)), double>);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    static_assert(std::is_same_v<decltype(remainder(6.l, 6.l)), long double>);
#endif
  }
};

TEST(TEST_CATEGORY, mathematical_functions_ieee_remainder_function) {
  TestIEEEFloatingPointRemainderFunction<TEST_EXECSPACE>();
}

// TODO: TestFpClassify, see https://github.com/kokkos/kokkos/issues/6279

#ifndef KOKKOS_MATHEMATICAL_FUNCTIONS_SKIP_2

// Known to fail with
// * CUDA 12.4 and GCC 13.2
// * CUDA 12.8 and GCC 13.3, 14.2
// * CUDA 13.0 and GCC 15.1
#if defined(KOKKOS_COMPILER_NVCC) && \
    (defined(KOKKOS_COMPILER_GNU) && \
     (KOKKOS_COMPILER_GNU >= 1300 && KOKKOS_COMPILER_GNU < 1600))
#define KOKKOS_TEST_WORKAROUND_DEPRECATED_STD_ITERATOR_WARNINGS_PUSH() \
  KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_PUSH()

#define KOKKOS_TEST_WORKAROUND_DEPRECATED_STD_ITERATOR_WARNINGS_POP() \
  KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_POP()
#else
#define KOKKOS_TEST_WORKAROUND_DEPRECATED_STD_ITERATOR_WARNINGS_PUSH()
#define KOKKOS_TEST_WORKAROUND_DEPRECATED_STD_ITERATOR_WARNINGS_POP()
#endif

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
    if (!isfinite(static_cast<KE::half_t>(2.f)) ||
        isfinite(quiet_NaN<KE::half_t>::value) ||
        isfinite(signaling_NaN<KE::half_t>::value) ||
        isfinite(infinity<KE::half_t>::value)) {
      ++e;
      Kokkos::printf("failed isfinite(KE::half_t)\n");
    }
    if (!isfinite(static_cast<KE::bhalf_t>(2.f)) ||
        isfinite(quiet_NaN<KE::bhalf_t>::value) ||
        isfinite(signaling_NaN<KE::bhalf_t>::value) ||
        isfinite(infinity<KE::bhalf_t>::value)) {
      ++e;
      Kokkos::printf("failed isfinite(KE::bhalf_t)\n");
    }
#endif
    if (!isfinite(3.) || isfinite(quiet_NaN<double>::value) ||
        isfinite(signaling_NaN<double>::value) ||
        isfinite(infinity<double>::value)) {
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

    KOKKOS_TEST_WORKAROUND_DEPRECATED_STD_ITERATOR_WARNINGS_PUSH()
    static_assert(std::is_same_v<decltype(isfinite(1)), bool>);
    static_assert(std::is_same_v<decltype(isfinite(2.f)), bool>);
    static_assert(std::is_same_v<decltype(isfinite(3.)), bool>);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    static_assert(std::is_same_v<decltype(isfinite(4.l)), bool>);
#endif

    KOKKOS_TEST_STATIC_ASSERT_UNARY_PREDICATE(Kokkos::isfinite, float, bool);
    KOKKOS_TEST_STATIC_ASSERT_UNARY_PREDICATE(Kokkos::isfinite, double, bool);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    KOKKOS_TEST_STATIC_ASSERT_UNARY_PREDICATE(Kokkos::isfinite, long double,
                                              bool);
#endif
    KOKKOS_TEST_WORKAROUND_DEPRECATED_STD_ITERATOR_WARNINGS_POP()
  }
};

TEST(TEST_CATEGORY, mathematical_functions_isfinite) {
#if __FINITE_MATH_ONLY__
  GTEST_SKIP() << "skipping when compiling with -ffinite-math-only";
#endif
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
    if (isinf(static_cast<KE::half_t>(2.f)) ||
        isinf(quiet_NaN<KE::half_t>::value) ||
        isinf(signaling_NaN<KE::half_t>::value) ||
        !isinf(infinity<KE::half_t>::value)) {
      ++e;
      Kokkos::printf("failed isinf(KE::half_t)\n");
    }
    if (isinf(static_cast<KE::bhalf_t>(2.f)) ||
        isinf(quiet_NaN<KE::bhalf_t>::value) ||
        isinf(signaling_NaN<KE::bhalf_t>::value) ||
        !isinf(infinity<KE::bhalf_t>::value)) {
      ++e;
      Kokkos::printf("failed isinf(KE::bhalf_t)\n");
    }
#endif
    if (isinf(3.) || isinf(quiet_NaN<double>::value) ||
        isinf(signaling_NaN<double>::value) ||
        !isinf(infinity<double>::value)) {
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

    KOKKOS_TEST_WORKAROUND_DEPRECATED_STD_ITERATOR_WARNINGS_PUSH()
    static_assert(std::is_same_v<decltype(isinf(1)), bool>);
    static_assert(std::is_same_v<decltype(isinf(2.f)), bool>);
    static_assert(std::is_same_v<decltype(isinf(3.)), bool>);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    static_assert(std::is_same_v<decltype(isinf(4.l)), bool>);
#endif

    KOKKOS_TEST_STATIC_ASSERT_UNARY_PREDICATE(Kokkos::isinf, float, bool);
    KOKKOS_TEST_STATIC_ASSERT_UNARY_PREDICATE(Kokkos::isinf, double, bool);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    KOKKOS_TEST_STATIC_ASSERT_UNARY_PREDICATE(Kokkos::isinf, long double, bool);
#endif
    KOKKOS_TEST_WORKAROUND_DEPRECATED_STD_ITERATOR_WARNINGS_POP()
  }
};

TEST(TEST_CATEGORY, mathematical_functions_isinf) {
#if __FINITE_MATH_ONLY__
  GTEST_SKIP() << "skipping when compiling with -ffinite-math-only";
#endif
  TestIsInf<TEST_EXECSPACE>();
}

template <class Space>
struct TestFpClassify {
  TestFpClassify() { run(); }
  void run() const {
    int errors = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Space>(0, 1), *this, errors);
    ASSERT_EQ(errors, 0);
  }
  KOKKOS_FUNCTION void operator()(int, int& e) const {
    using KE::denorm_min;
    using KE::infinity;
    using KE::quiet_NaN;
    using KE::signaling_NaN;
    using Kokkos::fpclassify;

    if (fpclassify(0) != FP_ZERO || fpclassify(1) != FP_NORMAL) {
      ++e;
      Kokkos::printf("failed fpclassify(integral)\n");
    }

    if (fpclassify(0.f) != FP_ZERO || fpclassify(-0.f) != FP_ZERO ||
        fpclassify(1.f) != FP_NORMAL
#if !__FINITE_MATH_ONLY__
        || fpclassify(signaling_NaN<float>::value) != FP_NAN ||
        fpclassify(quiet_NaN<float>::value) != FP_NAN ||
        fpclassify(infinity<float>::value) != FP_INFINITE ||
        fpclassify(denorm_min<float>::value) != FP_SUBNORMAL
#endif
    ) {
      ++e;
      Kokkos::printf("failed fpclassify(float)\n");
    }

    if (fpclassify(0.) != FP_ZERO || fpclassify(-0.) != FP_ZERO ||
        fpclassify(1.) != FP_NORMAL
#if !__FINITE_MATH_ONLY__
        || fpclassify(signaling_NaN<double>::value) != FP_NAN ||
        fpclassify(quiet_NaN<double>::value) != FP_NAN ||
        fpclassify(infinity<double>::value) != FP_INFINITE ||
        fpclassify(denorm_min<double>::value) != FP_SUBNORMAL
#endif
    ) {
      ++e;
      Kokkos::printf("failed fpclassify(double)\n");
    }

#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    if (fpclassify(0.l) != FP_ZERO || fpclassify(-0.l) != FP_ZERO ||
        fpclassify(1.l) != FP_NORMAL
#if !__FINITE_MATH_ONLY__
        || fpclassify(signaling_NaN<long double>::value) != FP_NAN ||
        fpclassify(quiet_NaN<long double>::value) != FP_NAN ||
        fpclassify(infinity<long double>::value) != FP_INFINITE ||
        fpclassify(denorm_min<long double>::value) != FP_SUBNORMAL
#endif
    ) {
      ++e;
      Kokkos::printf("failed fpclassify(long double)\n");
    }
#endif

    if (fpclassify(static_cast<KE::half_t>(0.f)) != FP_ZERO ||
        fpclassify(static_cast<KE::half_t>(-0.f)) != FP_ZERO ||
        fpclassify(static_cast<KE::half_t>(1.f)) != FP_NORMAL
#if !__FINITE_MATH_ONLY__
#if !(defined(KOKKOS_ENABLE_CUDA) &&                         \
      defined(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE) && \
      defined(KOKKOS_COMPILER_CLANG))
        // FIXME internal compiler error for Clang+Cuda and RDC
        || fpclassify(signaling_NaN<KE::half_t>::value) != FP_NAN ||
        fpclassify(quiet_NaN<KE::half_t>::value) != FP_NAN ||
        fpclassify(infinity<KE::half_t>::value) != FP_INFINITE ||
        fpclassify(denorm_min<KE::half_t>::value) != FP_SUBNORMAL
#endif
#endif
    ) {
      ++e;
      Kokkos::printf("failed fpclassify(Kokkos::Experimental::half_t)\n");
    }

    if (fpclassify(static_cast<KE::bhalf_t>(0.f)) != FP_ZERO ||
        fpclassify(static_cast<KE::bhalf_t>(-0.f)) != FP_ZERO ||
        fpclassify(static_cast<KE::bhalf_t>(1.f)) != FP_NORMAL
#if !__FINITE_MATH_ONLY__
        || fpclassify(signaling_NaN<KE::bhalf_t>::value) != FP_NAN ||
        fpclassify(quiet_NaN<KE::bhalf_t>::value) != FP_NAN ||
        fpclassify(infinity<KE::bhalf_t>::value) != FP_INFINITE ||
        fpclassify(denorm_min<KE::bhalf_t>::value) != FP_SUBNORMAL
#endif
    ) {
      ++e;
      Kokkos::printf("failed fpclassify(Kokkos::Experimental::bhalf_t)\n");
    }
  }
};

TEST(TEST_CATEGORY, mathematical_functions_fpclassify) {
  TestFpClassify<TEST_EXECSPACE>();
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
    if (isnan(static_cast<KE::half_t>(2.f)) ||
        !isnan(quiet_NaN<KE::half_t>::value) ||
        !isnan(signaling_NaN<KE::half_t>::value) ||
        isnan(infinity<KE::half_t>::value)) {
      ++e;
      Kokkos::printf("failed isnan(KE::half_t)\n");
    }
    if (isnan(static_cast<KE::bhalf_t>(2.f)) ||
        !isnan(quiet_NaN<KE::bhalf_t>::value) ||
        !isnan(signaling_NaN<KE::bhalf_t>::value) ||
        isnan(infinity<KE::bhalf_t>::value)) {
      ++e;
      Kokkos::printf("failed isnan(KE::bhalf_t)\n");
    }
    if (isnan(3.) || !isnan(quiet_NaN<double>::value) ||
        !isnan(signaling_NaN<double>::value) ||
        isnan(infinity<double>::value)) {
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

    KOKKOS_TEST_WORKAROUND_DEPRECATED_STD_ITERATOR_WARNINGS_PUSH()
    static_assert(std::is_same_v<decltype(isnan(1)), bool>);
    static_assert(std::is_same_v<decltype(isnan(2.f)), bool>);
    static_assert(std::is_same_v<decltype(isnan(3.)), bool>);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    static_assert(std::is_same_v<decltype(isnan(4.l)), bool>);
#endif

    KOKKOS_TEST_STATIC_ASSERT_UNARY_PREDICATE(Kokkos::isnan, float, bool);
    KOKKOS_TEST_STATIC_ASSERT_UNARY_PREDICATE(Kokkos::isnan, double, bool);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    KOKKOS_TEST_STATIC_ASSERT_UNARY_PREDICATE(Kokkos::isnan, long double, bool);
#endif
    KOKKOS_TEST_WORKAROUND_DEPRECATED_STD_ITERATOR_WARNINGS_POP()
  }
};

TEST(TEST_CATEGORY, mathematical_functions_isnan) {
#if __FINITE_MATH_ONLY__
  GTEST_SKIP() << "skipping when compiling with -ffinite-math-only";
#endif
  TestIsNaN<TEST_EXECSPACE>();
}

#define DEVICE_ASSERT(CALL)                          \
  if (!(CALL)) {                                     \
    printf(KOKKOS_IMPL_STRINGIFY(CALL) " failed\n"); \
    ++e;                                             \
  }

template <class Space>
struct TestIsNormal {
  TestIsNormal() { run(); }
  void run() const {
    int errors = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Space>(0, 1), *this, errors);
    ASSERT_EQ(errors, 0);
  }
  KOKKOS_FUNCTION void operator()(int, int& e) const {
    using KE::denorm_min;
    using KE::infinity;
    using KE::norm_min;
    using KE::quiet_NaN;
    using KE::signaling_NaN;
    using Kokkos::isnormal;
    if (isnormal(0) || !isnormal(1) || !isnormal(INT_MAX)) {
      ++e;
      Kokkos::printf("failed isnormal(integral)\n");
    }
    if (isnormal(0.f) || !isnormal(2.f) || !isnormal(-3.f) ||
        isnormal(quiet_NaN<float>::value) ||
        isnormal(signaling_NaN<float>::value) ||
        isnormal(infinity<float>::value) ||
        isnormal(denorm_min<float>::value) ||
        !isnormal(norm_min<float>::value)) {
      ++e;
      Kokkos::printf("failed isnormal(float)\n");
    }
#if !(defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_COMPILER_MSVC))
    if (isnormal(static_cast<KE::half_t>(0.f)) ||
        !isnormal(static_cast<KE::half_t>(2.f)) ||
        !isnormal(static_cast<KE::half_t>(-2.f))
#if !(defined(KOKKOS_ENABLE_CUDA) &&                         \
      defined(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE) && \
      defined(KOKKOS_COMPILER_CLANG))
        // FIXME internal compiler error for Clang+Cuda and RDC
        || isnormal(quiet_NaN<KE::half_t>::value) ||
        isnormal(signaling_NaN<KE::half_t>::value) ||
        isnormal(infinity<KE::half_t>::value) ||
        isnormal(denorm_min<KE::half_t>::value) ||
        !isnormal(norm_min<KE::half_t>::value)
#endif
    ) {
      ++e;
      Kokkos::printf("failed isnormal(KE::half_t)\n");
    }
    if (isnormal(static_cast<KE::bhalf_t>(0.f)) ||
        !isnormal(static_cast<KE::bhalf_t>(2.f)) ||
        !isnormal(static_cast<KE::bhalf_t>(-2.f)) ||
        isnormal(quiet_NaN<KE::bhalf_t>::value) ||
        isnormal(signaling_NaN<KE::bhalf_t>::value) ||
        isnormal(infinity<KE::bhalf_t>::value) ||
        isnormal(denorm_min<KE::bhalf_t>::value) ||
        !isnormal(norm_min<KE::bhalf_t>::value)) {
      ++e;
      Kokkos::printf("failed isnormal(KE::bhalf_t)\n");
    }
#endif
    if (isnormal(0.) || !isnormal(3.) || !isnormal(-3.) ||
        isnormal(quiet_NaN<double>::value) ||
        isnormal(signaling_NaN<double>::value) ||
        isnormal(infinity<double>::value) ||
        isnormal(denorm_min<double>::value) ||
        !isnormal(norm_min<double>::value)) {
      ++e;
      Kokkos::printf("failed isnormal(double)\n");
    }
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    if (isnormal(0.l) || !isnormal(4.l) || !isnormal(-4.l) ||
        isnormal(quiet_NaN<long double>::value) ||
        isnormal(signaling_NaN<long double>::value) ||
        isnormal(infinity<long double>::value) ||
        isnormal(denorm_min<long double>::value) ||
        !isnormal(norm_min<long double>::value)) {
      ++e;
      Kokkos::printf("failed isnormal(long double)\n");
    }
#endif
    // special values
    if (isnormal(INFINITY) || isnormal(NAN) ||
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
        !isnormal(LDBL_MAX) || !isnormal(LDBL_MIN) || isnormal(LDBL_TRUE_MIN) ||
#endif
        !isnormal(FLT_MAX) || !isnormal(FLT_MIN) || isnormal(FLT_TRUE_MIN) ||
        !isnormal(DBL_MAX) || !isnormal(DBL_MIN) || isnormal(DBL_TRUE_MIN)) {
      ++e;
      Kokkos::printf("failed isnormal(floating_point) special values\n");
    }

    KOKKOS_TEST_WORKAROUND_DEPRECATED_STD_ITERATOR_WARNINGS_PUSH()
    static_assert(std::is_same_v<decltype(isnormal(1)), bool>);
    static_assert(std::is_same_v<decltype(isnormal(2.f)), bool>);
    static_assert(std::is_same_v<decltype(isnormal(3.)), bool>);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    static_assert(std::is_same_v<decltype(isnormal(4.l)), bool>);
#endif

    KOKKOS_TEST_STATIC_ASSERT_UNARY_PREDICATE(Kokkos::isnormal, float, bool);
    KOKKOS_TEST_STATIC_ASSERT_UNARY_PREDICATE(Kokkos::isnormal, double, bool);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    KOKKOS_TEST_STATIC_ASSERT_UNARY_PREDICATE(Kokkos::isnormal, long double,
                                              bool);
#endif
    KOKKOS_TEST_WORKAROUND_DEPRECATED_STD_ITERATOR_WARNINGS_POP()
  }
};

TEST(TEST_CATEGORY, mathematical_functions_isnormal) {
#if __FINITE_MATH_ONLY__
  GTEST_SKIP() << "skipping when compiling with -ffinite-math-only";
#endif
  TestIsNormal<TEST_EXECSPACE>();
}

template <class Space>
struct TestSignbit {
  TestSignbit() { run(); }
  void run() const {
    int errors = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Space>(0, 1), *this, errors);
    ASSERT_EQ(errors, 0);
  }
  KOKKOS_FUNCTION void operator()(int, int& e) const {
    using KE::denorm_min;
    using KE::finite_max;
    using KE::finite_min;
    using KE::infinity;
    using KE::quiet_NaN;
    using KE::signaling_NaN;
    using Kokkos::signbit;
    if (signbit(1) || signbit(INT_MAX) || !signbit(-2) || !signbit(INT_MIN) ||
        signbit(0)) {
      ++e;
      Kokkos::printf("failed signbit(integral)\n");
    }
    if (signbit(3.f) || signbit(finite_max<float>::value) ||
        signbit(infinity<float>::value) || signbit(denorm_min<float>::value) ||
        signbit(quiet_NaN<float>::value) ||
        signbit(signaling_NaN<float>::value) || signbit(0.f) ||
        !signbit(-0.4f) || !signbit(finite_min<float>::value) ||
        !signbit(-infinity<float>::value) ||
        !signbit(-denorm_min<float>::value) ||
        !signbit(-quiet_NaN<float>::value) ||
        !signbit(-signaling_NaN<float>::value) || !signbit(-0.f)) {
      ++e;
      Kokkos::printf("failed signbit(float)\n");
    }
#if !(defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_COMPILER_MSVC))
    if (signbit(static_cast<KE::half_t>(0.f)) ||
        !signbit(static_cast<KE::half_t>(-0.f))
#if !(defined(KOKKOS_ENABLE_CUDA) &&                         \
      defined(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE) && \
      defined(KOKKOS_COMPILER_CLANG))
        // FIXME internal compiler error for Clang+Cuda and RDC
        || signbit(finite_max<KE::half_t>::value) ||
        signbit(infinity<KE::half_t>::value) ||
        signbit(denorm_min<KE::half_t>::value) ||
        signbit(quiet_NaN<KE::half_t>::value) ||
        signbit(signaling_NaN<KE::half_t>::value) ||
        !signbit(finite_min<KE::half_t>::value) ||
        !signbit(-static_cast<KE::half_t>(infinity<KE::half_t>::value)) ||
        !signbit(-static_cast<KE::half_t>(denorm_min<KE::half_t>::value))
    // https://docs.nvidia.com/cuda/cuda-programming-guide/05-appendices/mathematical-functions.html#cuda-and-ieee-754-compliance:
    // "[...] result in the sign of a NaN being updated in an
    // implementation-defined manner."
#ifndef KOKKOS_ENABLE_CUDA
        || !signbit(-static_cast<KE::half_t>(quiet_NaN<KE::half_t>::value)) ||
        !signbit(-static_cast<KE::half_t>(signaling_NaN<KE::half_t>::value))
#endif
#endif
    ) {
      ++e;
      Kokkos::printf("failed signbit(KE::half_t)\n");
    }
    if (signbit(static_cast<KE::bhalf_t>(0.f)) ||
        signbit(finite_max<KE::bhalf_t>::value) ||
        signbit(infinity<KE::bhalf_t>::value) ||
        signbit(denorm_min<KE::bhalf_t>::value) ||
        signbit(quiet_NaN<KE::bhalf_t>::value) ||
        signbit(signaling_NaN<KE::bhalf_t>::value) ||
        !signbit(static_cast<KE::bhalf_t>(-0.f)) ||
        !signbit(finite_min<KE::bhalf_t>::value) ||
        !signbit(-static_cast<KE::bhalf_t>(infinity<KE::bhalf_t>::value)) ||
        !signbit(-static_cast<KE::bhalf_t>(denorm_min<KE::bhalf_t>::value))
// the bhalf test also fails for SYCL+Cuda
#ifndef KOKKOS_IMPL_ARCH_NVIDIA_GPU
        || !signbit(-static_cast<KE::bhalf_t>(quiet_NaN<KE::bhalf_t>::value)) ||
        !signbit(-static_cast<KE::bhalf_t>(signaling_NaN<KE::bhalf_t>::value))
#endif
    ) {
      ++e;
      Kokkos::printf("failed signbit(KE::bhalf_t)\n");
    }
#endif
    if (signbit(.5) || signbit(finite_max<double>::value) ||
        signbit(infinity<double>::value) ||
        signbit(denorm_min<double>::value) ||
        signbit(quiet_NaN<double>::value) ||
        signbit(signaling_NaN<double>::value) || signbit(0.) || !signbit(-6.) ||
        !signbit(finite_min<double>::value) ||
        !signbit(-infinity<double>::value) ||
        !signbit(-denorm_min<double>::value) ||
        !signbit(-quiet_NaN<double>::value) ||
        !signbit(-signaling_NaN<double>::value) || !signbit(-0.)) {
      ++e;
      Kokkos::printf("failed signbit(double)\n");
    }
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    if (signbit(7.l) || signbit(finite_max<long double>::value) ||
        signbit(infinity<long double>::value) ||
        signbit(denorm_min<long double>::value) ||
        signbit(quiet_NaN<long double>::value) ||
        signbit(signaling_NaN<long double>::value) || signbit(0.l) ||
        !signbit(-.8l) || !signbit(finite_min<long double>::value) ||
        !signbit(-infinity<long double>::value) ||
        !signbit(-denorm_min<long double>::value) ||
        !signbit(-quiet_NaN<long double>::value) ||
        !signbit(-signaling_NaN<long double>::value) || !signbit(-0.l)) {
      ++e;
      Kokkos::printf("failed signbit(long double)\n");
    }
#endif
    // special values
    if (signbit(INFINITY) || signbit(NAN)) {
      ++e;
      Kokkos::printf("failed signbit(floating_point) special values\n");
    }

    KOKKOS_TEST_WORKAROUND_DEPRECATED_STD_ITERATOR_WARNINGS_PUSH()
    static_assert(std::is_same_v<decltype(signbit(1)), bool>);
    static_assert(std::is_same_v<decltype(signbit(2.f)), bool>);
    static_assert(std::is_same_v<decltype(signbit(3.)), bool>);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    static_assert(std::is_same_v<decltype(signbit(4.l)), bool>);
#endif

    KOKKOS_TEST_STATIC_ASSERT_UNARY_PREDICATE(Kokkos::signbit, float, bool);
    KOKKOS_TEST_STATIC_ASSERT_UNARY_PREDICATE(Kokkos::signbit, double, bool);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    KOKKOS_TEST_STATIC_ASSERT_UNARY_PREDICATE(Kokkos::signbit, long double,
                                              bool);
#endif
    KOKKOS_TEST_WORKAROUND_DEPRECATED_STD_ITERATOR_WARNINGS_POP()
  }
};

TEST(TEST_CATEGORY, mathematical_functions_signbit) {
#if __FINITE_MATH_ONLY__
  GTEST_SKIP() << "skipping when compiling with -ffinite-math-only";
#endif
  TestSignbit<TEST_EXECSPACE>();
}

TEST(TEST_CATEGORY, mathematical_functions_binary_predicates) {
  auto test_all_predicates = [](auto x, auto y) {
    do_test_math_binary_predicate<TEST_EXECSPACE, kk_isgreater>(x, y);
    do_test_math_binary_predicate<TEST_EXECSPACE, kk_isgreaterequal>(x, y);
    do_test_math_binary_predicate<TEST_EXECSPACE, kk_isless>(x, y);
    do_test_math_binary_predicate<TEST_EXECSPACE, kk_islessequal>(x, y);
    do_test_math_binary_predicate<TEST_EXECSPACE, kk_islessgreater>(x, y);
    do_test_math_binary_predicate<TEST_EXECSPACE, kk_isunordered>(x, y);
  };

  test_all_predicates(2.f, 3.f);
  test_all_predicates(2., 3.);
  test_all_predicates(2, 3.f);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  test_all_predicates(2.l, 3.l);
#endif
}

KE::half_t ref_test_fallback_half(KE::half_t) {
#if defined(KOKKOS_ENABLE_SYCL) && defined(KOKKOS_IMPL_SYCL_HALF_TYPE_DEFINED)
  // When SYCL is enabled, half_t is available on both the GPU and the CPU.
  return KE::half_t(0.f);
#elif defined(KOKKOS_ENABLE_CUDA)
  if constexpr (std::is_same_v<TEST_EXECSPACE, Kokkos::Cuda>) {
    return KE::half_t(0.f);
  } else {
    return KE::half_t(1.f);
  }
#elif defined(KOKKOS_ENABLE_HIP)
  if constexpr (std::is_same_v<TEST_EXECSPACE, Kokkos::HIP>) {
    return KE::half_t(0.f);
  } else {
    return KE::half_t(1.f);
  }
#else
  return KE::half_t(1.f);
#endif
}

KE::bhalf_t ref_test_fallback_bhalf(KE::bhalf_t) {
#if defined(KOKKOS_ENABLE_SYCL) && defined(KOKKOS_IMPL_SYCL_BHALF_TYPE_DEFINED)
  // When SYCL is enabled, bhalf_t is available on both the GPU and the CPU.
  return KE::bhalf_t(0.f);
#elif defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_IMPL_ARCH_NVIDIA_GPU) && \
    (KOKKOS_IMPL_ARCH_NVIDIA_GPU >= 80)
  // bhalf_t support for CUDA is only available starting with Ampere (80)
  if constexpr (std::is_same_v<TEST_EXECSPACE, Kokkos::Cuda>) {
    return KE::bhalf_t(0.f);
  } else {
    return KE::bhalf_t(1.f);
  }
#elif defined(KOKKOS_ENABLE_HIP)
  // bhalf_t is supported on host and device for HIP but the mathematical
  // functions only have have a native implementation on the device
  if constexpr (std::is_same_v<TEST_EXECSPACE, Kokkos::HIP>) {
    return KE::bhalf_t(0.f);
  } else {
    return KE::bhalf_t(1.f);
  }
#else
  return KE::bhalf_t(1.f);
#endif
}

DEFINE_UNARY_FUNCTION_EVAL_CUSTOM(test_fallback_half, 0,
                                  ref_test_fallback_half(x));
DEFINE_UNARY_FUNCTION_EVAL_CUSTOM(test_fallback_bhalf, 0,
                                  ref_test_fallback_bhalf(x));

TEST(TEST_CATEGORY, mathematical_functions_impl_half_fallback) {
  TestMathUnaryFunction<TEST_EXECSPACE, MathUnaryFunction_test_fallback_half,
                        KE::half_t, 1>({KE::half_t(1.f)});
  TestMathUnaryFunction<TEST_EXECSPACE, MathUnaryFunction_test_fallback_bhalf,
                        KE::bhalf_t, 1>({KE::bhalf_t(1.f)});
}

template <class Space, class FP16Type>
struct TestNextAfterHalf {
  TestNextAfterHalf() { run(); }
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
    using Kokkos::nextafter;

    // Define useful constants
    const std::uint16_t FP16_POS_ZERO     = 0x0000;
    const std::uint16_t FP16_NEG_ZERO     = 0x8000;
    const std::uint16_t FP16_SMALLEST_POS = 0x0001;
    const std::uint16_t FP16_SMALLEST_NEG = 0x8001;

    const FP16Type pos_one{1.0f}, pos_two{2.0f};
    const FP16Type neg_one{-1.0f}, neg_two{-2.0f};
    const FP16Type pos_zero     = Kokkos::bit_cast<FP16Type>(FP16_POS_ZERO);
    const FP16Type neg_zero     = Kokkos::bit_cast<FP16Type>(FP16_NEG_ZERO);
    const FP16Type pos_smallest = Kokkos::bit_cast<FP16Type>(FP16_SMALLEST_POS);
    const FP16Type neg_smallest = Kokkos::bit_cast<FP16Type>(FP16_SMALLEST_NEG);
    const FP16Type pos_max = Kokkos::Experimental::finite_max<FP16Type>::value;
    const FP16Type neg_max = Kokkos::Experimental::finite_min<FP16Type>::value;
    const FP16Type pos_inf = Kokkos::Experimental::infinity<FP16Type>::value;
    const FP16Type neg_inf =
        -static_cast<FP16Type>(Kokkos::Experimental::infinity<FP16Type>::value);

    // NaN Handling
    if (!isnan(nextafter(quiet_NaN<FP16Type>::value, pos_one)) ||
        !isnan(nextafter(signaling_NaN<FP16Type>::value, pos_one)) ||
        !isnan(nextafter(pos_one, quiet_NaN<FP16Type>::value)) ||
        !isnan(nextafter(pos_one, signaling_NaN<FP16Type>::value)) ||
        !isnan(nextafter(quiet_NaN<FP16Type>::value,
                         quiet_NaN<FP16Type>::value)) ||
        !isnan(nextafter(quiet_NaN<FP16Type>::value,
                         signaling_NaN<FP16Type>::value)) ||
        !isnan(nextafter(signaling_NaN<FP16Type>::value,
                         quiet_NaN<FP16Type>::value)) ||
        !isnan(nextafter(signaling_NaN<FP16Type>::value,
                         signaling_NaN<FP16Type>::value))) {
      ++e;
      Kokkos::printf("failed half precision nextafter(NaN)\n");
    }

    // Equality (from==toward) Handling
    if (nextafter(pos_one, pos_one) != pos_one ||
        nextafter(pos_zero, pos_zero) != pos_zero ||
        nextafter(neg_zero, neg_zero) != neg_zero ||
        nextafter(pos_inf, pos_inf) != pos_inf ||
        nextafter(neg_inf, neg_inf) != neg_inf) {
      ++e;
      Kokkos::printf("failed half precision nextafter(equality)\n");
    }

    // Zero Handling
    if (nextafter(pos_zero, pos_one) != pos_smallest ||
        nextafter(pos_zero, neg_one) != neg_smallest ||
        nextafter(pos_zero, neg_zero) != neg_zero ||
        nextafter(neg_zero, pos_one) != pos_smallest ||
        nextafter(neg_zero, neg_one) != neg_smallest ||
        nextafter(neg_zero, pos_zero) != pos_zero) {
      ++e;
      Kokkos::printf("failed half precision nextafter(zero)\n");
    }

    // From Negative Non Zero Handling
    const FP16Type after_neg_one = Kokkos::bit_cast<FP16Type>(
        std::uint16_t(Kokkos::bit_cast<std::uint16_t>(neg_one) - 1));
    const FP16Type before_neg_one = Kokkos::bit_cast<FP16Type>(
        std::uint16_t(Kokkos::bit_cast<std::uint16_t>(neg_one) + 1));
    if (nextafter(neg_smallest, pos_zero) != neg_zero ||
        nextafter(neg_one, pos_one) != after_neg_one ||
        nextafter(neg_one, neg_two) != before_neg_one ||
        nextafter(neg_max, neg_inf) != neg_inf) {
      ++e;
      Kokkos::printf("failed half precision nextafter(negative)\n");
    }

    // From Positive Non Zero Handling
    const FP16Type after_pos_one = Kokkos::bit_cast<FP16Type>(
        std::uint16_t(Kokkos::bit_cast<std::uint16_t>(pos_one) + 1));
    const FP16Type before_pos_one = Kokkos::bit_cast<FP16Type>(
        std::uint16_t(Kokkos::bit_cast<std::uint16_t>(pos_one) - 1));
    if (nextafter(pos_smallest, neg_zero) != pos_zero ||
        nextafter(pos_one, neg_one) != before_pos_one ||
        nextafter(pos_one, pos_two) != after_pos_one ||
        nextafter(pos_max, pos_inf) != pos_inf) {
      ++e;
      Kokkos::printf("failed half precision nextafter(positive)\n");
    }

    // From Inf Handling
    // Note: The behavior of nextafter with infinities is
    // implementation-defined, but in Kokkos it returns the maximum
    // finite value when moving towards a finite value.
    if (nextafter(pos_inf, pos_one) != pos_max ||
        nextafter(neg_inf, neg_one) != neg_max ||
        nextafter(pos_inf, pos_inf) != pos_inf ||
        nextafter(neg_inf, neg_inf) != neg_inf) {
      ++e;
      Kokkos::printf("failed half precision nextafter(inf)\n");
    }
  }
};

TEST(TEST_CATEGORY, mathematical_functions_nextafter_fp16) {
#if defined(KOKKOS_ENABLE_CUDA) &&                         \
    defined(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE) && \
    defined(KOKKOS_COMPILER_CLANG)
  GTEST_SKIP() << "FIXME internal compiler error for Clang+Cuda and RDC";
#else
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_COMPILER_MSVC)
  GTEST_SKIP() << "FIXME MSVC nextafter for half precision "
                  "not implemented yet";
#else
  bool skipped = true;
#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
  skipped      = false;
  TestNextAfterHalf<TEST_EXECSPACE, Kokkos::Experimental::half_t>();
#endif
#if defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT
  skipped = false;
  TestNextAfterHalf<TEST_EXECSPACE, Kokkos::Experimental::bhalf_t>();
#endif
  if (skipped) GTEST_SKIP() << "no 16-bit floating-point precision support";
#endif
#endif
}
#endif

template <class T>
struct TestNextToward {
  TestNextToward() {
#if __FINITE_MATH_ONLY__
    constexpr bool finite_math_only = true;
#else
    constexpr bool finite_math_only = false;
#endif
    if constexpr (std::is_integral_v<T>) {
      test_integral<finite_math_only>();
    } else {
      test_float<finite_math_only>();
    }
  }

  template <bool finite_math_only>
  void test_float() const {
    using Kokkos::isnan;
    using Kokkos::nexttoward;

    const T pos_one{1.0f}, neg_one{-1.0f};
    const T pos_zero{0.0f}, neg_zero{-0.0f};
    const T pos_smallest = KE::denorm_min<T>::value;
    const T neg_smallest = -static_cast<T>(KE::denorm_min<T>::value);
    const T pos_max      = KE::finite_max<T>::value;
    const T neg_max      = KE::finite_min<T>::value;
    const T pos_inf      = KE::infinity<T>::value;
    const T neg_inf      = -static_cast<T>(KE::infinity<T>::value);

    long double target_pos_one{1.0l}, target_pos_two{2.0l};
    long double target_neg_one{-1.0l}, target_neg_two{-2.0l};
    long double target_pos_zero{0.0l}, target_neg_zero{-0.0l};
    long double target_pos_inf{KE::infinity<long double>::value};
    long double target_neg_inf{-KE::infinity<long double>::value};

    // Check return type
    testing::StaticAssertTypeEq<decltype(nexttoward(pos_one, target_pos_one)),
                                T>();

    // NaN Handling. If finite-math is enabled, skip this
    if constexpr (!finite_math_only) {
      std::string nan_err_msg =
          "failed " +
          static_cast<std::string>(Kokkos::Impl::TypeInfo<T>::name()) +
          " type nexttoward(NaN)";

      EXPECT_TRUE(isnan(nexttoward(KE::quiet_NaN<T>::value, target_pos_one)))
          << nan_err_msg;
      EXPECT_TRUE(
          isnan(nexttoward(KE::signaling_NaN<T>::value, target_pos_one)))
          << nan_err_msg;
      EXPECT_TRUE(isnan(nexttoward(pos_one, KE::quiet_NaN<long double>::value)))
          << nan_err_msg;
      EXPECT_TRUE(
          isnan(nexttoward(pos_one, KE::signaling_NaN<long double>::value)))
          << nan_err_msg;
      EXPECT_TRUE(isnan(nexttoward(KE::quiet_NaN<T>::value,
                                   KE::quiet_NaN<long double>::value)))
          << nan_err_msg;
      EXPECT_TRUE(isnan(nexttoward(KE::quiet_NaN<T>::value,
                                   KE::signaling_NaN<long double>::value)))
          << nan_err_msg;
      EXPECT_TRUE(isnan(nexttoward(KE::signaling_NaN<T>::value,
                                   KE::quiet_NaN<long double>::value)))
          << nan_err_msg;
      EXPECT_TRUE(isnan(nexttoward(KE::signaling_NaN<T>::value,
                                   KE::signaling_NaN<long double>::value)))
          << nan_err_msg;
    }

    // Zero Handling
    std::string zero_err_msg =
        "failed " +
        static_cast<std::string>(Kokkos::Impl::TypeInfo<T>::name()) +
        " type nexttoward(zero)";
    EXPECT_EQ(nexttoward(pos_zero, target_pos_one), pos_smallest)
        << zero_err_msg;
    EXPECT_EQ(nexttoward(pos_zero, target_neg_one), neg_smallest)
        << zero_err_msg;
    EXPECT_EQ(nexttoward(pos_zero, target_neg_zero), neg_zero) << zero_err_msg;
    EXPECT_EQ(nexttoward(neg_zero, target_pos_one), pos_smallest)
        << zero_err_msg;
    EXPECT_EQ(nexttoward(neg_zero, target_neg_one), neg_smallest)
        << zero_err_msg;
    EXPECT_EQ(nexttoward(neg_zero, target_pos_zero), pos_zero) << zero_err_msg;

    // From Negative Non Zero Handling
    std::string negative_err_msg =
        "failed " +
        static_cast<std::string>(Kokkos::Impl::TypeInfo<T>::name()) +
        " type nexttoward(negative)";
    const T after_neg_one  = Kokkos::nextafter(neg_one, pos_inf);
    const T before_neg_one = Kokkos::nextafter(neg_one, neg_inf);
    EXPECT_EQ(nexttoward(neg_smallest, target_pos_zero), neg_zero)
        << negative_err_msg;
    EXPECT_EQ(nexttoward(neg_one, target_pos_one), after_neg_one)
        << negative_err_msg;
    EXPECT_EQ(nexttoward(neg_one, target_neg_two), before_neg_one)
        << negative_err_msg;

    // From Positive Non Zero Handling
    std::string positive_err_msg =
        "failed " +
        static_cast<std::string>(Kokkos::Impl::TypeInfo<T>::name()) +
        " type nexttoward(positive)";
    const T after_pos_one  = Kokkos::nextafter(pos_one, pos_inf);
    const T before_pos_one = Kokkos::nextafter(pos_one, neg_inf);
    EXPECT_EQ(nexttoward(pos_smallest, target_neg_zero), pos_zero)
        << positive_err_msg;
    EXPECT_EQ(nexttoward(pos_one, target_neg_one), before_pos_one)
        << positive_err_msg;
    EXPECT_EQ(nexttoward(pos_one, target_pos_two), after_pos_one)
        << positive_err_msg;

    // From Inf Handling. If finite-math is enabled, skip this
    if constexpr (!finite_math_only) {
      std::string inf_err_msg =
          "failed " +
          static_cast<std::string>(Kokkos::Impl::TypeInfo<T>::name()) +
          " type nexttoward(inf)";
      EXPECT_EQ(nexttoward(neg_max, target_neg_inf), neg_inf) << inf_err_msg;
      EXPECT_EQ(nexttoward(pos_max, target_pos_inf), pos_inf) << inf_err_msg;
      EXPECT_EQ(nexttoward(pos_inf, target_pos_one), pos_max) << inf_err_msg;
      EXPECT_EQ(nexttoward(neg_inf, target_neg_one), neg_max) << inf_err_msg;
      EXPECT_EQ(nexttoward(pos_inf, target_pos_inf), pos_inf) << inf_err_msg;
      EXPECT_EQ(nexttoward(neg_inf, target_neg_inf), neg_inf) << inf_err_msg;
    }
  }

  // As integral types do not support NaNs/Infinity/Denorms,
  // we skip tests related to these values
  template <bool finite_math_only>
  void test_integral() const {
    using Kokkos::nexttoward;

    // Since T is an integral type, we need to declare input constants in
    // integral type, but reference constants are in double
    const T pos_one{1}, neg_one{-1};
    const T pos_zero{0}, neg_zero{-0};

    const double ref_pos_one{1.0}, ref_pos_two{2.0};
    const double ref_neg_one{-1.0}, ref_neg_two{-2.0};
    const double ref_pos_zero{0.0}, ref_neg_zero{-0.0};
    const double ref_pos_smallest{std::numeric_limits<double>::denorm_min()};
    const double ref_neg_smallest{-std::numeric_limits<double>::denorm_min()};

    const long double target_pos_one{1.0l}, target_pos_two{2.0l};
    const long double target_neg_one{-1.0l}, target_neg_two{-2.0l};
    const long double target_pos_zero{0.0l}, target_neg_zero{-0.0l};

    // Check return type. For integral types, the return type is double
    testing::StaticAssertTypeEq<decltype(nexttoward(pos_one, target_pos_one)),
                                double>();

    // NaN Handling
    if constexpr (!finite_math_only) {
      std::string nan_err_msg =
          "failed " +
          static_cast<std::string>(Kokkos::Impl::TypeInfo<T>::name()) +
          " type nexttoward(NaN)";
      EXPECT_TRUE(std::isnan(
          nexttoward(pos_one, std::numeric_limits<long double>::quiet_NaN())))
          << nan_err_msg;
      EXPECT_TRUE(std::isnan(nexttoward(
          pos_one, std::numeric_limits<long double>::signaling_NaN())))
          << nan_err_msg;
    }

    // Zero Handling
    std::string zero_err_msg =
        "failed " +
        static_cast<std::string>(Kokkos::Impl::TypeInfo<T>::name()) +
        " type nexttoward(zero)";
    EXPECT_EQ(nexttoward(pos_zero, target_pos_one), ref_pos_smallest)
        << zero_err_msg;
    EXPECT_EQ(nexttoward(pos_zero, target_neg_one), ref_neg_smallest)
        << zero_err_msg;
    EXPECT_EQ(nexttoward(pos_zero, target_neg_zero), ref_neg_zero)
        << zero_err_msg;
    EXPECT_EQ(nexttoward(neg_zero, target_pos_one), ref_pos_smallest)
        << zero_err_msg;
    EXPECT_EQ(nexttoward(neg_zero, target_neg_one), ref_neg_smallest)
        << zero_err_msg;
    EXPECT_EQ(nexttoward(neg_zero, target_pos_zero), ref_pos_zero)
        << zero_err_msg;

    // From Negative Non Zero Handling
    std::string negative_err_msg =
        "failed " +
        static_cast<std::string>(Kokkos::Impl::TypeInfo<T>::name()) +
        " type nexttoward(negative)";
    const double after_neg_one  = std::nextafter(ref_neg_one, ref_pos_one);
    const double before_neg_one = std::nextafter(ref_neg_one, ref_neg_two);
    EXPECT_EQ(nexttoward(neg_one, target_pos_one), after_neg_one)
        << negative_err_msg;
    EXPECT_EQ(nexttoward(neg_one, target_neg_two), before_neg_one)
        << negative_err_msg;

    // From Positive Non Zero Handling
    std::string positive_err_msg =
        "failed " +
        static_cast<std::string>(Kokkos::Impl::TypeInfo<T>::name()) +
        " type nexttoward(positive)";
    const double after_pos_one  = std::nextafter(ref_pos_one, ref_pos_two);
    const double before_pos_one = std::nextafter(ref_pos_one, ref_neg_one);
    EXPECT_EQ(nexttoward(pos_one, target_neg_one), before_pos_one)
        << positive_err_msg;
    EXPECT_EQ(nexttoward(pos_one, target_pos_two), after_pos_one)
        << positive_err_msg;

    // If finite-math is enabled, skip this
    // Special treatments for max (integer)
    // Max (integer) corresponds to a big but finite number in double
    // Thus, nexttoward from max (integer) to inf should return a slightly
    // bigger number in double
    if constexpr (!finite_math_only) {
      const T pos_max{std::numeric_limits<T>::max()};
      const T neg_max{-std::numeric_limits<T>::max()};
      double pos_max_as_double = static_cast<double>(pos_max);
      double neg_max_as_double = static_cast<double>(neg_max);
      const long double target_pos_inf{
          std::numeric_limits<long double>::infinity()};
      const long double target_neg_inf{
          -std::numeric_limits<long double>::infinity()};
      const double before_neg_max = std::nextafter(
          neg_max_as_double, static_cast<double>(target_neg_inf));
      const double after_pos_max = std::nextafter(
          pos_max_as_double, static_cast<double>(target_pos_inf));
      std::string max_err_msg =
          "failed " +
          static_cast<std::string>(Kokkos::Impl::TypeInfo<T>::name()) +
          " type nexttoward(max)";
      EXPECT_EQ(nexttoward(neg_max, target_neg_inf), before_neg_max)
          << max_err_msg;
      EXPECT_EQ(nexttoward(pos_max, target_pos_inf), after_pos_max)
          << max_err_msg;
    }
  }
};

TEST(TEST_CATEGORY, mathematical_functions_nexttoward) {
  if constexpr (std::is_same_v<TEST_EXECSPACE,
                               Kokkos::DefaultHostExecutionSpace>) {
    TestNextToward<int>();
    TestNextToward<float>();
    TestNextToward<double>();
    TestNextToward<long double>();
  } else {
    GTEST_SKIP() << "nexttoward only defined for host";
  }
}

#endif
