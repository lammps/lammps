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
#include <algorithm>
#include <initializer_list>
#include <type_traits>
#include "Kokkos_ExecPolicy.hpp"
#include "Kokkos_Parallel_Reduce.hpp"

#include <cfloat>

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || \
    defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ENABLE_OPENMPTARGET)
#else
#define MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
#endif

// WORKAROUND icpx changing default FP model when optimization level is >= 1
// using -fp-model=precise works too
#if defined(__INTEL_LLVM_COMPILER)
#define KOKKOS_IMPL_WORKAROUND_INTEL_LLVM_DEFAULT_FLOATING_POINT_MODEL
#endif

// clang-format off
template <class>
struct math_unary_function_return_type;
// Floating-point types
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

struct FloatingPointComparison {
 private:
  template <class T>
  KOKKOS_FUNCTION double eps(T) const {
    return DBL_EPSILON;
  }
  KOKKOS_FUNCTION
  double eps(float) const { return FLT_EPSILON; }
  KOKKOS_FUNCTION
  double eps(long double) const { return LDBL_EPSILON; }

  // Using absolute here instead of abs, since we actually test abs ...
  template <class T>
  KOKKOS_FUNCTION typename std::enable_if<std::is_signed<T>::value, T>::type
  absolute(T val) const {
    return val < T(0) ? -val : val;
  }

  template <class T>
  KOKKOS_FUNCTION typename std::enable_if<!std::is_signed<T>::value, T>::type
  absolute(T val) const {
    return val;
  }

 public:
  template <class FPT>
  KOKKOS_FUNCTION bool compare_near_zero(FPT const& fpv, double ulp) const {
    auto abs_tol = eps(fpv) * ulp;

    bool ar = absolute(fpv) < abs_tol;
    if (!ar) {
#if !defined(KOKKOS_ENABLE_SYCL) && !defined(KOKKOS_ENABLE_HIP)
      printf("absolute value exceeds tolerance [|%e| > %e]\n", (double)fpv,
             abs_tol);
#endif
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
      bool ar         = rel_diff < rel_tol;
      if (!ar) {
#if !defined(KOKKOS_ENABLE_SYCL) && !defined(KOKKOS_ENABLE_HIP)
        printf("relative difference exceeds tolerance [%e > %e]\n",
               (double)rel_diff, rel_tol);
#endif
      }

      return ar;
    }
  }
};

template <class>
struct math_function_name;

#define DEFINE_UNARY_FUNCTION_EVAL(FUNC, ULP_FACTOR)                           \
  struct MathUnaryFunction_##FUNC {                                            \
    template <typename T>                                                      \
    static KOKKOS_FUNCTION auto eval(T x) {                                    \
      static_assert(std::is_same<decltype(Kokkos::Experimental::FUNC((T)0)),   \
                                 math_unary_function_return_type_t<T>>::value, \
                    "");                                                       \
      return Kokkos::Experimental::FUNC(x);                                    \
    }                                                                          \
    template <typename T>                                                      \
    static auto eval_std(T x) {                                                \
      static_assert(std::is_same<decltype(std::FUNC((T)0)),                    \
                                 math_unary_function_return_type_t<T>>::value, \
                    "");                                                       \
      return std::FUNC(x);                                                     \
    }                                                                          \
    static KOKKOS_FUNCTION double ulp_factor() { return ULP_FACTOR; }          \
  };                                                                           \
  using kk_##FUNC = MathUnaryFunction_##FUNC;                                  \
  template <>                                                                  \
  struct math_function_name<MathUnaryFunction_##FUNC> {                        \
    static constexpr char name[] = #FUNC;                                      \
  };                                                                           \
  constexpr char math_function_name<MathUnaryFunction_##FUNC>::name[]

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

DEFINE_UNARY_FUNCTION_EVAL(sqrt, 2);
DEFINE_UNARY_FUNCTION_EVAL(cbrt, 2);

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
#ifndef KOKKOS_ENABLE_SYCL
DEFINE_UNARY_FUNCTION_EVAL(nearbyint, 2);
#endif

#undef DEFINE_UNARY_FUNCTION_EVAL

#define DEFINE_BINARY_FUNCTION_EVAL(FUNC, ULP_FACTOR)                    \
  struct MathBinaryFunction_##FUNC {                                     \
    template <typename T, typename U>                                    \
    static KOKKOS_FUNCTION auto eval(T x, U y) {                         \
      static_assert(                                                     \
          std::is_same<decltype(Kokkos::Experimental::FUNC((T)0, (U)0)), \
                       math_binary_function_return_type_t<T, U>>::value, \
          "");                                                           \
      return Kokkos::Experimental::FUNC(x, y);                           \
    }                                                                    \
    template <typename T, typename U>                                    \
    static auto eval_std(T x, U y) {                                     \
      static_assert(                                                     \
          std::is_same<decltype(std::FUNC((T)0, (U)0)),                  \
                       math_binary_function_return_type_t<T, U>>::value, \
          "");                                                           \
      return std::FUNC(x, y);                                            \
    }                                                                    \
    static KOKKOS_FUNCTION double ulp_factor() { return ULP_FACTOR; }    \
  };                                                                     \
  using kk_##FUNC = MathBinaryFunction_##FUNC;                           \
  template <>                                                            \
  struct math_function_name<MathBinaryFunction_##FUNC> {                 \
    static constexpr char name[] = #FUNC;                                \
  };                                                                     \
  constexpr char math_function_name<MathBinaryFunction_##FUNC>::name[]

DEFINE_BINARY_FUNCTION_EVAL(pow, 2);
DEFINE_BINARY_FUNCTION_EVAL(hypot, 2);

#undef DEFINE_BINARY_FUNCTION_EVAL

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
#if !defined(KOKKOS_ENABLE_SYCL) && !defined(KOKKOS_ENABLE_HIP)
      printf("value at %f which is %f was expected to be %f\n", (double)val_[i],
             (double)Func::eval(val_[i]), (double)res_[i]);
#endif
    }
  }
};

template <class Space, class... Func, class Arg, std::size_t N>
void do_test_math_unary_function(const Arg (&x)[N]) {
  (void)std::initializer_list<int>{
      (TestMathUnaryFunction<Space, Func, Arg, N>(x), 0)...};
}

#define TEST_MATH_FUNCTION(FUNC) \
  do_test_math_unary_function<TEST_EXECSPACE, MathUnaryFunction_##FUNC>

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
#if !defined(KOKKOS_ENABLE_SYCL) && !defined(KOKKOS_ENABLE_HIP)
      printf("value at %f, %f which is %f was expected to be %f\n",
             (double)val1_, (double)val2_, (double)Func::eval(val1_, val2_),
             (double)res_);
#endif
    }
  }
};

template <class Space, class... Func, class Arg1, class Arg2>
void do_test_math_binary_function(Arg1 arg1, Arg2 arg2) {
  (void)std::initializer_list<int>{
      (TestMathBinaryFunction<Space, Func, Arg1, Arg2>(arg1, arg2), 0)...};
}

TEST(TEST_CATEGORY, mathematical_functions_trigonometric_functions) {
  TEST_MATH_FUNCTION(sin)({true, false});
  TEST_MATH_FUNCTION(sin)({-3, -2, -1, 0, 1});
  TEST_MATH_FUNCTION(sin)({-3l, -2l, -1l, 0l, 1l});
  TEST_MATH_FUNCTION(sin)({-3ll, -2ll, -1ll, 0ll, 1ll});
  TEST_MATH_FUNCTION(sin)({2u, 3u, 4u, 5u, 6u});
  TEST_MATH_FUNCTION(sin)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_MATH_FUNCTION(sin)({2ull, 3ull, 4ull, 5ull, 6ull});
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
  TEST_MATH_FUNCTION(atan)({-1.5f, 1.3f, -1.1f, .9f, -.7f, .5f});
  TEST_MATH_FUNCTION(atan)({1.4, -1.2, 1., -.8, .6, -.4, .2, -0.});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(atan)({-.98l, .67l, -54.l, .34l, -.21l});
#endif

  // TODO atan2
}

TEST(TEST_CATEGORY, mathematical_functions_power_functions) {
  TEST_MATH_FUNCTION(sqrt)({0, 1, 2, 3, 5, 7, 11});
  TEST_MATH_FUNCTION(sqrt)({0l, 1l, 2l, 3l, 5l, 7l, 11l});
  TEST_MATH_FUNCTION(sqrt)({0ll, 1ll, 2ll, 3ll, 5ll, 7ll, 11ll});
  TEST_MATH_FUNCTION(sqrt)({0u, 1u, 2u, 3u, 5u, 7u});
  TEST_MATH_FUNCTION(sqrt)({0ul, 1ul, 2ul, 3ul, 5ul, 7ul});
  TEST_MATH_FUNCTION(sqrt)({0ull, 1ull, 2ull, 3ull, 5ull, 7ull});
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
  TEST_MATH_FUNCTION(cbrt)({-1.f, .2f, -3.f, .4f, -5.f});
  TEST_MATH_FUNCTION(cbrt)({11.1, -2.2, 33.3, -4.4, 55.5});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(cbrt)({-10.l, 20.l, -30.l, 40.l, -50.l});
#endif

  do_test_math_binary_function<TEST_EXECSPACE, kk_pow>(2.f, 3.f);
  do_test_math_binary_function<TEST_EXECSPACE, kk_pow>(2., 3.);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  do_test_math_binary_function<TEST_EXECSPACE, kk_pow>(2.l, 3.l);
#endif

  do_test_math_binary_function<TEST_EXECSPACE, kk_hypot>(2.f, 3.f);
  do_test_math_binary_function<TEST_EXECSPACE, kk_hypot>(2., 3.);
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
// FIXME: fails with gcc on Power platforms
#if !(defined(KOKKOS_ARCH_POWER8) || defined(KOKKOS_ARCH_POWER9))
  do_test_math_binary_function<TEST_EXECSPACE, kk_hypot>(2.l, 3.l);
#endif
#endif
}

TEST(TEST_CATEGORY, mathematical_functions_exponential_functions) {
  TEST_MATH_FUNCTION(exp)({-9, -8, -7, -6, -5, 4, 3, 2, 1, 0});
  TEST_MATH_FUNCTION(exp)({-9l, -8l, -7l, -6l, -5l, 4l, 3l, 2l, 1l, 0l});
  TEST_MATH_FUNCTION(exp)({-9ll, -8ll, -7ll, -6ll, -5ll, 4ll, 3ll, 2ll, 1ll});
  TEST_MATH_FUNCTION(exp)({0u, 1u, 2u, 3u, 4u, 5u});
  TEST_MATH_FUNCTION(exp)({0ul, 1ul, 2ul, 3ul, 4ul, 5ul});
  TEST_MATH_FUNCTION(exp)({0ull, 1ull, 2ull, 3ull, 4ull, 5ull});
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
  TEST_MATH_FUNCTION(log10)({1234.f, 567.f, 89.f, .1f});
  TEST_MATH_FUNCTION(log10)({1234., 567., 89., .02});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(log10)({1234.l, 567.l, 89.l, .003l});
#endif

// FIXME_OPENMPTARGET FIXME_AMD
#if defined(KOKKOS_ENABLE_OPENMPTARGET) &&                           \
    (defined(KOKKOS_ARCH_VEGA906) || defined(KOKKOS_ARCH_VEGA908) || \
     defined(KOKKOS_ARCH_VEGA90A))

  TEST_MATH_FUNCTION(log2)({1, 23, 456, 7890});
#endif
  TEST_MATH_FUNCTION(log2)({1l, 23l, 456l, 7890l});
  TEST_MATH_FUNCTION(log2)({1ll, 23ll, 456ll, 7890ll});
  TEST_MATH_FUNCTION(log2)({1u, 23u, 456u, 7890u});
  TEST_MATH_FUNCTION(log2)({1ul, 23ul, 456ul, 7890ul});
  TEST_MATH_FUNCTION(log2)({1ull, 23ull, 456ull, 7890ull});
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
  TEST_MATH_FUNCTION(log1p)({1234.f, 567.f, 89.f, -.9f});
  TEST_MATH_FUNCTION(log1p)({1234., 567., 89., -.08});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(log1p)({1234.l, 567.l, 89.l, -.007l});
#endif
}

TEST(TEST_CATEGORY, mathematical_functions_hyperbolic_functions) {
  TEST_MATH_FUNCTION(sinh)({-3, -2, -1, 0, 1});
  TEST_MATH_FUNCTION(sinh)({-3l, -2l, -1l, 0l, 1l});
  TEST_MATH_FUNCTION(sinh)({-3ll, -2ll, -1ll, 0ll, 1ll});
  TEST_MATH_FUNCTION(sinh)({2u, 3u, 4u, 5u, 6u});
  TEST_MATH_FUNCTION(sinh)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_MATH_FUNCTION(sinh)({2ull, 3ull, 4ull, 5ull, 6ull});
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
  TEST_MATH_FUNCTION(atanh)({-.97f, .86f, -.53f, .42f, -.1f, 0.f});
  TEST_MATH_FUNCTION(atanh)({-.97, .86, -.53, .42, -.1, 0.});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(atanh)({-.97l, .86l, -.53l, .42l, -.1l, 0.l});
#endif
}

TEST(TEST_CATEGORY, mathematical_functions_error_and_gamma_functions) {
  TEST_MATH_FUNCTION(erf)({-3, -2, -1, 0, 1});
  TEST_MATH_FUNCTION(erf)({-3l, -2l, -1l, 0l, 1l});
  TEST_MATH_FUNCTION(erf)({-3ll, -2ll, -1ll, 0ll, 1ll});
  TEST_MATH_FUNCTION(erf)({2u, 3u, 4u, 5u, 6u});
  TEST_MATH_FUNCTION(erf)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_MATH_FUNCTION(erf)({2ull, 3ull, 4ull, 5ull, 6ull});
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
  TEST_MATH_FUNCTION(trunc)({-1.1f, 2.2f, -3.3f, 4.4f, -5.5f});
  TEST_MATH_FUNCTION(trunc)({-6.6, 7.7, -8.8, 9.9});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(trunc)({12.3l, 4.56l, 789.l});
#endif

#ifndef KOKKOS_ENABLE_SYCL
  TEST_MATH_FUNCTION(nearbyint)({-3, -2, -1, 0, 1});
  TEST_MATH_FUNCTION(nearbyint)({-3l, -2l, -1l, 0l, 1l});
  TEST_MATH_FUNCTION(nearbyint)({-3ll, -2ll, -1ll, 0ll, 1ll});
  TEST_MATH_FUNCTION(nearbyint)({2u, 3u, 4u, 5u, 6u});
  TEST_MATH_FUNCTION(nearbyint)({2ul, 3ul, 4ul, 5ul, 6ul});
  TEST_MATH_FUNCTION(nearbyint)({2ull, 3ull, 4ull, 5ull, 6ull});
  TEST_MATH_FUNCTION(nearbyint)({-1.1f, 2.2f, -3.3f, 4.4f, -5.5f});
  TEST_MATH_FUNCTION(nearbyint)({-6.6, 7.7, -8.8, 9.9});
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
  TEST_MATH_FUNCTION(nearbyint)({12.3l, 4.56l, 789.l});
#endif
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
    using Kokkos::Experimental::abs;
    if (abs(1) != 1 || abs(-1) != 1) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed abs(int)\n");
    }
    if (abs(2l) != 2l || abs(-2l) != 2l) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed abs(long int)\n");
    }
    if (abs(3ll) != 3ll || abs(-3ll) != 3ll) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed abs(long long int)\n");
    }
    if (abs(4.f) != 4.f || abs(-4.f) != 4.f) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed abs(float)\n");
    }
    if (abs(5.) != 5. || abs(-5.) != 5.) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed abs(double)\n");
    }
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    if (abs(6.l) != 6.l || abs(-6.l) != 6.l) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed abs(long double)\n");
    }
#endif
    // special values
    using Kokkos::Experimental::isinf;
    using Kokkos::Experimental::isnan;
    if (abs(-0.) != 0.
#ifndef KOKKOS_IMPL_WORKAROUND_INTEL_LLVM_DEFAULT_FLOATING_POINT_MODEL
        || !isinf(abs(-INFINITY)) || !isnan(abs(-NAN))
#endif
    ) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF(
          "failed abs(floating_point) special values\n");
    }

    static_assert(std::is_same<decltype(abs(1)), int>::value, "");
    static_assert(std::is_same<decltype(abs(2l)), long>::value, "");
    static_assert(std::is_same<decltype(abs(3ll)), long long>::value, "");
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
struct TestIsNaN {
  TestIsNaN() { run(); }
  void run() const {
    int errors = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Space>(0, 1), *this, errors);
    ASSERT_EQ(errors, 0);
  }
  KOKKOS_FUNCTION void operator()(int, int& e) const {
    using Kokkos::Experimental::isnan;
    using Kokkos::Experimental::quiet_NaN;
    using Kokkos::Experimental::signaling_NaN;
    if (isnan(1) || isnan(INT_MAX)) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed isnan(integral)\n");
    }
    if (isnan(2.f)
#ifndef KOKKOS_IMPL_WORKAROUND_INTEL_LLVM_DEFAULT_FLOATING_POINT_MODEL
        || !isnan(quiet_NaN<float>::value) ||
        !isnan(signaling_NaN<float>::value)
#endif

    ) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed isnan(float)\n");
    }
    if (isnan(3.)
#ifndef KOKKOS_IMPL_WORKAROUND_INTEL_LLVM_DEFAULT_FLOATING_POINT_MODEL
        || !isnan(quiet_NaN<double>::value) ||
        !isnan(signaling_NaN<double>::value)
#endif
    ) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed isnan(double)\n");
    }
#ifdef MATHEMATICAL_FUNCTIONS_HAVE_LONG_DOUBLE_OVERLOADS
    if (isnan(4.l)
#ifndef KOKKOS_IMPL_WORKAROUND_INTEL_LLVM_DEFAULT_FLOATING_POINT_MODEL
        || !isnan(quiet_NaN<long double>::value) ||
        !isnan(signaling_NaN<long double>::value)
#endif
    ) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed isnan(long double)\n");
    }
#endif
    // special values
    if (isnan(INFINITY)
#ifndef KOKKOS_IMPL_WORKAROUND_INTEL_LLVM_DEFAULT_FLOATING_POINT_MODEL
        || !isnan(NAN)
#endif
    ) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF(
          "failed isnan(floating_point) special values\n");
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
