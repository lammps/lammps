// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_QUAD_PRECISION_MATH_HPP
#define KOKKOS_QUAD_PRECISION_MATH_HPP

#include <Kokkos_Macros.hpp>

#if defined(KOKKOS_ENABLE_LIBQUADMATH)

#include <Kokkos_MathematicalConstants.hpp>
#include <Kokkos_MathematicalFunctions.hpp>

#include <quadmath.h>

#if !(defined(__FLOAT128__) || defined(__SIZEOF_FLOAT128__))
#error __float128 not supported on this host
#endif

//<editor-fold desc="Common mathematical functions __float128 overloads">
namespace Kokkos {
// clang-format off
namespace Impl {
template <> struct promote<__float128> { using type = __float128; };
}
// Basic operations
inline __float128 abs(__float128 x) { return ::fabsq(x); }
inline __float128 fabs(__float128 x) { return ::fabsq(x); }
inline __float128 fmod(__float128 x, __float128 y) { return ::fmodq(x, y); }
inline __float128 remainder(__float128 x, __float128 y) { return ::remainderq(x, y); }
// remquo
inline __float128 fma(__float128 x, __float128 y, __float128 z) { return ::fmaq(x, y, z); }
inline __float128 fmax(__float128 x, __float128 y) { return ::fmaxq(x, y); }
inline __float128 fmin(__float128 x, __float128 y) { return ::fminq(x, y); }
inline __float128 fdim(__float128 x, __float128 y) { return ::fdimq(x, y); }
inline __float128 nanq(char const* arg) { return ::nanq(arg); }
// Exponential functions
inline __float128 exp(__float128 x) { return ::expq(x); }
inline __float128 exp2(__float128 x) { return ::exp2q(x); }
inline __float128 expm1(__float128 x) { return ::expm1q(x); }
inline __float128 log(__float128 x) { return ::logq(x); }
inline __float128 log10(__float128 x) { return ::log10q(x); }
inline __float128 log2(__float128 x) { return ::log2q(x); }
inline __float128 log1p(__float128 x) { return ::log1pq(x); }
// Power functions
inline __float128 pow(__float128 x, __float128 y) { return ::powq(x, y); }
inline __float128 sqrt(__float128 x) { return ::sqrtq(x); }
inline __float128 cbrt(__float128 x) { return ::cbrtq(x); }
inline __float128 hypot(__float128 x, __float128 y) { return ::hypotq(x, y); }
// Trigonometric functions
inline __float128 sin(__float128 x) { return ::sinq(x); }
inline __float128 cos(__float128 x) { return ::cosq(x); }
inline __float128 tan(__float128 x) { return ::tanq(x); }
inline __float128 asin(__float128 x) { return ::asinq(x); }
inline __float128 acos(__float128 x) { return ::acosq(x); }
inline __float128 atan(__float128 x) { return ::atanq(x); }
inline __float128 atan2(__float128 x, __float128 y) { return ::atan2q(x, y); }
// Hyperbolic functions
inline __float128 sinh(__float128 x) { return ::sinhq(x); }
inline __float128 cosh(__float128 x) { return ::coshq(x); }
inline __float128 tanh(__float128 x) { return ::tanhq(x); }
inline __float128 asinh(__float128 x) { return ::asinhq(x); }
inline __float128 acosh(__float128 x) { return ::acoshq(x); }
inline __float128 atanh(__float128 x) { return ::atanhq(x); }
// Error and gamma functions
inline __float128 erf(__float128 x) { return ::erfq(x); }
inline __float128 erfc(__float128 x) { return ::erfcq(x); }
inline __float128 tgamma(__float128 x) { return ::tgammaq(x); }
inline __float128 lgamma(__float128 x) { return ::lgammaq(x); }
// Nearest integer floating point operations
inline __float128 ceil(__float128 x) { return ::ceilq(x); }
inline __float128 floor(__float128 x) { return ::floorq(x); }
inline __float128 trunc(__float128 x) { return ::truncq(x); }
inline __float128 round(__float128 x) { return ::roundq(x); }
// lround
// llround
inline __float128 nearbyint(__float128 x) { return ::nearbyintq(x); }
// rint
// lrint
// llrint
// Floating point manipulation functions
// frexp
// ldexp
// modf
// scalbn
// scalbln
// ilog
inline __float128 logb(__float128 x) { return ::logbq(x); }
inline __float128 nextafter(__float128 x, __float128 y) { return ::nextafterq(x, y); }
// nexttoward
inline __float128 copysign(__float128 x, __float128 y) { return ::copysignq(x, y); }
// Classification and comparison
// fpclassify
inline bool isfinite(__float128 x) { return !::isinfq(x); }  // isfiniteq not provided
inline bool isinf(__float128 x) { return ::isinfq(x); }
inline bool isnan(__float128 x) { return ::isnanq(x); }
// isnormal
inline bool signbit(__float128 x) { return ::signbitq(x); }
// isgreater
// isgreaterequal
// isless
// islessequal
// islessgreater
// isunordered
// clang-format on
}  // namespace Kokkos
//</editor-fold>

//<editor-fold desc="Mathematical constants __float128 specializations">
namespace Kokkos::numbers {
// clang-format off
template <> constexpr __float128 e_v         <__float128> = 2.718281828459045235360287471352662498Q;
template <> constexpr __float128 log2e_v     <__float128> = 1.442695040888963407359924681001892137Q;
template <> constexpr __float128 log10e_v    <__float128> = 0.434294481903251827651128918916605082Q;
template <> constexpr __float128 pi_v        <__float128> = 3.141592653589793238462643383279502884Q;
template <> constexpr __float128 inv_pi_v    <__float128> = 0.318309886183790671537767526745028724Q;
template <> constexpr __float128 inv_sqrtpi_v<__float128> = 0.564189583547756286948079451560772586Q;
template <> constexpr __float128 ln2_v       <__float128> = 0.693147180559945309417232121458176568Q;
template <> constexpr __float128 ln10_v      <__float128> = 2.302585092994045684017991454684364208Q;
template <> constexpr __float128 sqrt2_v     <__float128> = 1.414213562373095048801688724209698079Q;
template <> constexpr __float128 sqrt3_v     <__float128> = 1.732050807568877293527446341505872367Q;
template <> constexpr __float128 inv_sqrt3_v <__float128> = 0.577350269189625764509148780501957456Q;
template <> constexpr __float128 egamma_v    <__float128> = 0.577215664901532860606512090082402431Q;
template <> constexpr __float128 phi_v       <__float128> = 1.618033988749894848204586834365638118Q;
// clang-format on
}  // namespace Kokkos::numbers
//</editor-fold>

#endif

#endif
