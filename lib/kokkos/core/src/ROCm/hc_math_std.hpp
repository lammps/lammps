#pragma once

#include "hc.hpp"
#include <cmath>

// Math functions with integer overloads will be converted to
// this floating point type.
#define HC_IMPLICIT_FLOAT_CONV double

#ifdef __KALMAR_ACCELERATOR__

#define HC_MATH_WRAPPER_1(function, arg1)              \
  template <typename T>                                \
  inline T function(T arg1) __attribute__((hc, cpu)) { \
    return hc::precise_math::function(arg1);           \
  }

#define KALMAR_MATH_WRAPPER_1(function, arg1) HC_MATH_WRAPPER_1(function, arg1)

#define HC_MATH_WRAPPER_FP_OVERLOAD_1(function, arg1)                       \
  template <typename T>                                                     \
  inline typename std::enable_if<std::is_integral<T>::value,                \
                                 HC_IMPLICIT_FLOAT_CONV>::type              \
  function(T arg1) __attribute__((hc, cpu)) {                               \
    return hc::precise_math::function(                                      \
        static_cast<HC_IMPLICIT_FLOAT_CONV>(arg1));                         \
  }                                                                         \
  template <typename T>                                                     \
  inline typename std::enable_if<std::is_floating_point<T>::value, T>::type \
  function(T arg1) __attribute__((hc, cpu)) {                               \
    return hc::precise_math::function(arg1);                                \
  }

#define KALMAR_MATH_WRAPPER_FP_OVERLOAD_1(function, arg1) \
  HC_MATH_WRAPPER_FP_OVERLOAD_1(function, arg1)

#define HC_MATH_WRAPPER_2(function, arg1, arg2)                \
  template <typename T>                                        \
  inline T function(T arg1, T arg2) __attribute__((hc, cpu)) { \
    return hc::precise_math::function(arg1, arg2);             \
  }

#define HC_MATH_ALIAS_2(alias, function, arg1, arg2)        \
  template <typename T>                                     \
  inline T alias(T arg1, T arg2) __attribute__((hc, cpu)) { \
    return hc::precise_math::function(arg1, arg2);          \
  }

#define HC_MATH_WRAPPER_3(function, arg1, arg2, arg3)                  \
  template <typename T>                                                \
  inline T function(T arg1, T arg2, T arg3) __attribute__((hc, cpu)) { \
    return hc::precise_math::function(arg1, arg2, arg3);               \
  }

#define HC_MATH_WRAPPER_TQ(function, arg1)             \
  template <typename T, typename Q>                    \
  inline T function(Q arg1) __attribute__((hc, cpu)) { \
    return hc::precise_math::function(arg1);           \
  }

#define HC_MATH_WRAPPER_FP_OVERLOAD_TQ(function, T, arg1)                   \
  template <typename Q>                                                     \
  inline typename std::enable_if<std::is_integral<Q>::value, T>::type       \
  function(Q arg1) __attribute__((hc, cpu)) {                               \
    return hc::precise_math::function(                                      \
        static_cast<HC_IMPLICIT_FLOAT_CONV>(arg1));                         \
  }                                                                         \
  template <typename Q>                                                     \
  inline typename std::enable_if<std::is_floating_point<Q>::value, T>::type \
  function(Q arg1) __attribute__((hc, cpu)) {                               \
    return hc::precise_math::function(arg1);                                \
  }

#define HC_MATH_WRAPPER_TTQ(function, arg1, arg2)              \
  template <typename T, typename Q>                            \
  inline T function(T arg1, Q arg2) __attribute__((hc, cpu)) { \
    return hc::precise_math::function(arg1, arg2);             \
  }

#define HC_MATH_WRAPPER_FP_OVERLOAD_TTQ(function, arg1, arg2)          \
  template <typename T, typename Q>                                    \
  inline typename std::enable_if<std::is_integral<T>::value ||         \
                                     std::is_integral<Q>::value,       \
                                 HC_IMPLICIT_FLOAT_CONV>::type         \
  function(T arg1, Q arg2) __attribute__((hc, cpu)) {                  \
    return hc::precise_math::function(                                 \
        static_cast<HC_IMPLICIT_FLOAT_CONV>(arg1),                     \
        static_cast<HC_IMPLICIT_FLOAT_CONV>(arg2));                    \
  }                                                                    \
  template <typename T, typename Q>                                    \
  inline typename std::enable_if<std::is_floating_point<T>::value &&   \
                                     std::is_floating_point<Q>::value, \
                                 T>::type                              \
  function(T arg1, Q arg2) __attribute__((hc, cpu)) {                  \
    return hc::precise_math::function(arg1, arg2);                     \
  }

#define HC_MATH_WRAPPER_TTTQ(function, arg1, arg2, arg3)               \
  template <typename T, typename Q>                                    \
  inline T function(T arg1, T arg2, Q arg3) __attribute__((hc, cpu)) { \
    return hc::precise_math::function(arg1, arg2, arg3);               \
  }

#define HC_MATH_WRAPPER_VTQQ(function, arg1, arg2, arg3)                  \
  template <typename T, typename Q>                                       \
  inline void function(T arg1, Q arg2, Q arg3) __attribute__((hc, cpu)) { \
    hc::precise_math::function(arg1, arg2, arg3);                         \
  }

#else

#define HC_MATH_WRAPPER_1(function, arg1)              \
  template <typename T>                                \
  inline T function(T arg1) __attribute__((hc, cpu)) { \
    return std::function(arg1);                        \
  }

#define KALMAR_MATH_WRAPPER_1(function, arg1)          \
  template <typename T>                                \
  inline T function(T arg1) __attribute__((hc, cpu)) { \
    return hc::precise_math::function(arg1);           \
  }

#define HC_MATH_WRAPPER_FP_OVERLOAD_1(function, arg1)                       \
  template <typename T>                                                     \
  inline typename std::enable_if<std::is_integral<T>::value,                \
                                 HC_IMPLICIT_FLOAT_CONV>::type              \
  function(T arg1) __attribute__((hc, cpu)) {                               \
    return ::function(static_cast<HC_IMPLICIT_FLOAT_CONV>(arg1));           \
  }                                                                         \
  template <typename T>                                                     \
  inline typename std::enable_if<std::is_floating_point<T>::value, T>::type \
  function(T arg1) __attribute__((hc, cpu)) {                               \
    return std::function(arg1);                                             \
  }

#define KALMAR_MATH_WRAPPER_FP_OVERLOAD_1(function, arg1)                   \
  template <typename T>                                                     \
  inline typename std::enable_if<std::is_integral<T>::value,                \
                                 HC_IMPLICIT_FLOAT_CONV>::type              \
  function(T arg1) __attribute__((hc, cpu)) {                               \
    return hc::precise_math::function(                                      \
        static_cast<HC_IMPLICIT_FLOAT_CONV>(arg1));                         \
  }                                                                         \
  template <typename T>                                                     \
  inline typename std::enable_if<std::is_floating_point<T>::value, T>::type \
  function(T arg1) __attribute__((hc, cpu)) {                               \
    return hc::precise_math::function(arg1);                                \
  }

#define HC_MATH_WRAPPER_2(function, arg1, arg2)                \
  template <typename T>                                        \
  inline T function(T arg1, T arg2) __attribute__((hc, cpu)) { \
    return std::function(arg1, arg2);                          \
  }

#define HC_MATH_ALIAS_2(alias, function, arg1, arg2)        \
  template <typename T>                                     \
  inline T alias(T arg1, T arg2) __attribute__((hc, cpu)) { \
    return std::function(arg1, arg2);                       \
  }

#define HC_MATH_WRAPPER_3(function, arg1, arg2, arg3)                  \
  template <typename T>                                                \
  inline T function(T arg1, T arg2, T arg3) __attribute__((hc, cpu)) { \
    return std::function(arg1, arg2, arg3);                            \
  }

#define HC_MATH_WRAPPER_TQ(function, arg1)             \
  template <typename T, typename Q>                    \
  inline T function(Q arg1) __attribute__((hc, cpu)) { \
    return std::function(arg1);                        \
  }

#define HC_MATH_WRAPPER_FP_OVERLOAD_TQ(function, T, arg1)                   \
  template <typename Q>                                                     \
  inline typename std::enable_if<std::is_integral<Q>::value, T>::type       \
  function(Q arg1) __attribute__((hc)) {                                    \
    return std::function(static_cast<HC_IMPLICIT_FLOAT_CONV>(arg1));        \
  }                                                                         \
  template <typename Q>                                                     \
  inline typename std::enable_if<std::is_floating_point<Q>::value, T>::type \
  function(Q arg1) __attribute__((hc)) {                                    \
    return std::function(arg1);                                             \
  }

#define HC_MATH_WRAPPER_TTQ(function, arg1, arg2)              \
  template <typename T, typename Q>                            \
  inline T function(T arg1, Q arg2) __attribute__((hc, cpu)) { \
    return std::function(arg1, arg2);                          \
  }

#define HC_MATH_WRAPPER_FP_OVERLOAD_TTQ(function, arg1, arg2)          \
  template <typename T, typename Q>                                    \
  inline typename std::enable_if<std::is_integral<T>::value ||         \
                                     std::is_integral<Q>::value,       \
                                 HC_IMPLICIT_FLOAT_CONV>::type         \
  function(T arg1, Q arg2) __attribute__((hc, cpu)) {                  \
    return std::function(static_cast<HC_IMPLICIT_FLOAT_CONV>(arg1),    \
                         static_cast<HC_IMPLICIT_FLOAT_CONV>(arg2));   \
  }                                                                    \
  template <typename T, typename Q>                                    \
  inline typename std::enable_if<std::is_floating_point<T>::value &&   \
                                     std::is_floating_point<Q>::value, \
                                 T>::type                              \
  function(T arg1, Q arg2) __attribute__((hc, cpu)) {                  \
    return std::function(arg1, arg2);                                  \
  }

#define HC_MATH_WRAPPER_TTTQ(function, arg1, arg2, arg3)               \
  template <typename T, typename Q>                                    \
  inline T function(T arg1, T arg2, Q arg3) __attribute__((hc, cpu)) { \
    return std::function(arg1, arg2, arg3);                            \
  }

#define HC_MATH_WRAPPER_VTQQ(function, arg1, arg2, arg3)                  \
  template <typename T, typename Q>                                       \
  inline void function(T arg1, Q arg2, Q arg3) __attribute__((hc, cpu)) { \
    std::function(arg1, arg2, arg3);                                      \
  }

#endif

// override global math functions
namespace std {

// following math functions are NOT available because they don't have a GPU
// implementation
//
// erfinv
// erfcinv
// fpclassify
//
// following math functions are NOT available because they don't have a CPU
// implementation
//
// cospif
// cospi
// rsqrtf
// rsqrt
// sinpif
// sinpi
// tanpi
//

HC_MATH_WRAPPER_TQ(ilogbf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_TQ(ilogb, int, x)
HC_MATH_WRAPPER_FP_OVERLOAD_TQ(isfinite, bool, x)
HC_MATH_WRAPPER_FP_OVERLOAD_TQ(isinf, bool, x)
HC_MATH_WRAPPER_FP_OVERLOAD_TQ(isnan, bool, x)
HC_MATH_WRAPPER_FP_OVERLOAD_TQ(isnormal, bool, x)
HC_MATH_WRAPPER_TQ(nanf, tagp)
HC_MATH_WRAPPER_TQ(nan, tagp)
// HC_MATH_WRAPPER_TQ(signbitf, x)
HC_MATH_WRAPPER_TQ(signbit, x)
HC_MATH_WRAPPER_TTQ(frexpf, x, exp)
HC_MATH_WRAPPER_TTQ(frexp, x, exp)
HC_MATH_WRAPPER_TTQ(ldexpf, x, exp)
HC_MATH_WRAPPER_TTQ(ldexp, x, exp)
HC_MATH_WRAPPER_TTQ(lgammaf, x, exp)
HC_MATH_WRAPPER_TTQ(lgamma, x, exp)
HC_MATH_WRAPPER_TTQ(modff, x, exp)
HC_MATH_WRAPPER_TTQ(modf, x, exp)
HC_MATH_WRAPPER_TTQ(scalbnf, x, exp)
HC_MATH_WRAPPER_TTQ(scalbn, x, exp)
HC_MATH_WRAPPER_TTTQ(remquof, x, y, quo)
HC_MATH_WRAPPER_TTTQ(remquo, x, y, quo)
HC_MATH_WRAPPER_VTQQ(sincosf, x, s, c)
HC_MATH_WRAPPER_VTQQ(sincos, x, s, c)

HC_MATH_WRAPPER_1(acosf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(acos, x)
HC_MATH_WRAPPER_1(acoshf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(acosh, x)
HC_MATH_WRAPPER_1(asinf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(asin, x)
HC_MATH_WRAPPER_1(asinhf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(asinh, x)
HC_MATH_WRAPPER_1(atanf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(atan, x)
HC_MATH_WRAPPER_1(atanhf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(atanh, x)
HC_MATH_WRAPPER_2(atan2f, x, y)
HC_MATH_WRAPPER_2(atan2, x, y)
HC_MATH_WRAPPER_1(cbrtf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(cbrt, x)
HC_MATH_WRAPPER_1(ceilf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(ceil, x)
HC_MATH_WRAPPER_2(copysignf, x, y)
HC_MATH_WRAPPER_2(copysign, x, y)
HC_MATH_WRAPPER_1(cosf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(cos, x)
HC_MATH_WRAPPER_1(coshf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(cosh, x)
KALMAR_MATH_WRAPPER_1(cospif, x)
KALMAR_MATH_WRAPPER_FP_OVERLOAD_1(cospi, x)
HC_MATH_WRAPPER_1(erff, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(erf, x)
HC_MATH_WRAPPER_1(erfcf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(erfc, x)
HC_MATH_WRAPPER_1(expf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(exp, x)
HC_MATH_WRAPPER_1(exp2f, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(exp2, x)
HC_MATH_WRAPPER_1(exp10f, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(exp10, x)
HC_MATH_WRAPPER_1(expm1f, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(expm1, x)
HC_MATH_WRAPPER_1(fabsf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(fabs, x)
HC_MATH_WRAPPER_2(fdimf, x, y)
HC_MATH_WRAPPER_2(fdim, x, y)
HC_MATH_WRAPPER_1(floorf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(floor, x)
HC_MATH_WRAPPER_3(fmaf, x, y, z)
HC_MATH_WRAPPER_3(fma, x, y, z)
HC_MATH_WRAPPER_2(fmaxf, x, y)
HC_MATH_WRAPPER_2(fmax, x, y)
HC_MATH_WRAPPER_2(fminf, x, y)
HC_MATH_WRAPPER_2(fmin, x, y)
HC_MATH_WRAPPER_2(fmodf, x, y)
HC_MATH_WRAPPER_2(fmod, x, y)
HC_MATH_WRAPPER_2(hypotf, x, y)
HC_MATH_WRAPPER_2(hypot, x, y)
HC_MATH_WRAPPER_1(logf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(log, x)
HC_MATH_WRAPPER_1(log10f, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(log10, x)
HC_MATH_WRAPPER_1(log2f, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(log2, x)
HC_MATH_WRAPPER_1(log1pf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(log1p, x)
HC_MATH_WRAPPER_1(logbf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(logb, x)
HC_MATH_WRAPPER_1(nearbyintf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(nearbyint, x)
HC_MATH_WRAPPER_2(nextafterf, x, y)
HC_MATH_WRAPPER_2(nextafter, x, y)
HC_MATH_WRAPPER_2(powf, x, y)
HC_MATH_WRAPPER_FP_OVERLOAD_TTQ(pow, x, y)
// HC_MATH_WRAPPER_1(rcbrtf, x)
// HC_MATH_WRAPPER_1(rcbrt, x)
HC_MATH_WRAPPER_2(remainderf, x, y)
HC_MATH_WRAPPER_2(remainder, x, y)
HC_MATH_WRAPPER_1(roundf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(round, x)
KALMAR_MATH_WRAPPER_1(rsqrtf, x)
KALMAR_MATH_WRAPPER_FP_OVERLOAD_1(rsqrt, x)
HC_MATH_WRAPPER_2(scalbf, x, exp)
HC_MATH_WRAPPER_2(scalb, x, exp)
HC_MATH_WRAPPER_1(sinf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(sin, x)
HC_MATH_WRAPPER_1(sinhf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(sinh, x)
KALMAR_MATH_WRAPPER_1(sinpif, x)
KALMAR_MATH_WRAPPER_FP_OVERLOAD_1(sinpi, x)
HC_MATH_WRAPPER_1(sqrtf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(sqrt, x)
HC_MATH_WRAPPER_1(tgammaf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(tgamma, x)
HC_MATH_WRAPPER_1(tanf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(tan, x)
HC_MATH_WRAPPER_1(tanhf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(tanh, x)
HC_MATH_WRAPPER_1(truncf, x)
HC_MATH_WRAPPER_FP_OVERLOAD_1(trunc, x)

// HC_MATH_ALIAS_2(min, fmin, x, y)
// HC_MATH_ALIAS_2(max, fmax, x, y)

}  // namespace std
