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

#ifndef KOKKOS_SIMD_TESTING_OPS_HPP
#define KOKKOS_SIMD_TESTING_OPS_HPP

#include <Kokkos_SIMD.hpp>

class plus {
 public:
  template <class T>
  auto on_host(T const& a, T const& b) const {
    return a + b;
  }
  template <class T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, T const& b) const {
    return a + b;
  }
};

class plus_eq {
 public:
  template <class T>
  auto on_host(T&& a, T&& b) const {
    return a += b;
  }
  template <class T>
  KOKKOS_INLINE_FUNCTION auto on_device(T&& a, T&& b) const {
    return a += b;
  }
};

class minus {
 public:
  template <class T>
  auto on_host(T const& a, T const& b) const {
    return a - b;
  }
  template <class T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, T const& b) const {
    return a - b;
  }
};

class minus_eq {
 public:
  template <class T>
  auto on_host(T&& a, T&& b) const {
    return a -= b;
  }
  template <class T>
  KOKKOS_INLINE_FUNCTION auto on_device(T&& a, T&& b) const {
    return a -= b;
  }
};

class multiplies {
 public:
  template <class T>
  auto on_host(T const& a, T const& b) const {
    return a * b;
  }
  template <class T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, T const& b) const {
    return a * b;
  }
};

class multiplies_eq {
 public:
  template <class T>
  auto on_host(T&& a, T&& b) const {
    return a *= b;
  }
  template <class T>
  KOKKOS_INLINE_FUNCTION auto on_device(T&& a, T&& b) const {
    return a *= b;
  }
};

class divides {
 public:
  template <class T>
  auto on_host(T const& a, T const& b) const {
    return a / b;
  }
  template <class T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, T const& b) const {
    return a / b;
  }
};

class divides_eq {
 public:
  template <class T>
  auto on_host(T&& a, T&& b) const {
    return a /= b;
  }
  template <class T>
  KOKKOS_INLINE_FUNCTION auto on_device(T&& a, T&& b) const {
    return a /= b;
  }
};

class absolutes {
  template <typename T>
  static KOKKOS_FUNCTION auto abs_impl(T const& x) {
    if constexpr (std::is_signed_v<T>) {
      return Kokkos::abs(x);
    }
    return x;
  }

 public:
  template <typename T>
  auto on_host(T const& a) const {
    if constexpr (std::is_signed_v<typename T::value_type>) {
#if defined(KOKKOS_ENABLE_DEPRECATED_CODE_4)
#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
      KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_PUSH()
#endif
      return Kokkos::Experimental::abs(a);
#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
      KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_POP()
#endif
#else
      return Kokkos::abs(a);
#endif
    }
    return a;
  }
  template <typename T>
  auto on_host_serial(T const& a) const {
    return abs_impl(a);
  }
  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a) const {
    if constexpr (std::is_signed_v<typename T::value_type>) {
      return Kokkos::abs(a);
    }
    return a;
  }
  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device_serial(T const& a) const {
    return abs_impl(a);
  }
};

class floors {
 public:
  template <typename T>
  auto on_host(T const& a) const {
    return Kokkos::floor(a);
  }
  template <typename T>
  auto on_host_serial(T const& a) const {
    return Kokkos::floor(a);
  }
  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a) const {
    return Kokkos::floor(a);
  }
  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device_serial(T const& a) const {
    return Kokkos::floor(a);
  }
};

class ceils {
 public:
  template <typename T>
  auto on_host(T const& a) const {
    return Kokkos::ceil(a);
  }
  template <typename T>
  auto on_host_serial(T const& a) const {
    return Kokkos::ceil(a);
  }
  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a) const {
    return Kokkos::ceil(a);
  }
  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device_serial(T const& a) const {
    return Kokkos::ceil(a);
  }
};

class rounds {
 public:
  template <typename T>
  auto on_host(T const& a) const {
    return Kokkos::round(a);
  }
  template <typename T>
  auto on_host_serial(T const& a) const {
    return std::rint(a);
  }
  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a) const {
    return Kokkos::round(a);
  }
  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device_serial(T const& a) const {
    return Kokkos::Experimental::round_half_to_nearest_even(a);
  }
};

class truncates {
 public:
  template <typename T>
  auto on_host(T const& a) const {
    return Kokkos::trunc(a);
  }
  template <typename T>
  auto on_host_serial(T const& a) const {
    return Kokkos::trunc(a);
  }
  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a) const {
    return Kokkos::trunc(a);
  }
  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device_serial(T const& a) const {
    return Kokkos::trunc(a);
  }
};

class shift_right {
 public:
  template <typename T, typename U>
  auto on_host(T&& a, U&& b) const {
    return a >> b;
  }
  template <typename T, typename U>
  KOKKOS_INLINE_FUNCTION auto on_device(T&& a, U&& b) const {
    return a >> b;
  }
};

class shift_right_eq {
 public:
  template <typename T, typename U>
  auto on_host(T&& a, U&& b) const {
    return a >>= b;
  }
  template <typename T, typename U>
  KOKKOS_INLINE_FUNCTION auto on_device(T&& a, U&& b) const {
    return a >>= b;
  }
};

class shift_left {
 public:
  template <typename T, typename U>
  auto on_host(T&& a, U&& b) const {
    return a << b;
  }
  template <typename T, typename U>
  KOKKOS_INLINE_FUNCTION auto on_device(T&& a, U&& b) const {
    return a << b;
  }
};

class shift_left_eq {
 public:
  template <typename T, typename U>
  auto on_host(T&& a, U&& b) const {
    return a <<= b;
  }
  template <typename T, typename U>
  KOKKOS_INLINE_FUNCTION auto on_device(T&& a, U&& b) const {
    return a <<= b;
  }
};

class minimum {
 public:
  template <typename T>
  auto on_host(T const& a, T const& b) const {
#if defined(KOKKOS_ENABLE_DEPRECATED_CODE_4)
    if constexpr (std::is_arithmetic_v<T>) {
      return Kokkos::min(a, b);
    } else {
#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
      KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_PUSH()
#endif
      return Kokkos::Experimental::min(a, b);
#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
      KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_POP()
#endif
    }
#else
    return Kokkos::min(a, b);
#endif
  }
  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, T const& b) const {
    return Kokkos::min(a, b);
  }

  template <typename T>
  KOKKOS_INLINE_FUNCTION T operator()(T const& a, T const& b) const {
    return Kokkos::min(a, b);
  }
};

class maximum {
 public:
  template <typename T>
  auto on_host(T const& a, T const& b) const {
#if defined(KOKKOS_ENABLE_DEPRECATED_CODE_4)
    if constexpr (std::is_arithmetic_v<T>) {
      return Kokkos::max(a, b);
    } else {
#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
      KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_PUSH()
#endif
      return Kokkos::Experimental::max(a, b);
#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
      KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_POP()
#endif
    }
#else
    return Kokkos::max(a, b);
#endif
  }
  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, T const& b) const {
    return Kokkos::max(a, b);
  }

  template <typename T>
  KOKKOS_INLINE_FUNCTION T operator()(T const& a, T const& b) const {
    return Kokkos::max(a, b);
  }
};

class reduce_min {
 public:
  template <typename T, typename U, typename MaskType>
  auto on_host(T const& a, U, MaskType) const {
    return Kokkos::Experimental::reduce_min(a);
  }
  template <typename T, typename U, typename MaskType>
  auto on_host_serial(T const& a, U, MaskType) const {
    auto result = Kokkos::reduction_identity<U>::min();
    for (std::size_t i = 0; i < a.size(); ++i) {
      result = Kokkos::min(result, a[i]);
    }
    return result;
  }

  template <typename T, typename U, typename MaskType>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, U, MaskType) const {
    return Kokkos::Experimental::reduce_min(a);
  }
  template <typename T, typename U, typename MaskType>
  KOKKOS_INLINE_FUNCTION auto on_device_serial(T const& a, U, MaskType) const {
    auto result = Kokkos::reduction_identity<U>::min();
    for (std::size_t i = 0; i < a.size(); ++i) {
      result = Kokkos::min(result, a[i]);
    }
    return result;
  }
};

class reduce_max {
 public:
  template <typename T, typename U, typename MaskType>
  auto on_host(T const& a, U, MaskType) const {
    return Kokkos::Experimental::reduce_max(a);
  }
  template <typename T, typename U, typename MaskType>
  auto on_host_serial(T const& a, U, MaskType) const {
    auto result = Kokkos::reduction_identity<U>::max();
    for (std::size_t i = 0; i < a.size(); ++i) {
      result = Kokkos::max(result, a[i]);
    }
    return result;
  }

  template <typename T, typename U, typename MaskType>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, U, MaskType) const {
    return Kokkos::Experimental::reduce_max(a);
  }
  template <typename T, typename U, typename MaskType>
  KOKKOS_INLINE_FUNCTION auto on_device_serial(T const& a, U, MaskType) const {
    auto result = Kokkos::reduction_identity<U>::max();
    for (std::size_t i = 0; i < a.size(); ++i) {
      result = Kokkos::max(result, a[i]);
    }
    return result;
  }
};

template <typename BinaryOperation = std::plus<>>
class reduce {
 public:
  template <typename T, typename U, typename MaskType>
  auto on_host(T const& a, U, MaskType) const {
    return Kokkos::Experimental::reduce(a, BinaryOperation());
  }
  template <typename T, typename U, typename MaskType>
  auto on_host_serial(T const& a, U, MaskType) const {
    U result = a[0];
    for (std::size_t i = 1; i < a.size(); ++i) {
      result = BinaryOperation()(result, a[i]);
    }
    return result;
  }

  template <typename T, typename U, typename MaskType>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, U, MaskType) const {
    return Kokkos::Experimental::reduce(a, BinaryOperation());
  }
  template <typename T, typename U, typename MaskType>
  KOKKOS_INLINE_FUNCTION auto on_device_serial(T const& a, U, MaskType) const {
    U result = a[0];
    for (std::size_t i = 1; i < a.size(); ++i) {
      if constexpr (std::is_same_v<BinaryOperation, std::plus<>>) {
        result = result + a[i];
      } else if constexpr (std::is_same_v<BinaryOperation, std::multiplies<>>) {
        result = result * a[i];
      } else if constexpr (std::is_same_v<BinaryOperation, std::bit_and<>>) {
        result = result & a[i];
      } else if constexpr (std::is_same_v<BinaryOperation, std::bit_or<>>) {
        result = result | a[i];
      } else if constexpr (std::is_same_v<BinaryOperation, std::bit_xor<>>) {
        result = result ^ a[i];
      } else {
        result = BinaryOperation()(result, a[i]);
      }
    }
    return result;
  }
};

class masked_reduce_min {
 public:
  template <typename T, typename U, typename MaskType>
  auto on_host(T const& a, U, MaskType mask) const {
    return Kokkos::Experimental::reduce_min(a, mask);
  }
  template <typename T, typename U, typename MaskType>
  auto on_host_serial(T const& a, U, MaskType mask) const {
    if (Kokkos::Experimental::none_of(mask))
      return Kokkos::reduction_identity<U>::min();
    auto w        = Kokkos::Experimental::where(mask, a);
    auto const& v = w.impl_get_value();
    auto const& m = w.impl_get_mask();
    auto result   = Kokkos::reduction_identity<U>::min();
    for (std::size_t i = 0; i < v.size(); ++i) {
      if (m[i]) result = Kokkos::min(result, v[i]);
    }
    return result;
  }

  template <typename T, typename U, typename MaskType>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, U, MaskType mask) const {
    return Kokkos::Experimental::reduce_min(a, mask);
  }
  template <typename T, typename U, typename MaskType>
  KOKKOS_INLINE_FUNCTION auto on_device_serial(T const& a, U,
                                               MaskType mask) const {
    if (Kokkos::Experimental::none_of(mask))
      return Kokkos::reduction_identity<U>::min();
    auto w        = Kokkos::Experimental::where(mask, a);
    auto const& v = w.impl_get_value();
    auto const& m = w.impl_get_mask();
    auto result   = Kokkos::reduction_identity<U>::min();
    for (std::size_t i = 0; i < v.size(); ++i) {
      if (m[i]) result = Kokkos::min(result, v[i]);
    }
    return result;
  }
};

class masked_reduce_max {
 public:
  template <typename T, typename U, typename MaskType>
  auto on_host(T const& a, U, MaskType mask) const {
    return Kokkos::Experimental::reduce_max(a, mask);
  }
  template <typename T, typename U, typename MaskType>
  auto on_host_serial(T const& a, U, MaskType mask) const {
    if (Kokkos::Experimental::none_of(mask))
      return Kokkos::reduction_identity<U>::max();
    auto w        = Kokkos::Experimental::where(mask, a);
    auto const& v = w.impl_get_value();
    auto const& m = w.impl_get_mask();
    auto result   = Kokkos::reduction_identity<U>::max();
    for (std::size_t i = 0; i < v.size(); ++i) {
      if (m[i]) result = Kokkos::max(result, v[i]);
    }
    return result;
  }

  template <typename T, typename U, typename MaskType>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, U, MaskType mask) const {
    return Kokkos::Experimental::reduce_max(a, mask);
  }
  template <typename T, typename U, typename MaskType>
  KOKKOS_INLINE_FUNCTION auto on_device_serial(T const& a, U,
                                               MaskType mask) const {
    if (Kokkos::Experimental::none_of(mask))
      return Kokkos::reduction_identity<U>::max();
    auto w        = Kokkos::Experimental::where(mask, a);
    auto const& v = w.impl_get_value();
    auto const& m = w.impl_get_mask();
    auto result   = Kokkos::reduction_identity<U>::max();
    for (std::size_t i = 0; i < v.size(); ++i) {
      if (m[i]) result = Kokkos::max(result, v[i]);
    }
    return result;
  }
};

template <typename BinaryOperation = std::plus<>>
class masked_reduce {
 public:
  template <typename T, typename U, typename MaskType>
  auto on_host(T const& a, U const& identity, MaskType mask) const {
    return Kokkos::Experimental::reduce(a, mask, BinaryOperation(), identity);
  }
  template <typename T, typename U, typename MaskType>
  auto on_host_serial(T const& a, U const& identity, MaskType mask) const {
    if (Kokkos::Experimental::none_of(mask)) return identity;
    auto w        = Kokkos::Experimental::where(mask, a);
    auto const& v = w.impl_get_value();
    auto const& m = w.impl_get_mask();
    U result      = identity;
    for (std::size_t i = 0; i < v.size(); ++i) {
      if (m[i]) result = BinaryOperation()(result, v[i]);
    }
    return result;
  }

  template <typename T, typename U, typename MaskType>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, U const& identity,
                                        MaskType mask) const {
    return Kokkos::Experimental::reduce(a, mask, BinaryOperation(), identity);
  }
  template <typename T, typename U, typename MaskType>
  KOKKOS_INLINE_FUNCTION auto on_device_serial(T const& a, U const& identity,
                                               MaskType mask) const {
    if (Kokkos::Experimental::none_of(mask)) return identity;
    auto w        = Kokkos::Experimental::where(mask, a);
    auto const& v = w.impl_get_value();
    auto const& m = w.impl_get_mask();
    U result      = identity;
    for (std::size_t i = 0; i < v.size(); ++i) {
      if constexpr (std::is_same_v<BinaryOperation, std::plus<>>) {
        if (m[i]) result = result + v[i];
      } else if constexpr (std::is_same_v<BinaryOperation, std::multiplies<>>) {
        if (m[i]) result = result * v[i];
      } else if constexpr (std::is_same_v<BinaryOperation, std::bit_and<>>) {
        if (m[i]) result = result & v[i];
      } else if constexpr (std::is_same_v<BinaryOperation, std::bit_or<>>) {
        if (m[i]) result = result | v[i];
      } else if constexpr (std::is_same_v<BinaryOperation, std::bit_xor<>>) {
        if (m[i]) result = result ^ v[i];
      } else {
        if (m[i]) result = BinaryOperation()(result, v[i]);
      }
    }
    return result;
  }
};

#if defined(KOKKOS_ENABLE_DEPRECATED_CODE_4)
#define KOKKOS_IMPL_SIMD_UNARY_TEST_OP(FUNC)                         \
  class FUNC##_op {                                                  \
   public:                                                           \
    template <typename T>                                            \
    auto on_host(T const& a) const {                                 \
      return Kokkos::Experimental::FUNC(a);                          \
    }                                                                \
    template <typename T>                                            \
    auto on_host_serial(T const& a) const {                          \
      return Kokkos::FUNC(a);                                        \
    }                                                                \
    template <typename T>                                            \
    KOKKOS_INLINE_FUNCTION auto on_device(T const& a) const {        \
      return Kokkos::FUNC(a);                                        \
    }                                                                \
    template <typename T>                                            \
    KOKKOS_INLINE_FUNCTION auto on_device_serial(T const& a) const { \
      return Kokkos::FUNC(a);                                        \
    }                                                                \
  };

#define KOKKOS_IMPL_SIMD_BINARY_TEST_OP(FUNC)                             \
  class FUNC##_op {                                                       \
   public:                                                                \
    template <typename T>                                                 \
    auto on_host(T const& a, T const& b) const {                          \
      if constexpr (std::is_arithmetic_v<T>) {                            \
        return Kokkos::FUNC(a, b);                                        \
      } else {                                                            \
        return Kokkos::Experimental::FUNC(a, b);                          \
      }                                                                   \
    }                                                                     \
    template <typename T>                                                 \
    KOKKOS_INLINE_FUNCTION auto on_device(T const& a, T const& b) const { \
      return Kokkos::FUNC(a, b);                                          \
    }                                                                     \
  };

#define KOKKOS_IMPL_SIMD_TERNARY_TEST_OP(FUNC)                    \
  class FUNC##_op {                                               \
   public:                                                        \
    template <typename T>                                         \
    auto on_host(T const& a, T const& b, T const& c) const {      \
      if constexpr (std::is_arithmetic_v<T>) {                    \
        return Kokkos::FUNC(a, b, c);                             \
      } else {                                                    \
        return Kokkos::Experimental::FUNC(a, b, c);               \
      }                                                           \
    }                                                             \
    template <typename T>                                         \
    KOKKOS_INLINE_FUNCTION auto on_device(T const& a, T const& b, \
                                          T const& c) const {     \
      return Kokkos::FUNC(a, b, c);                               \
    }                                                             \
  };
#else
#define KOKKOS_IMPL_SIMD_UNARY_TEST_OP(FUNC)                         \
  class FUNC##_op {                                                  \
   public:                                                           \
    template <typename T>                                            \
    auto on_host(T const& a) const {                                 \
      return Kokkos::FUNC(a);                                        \
    }                                                                \
    template <typename T>                                            \
    auto on_host_serial(T const& a) const {                          \
      return Kokkos::FUNC(a);                                        \
    }                                                                \
    template <typename T>                                            \
    KOKKOS_INLINE_FUNCTION auto on_device(T const& a) const {        \
      return Kokkos::FUNC(a);                                        \
    }                                                                \
    template <typename T>                                            \
    KOKKOS_INLINE_FUNCTION auto on_device_serial(T const& a) const { \
      return Kokkos::FUNC(a);                                        \
    }                                                                \
  };

#define KOKKOS_IMPL_SIMD_BINARY_TEST_OP(FUNC)                             \
  class FUNC##_op {                                                       \
   public:                                                                \
    template <typename T>                                                 \
    auto on_host(T const& a, T const& b) const {                          \
      return Kokkos::FUNC(a, b);                                          \
    }                                                                     \
    template <typename T>                                                 \
    KOKKOS_INLINE_FUNCTION auto on_device(T const& a, T const& b) const { \
      return Kokkos::FUNC(a, b);                                          \
    }                                                                     \
  };

#define KOKKOS_IMPL_SIMD_TERNARY_TEST_OP(FUNC)                    \
  class FUNC##_op {                                               \
   public:                                                        \
    template <typename T>                                         \
    auto on_host(T const& a, T const& b, T const& c) const {      \
      return Kokkos::FUNC(a, b, c);                               \
    }                                                             \
    template <typename T>                                         \
    KOKKOS_INLINE_FUNCTION auto on_device(T const& a, T const& b, \
                                          T const& c) const {     \
      return Kokkos::FUNC(a, b, c);                               \
    }                                                             \
  };
#endif

#if defined(KOKKOS_ENABLE_DEPRECATED_CODE_4)
KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_PUSH()
#endif

KOKKOS_IMPL_SIMD_UNARY_TEST_OP(abs)
KOKKOS_IMPL_SIMD_UNARY_TEST_OP(exp)
KOKKOS_IMPL_SIMD_UNARY_TEST_OP(exp2)
KOKKOS_IMPL_SIMD_UNARY_TEST_OP(log)
KOKKOS_IMPL_SIMD_UNARY_TEST_OP(log10)
KOKKOS_IMPL_SIMD_UNARY_TEST_OP(log2)
KOKKOS_IMPL_SIMD_UNARY_TEST_OP(sqrt)
KOKKOS_IMPL_SIMD_UNARY_TEST_OP(cbrt)
KOKKOS_IMPL_SIMD_UNARY_TEST_OP(sin)
KOKKOS_IMPL_SIMD_UNARY_TEST_OP(cos)
KOKKOS_IMPL_SIMD_UNARY_TEST_OP(tan)
KOKKOS_IMPL_SIMD_UNARY_TEST_OP(asin)
KOKKOS_IMPL_SIMD_UNARY_TEST_OP(acos)
KOKKOS_IMPL_SIMD_UNARY_TEST_OP(atan)
KOKKOS_IMPL_SIMD_UNARY_TEST_OP(sinh)
KOKKOS_IMPL_SIMD_UNARY_TEST_OP(cosh)
KOKKOS_IMPL_SIMD_UNARY_TEST_OP(tanh)
KOKKOS_IMPL_SIMD_UNARY_TEST_OP(asinh)
KOKKOS_IMPL_SIMD_UNARY_TEST_OP(acosh)
KOKKOS_IMPL_SIMD_UNARY_TEST_OP(atanh)
KOKKOS_IMPL_SIMD_UNARY_TEST_OP(erf)
KOKKOS_IMPL_SIMD_UNARY_TEST_OP(erfc)
KOKKOS_IMPL_SIMD_UNARY_TEST_OP(tgamma)
KOKKOS_IMPL_SIMD_UNARY_TEST_OP(lgamma)

KOKKOS_IMPL_SIMD_BINARY_TEST_OP(pow)
KOKKOS_IMPL_SIMD_BINARY_TEST_OP(hypot)
KOKKOS_IMPL_SIMD_BINARY_TEST_OP(atan2)
KOKKOS_IMPL_SIMD_BINARY_TEST_OP(copysign)

KOKKOS_IMPL_SIMD_TERNARY_TEST_OP(fma)

class ternary_hypot_op {
 public:
  template <typename T>
  auto on_host(T const& a, T const& b, T const& c) const {
#if defined(KOKKOS_ENABLE_DEPRECATED_CODE_4)
    if constexpr (std::is_arithmetic_v<T>) {
      return Kokkos::hypot(a, b, c);
    } else {
      return Kokkos::Experimental::hypot(a, b, c);
    }
#else
    return Kokkos::hypot(a, b, c);
#endif
  }
  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, T const& b,
                                        T const& c) const {
    return Kokkos::hypot(a, b, c);
  }
};

#if defined(KOKKOS_ENABLE_DEPRECATED_CODE_4)
KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_POP()
#endif

#endif
