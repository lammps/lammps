// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_SIMD_PERFTEST_OPERATORS_HPP
#define KOKKOS_SIMD_PERFTEST_OPERATORS_HPP

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

class absolutes {
 public:
  template <typename T>
  auto on_host(T const& a) const {
    if constexpr (std::is_signed_v<typename T::value_type>) {
      return Kokkos::abs(a);
    }
    return a;
  }
  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a) const {
    if constexpr (std::is_signed_v<typename T::value_type>) {
      return Kokkos::abs(a);
    }
    return a;
  }
};

class floors {
 public:
  template <typename T>
  auto on_host(T const& a) const {
    return Kokkos::floor(a);
  }
  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a) const {
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
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a) const {
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
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a) const {
    return Kokkos::round(a);
  }
};

class truncates {
 public:
  template <typename T>
  auto on_host(T const& a) const {
    return Kokkos::trunc(a);
  }
  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a) const {
    return Kokkos::trunc(a);
  }
};

class shift_right {
 public:
  template <typename T, typename U>
  auto on_host(T const& a, U const& b) const {
    return a >> b;
  }
  template <typename T, typename U>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, U const& b) const {
    return a >> b;
  }
};

class shift_left {
 public:
  template <typename T, typename U>
  auto on_host(T const& a, U const& b) const {
    return a << b;
  }
  template <typename T, typename U>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, U const& b) const {
    return a << b;
  }
};

class minimum {
 public:
  template <typename T>
  auto on_host(T const& a, T const& b) const {
    return Kokkos::min(a, b);
  }
  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, T const& b) const {
    return Kokkos::min(a, b);
  }
};

class maximum {
 public:
  template <typename T>
  auto on_host(T const& a, T const& b) const {
    return Kokkos::max(a, b);
  }
  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, T const& b) const {
    return Kokkos::max(a, b);
  }
};

class reduce {
 public:
  template <typename T, typename MaskType>
  auto on_host(T const& a, MaskType) const {
    return Kokkos::Experimental::reduce(a, std::plus<>());
  }
  template <typename T, typename MaskType>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, MaskType) const {
    return Kokkos::Experimental::reduce(a, std::plus<>());
  }
};

class reduce_min {
 public:
  template <typename T, typename MaskType>
  auto on_host(T const& a, MaskType) const {
    return Kokkos::Experimental::reduce_min(a);
  }
  template <typename T, typename MaskType>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, MaskType) const {
    return Kokkos::Experimental::reduce_min(a);
  }
};

class reduce_max {
 public:
  template <typename T, typename MaskType>
  auto on_host(T const& a, MaskType) const {
    return Kokkos::Experimental::reduce_max(a);
  }
  template <typename T, typename MaskType>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, MaskType) const {
    return Kokkos::Experimental::reduce_max(a);
  }
};

class masked_reduce {
 public:
  template <typename T, typename MaskType>
  auto on_host(T const& a, MaskType const& mask) const {
    using DataType = typename T::value_type;
    return Kokkos::Experimental::reduce(
        a, mask, std::plus<>(),
        DataType(
            Kokkos::Experimental::Impl::Identity<DataType, std::plus<>>()));
  }
  template <typename T, typename MaskType>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a,
                                        MaskType const& mask) const {
    using DataType = typename T::value_type;
    return Kokkos::Experimental::reduce(
        a, mask, std::plus<>(),
        DataType(
            Kokkos::Experimental::Impl::Identity<DataType, std::plus<>>()));
  }
};

class masked_reduce_min {
 public:
  template <typename T, typename MaskType>
  auto on_host(T const& a, MaskType const& mask) const {
    return Kokkos::Experimental::reduce_min(a, mask);
  }
  template <typename T, typename MaskType>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a,
                                        MaskType const& mask) const {
    return Kokkos::Experimental::reduce_min(a, mask);
  }
};

class masked_reduce_max {
 public:
  template <typename T, typename MaskType>
  auto on_host(T const& a, MaskType const& mask) const {
    return Kokkos::Experimental::reduce_max(a, mask);
  }
  template <typename T, typename MaskType>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a,
                                        MaskType const& mask) const {
    return Kokkos::Experimental::reduce_max(a, mask);
  }
};

#define KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(name)        \
  class name##_op {                                           \
   public:                                                    \
    template <typename T>                                     \
    auto on_host(T const& a) const {                          \
      return Kokkos::name(a);                                 \
    }                                                         \
    template <typename T>                                     \
    KOKKOS_INLINE_FUNCTION auto on_device(T const& a) const { \
      return Kokkos::name(a);                                 \
    }                                                         \
  };

#define KOKKOS_IMPL_SIMD_PERFTEST_BINARY_OPERATOR(name)                   \
  class name##_op {                                                       \
   public:                                                                \
    template <typename T>                                                 \
    auto on_host(T const& a, T const& b) const {                          \
      return Kokkos::name(a, b);                                          \
    }                                                                     \
    template <typename T>                                                 \
    KOKKOS_INLINE_FUNCTION auto on_device(T const& a, T const& b) const { \
      return Kokkos::name(a, b);                                          \
    }                                                                     \
  };

#define KOKKOS_IMPL_SIMD_PERFTEST_TERNARY_OPERATOR(name)          \
  class name##_op {                                               \
   public:                                                        \
    template <typename T>                                         \
    auto on_host(T const& a, T const& b, T const& c) const {      \
      return Kokkos::name(a, b, c);                               \
    }                                                             \
    template <typename T>                                         \
    KOKKOS_INLINE_FUNCTION auto on_device(T const& a, T const& b, \
                                          T const& c) const {     \
      return Kokkos::name(a, b, c);                               \
    }                                                             \
  };

KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(exp)
KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(exp2)
KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(log)
KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(log10)
KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(log2)
KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(sqrt)
KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(cbrt)
KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(sin)
KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(cos)
KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(tan)
KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(asin)
KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(acos)
KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(atan)
KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(sinh)
KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(cosh)
KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(tanh)
KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(asinh)
KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(acosh)
KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(atanh)
KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(erf)
KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(erfc)
KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(tgamma)
KOKKOS_IMPL_SIMD_PERFTEST_UNARY_OPERATOR(lgamma)

KOKKOS_IMPL_SIMD_PERFTEST_BINARY_OPERATOR(pow)
KOKKOS_IMPL_SIMD_PERFTEST_BINARY_OPERATOR(hypot)
KOKKOS_IMPL_SIMD_PERFTEST_BINARY_OPERATOR(atan2)
KOKKOS_IMPL_SIMD_PERFTEST_BINARY_OPERATOR(copysign)

KOKKOS_IMPL_SIMD_PERFTEST_TERNARY_OPERATOR(fma)

class ternary_hypot_op {
 public:
  template <typename T>
  auto on_host(T const& a, T const& b, T const& c) const {
    return Kokkos::hypot(a, b, c);
  }
  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, T const& b,
                                        T const& c) const {
    return Kokkos::hypot(a, b, c);
  }
};

#endif
