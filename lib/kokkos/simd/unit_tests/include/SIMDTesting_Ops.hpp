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
      KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_PUSH()
      return Kokkos::Experimental::abs(a);
      KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_POP()
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

class cbrt_op {
 public:
  template <typename T>
  auto on_host(T const& a) const {
#if defined(KOKKOS_ENABLE_DEPRECATED_CODE_4)
    return Kokkos::Experimental::cbrt(a);
#else
    return Kokkos::cbrt(a);
#endif
  }
  template <typename T>
  auto on_host_serial(T const& a) const {
    return Kokkos::cbrt(a);
  }
};

class exp_op {
 public:
  template <typename T>
  auto on_host(T const& a) const {
#if defined(KOKKOS_ENABLE_DEPRECATED_CODE_4)
    return Kokkos::Experimental::exp(a);
#else
    return Kokkos::exp(a);
#endif
  }
  template <typename T>
  auto on_host_serial(T const& a) const {
    return Kokkos::exp(a);
  }
};

class log_op {
 public:
  template <typename T>
  auto on_host(T const& a) const {
#if defined(KOKKOS_ENABLE_DEPRECATED_CODE_4)
    return Kokkos::Experimental::log(a);
#else
    return Kokkos::log(a);
#endif
  }
  template <typename T>
  auto on_host_serial(T const& a) const {
    return Kokkos::log(a);
  }
};

class hmin {
 public:
  template <typename T>
  auto on_host(T const& a) const {
    return Kokkos::Experimental::hmin(a);
  }
  template <typename T>
  auto on_host_serial(T const& a) const {
    using DataType = typename T::value_type::value_type;

    auto const& v = a.impl_get_value();
    auto const& m = a.impl_get_mask();
    auto result   = Kokkos::reduction_identity<DataType>::min();
    for (std::size_t i = 0; i < v.size(); ++i) {
      if (m[i]) result = Kokkos::min(result, v[i]);
    }
    return result;
  }

  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a) const {
    return Kokkos::Experimental::hmin(a);
  }
  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device_serial(T const& a) const {
    using DataType = typename T::value_type::value_type;

    auto const& v = a.impl_get_value();
    auto const& m = a.impl_get_mask();
    auto result   = Kokkos::reduction_identity<DataType>::min();
    for (std::size_t i = 0; i < v.size(); ++i) {
      if (m[i]) result = Kokkos::min(result, v[i]);
    }
    return result;
  }
};

class hmax {
 public:
  template <typename T>
  auto on_host(T const& a) const {
    return Kokkos::Experimental::hmax(a);
  }
  template <typename T>
  auto on_host_serial(T const& a) const {
    using DataType = typename T::value_type::value_type;

    auto const& v = a.impl_get_value();
    auto const& m = a.impl_get_mask();
    auto result   = Kokkos::reduction_identity<DataType>::max();
    for (std::size_t i = 0; i < v.size(); ++i) {
      if (m[i]) result = Kokkos::max(result, v[i]);
    }
    return result;
  }

  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a) const {
    return Kokkos::Experimental::hmax(a);
  }
  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device_serial(T const& a) const {
    using DataType = typename T::value_type::value_type;

    auto const& v = a.impl_get_value();
    auto const& m = a.impl_get_mask();
    auto result   = Kokkos::reduction_identity<DataType>::max();
    for (std::size_t i = 0; i < v.size(); ++i) {
      if (m[i]) result = Kokkos::max(result, v[i]);
    }
    return result;
  }
};

class reduce {
 public:
  template <typename T>
  auto on_host(T const& a) const {
    using DataType = typename T::value_type::value_type;
    return Kokkos::Experimental::reduce(a, DataType(0), std::plus<>());
  }
  template <typename T>
  auto on_host_serial(T const& a) const {
    using DataType = typename T::value_type::value_type;

    auto const& v = a.impl_get_value();
    auto const& m = a.impl_get_mask();
    auto result   = Kokkos::reduction_identity<DataType>::sum();
    for (std::size_t i = 0; i < v.size(); ++i) {
      if (m[i]) result += v[i];
    }
    return result;
  }

  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a) const {
    using DataType = typename T::value_type::value_type;
    return Kokkos::Experimental::reduce(a, DataType(0), std::plus<>());
  }
  template <typename T>
  KOKKOS_INLINE_FUNCTION auto on_device_serial(T const& a) const {
    using DataType = typename T::value_type::value_type;

    auto const& v = a.impl_get_value();
    auto const& m = a.impl_get_mask();
    auto result   = Kokkos::reduction_identity<DataType>::sum();
    for (std::size_t i = 0; i < v.size(); ++i) {
      if (m[i]) result += v[i];
    }
    return result;
  }
};

#endif
