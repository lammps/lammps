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

#ifndef KOKKOS_HIP_HALF_HPP_
#define KOKKOS_HIP_HALF_HPP_

#ifdef KOKKOS_IMPL_HALF_TYPE_DEFINED

#include <Kokkos_Half.hpp>
#include <Kokkos_ReductionIdentity.hpp>

namespace Kokkos {
namespace Experimental {

/************************** half conversions **********************************/
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(half_t val) { return val; }

KOKKOS_INLINE_FUNCTION
half_t cast_to_half(float val) { return half_t(__float2half(val)); }

KOKKOS_INLINE_FUNCTION
half_t cast_to_half(bool val) { return cast_to_half(static_cast<float>(val)); }

KOKKOS_INLINE_FUNCTION
half_t cast_to_half(double val) {
  return half_t(__float2half(static_cast<float>(val)));
}

KOKKOS_INLINE_FUNCTION
half_t cast_to_half(short val) {
#ifdef __HIP_DEVICE_COMPILE__
  return half_t(__short2half_rn(val));
#else
  return half_t(__float2half(static_cast<float>(val)));
#endif
}

KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned short val) {
#ifdef __HIP_DEVICE_COMPILE__
  return half_t(__ushort2half_rn(val));
#else
  return half_t(__float2half(static_cast<float>(val)));
#endif
}

KOKKOS_INLINE_FUNCTION
half_t cast_to_half(int val) {
#ifdef __HIP_DEVICE_COMPILE__
  return half_t(__int2half_rn(val));
#else
  return half_t(__float2half(static_cast<float>(val)));
#endif
}

KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned int val) {
#ifdef __HIP_DEVICE_COMPILE__
  return half_t(__uint2half_rn(val));
#else
  return half_t(__float2half(static_cast<float>(val)));
#endif
}

KOKKOS_INLINE_FUNCTION
half_t cast_to_half(long long val) {
#ifdef __HIP_DEVICE_COMPILE__
  return half_t(__ll2half_rn(val));
#else
  return half_t(__float2half(static_cast<float>(val)));
#endif
}

KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned long long val) {
#ifdef __HIP_DEVICE_COMPILE__
  return half_t(__ull2half_rn(val));
#else
  return half_t(__float2half(static_cast<float>(val)));
#endif
}

KOKKOS_INLINE_FUNCTION
half_t cast_to_half(long val) {
  return cast_to_half(static_cast<long long>(val));
}

KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned long val) {
  return cast_to_half(static_cast<unsigned long long>(val));
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, float>::value, T>
cast_from_half(half_t val) {
  return __half2float(half_t::impl_type(val));
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, bool>::value, T>
cast_from_half(half_t val) {
  return static_cast<T>(cast_from_half<float>(val));
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, double>::value, T>
cast_from_half(half_t val) {
  return static_cast<T>(__half2float(half_t::impl_type(val)));
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, short>::value, T>
cast_from_half(half_t val) {
#ifdef __HIP_DEVICE_COMPILE__
  return __half2short_rz(half_t::impl_type(val));
#else
  return static_cast<T>(__half2float(half_t::impl_type(val)));
#endif
}

template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same<T, unsigned short>::value, T>
    cast_from_half(half_t val) {
#ifdef __HIP_DEVICE_COMPILE__
  return __half2ushort_rz(half_t::impl_type(val));
#else
  return static_cast<T>(__half2float(half_t::impl_type(val)));
#endif
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, int>::value, T>
cast_from_half(half_t val) {
#ifdef __HIP_DEVICE_COMPILE__
  return __half2int_rz(half_t::impl_type(val));
#else
  return static_cast<T>(__half2float(half_t::impl_type(val)));
#endif
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, unsigned>::value, T>
cast_from_half(half_t val) {
#ifdef __HIP_DEVICE_COMPILE__
  return __half2uint_rz(half_t::impl_type(val));
#else
  return static_cast<T>(__half2float(half_t::impl_type(val)));
#endif
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, long long>::value, T>
cast_from_half(half_t val) {
#ifdef __HIP_DEVICE_COMPILE__
  return __half2ll_rz(half_t::impl_type(val));
#else
  return static_cast<T>(__half2float(half_t::impl_type(val)));
#endif
}

template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same<T, unsigned long long>::value, T>
    cast_from_half(half_t val) {
#ifdef __HIP_DEVICE_COMPILE__
  return __half2ull_rz(half_t::impl_type(val));
#else
  return static_cast<T>(__half2float(half_t::impl_type(val)));
#endif
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, long>::value, T>
cast_from_half(half_t val) {
  return static_cast<T>(cast_from_half<long long>(val));
}

template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same<T, unsigned long>::value, T>
    cast_from_half(half_t val) {
  return static_cast<T>(cast_from_half<unsigned long long>(val));
}
}  // namespace Experimental

// use float as the return type for sum and prod since hip_fp16.h
// has no constexpr functions for casting to __half
template <>
struct reduction_identity<Kokkos::Experimental::half_t> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static float sum() noexcept {
    return 0.0F;
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static float prod() noexcept {
    return 1.0F;
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static float max() noexcept {
    return -65504.0F;
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static float min() noexcept {
    return 65504.0F;
  }
};

}  // namespace Kokkos
#endif
#endif
