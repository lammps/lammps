// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_CUDA_HALF_HPP_
#define KOKKOS_CUDA_HALF_HPP_

#include <Kokkos_Half.hpp>
#include <impl/Kokkos_NvidiaGpuArchitectures.hpp>

#include <cuda_bf16.h>

namespace Kokkos::Experimental {

/************************** half conversions **********************************/
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(half_t val) { return val; }

KOKKOS_INLINE_FUNCTION
half_t cast_to_half(float val) { return __float2half(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(double val) { return __double2half(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(short val) { return __short2half_rn(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned short val) { return __ushort2half_rn(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(int val) { return __int2half_rn(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned int val) { return __uint2half_rn(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(long long val) { return __ll2half_rn(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned long long val) { return __ull2half_rn(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(long val) {
  return cast_to_half(static_cast<long long>(val));
}
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned long val) {
  return cast_to_half(static_cast<unsigned long long>(val));
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, float>, T>
cast_from_half(half_t val) {
  return __half2float(__half(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, double>, T>
cast_from_half(half_t val) {
  return static_cast<double>(__half2float(__half(val)));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, short>, T>
cast_from_half(half_t val) {
  return __half2short_rz(__half(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, unsigned short>, T>
cast_from_half(half_t val) {
  return __half2ushort_rz(__half(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, int>, T>
cast_from_half(half_t val) {
  return __half2int_rz(__half(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, unsigned int>, T>
cast_from_half(half_t val) {
  return __half2uint_rz(__half(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, long long>, T>
cast_from_half(half_t val) {
  return __half2ll_rz(__half(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same_v<T, unsigned long long>, T>
    cast_from_half(half_t val) {
  return __half2ull_rz(__half(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, long>, T>
cast_from_half(half_t val) {
  return static_cast<T>(cast_from_half<long long>(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, unsigned long>, T>
cast_from_half(half_t val) {
  return static_cast<T>(cast_from_half<unsigned long long>(val));
}

/************************** bhalf conversions *********************************/
// if architecture is older than Ampere
#if KOKKOS_IMPL_ARCH_NVIDIA_GPU < 80
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(bhalf_t val) { return val; }

KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(float val) { return bhalf_t(__float2bfloat16(val)); }

KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(bool val) {
  return cast_to_bhalf(static_cast<float>(val));
}

KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(double val) {
  // double2bfloat16 was only introduced in CUDA 11 too
  return bhalf_t(__float2bfloat16(static_cast<float>(val)));
}

KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(short val) {
  return bhalf_t(__float2bfloat16(static_cast<float>(val)));
}

KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(unsigned short val) {
  return bhalf_t(__float2bfloat16(static_cast<float>(val)));
}

KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(int val) {
  return bhalf_t(__float2bfloat16(static_cast<float>(val)));
}

KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(unsigned int val) {
  return bhalf_t(__float2bfloat16(static_cast<float>(val)));
}

KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(long long val) {
  return bhalf_t(__float2bfloat16(static_cast<float>(val)));
}

KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(unsigned long long val) {
  return bhalf_t(__float2bfloat16(static_cast<float>(val)));
}

KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(long val) {
  return cast_to_bhalf(static_cast<long long>(val));
}

KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(unsigned long val) {
  return cast_to_bhalf(static_cast<unsigned long long>(val));
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, float>, T>
cast_from_bhalf(bhalf_t val) {
  return __bfloat162float(bhalf_t::impl_type(val));
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, bool>, T>
cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(cast_from_bhalf<float>(val));
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, double>, T>
cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(__bfloat162float(bhalf_t::impl_type(val)));
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, short>, T>
cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(__bfloat162float(bhalf_t::impl_type(val)));
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, unsigned short>, T>
cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(__bfloat162float(bhalf_t::impl_type(val)));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, int>, T>
cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(__bfloat162float(bhalf_t::impl_type(val)));
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, unsigned>, T>
cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(__bfloat162float(bhalf_t::impl_type(val)));
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, long long>, T>
cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(__bfloat162float(bhalf_t::impl_type(val)));
}

template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same_v<T, unsigned long long>, T>
    cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(__bfloat162float(bhalf_t::impl_type(val)));
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, long>, T>
cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(cast_from_bhalf<long long>(val));
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, unsigned long>, T>
cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(cast_from_bhalf<unsigned long long>(val));
}
#else
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(bhalf_t val) { return val; }
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(float val) { return __float2bfloat16(val); }
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(double val) { return __double2bfloat16(val); }
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(short val) { return __short2bfloat16_rn(val); }
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(unsigned short val) { return __ushort2bfloat16_rn(val); }
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(int val) { return __int2bfloat16_rn(val); }
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(unsigned int val) { return __uint2bfloat16_rn(val); }
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(long long val) { return __ll2bfloat16_rn(val); }
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(unsigned long long val) { return __ull2bfloat16_rn(val); }
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(long val) {
  return cast_to_bhalf(static_cast<long long>(val));
}
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(unsigned long val) {
  return cast_to_bhalf(static_cast<unsigned long long>(val));
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, float>, T>
cast_from_bhalf(bhalf_t val) {
  return __bfloat162float(__nv_bfloat16(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, double>, T>
cast_from_bhalf(bhalf_t val) {
  return static_cast<double>(__bfloat162float(__nv_bfloat16(val)));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, short>, T>
cast_from_bhalf(bhalf_t val) {
  return __bfloat162short_rz(__nv_bfloat16(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, unsigned short>, T>
cast_from_bhalf(bhalf_t val) {
  return __bfloat162ushort_rz(__nv_bfloat16(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, int>, T>
cast_from_bhalf(bhalf_t val) {
  return __bfloat162int_rz(__nv_bfloat16(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, unsigned int>, T>
cast_from_bhalf(bhalf_t val) {
  return __bfloat162uint_rz(__nv_bfloat16(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, long long>, T>
cast_from_bhalf(bhalf_t val) {
  return __bfloat162ll_rz(__nv_bfloat16(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same_v<T, unsigned long long>, T>
    cast_from_bhalf(bhalf_t val) {
  return __bfloat162ull_rz(__nv_bfloat16(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, long>, T>
cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(cast_from_bhalf<long long>(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same_v<T, unsigned long>, T>
cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(cast_from_bhalf<unsigned long long>(val));
}
#endif

}  // namespace Kokkos::Experimental

#endif
