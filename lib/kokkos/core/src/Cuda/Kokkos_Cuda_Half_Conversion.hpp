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

#ifndef KOKKOS_CUDA_HALF_HPP_
#define KOKKOS_CUDA_HALF_HPP_

#ifdef KOKKOS_IMPL_CUDA_HALF_TYPE_DEFINED

#include <Kokkos_Half.hpp>
#include <Kokkos_ReductionIdentity.hpp>

#if CUDA_VERSION >= 11000
#include <cuda_bf16.h>
#endif

namespace Kokkos {
namespace Experimental {

/************************** half conversions **********************************/
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(half_t val) { return val; }

// CUDA before 11.1 only has the half <-> float conversions marked host device
// So we will largely convert to float on the host for conversion
// But still call the correct functions on the device
#if (CUDA_VERSION < 11010)

KOKKOS_INLINE_FUNCTION
half_t cast_to_half(float val) { return half_t(__float2half(val)); }

KOKKOS_INLINE_FUNCTION
half_t cast_to_half(bool val) { return cast_to_half(static_cast<float>(val)); }

KOKKOS_INLINE_FUNCTION
half_t cast_to_half(double val) {
  // double2half was only introduced in CUDA 11 too
  return half_t(__float2half(static_cast<float>(val)));
}

KOKKOS_INLINE_FUNCTION
half_t cast_to_half(short val) {
  KOKKOS_IF_ON_DEVICE((return half_t(__short2half_rn(val));))
  KOKKOS_IF_ON_HOST((return half_t(__float2half(static_cast<float>(val)));))
}

KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned short val) {
  KOKKOS_IF_ON_DEVICE((return half_t(__ushort2half_rn(val));))
  KOKKOS_IF_ON_HOST((return half_t(__float2half(static_cast<float>(val)));))
}

KOKKOS_INLINE_FUNCTION
half_t cast_to_half(int val) {
  KOKKOS_IF_ON_DEVICE((return half_t(__int2half_rn(val));))
  KOKKOS_IF_ON_HOST((return half_t(__float2half(static_cast<float>(val)));))
}

KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned int val) {
  KOKKOS_IF_ON_DEVICE((return half_t(__uint2half_rn(val));))
  KOKKOS_IF_ON_HOST((return half_t(__float2half(static_cast<float>(val)));))
}

KOKKOS_INLINE_FUNCTION
half_t cast_to_half(long long val) {
  KOKKOS_IF_ON_DEVICE((return half_t(__ll2half_rn(val));))
  KOKKOS_IF_ON_HOST((return half_t(__float2half(static_cast<float>(val)));))
}

KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned long long val) {
  KOKKOS_IF_ON_DEVICE((return half_t(__ull2half_rn(val));))
  KOKKOS_IF_ON_HOST((return half_t(__float2half(static_cast<float>(val)));))
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
  KOKKOS_IF_ON_DEVICE((return __half2short_rz(half_t::impl_type(val));))
  KOKKOS_IF_ON_HOST(
      (return static_cast<T>(__half2float(half_t::impl_type(val)));))
}

template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same<T, unsigned short>::value, T>
    cast_from_half(half_t val) {
  KOKKOS_IF_ON_DEVICE((return __half2ushort_rz(half_t::impl_type(val));))
  KOKKOS_IF_ON_HOST(
      (return static_cast<T>(__half2float(half_t::impl_type(val)));))
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, int>::value, T>
cast_from_half(half_t val) {
  KOKKOS_IF_ON_DEVICE((return __half2int_rz(half_t::impl_type(val));))
  KOKKOS_IF_ON_HOST(
      (return static_cast<T>(__half2float(half_t::impl_type(val)));))
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, unsigned>::value, T>
cast_from_half(half_t val) {
  KOKKOS_IF_ON_DEVICE((return __half2uint_rz(half_t::impl_type(val));))
  KOKKOS_IF_ON_HOST(
      (return static_cast<T>(__half2float(half_t::impl_type(val)));))
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, long long>::value, T>
cast_from_half(half_t val) {
  KOKKOS_IF_ON_DEVICE((return __half2ll_rz(half_t::impl_type(val));))
  KOKKOS_IF_ON_HOST(
      (return static_cast<T>(__half2float(half_t::impl_type(val)));))
}

template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same<T, unsigned long long>::value, T>
    cast_from_half(half_t val) {
  KOKKOS_IF_ON_DEVICE((return __half2ull_rz(half_t::impl_type(val));))
  KOKKOS_IF_ON_HOST(
      (return static_cast<T>(__half2float(half_t::impl_type(val)));))
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

#else  // CUDA 11.1 versions follow

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
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, float>::value, T>
cast_from_half(half_t val) {
  return __half2float(__half(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, double>::value, T>
cast_from_half(half_t val) {
  return static_cast<double>(__half2float(__half(val)));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, short>::value, T>
cast_from_half(half_t val) {
  return __half2short_rz(__half(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same<T, unsigned short>::value, T>
    cast_from_half(half_t val) {
  return __half2ushort_rz(__half(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, int>::value, T>
cast_from_half(half_t val) {
  return __half2int_rz(__half(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, unsigned int>::value, T>
cast_from_half(half_t val) {
  return __half2uint_rz(__half(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, long long>::value, T>
cast_from_half(half_t val) {
  return __half2ll_rz(__half(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same<T, unsigned long long>::value, T>
    cast_from_half(half_t val) {
  return __half2ull_rz(__half(val));
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
#endif

/************************** bhalf conversions *********************************/
// Go in this branch if CUDA version is >= 11.0.0 and less than 11.1.0 or if the
// architecture is older than Ampere
#if !defined(KOKKOS_ARCH_KEPLER) && !defined(KOKKOS_ARCH_MAXWELL) && \
    !defined(KOKKOS_ARCH_PASCAL) && !defined(KOKKOS_ARCH_VOLTA) &&   \
    !defined(KOKKOS_ARCH_TURING75)
#define KOKKOS_IMPL_NVIDIA_GPU_ARCH_SUPPORT_BHALF
#endif

#if CUDA_VERSION >= 11000 && \
    (CUDA_VERSION < 11010 || \
     !defined(KOKKOS_IMPL_NVIDIA_GPU_ARCH_SUPPORT_BHALF))
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
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, float>::value, T>
cast_from_bhalf(bhalf_t val) {
  return __bfloat162float(bhalf_t::impl_type(val));
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, bool>::value, T>
cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(cast_from_bhalf<float>(val));
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, double>::value, T>
cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(__bfloat162float(bhalf_t::impl_type(val)));
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, short>::value, T>
cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(__bfloat162float(bhalf_t::impl_type(val)));
}

template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same<T, unsigned short>::value, T>
    cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(__bfloat162float(bhalf_t::impl_type(val)));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, int>::value, T>
cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(__bfloat162float(bhalf_t::impl_type(val)));
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, unsigned>::value, T>
cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(__bfloat162float(bhalf_t::impl_type(val)));
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, long long>::value, T>
cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(__bfloat162float(bhalf_t::impl_type(val)));
}

template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same<T, unsigned long long>::value, T>
    cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(__bfloat162float(bhalf_t::impl_type(val)));
}

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, long>::value, T>
cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(cast_from_bhalf<long long>(val));
}

template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same<T, unsigned long>::value, T>
    cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(cast_from_bhalf<unsigned long long>(val));
}
#endif  // CUDA_VERSION >= 11000 && CUDA_VERSION < 11010

#if CUDA_VERSION >= 11010 && defined(KOKKOS_IMPL_NVIDIA_GPU_ARCH_SUPPORT_BHALF)
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
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, float>::value, T>
cast_from_bhalf(bhalf_t val) {
  return __bfloat162float(__nv_bfloat16(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, double>::value, T>
cast_from_bhalf(bhalf_t val) {
  return static_cast<double>(__bfloat162float(__nv_bfloat16(val)));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, short>::value, T>
cast_from_bhalf(bhalf_t val) {
  return __bfloat162short_rz(__nv_bfloat16(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same<T, unsigned short>::value, T>
    cast_from_bhalf(bhalf_t val) {
  return __bfloat162ushort_rz(__nv_bfloat16(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, int>::value, T>
cast_from_bhalf(bhalf_t val) {
  return __bfloat162int_rz(__nv_bfloat16(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, unsigned int>::value, T>
cast_from_bhalf(bhalf_t val) {
  return __bfloat162uint_rz(__nv_bfloat16(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, long long>::value, T>
cast_from_bhalf(bhalf_t val) {
  return __bfloat162ll_rz(__nv_bfloat16(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same<T, unsigned long long>::value, T>
    cast_from_bhalf(bhalf_t val) {
  return __bfloat162ull_rz(__nv_bfloat16(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, long>::value, T>
cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(cast_from_bhalf<long long>(val));
}
template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same<T, unsigned long>::value, T>
    cast_from_bhalf(bhalf_t val) {
  return static_cast<T>(cast_from_bhalf<unsigned long long>(val));
}
#endif  // CUDA_VERSION >= 11010

#undef KOKKOS_IMPL_NVIDIA_GPU_ARCH_SUPPORT_BHALF
}  // namespace Experimental

#if (CUDA_VERSION >= 11000)
template <>
struct reduction_identity<Kokkos::Experimental::bhalf_t> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static float sum() noexcept {
    return 0.0F;
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static float prod() noexcept {
    return 1.0F;
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static float max() noexcept {
    return -0x7f7f;
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static float min() noexcept {
    return 0x7f7f;
  }
};
#endif  // CUDA_VERSION >= 11000

// use float as the return type for sum and prod since cuda_fp16.h
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
#endif  // KOKKOS_IMPL_CUDA_HALF_TYPE_DEFINED
#endif
