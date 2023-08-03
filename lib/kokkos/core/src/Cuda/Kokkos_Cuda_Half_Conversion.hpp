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

#ifndef KOKKOS_CUDA_HALF_HPP_
#define KOKKOS_CUDA_HALF_HPP_

#ifdef KOKKOS_IMPL_CUDA_HALF_TYPE_DEFINED

#include <Kokkos_Half.hpp>
#include <Kokkos_NumericTraits.hpp>  // reduction_identity

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
// architecture is not Ampere
#if CUDA_VERSION >= 11000 && \
    (CUDA_VERSION < 11010 || \
     !(defined(KOKKOS_ARCH_AMPERE) || defined(KOKKOS_ARCH_HOPPER)))
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

#if CUDA_VERSION >= 11010 && \
    ((defined(KOKKOS_ARCH_AMPERE) || defined(KOKKOS_ARCH_HOPPER)))
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
