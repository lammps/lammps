/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.5
//       Copyright (2022) National Technology & Engineering
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

#ifndef KOKKOS_SYCL_HALF_HPP_
#define KOKKOS_SYCL_HALF_HPP_

#ifdef KOKKOS_IMPL_SYCL_HALF_TYPE_DEFINED

#include <Kokkos_Half.hpp>
#include <Kokkos_NumericTraits.hpp>  // reduction_identity

namespace Kokkos {
namespace Experimental {

/************************** half conversions **********************************/
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(half_t val) { return val; }

KOKKOS_INLINE_FUNCTION
half_t cast_to_half(float val) { return half_t::impl_type(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(double val) { return half_t::impl_type(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(short val) { return half_t::impl_type(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned short val) { return half_t::impl_type(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(int val) { return half_t::impl_type(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned int val) { return half_t::impl_type(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(long long val) { return half_t::impl_type(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned long long val) { return half_t::impl_type(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(long val) { return half_t::impl_type(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned long val) { return half_t::impl_type(val); }

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, float>::value, T>
cast_from_half(half_t val) {
  return half_t::impl_type(val);
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, double>::value, T>
cast_from_half(half_t val) {
  return half_t::impl_type(val);
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, short>::value, T>
cast_from_half(half_t val) {
  return half_t::impl_type(val);
}
template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same<T, unsigned short>::value, T>
    cast_from_half(half_t val) {
  return half_t::impl_type(val);
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, int>::value, T>
cast_from_half(half_t val) {
  return half_t::impl_type(val);
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, unsigned int>::value, T>
cast_from_half(half_t val) {
  return half_t::impl_type(val);
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, long long>::value, T>
cast_from_half(half_t val) {
  return half_t::impl_type(val);
}
template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same<T, unsigned long long>::value, T>
    cast_from_half(half_t val) {
  return half_t::impl_type(val);
}
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, long>::value, T>
cast_from_half(half_t val) {
  return half_t::impl_type(val);
}
template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same<T, unsigned long>::value, T>
    cast_from_half(half_t val) {
  return half_t::impl_type(val);
}
}  // namespace Experimental

template <>
struct reduction_identity<Kokkos::Experimental::half_t> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static Kokkos::Experimental::half_t
  sum() noexcept {
    return Kokkos::Experimental::half_t::impl_type(0.0F);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static Kokkos::Experimental::half_t
  prod() noexcept {
    return Kokkos::Experimental::half_t::impl_type(1.0F);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static Kokkos::Experimental::half_t
  max() noexcept {
    return std::numeric_limits<
        Kokkos::Experimental::half_t::impl_type>::lowest();
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static Kokkos::Experimental::half_t
  min() noexcept {
    return std::numeric_limits<Kokkos::Experimental::half_t::impl_type>::max();
  }
};

}  // namespace Kokkos
#endif  // KOKKOS_IMPL_SYCL_HALF_TYPE_DEFINED
#endif
