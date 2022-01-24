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

#ifndef KOKKOS_HALF_HPP_
#define KOKKOS_HALF_HPP_

#include <type_traits>
#include <Kokkos_Macros.hpp>

// Include special backend specific versions here
#include <Cuda/Kokkos_Cuda_Half.hpp>

// Potentially include special compiler specific versions here
// e.g. for Intel

// If none of the above actually did anything and defined a half precision type
// define a fallback implementation here using float
#ifndef KOKKOS_IMPL_HALF_TYPE_DEFINED
#define KOKKOS_IMPL_HALF_TYPE_DEFINED
#define KOKKOS_HALF_T_IS_FLOAT true
namespace Kokkos {
namespace Impl {
struct half_impl_t {
  using type = float;
};
}  // namespace Impl
namespace Experimental {

using half_t = Kokkos::Impl::half_impl_t::type;

// cast_to_half
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(float val) { return half_t(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(bool val) { return half_t(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(double val) { return half_t(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(short val) { return half_t(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned short val) { return half_t(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(int val) { return half_t(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned int val) { return half_t(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(long val) { return half_t(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned long val) { return half_t(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(long long val) { return half_t(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned long long val) { return half_t(val); }

// cast_from_half
// Using an explicit list here too, since the other ones are explicit and for
// example don't include char
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<
    std::is_same<T, float>::value || std::is_same<T, bool>::value ||
        std::is_same<T, double>::value || std::is_same<T, short>::value ||
        std::is_same<T, unsigned short>::value || std::is_same<T, int>::value ||
        std::is_same<T, unsigned int>::value || std::is_same<T, long>::value ||
        std::is_same<T, unsigned long>::value ||
        std::is_same<T, long long>::value ||
        std::is_same<T, unsigned long long>::value,
    T>
cast_from_half(half_t val) {
  return T(val);
}

}  // namespace Experimental
}  // namespace Kokkos

#else
#define KOKKOS_HALF_T_IS_FLOAT false
#endif  // KOKKOS_IMPL_HALF_TYPE_DEFINED
#endif  // KOKKOS_HALF_HPP_
