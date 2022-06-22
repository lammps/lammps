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
#ifndef KOKKOS_MATHEMATICAL_CONSTANTS_HPP
#define KOKKOS_MATHEMATICAL_CONSTANTS_HPP

#include <Kokkos_Macros.hpp>
#include <type_traits>

namespace Kokkos {
namespace Experimental {

#if defined(KOKKOS_ENABLE_CXX17)
#define KOKKOS_IMPL_MATH_CONSTANT(TRAIT, VALUE) \
  template <class T>                            \
  inline constexpr auto TRAIT##_v =             \
      std::enable_if_t<std::is_floating_point_v<T>, T>(VALUE)
#else
#define KOKKOS_IMPL_MATH_CONSTANT(TRAIT, VALUE) \
  template <class T>                            \
  constexpr auto TRAIT##_v =                    \
      std::enable_if_t<std::is_floating_point<T>::value, T>(VALUE)
#endif

// clang-format off
KOKKOS_IMPL_MATH_CONSTANT(e,          2.718281828459045235360287471352662498L);
KOKKOS_IMPL_MATH_CONSTANT(log2e,      1.442695040888963407359924681001892137L);
KOKKOS_IMPL_MATH_CONSTANT(log10e,     0.434294481903251827651128918916605082L);
KOKKOS_IMPL_MATH_CONSTANT(pi,         3.141592653589793238462643383279502884L);
KOKKOS_IMPL_MATH_CONSTANT(inv_pi,     0.318309886183790671537767526745028724L);
KOKKOS_IMPL_MATH_CONSTANT(inv_sqrtpi, 0.564189583547756286948079451560772586L);
KOKKOS_IMPL_MATH_CONSTANT(ln2,        0.693147180559945309417232121458176568L);
KOKKOS_IMPL_MATH_CONSTANT(ln10,       2.302585092994045684017991454684364208L);
KOKKOS_IMPL_MATH_CONSTANT(sqrt2,      1.414213562373095048801688724209698079L);
KOKKOS_IMPL_MATH_CONSTANT(sqrt3,      1.732050807568877293527446341505872367L);
KOKKOS_IMPL_MATH_CONSTANT(inv_sqrt3,  0.577350269189625764509148780501957456L);
KOKKOS_IMPL_MATH_CONSTANT(egamma,     0.577215664901532860606512090082402431L);
KOKKOS_IMPL_MATH_CONSTANT(phi,        1.618033988749894848204586834365638118L);
// clang-format on

#undef KOKKOS_IMPL_MATH_CONSTANT

}  // namespace Experimental
}  // namespace Kokkos
#endif
