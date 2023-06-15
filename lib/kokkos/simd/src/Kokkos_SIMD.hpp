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

#ifndef KOKKOS_SIMD_HPP
#define KOKKOS_SIMD_HPP

#include <Kokkos_SIMD_Common.hpp>

#include <Kokkos_SIMD_Scalar.hpp>

#ifdef KOKKOS_ARCH_AVX512XEON
#include <Kokkos_SIMD_AVX512.hpp>
#endif

namespace Kokkos {
namespace Experimental {

namespace simd_abi {

namespace Impl {

#if defined(KOKKOS_ARCH_AVX512XEON)
using host_native = avx512_fixed_size<8>;
#else
using host_native  = scalar;
#endif

template <class T>
struct ForSpace;

#ifdef KOKKOS_ENABLE_SERIAL
template <>
struct ForSpace<Kokkos::Serial> {
  using type = host_native;
};
#endif

#ifdef KOKKOS_ENABLE_CUDA
template <>
struct ForSpace<Kokkos::Cuda> {
  using type = scalar;
};
#endif

#ifdef KOKKOS_ENABLE_THREADS
template <>
struct ForSpace<Kokkos::Threads> {
  using type = host_native;
};
#endif

#ifdef KOKKOS_ENABLE_HPX
template <>
struct ForSpace<Kokkos::Experimental::HPX> {
  using type = scalar;
};
#endif

#ifdef KOKKOS_ENABLE_OPENMP
template <>
struct ForSpace<Kokkos::OpenMP> {
  using type = host_native;
};
#endif

#ifdef KOKKOS_ENABLE_OPENMPTARGET
template <>
struct ForSpace<Kokkos::Experimental::OpenMPTarget> {
  using type = scalar;
};
#endif

#ifdef KOKKOS_ENABLE_HIP
template <>
struct ForSpace<Kokkos::Experimental::HIP> {
  using type = scalar;
};
#endif

#ifdef KOKKOS_ENABLE_SYCL
template <>
struct ForSpace<Kokkos::Experimental::SYCL> {
  using type = scalar;
};
#endif

}  // namespace Impl

template <class Space>
using ForSpace = typename Impl::ForSpace<typename Space::execution_space>::type;

template <class T>
using native = ForSpace<Kokkos::DefaultExecutionSpace>;

}  // namespace simd_abi

template <class T>
using native_simd = simd<T, simd_abi::native<T>>;
template <class T>
using native_simd_mask = simd_mask<T, simd_abi::native<T>>;

namespace Impl {

template <class... Abis>
class abi_set {};

#ifdef KOKKOS_ARCH_AVX512XEON
using host_abi_set = abi_set<simd_abi::scalar, simd_abi::avx512_fixed_size<8>>;
#else
using host_abi_set = abi_set<simd_abi::scalar>;
#endif

using device_abi_set = abi_set<simd_abi::scalar>;

}  // namespace Impl

}  // namespace Experimental
}  // namespace Kokkos

#endif
