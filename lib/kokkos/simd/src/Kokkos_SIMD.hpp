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

#ifndef KOKKOS_SIMD_HPP
#define KOKKOS_SIMD_HPP

#include <Kokkos_SIMD_Common.hpp>

#include <Kokkos_SIMD_Scalar.hpp>

#ifdef KOKKOS_ARCH_AVX2
#include <Kokkos_SIMD_AVX2.hpp>
#endif

#ifdef KOKKOS_ARCH_AVX512XEON
#include <Kokkos_SIMD_AVX512.hpp>
#endif

#ifdef __ARM_NEON
#include <Kokkos_SIMD_NEON.hpp>
#endif

namespace Kokkos {
namespace Experimental {

namespace simd_abi {

namespace Impl {

#if defined(KOKKOS_ARCH_AVX512XEON)
using host_native = avx512_fixed_size<8>;
#elif defined(KOKKOS_ARCH_AVX2)
using host_native  = avx2_fixed_size<4>;
#elif defined(__ARM_NEON)
using host_native  = neon_fixed_size<2>;
#else
using host_native   = scalar;
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

#ifdef KOKKOS_ENABLE_OPENACC
template <>
struct ForSpace<Kokkos::Experimental::OpenACC> {
  using type = scalar;
};
#endif

#ifdef KOKKOS_ENABLE_HIP
template <>
struct ForSpace<Kokkos::HIP> {
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

template <typename... Ts>
class data_types {};

#if defined(KOKKOS_ARCH_AVX512XEON)
using host_abi_set  = abi_set<simd_abi::scalar, simd_abi::avx512_fixed_size<8>>;
using data_type_set = data_types<std::int32_t, std::uint32_t, std::int64_t,
                                 std::uint64_t, double>;
#elif defined(KOKKOS_ARCH_AVX2)
using host_abi_set = abi_set<simd_abi::scalar, simd_abi::avx2_fixed_size<4>>;
using data_type_set =
    data_types<std::int32_t, std::int64_t, std::uint64_t, double>;
#elif defined(__ARM_NEON)
using host_abi_set = abi_set<simd_abi::scalar, simd_abi::neon_fixed_size<2>>;
using data_type_set =
    data_types<std::int32_t, std::int64_t, std::uint64_t, double>;
#else
using host_abi_set  = abi_set<simd_abi::scalar>;
using data_type_set = data_types<std::int32_t, std::uint32_t, std::int64_t,
                                 std::uint64_t, double>;
#endif

using device_abi_set = abi_set<simd_abi::scalar>;

}  // namespace Impl

}  // namespace Experimental
}  // namespace Kokkos

#endif
