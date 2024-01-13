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

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

#if 1 != ((defined(KOKKOS_COMPILER_INTEL) ? 1 : 0) +      \
          (defined(KOKKOS_COMPILER_INTEL_LLVM) ? 1 : 0) + \
          (defined(KOKKOS_COMPILER_CRAYC) ? 1 : 0) +      \
          (defined(KOKKOS_COMPILER_CRAY_LLVM) ? 1 : 0) +  \
          (defined(KOKKOS_COMPILER_APPLECC) ? 1 : 0) +    \
          (defined(KOKKOS_COMPILER_CLANG) ? 1 : 0) +      \
          (defined(KOKKOS_COMPILER_GNU) ? 1 : 0) +        \
          (defined(KOKKOS_COMPILER_NVHPC) ? 1 : 0) +      \
          (defined(KOKKOS_COMPILER_MSVC) ? 1 : 0))
#error "Only one host compiler macro can be defined"
#endif

#if defined(KOKKOS_ENABLE_CUDA) && !defined(KOKKOS_ENABLE_CUDA_LAMBDA)
#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
#error "Macro bug: KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA shouldn't be defined"
#endif
#else
#if !defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
#error "Macro bug: KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA should be defined"
#endif
#endif

namespace TestCompilerMacros {

template <class DEVICE_TYPE>
struct AddFunctor {
  using execution_space = DEVICE_TYPE;
  using type            = typename Kokkos::View<int**, execution_space>;
  type a, b;
  int length;

  AddFunctor(type a_, type b_) : a(a_), b(b_), length(a.extent(1)) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_LOOPCOUNT
#pragma loop_count(128)
#endif
    for (int j = 0; j < length; j++) {
      a(i, j) += b(i, j);
    }
  }
};

template <class DeviceType>
bool Test() {
  using type = typename Kokkos::View<int**, DeviceType>;
  type a("A", 1024, 128);
  type b("B", 1024, 128);

  AddFunctor<DeviceType> f(a, b);
  Kokkos::parallel_for(1024, f);
  DeviceType().fence();

  return true;
}

}  // namespace TestCompilerMacros

namespace Test {
TEST(defaultdevicetype, compiler_macros) {
  ASSERT_TRUE((TestCompilerMacros::Test<Kokkos::DefaultHostExecutionSpace>()));
}
}  // namespace Test
