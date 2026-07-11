// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif

#if !defined(KOKKOS_ENABLE_CUDA) || defined(__CUDACC__)

#include <TestDefaultDeviceType_Category.hpp>
#include <TestReduceCombinatorical.hpp>

namespace Test {

TEST(defaultdevicetype, reduce_instantiation_c2) {
  TestReduceCombinatoricalInstantiation<>::execute_c2();
}

}  // namespace Test

#endif
