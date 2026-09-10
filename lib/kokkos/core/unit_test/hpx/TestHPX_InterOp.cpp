// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include <TestHPX_Category.hpp>

namespace Test {

// Test whether allocations survive Kokkos initialize/finalize if done via Raw
// HPX.
TEST(hpx, raw_hpx_interop) {
  // FIXME_HPX
  Kokkos::initialize();
  Kokkos::finalize();
}
}  // namespace Test
