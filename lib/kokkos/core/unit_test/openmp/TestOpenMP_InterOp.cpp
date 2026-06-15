// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include <TestOpenMP_Category.hpp>
#include <omp.h>

namespace Test {

// Test whether allocations survive Kokkos initialize/finalize if done via Raw
// Cuda.
TEST(openmp, raw_openmp_interop) {
  int count = 0;
  int num_threads, concurrency;
#pragma omp parallel
  {
#pragma omp atomic
    count++;
    if (omp_get_thread_num() == 0) num_threads = omp_get_num_threads();
  }

  ASSERT_EQ(count, num_threads);

  Kokkos::initialize();

  count = 0;
#pragma omp parallel
  {
#pragma omp atomic
    count++;
  }

  concurrency = Kokkos::OpenMP().concurrency();
  ASSERT_EQ(count, concurrency);

  Kokkos::finalize();

  count = 0;
#pragma omp parallel
  {
#pragma omp atomic
    count++;
  }

  ASSERT_EQ(count, concurrency);
}
}  // namespace Test
