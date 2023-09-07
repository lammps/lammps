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

#include <Kokkos_Core.hpp>
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
