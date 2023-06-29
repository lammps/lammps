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

#include <TestOpenMP_Category.hpp>
#include <Kokkos_Core.hpp>

#include <mutex>

namespace Test {

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
TEST(openmp, partition_master) {
  using Mutex = Kokkos::Experimental::MasterLock<Kokkos::OpenMP>;

  Mutex mtx;
  int errors = 0;

  auto master = [&errors, &mtx](int /*partition_id*/, int /*num_partitions*/) {
    const int pool_size = Kokkos::OpenMP().impl_thread_pool_size();

    {
      std::unique_lock<Mutex> lock(mtx);
      if (Kokkos::OpenMP::in_parallel()) {
        ++errors;
      }
      if (Kokkos::OpenMP::impl_thread_pool_rank() != 0) {
        ++errors;
      }
    }

    {
      int local_errors = 0;
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<Kokkos::OpenMP>(0, 1000),
          [pool_size](const int, int& errs) {
            if (Kokkos::OpenMP().impl_thread_pool_size() != pool_size) {
              ++errs;
            }
          },
          local_errors);
      Kokkos::atomic_add(&errors, local_errors);
    }

    Kokkos::Experimental::UniqueToken<Kokkos::OpenMP> token;

    Kokkos::View<int*, Kokkos::OpenMP> count("", token.size());

    Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::OpenMP>(0, 1000),
                         [=](const int) {
                           int i = token.acquire();
                           ++count[i];
                           token.release(i);
                         });

    Kokkos::View<int, Kokkos::OpenMP> sum("");
    Kokkos::parallel_for(
        Kokkos::RangePolicy<Kokkos::OpenMP>(0, token.size()),
        [=](const int i) { Kokkos::atomic_add(sum.data(), count[i]); });

    if (sum() != 1000) {
      Kokkos::atomic_add(&errors, 1);
    }
  };

  master(0, 1);

  ASSERT_EQ(errors, 0);

  Kokkos::OpenMP::partition_master(master);
  ASSERT_EQ(errors, 0);

  Kokkos::OpenMP::partition_master(master, 4, 0);
  ASSERT_EQ(errors, 0);

  Kokkos::OpenMP::partition_master(master, 0, 4);
  ASSERT_EQ(errors, 0);

  Kokkos::OpenMP::partition_master(master, 2, 2);
  ASSERT_EQ(errors, 0);

  Kokkos::OpenMP::partition_master(master, 8, 0);
  ASSERT_EQ(errors, 0);

  Kokkos::OpenMP::partition_master(master, 0, 8);
  ASSERT_EQ(errors, 0);

  Kokkos::OpenMP::partition_master(master, 8, 8);
  ASSERT_EQ(errors, 0);
}
#endif

}  // namespace Test
