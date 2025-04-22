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

#include <gtest/gtest.h>

#include <cstdlib>

#include "KokkosExecutionEnvironmentNeverInitializedFixture.hpp"

namespace {

using Legion_DeathTest = KokkosExecutionEnvironmentNeverInitialized;

struct ReductionFunctor {
  Kokkos::View<int*> d;

  KOKKOS_FUNCTION void operator()(int i, int& sum) const { sum += d(i); }
};

// The purpose of this test is to mimic Legion's use case of initializing and
// finalizing individual backends
TEST(Legion_DeathTest, individual_backend_initialization) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  EXPECT_EXIT(
      {
        bool success = true;
        Kokkos::InitializationSettings kokkos_init_settings;

        Kokkos::Impl::pre_initialize(kokkos_init_settings);

        // We need to have a host execution space initialized first.
        Kokkos::DefaultHostExecutionSpace::impl_initialize(
            kokkos_init_settings);

        if (!std::is_same_v<Kokkos::DefaultExecutionSpace,
                            Kokkos::DefaultHostExecutionSpace>)
          Kokkos::DefaultExecutionSpace::impl_initialize(kokkos_init_settings);

        Kokkos::Impl::post_initialize(kokkos_init_settings);

        success &= Kokkos::is_initialized();

        {
          Kokkos::View<int*> d("d", 1000);
          Kokkos::deep_copy(d, 1);
          int result;
          Kokkos::parallel_reduce("TestRed", d.extent(0), ReductionFunctor{d},
                                  result);
          success &= (result == d.extent_int(0));
        }

        Kokkos::Impl::pre_finalize();
        if (!std::is_same_v<Kokkos::DefaultExecutionSpace,
                            Kokkos::DefaultHostExecutionSpace>)
          Kokkos::DefaultExecutionSpace::impl_finalize();
        Kokkos::DefaultHostExecutionSpace::impl_finalize();
        Kokkos::Impl::post_finalize();

        success &= !Kokkos::is_initialized();
        success &= Kokkos::is_finalized();
        std::exit(success ? EXIT_SUCCESS : EXIT_FAILURE);
      },
      ::testing::ExitedWithCode(EXIT_SUCCESS), "");
}

}  // namespace
