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

namespace {

struct ReductionFunctor {
  Kokkos::View<int*> d;

  KOKKOS_FUNCTION void operator()(int i, int& sum) const { sum += d(i); }
};

// The purpose of this test is to mimic Legion's use case of initializing and
// finalizing individual backends
TEST(initialization, legion_initialization) {
  Kokkos::InitializationSettings kokkos_init_settings;

  Kokkos::Impl::pre_initialize(kokkos_init_settings);

  // We need to have a host execution space initialized first.
  Kokkos::DefaultHostExecutionSpace::impl_initialize(kokkos_init_settings);

  if (!std::is_same_v<Kokkos::DefaultExecutionSpace,
                      Kokkos::DefaultHostExecutionSpace>)
    Kokkos::DefaultExecutionSpace::impl_initialize(kokkos_init_settings);

  Kokkos::Impl::post_initialize(kokkos_init_settings);

  EXPECT_TRUE(Kokkos::is_initialized());

  {
    Kokkos::View<int*> d("d", 1000);
    Kokkos::deep_copy(d, 1);
    int result;
    Kokkos::parallel_reduce("TestRed", d.extent(0), ReductionFunctor{d},
                            result);
    EXPECT_EQ(result, d.extent_int(0));
  }

  Kokkos::Impl::pre_finalize();
  if (!std::is_same_v<Kokkos::DefaultExecutionSpace,
                      Kokkos::DefaultHostExecutionSpace>)
    Kokkos::DefaultExecutionSpace::impl_finalize();
  Kokkos::DefaultHostExecutionSpace::impl_finalize();
  Kokkos::Impl::post_finalize();

  EXPECT_FALSE(Kokkos::is_initialized());
  EXPECT_TRUE(Kokkos::is_finalized());
}
}  // namespace
