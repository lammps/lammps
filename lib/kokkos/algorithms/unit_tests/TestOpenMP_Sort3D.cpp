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

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_OPENMP

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

//----------------------------------------------------------------------------
#include <TestRandom.hpp>
#include <TestSort.hpp>
#include <iomanip>

namespace Test {

TEST(openmp, SortUnsigned3D) {
  Impl::test_3D_sort<Kokkos::OpenMP, unsigned>(171);
}

}  // namespace Test
#else
void KOKKOS_ALGORITHMS_UNITTESTS_TESTOPENMP_PREVENT_LINK_ERROR() {}
#endif
