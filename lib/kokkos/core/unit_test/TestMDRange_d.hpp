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

#include <TestMDRange.hpp>

namespace Test {

TEST(TEST_CATEGORY, mdrange_3d) {
  TestMDRange_3D<TEST_EXECSPACE>::test_for3(1, 10, 100);
  TestMDRange_3D<TEST_EXECSPACE>::test_for3(100, 10, 100);
#ifndef KOKKOS_ENABLE_OPENMPTARGET
  // FIXME_OPENMPTARGET requires MDRange parallel_reduce
  TestMDRange_3D<TEST_EXECSPACE>::test_reduce3(1, 10, 100);
  TestMDRange_3D<TEST_EXECSPACE>::test_reduce3(100, 10, 100);
#endif
}

#ifndef KOKKOS_ENABLE_OPENMPTARGET
TEST(TEST_CATEGORY, mdrange_neg_idx) {
  TestMDRange_2D_NegIdx<TEST_EXECSPACE>::test_2D_negidx(128, 32);
  TestMDRange_3D_NegIdx<TEST_EXECSPACE>::test_3D_negidx(128, 32, 8);
  TestMDRange_4D_NegIdx<TEST_EXECSPACE>::test_4D_negidx(128, 32, 8, 8);
  TestMDRange_5D_NegIdx<TEST_EXECSPACE>::test_5D_negidx(128, 32, 8, 8, 4);
  TestMDRange_6D_NegIdx<TEST_EXECSPACE>::test_6D_negidx(128, 32, 8, 8, 4, 2);
}
#endif

}  // namespace Test
