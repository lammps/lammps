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

TEST(TEST_CATEGORY, mdrange_6d) {
  TestMDRange_6D<TEST_EXECSPACE>::test_for6(10, 10, 10, 10, 5, 5);
#ifndef KOKKOS_ENABLE_OPENMPTARGET
  // FIXME_OPENMPTARGET requires MDRange parallel_reduce
  TestMDRange_6D<TEST_EXECSPACE>::test_reduce6(100, 10, 10, 10, 5, 5);
#endif
}

}  // namespace Test
