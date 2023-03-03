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

#ifndef KOKKOS_ALGORITHMS_UNITTESTS_TESTSORT_COMMON_HPP
#define KOKKOS_ALGORITHMS_UNITTESTS_TESTSORT_COMMON_HPP

#include <TestSort.hpp>

namespace Test {
TEST(TEST_CATEGORY, SortUnsigned) {
  Impl::test_sort<TEST_EXECSPACE, unsigned>(171);
}
}  // namespace Test
#endif
