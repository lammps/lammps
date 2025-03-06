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

#ifndef KOKKOS_TEST_SUBVIEW_B_HPP
#define KOKKOS_TEST_SUBVIEW_B_HPP
#include <TestViewSubview.hpp>

namespace Test {

TEST(TEST_CATEGORY, view_subview_layoutleft_to_layoutleft) {
  TestViewSubview::test_layoutleft_to_layoutleft<TEST_EXECSPACE>();
  TestViewSubview::test_layoutleft_to_layoutleft<
      TEST_EXECSPACE, Kokkos::MemoryTraits<Kokkos::Atomic> >();
  TestViewSubview::test_layoutleft_to_layoutleft<
      TEST_EXECSPACE, Kokkos::MemoryTraits<Kokkos::RandomAccess> >();
}

TEST(TEST_CATEGORY, view_subview_layoutright_to_layoutright) {
  TestViewSubview::test_layoutright_to_layoutright<TEST_EXECSPACE>();
  TestViewSubview::test_layoutright_to_layoutright<
      TEST_EXECSPACE, Kokkos::MemoryTraits<Kokkos::Atomic> >();
  TestViewSubview::test_layoutright_to_layoutright<
      TEST_EXECSPACE, Kokkos::MemoryTraits<Kokkos::RandomAccess> >();
}

}  // namespace Test
#endif
