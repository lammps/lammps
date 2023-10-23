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
#include <TestReducers.hpp>
#include <TestNonTrivialScalarTypes.hpp>

namespace Test {
TEST(TEST_CATEGORY, reducers_complex_double) {
  TestReducers<Kokkos::complex<double>, TEST_EXECSPACE>::execute_basic();
}

TEST(TEST_CATEGORY, reducers_struct) {
  TestReducers<array_reduce<float, 1>, TEST_EXECSPACE>::test_sum(1031);
  TestReducers<array_reduce<float, 2>, TEST_EXECSPACE>::test_sum(1031);
  TestReducers<array_reduce<float, 4>, TEST_EXECSPACE>::test_sum(1031);
  // FIXME_OPENMPTARGET - The size of data in array_reduce has to be a power of
  // 2 for OPENMPTARGET backend in Release and RelWithDebInfo builds.
#ifdef KOKKOS_ENABLE_OPENMPTARGET
  TestReducers<array_reduce<float, 8>, TEST_EXECSPACE>::test_sum(1031);
#else
  TestReducers<array_reduce<float, 3>, TEST_EXECSPACE>::test_sum(1031);
  TestReducers<array_reduce<float, 7>, TEST_EXECSPACE>::test_sum(1031);
#endif
}

TEST(TEST_CATEGORY, reducers_half_t) {
  using ThisTestType = Kokkos::Experimental::half_t;
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(2);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(101);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(202);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(303);

  TestReducers<ThisTestType, TEST_EXECSPACE>::test_prod(5);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_prod(10);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_prod(15);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_prod(20);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_prod(25);
}

TEST(TEST_CATEGORY, reducers_bhalf_t) {
  using ThisTestType = Kokkos::Experimental::bhalf_t;

  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(2);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(25);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(50);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(51);

  TestReducers<ThisTestType, TEST_EXECSPACE>::test_prod(1);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_prod(2);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_prod(3);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_prod(4);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_prod(5);
}

TEST(TEST_CATEGORY, reducers_int8_t) {
  using ThisTestType = int8_t;

  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(1);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(2);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(3);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(4);

  TestReducers<ThisTestType, TEST_EXECSPACE>::test_prod(1);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_prod(2);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_prod(3);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_prod(4);
}

#if !defined(KOKKOS_ENABLE_HIP) && !defined(KOKKOS_ENABLE_OPENMPTARGET)
// TODO - resolve: "Kokkos_HIP_Vectorization.hpp:80:15: error: call to
//                 implicitly-deleted default constructor of 'conv_type'
//                   conv_type tmp_in;"
//
// TODO - resolve:  4: [  FAILED  ] openmptarget.reducers_point_t (1 ms)
TEST(TEST_CATEGORY, reducers_point_t) {
  using ThisTestType = point_t;

  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(1);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(2);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(3);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(4);
}
#endif  // !KOKKOS_ENABLE_HIP && !KOKKOS_ENABLE_OPENMPTARGET

}  // namespace Test
