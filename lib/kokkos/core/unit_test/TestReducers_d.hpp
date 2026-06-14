// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
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
  TestReducers<array_reduce<float, 3>, TEST_EXECSPACE>::test_sum(1031);
  TestReducers<array_reduce<float, 7>, TEST_EXECSPACE>::test_sum(1031);
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

TEST(TEST_CATEGORY, reducers_int16_t) {
  using ThisTestType = int16_t;

  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(1);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(2);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(3);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(4);

  TestReducers<ThisTestType, TEST_EXECSPACE>::test_prod(1);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_prod(2);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_prod(3);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_prod(4);
}

#if !defined(KOKKOS_ENABLE_HIP)
// TODO - resolve: "Kokkos_HIP_Vectorization.hpp:80:15: error: call to
//                 implicitly-deleted default constructor of 'conv_type'
//                   conv_type tmp_in;"
//
TEST(TEST_CATEGORY, reducers_point_t) {
  using ThisTestType = point_t;

  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(1);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(2);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(3);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_sum(4);
}
#endif  // !KOKKOS_ENABLE_HIP

}  // namespace Test
