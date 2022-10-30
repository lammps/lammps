/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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

  TestReducers<ThisTestType, TEST_EXECSPACE>::test_prod(5);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_prod(10);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_prod(15);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_prod(20);
  TestReducers<ThisTestType, TEST_EXECSPACE>::test_prod(25);
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
