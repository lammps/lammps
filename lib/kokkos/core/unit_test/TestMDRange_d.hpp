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

#include <TestMDRange.hpp>

namespace Test {

TEST(TEST_CATEGORY, mdrange_3d) {
  TestMDRange_3D<TEST_EXECSPACE>::test_for3(1, 10, 100);
  TestMDRange_3D<TEST_EXECSPACE>::test_for3(100, 10, 100);
#if !defined(KOKKOS_ENABLE_ROCM)  // MDRange Reduced explicitly handled in its
                                  // own cpp file
  TestMDRange_3D<TEST_EXECSPACE>::test_reduce3(1, 10, 100);
  TestMDRange_3D<TEST_EXECSPACE>::test_reduce3(100, 10, 100);
#endif
}

TEST(TEST_CATEGORY, mdrange_neg_idx) {
  TestMDRange_2D_NegIdx<TEST_EXECSPACE>::test_2D_negidx(128, 32);
  TestMDRange_3D_NegIdx<TEST_EXECSPACE>::test_3D_negidx(128, 32, 8);
  TestMDRange_4D_NegIdx<TEST_EXECSPACE>::test_4D_negidx(128, 32, 8, 8);
  TestMDRange_5D_NegIdx<TEST_EXECSPACE>::test_5D_negidx(128, 32, 8, 8, 4);
  TestMDRange_6D_NegIdx<TEST_EXECSPACE>::test_6D_negidx(128, 32, 8, 8, 4, 2);
}

}  // namespace Test
