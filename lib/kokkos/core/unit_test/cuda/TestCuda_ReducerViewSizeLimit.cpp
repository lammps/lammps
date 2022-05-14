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

#include <TestCuda_Category.hpp>
#include <Kokkos_Core.hpp>

namespace Test {

using ValueType = double;
using MemSpace  = Kokkos::CudaSpace;
using Matrix2D  = Kokkos::View<ValueType**, MemSpace>;
using Matrix3D  = Kokkos::View<ValueType***, MemSpace>;
using Vector    = Kokkos::View<ValueType*, MemSpace>;

namespace Impl {

struct ArrayReduceFunctor {
  using value_type = ValueType[];

  int value_count;
  Matrix2D m;

  ArrayReduceFunctor(const Matrix2D& m_) : value_count(m_.extent(1)), m(m_) {}

  KOKKOS_INLINE_FUNCTION void operator()(const int i, value_type sum) const {
    const int numVecs = value_count;
    for (int j = 0; j < numVecs; ++j) {
      sum[j] += m(i, j);
    }
  }

  KOKKOS_INLINE_FUNCTION void init(value_type update) const {
    const int numVecs = value_count;
    for (int j = 0; j < numVecs; ++j) {
      update[j] = 0.0;
    }
  }

  KOKKOS_INLINE_FUNCTION void join(volatile value_type update,
                                   const volatile value_type source) const {
    const int numVecs = value_count;
    for (int j = 0; j < numVecs; ++j) {
      update[j] += source[j];
    }
  }

  KOKKOS_INLINE_FUNCTION void join(value_type update,
                                   const value_type source) const {
    const int numVecs = value_count;
    for (int j = 0; j < numVecs; ++j) {
      update[j] += source[j];
    }
  }

  KOKKOS_INLINE_FUNCTION void final(value_type) const {}
};

struct MDArrayReduceFunctor {
  using value_type = ValueType[];

  int value_count;
  Matrix3D m;

  MDArrayReduceFunctor(const Matrix3D& m_) : value_count(m_.extent(2)), m(m_) {}

  KOKKOS_INLINE_FUNCTION void operator()(const int i, const int j,
                                         value_type sum) const {
    const int numVecs = value_count;
    for (int k = 0; k < numVecs; ++k) {
      sum[k] += m(i, j, k);
    }
  }

  KOKKOS_INLINE_FUNCTION void init(value_type update) const {
    const int numVecs = value_count;
    for (int j = 0; j < numVecs; ++j) {
      update[j] = 0.0;
    }
  }

  KOKKOS_INLINE_FUNCTION void final(value_type) const {}
};

struct ReduceViewSizeLimitTester {
  const ValueType initValue           = 3;
  const size_t nGlobalEntries         = 100;
  const int testViewSize              = 200;
  const size_t expectedInitShmemLimit = 373584;
  const unsigned initBlockSize        = Kokkos::Impl::CudaTraits::WarpSize * 8;

  void run_test_range() {
    Matrix2D matrix;
    Vector sum;

    for (int i = 0; i < testViewSize; ++i) {
      size_t sumInitShmemSize = (initBlockSize + 2) * sizeof(ValueType) * i;

      Kokkos::resize(Kokkos::WithoutInitializing, sum, i);
      Kokkos::resize(Kokkos::WithoutInitializing, matrix, nGlobalEntries, i);
      Kokkos::deep_copy(matrix, initValue);

      auto policy  = Kokkos::RangePolicy<TEST_EXECSPACE>(0, nGlobalEntries);
      auto functor = ArrayReduceFunctor(matrix);

      if (sumInitShmemSize < expectedInitShmemLimit) {
        EXPECT_NO_THROW(Kokkos::parallel_reduce(policy, functor, sum));
      } else {
        EXPECT_THROW(Kokkos::parallel_reduce(policy, functor, sum),
                     std::runtime_error);
      }
    }
  }

  void run_test_md_range_2D() {
    Matrix3D matrix;
    Vector sum;

    for (int i = 0; i < testViewSize; ++i) {
      size_t sumInitShmemSize = (initBlockSize + 2) * sizeof(ValueType) * i;

      Kokkos::resize(Kokkos::WithoutInitializing, sum, i);
      Kokkos::resize(Kokkos::WithoutInitializing, matrix, nGlobalEntries,
                     nGlobalEntries, i);
      Kokkos::deep_copy(matrix, initValue);

      auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
          {0, 0}, {nGlobalEntries, nGlobalEntries});
      auto functor = MDArrayReduceFunctor(matrix);

      if (sumInitShmemSize < expectedInitShmemLimit) {
        EXPECT_NO_THROW(Kokkos::parallel_reduce(policy, functor, sum));
      } else {
        EXPECT_THROW(Kokkos::parallel_reduce(policy, functor, sum),
                     std::runtime_error);
      }
    }
  }
};

}  // namespace Impl

TEST(cuda, reduceRangePolicyViewSizeLimit) {
  Impl::ReduceViewSizeLimitTester reduceViewSizeLimitTester;

  reduceViewSizeLimitTester.run_test_range();
}

TEST(cuda, reduceMDRangePolicyViewSizeLimit) {
  Impl::ReduceViewSizeLimitTester reduceViewSizeLimitTester;

  reduceViewSizeLimitTester.run_test_md_range_2D();
}

}  // namespace Test
