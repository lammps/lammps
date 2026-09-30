// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include <cstdio>

// This test checks parallel_scan() calls which use RangePolicy and custom
// scalar variable.

namespace {

template <typename T, int N>
struct ArrayValueType {
  T v[N];

  ArrayValueType()                        = default;
  ArrayValueType(const ArrayValueType& b) = default;
  ArrayValueType(ArrayValueType&& b)      = default;
  KOKKOS_INLINE_FUNCTION
  ArrayValueType(const T b) {
    for (int i = 0; i < N; ++i) this->v[i] = b;
  }
  ArrayValueType& operator=(const ArrayValueType& b) = default;
  ArrayValueType& operator=(ArrayValueType&& b)      = default;
  KOKKOS_INLINE_FUNCTION
  ArrayValueType& operator=(const T& b) {
    for (int i = 0; i < N; ++i) this->v[i] = b;
    return *this;
  }
  ~ArrayValueType() = default;
};

template <typename T, int N>
static KOKKOS_INLINE_FUNCTION void operator+=(ArrayValueType<T, N>& a,
                                              const ArrayValueType<T, N>& b) {
  for (int i = 0; i < N; ++i) a.v[i] += b.v[i];
}

template <class ExecSpace>
void test_customscalar_parallel_scan() {
  using update_type = ArrayValueType<size_t, 2>;

  const int nrows = 4409601;
  Kokkos::View<size_t*, ExecSpace> rowptr("rowptr", nrows + 1);

  auto rowptr_h = Kokkos::create_mirror_view(rowptr);

  Kokkos::deep_copy(rowptr_h, 6);
  rowptr_h(nrows) = 0;
  Kokkos::deep_copy(rowptr, rowptr_h);

  // Calculate correct exclusive prefix sum
  Kokkos::View<size_t*, Kokkos::HostSpace> prefix_correct("prefix_correct",
                                                          nrows + 1);
  size_t accum = 0;
  for (int i = 0; i < nrows + 1; i++) {
    prefix_correct(i) = accum;
    // note: last element of rowptr_h is 0 so the final value of accum
    // is not used
    accum += rowptr_h(i);
  }

  Kokkos::RangePolicy<ExecSpace> policy(0, nrows + 1);
  Kokkos::parallel_scan(
      policy,
      KOKKOS_LAMBDA(const int& lr, update_type& update, const bool& final) {
        update_type val;
        val = rowptr(lr);

        if (final) {
          rowptr(lr) = update.v[0];
        }
        update += val;
      });
  Kokkos::fence();

  auto rowptr_final =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rowptr);
  for (int i = 0; i < nrows + 1; i++) {
    ASSERT_EQ(rowptr_final(i), prefix_correct(i));
  }
}

TEST(TEST_CATEGORY, customscalar_parallel_scan) {
  test_customscalar_parallel_scan<TEST_EXECSPACE>();
}

}  // namespace
