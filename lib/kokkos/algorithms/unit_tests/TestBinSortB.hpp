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

#ifndef KOKKOS_ALGORITHMS_UNITTESTS_TEST_BINSORTB_HPP
#define KOKKOS_ALGORITHMS_UNITTESTS_TEST_BINSORTB_HPP

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_Sort.hpp>
#include <Kokkos_StdAlgorithms.hpp>
#include <TestStdAlgorithmsCommon.hpp>
#include <random>
#include <numeric>  //needed for iota

namespace Test {
namespace BinSortSetB {

template <class KeyType, class ExecutionSpace>
auto create_rank1_dev_and_host_views_of_keys(const ExecutionSpace& exec,
                                             int N) {
  namespace KE = Kokkos::Experimental;
  Kokkos::DefaultHostExecutionSpace defaultHostExeSpace;

  using KeyViewType = Kokkos::View<KeyType*, ExecutionSpace>;
  KeyViewType keys("keys", N);
  auto keys_h = Kokkos::create_mirror_view(keys);
  std::iota(KE::begin(keys_h), KE::end(keys_h), KeyType(0));
  KE::reverse(defaultHostExeSpace, keys_h);
  // keys now is = [N-1,N-2,...,2,1,0], shuffle it for avoid trivial case
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(KE::begin(keys_h), KE::end(keys_h), g);
  Kokkos::deep_copy(exec, keys, keys_h);

  return std::make_pair(keys, keys_h);
}

template <class ExecutionSpace, class ValueType, int ValuesViewRank,
          std::enable_if_t<ValuesViewRank == 1, int> = 0>
auto create_strided_view(std::size_t numRows, std::size_t /*numCols*/) {
  Kokkos::LayoutStride layout{numRows, 2};
  using v_t = Kokkos::View<ValueType*, Kokkos::LayoutStride, ExecutionSpace>;
  v_t v("v", layout);
  return v;
}

template <class ExecutionSpace, class ValueType, int ValuesViewRank,
          std::enable_if_t<ValuesViewRank == 2, int> = 0>
auto create_strided_view(std::size_t numRows, std::size_t numCols) {
  Kokkos::LayoutStride layout{numRows, 2, numCols, numRows * 2};
  using v_t = Kokkos::View<ValueType**, Kokkos::LayoutStride, ExecutionSpace>;
  v_t v("v", layout);
  return v;
}

template <class ExecutionSpace, class KeyType, class ValueType,
          int ValuesViewRank>
void test_on_view_with_stride(std::size_t numRows, std::size_t indB,
                              std::size_t indE, std::size_t numCols = 1) {
  ExecutionSpace exec;
  Kokkos::DefaultHostExecutionSpace defaultHostExeSpace;
  namespace KE = Kokkos::Experimental;

  // 1. generate 1D view of keys
  auto [keys, keys_h] =
      create_rank1_dev_and_host_views_of_keys<KeyType>(exec, numRows);
  using KeyViewType = decltype(keys);

  // need this map key->row to use later for checking
  std::unordered_map<KeyType, std::size_t> keyToRowBeforeSort;
  for (std::size_t i = 0; i < numRows; ++i) {
    keyToRowBeforeSort[keys_h(i)] = i;
  }

  // 2. create binOp
  using BinOp = Kokkos::BinOp1D<KeyViewType>;
  auto itB    = KE::cbegin(keys_h) + indB;
  auto itE    = itB + indE - indB;
  auto it     = KE::minmax_element(defaultHostExeSpace, itB, itE);
  // seems like the behavior is odd when we use # buckets = # keys
  // so use +5 for using more buckets than keys.
  // This is something to investigate.
  BinOp binner(indE - indB + 5, *it.first, *it.second);

  // 3. create sorter
  Kokkos::BinSort<KeyViewType, BinOp> sorter(keys, indB, indE, binner, false);
  sorter.create_permute_vector(exec);
  sorter.sort(exec, keys, indB, indE);
  Kokkos::deep_copy(exec, keys_h, keys);

  auto v = create_strided_view<ExecutionSpace, ValueType, ValuesViewRank>(
      numRows, numCols);

  Kokkos::Random_XorShift64_Pool<ExecutionSpace> pool(73931);
  Kokkos::fill_random(v, pool, ValueType(545));
  auto v_before_sort_h = stdalgos::create_host_space_copy(v);
  sorter.sort(exec, v, indB, indE);
  auto v_after_sort_h = stdalgos::create_host_space_copy(v);

  for (size_t i = 0; i < v.extent(0); ++i) {
    // if i within [indB,indE), the sorting was done
    // so we need to do proper checking since rows have changed
    if (i >= size_t(indB) && i < size_t(indE)) {
      const KeyType key = keys_h(i);
      if constexpr (ValuesViewRank == 1) {
        ASSERT_TRUE(v_before_sort_h(keyToRowBeforeSort.at(key)) ==
                    v_after_sort_h(i));
      } else {
        for (size_t j = 0; j < v.extent(1); ++j) {
          ASSERT_TRUE(v_before_sort_h(keyToRowBeforeSort.at(key), j) ==
                      v_after_sort_h(i, j));
        }
      }
    }
    // outside the target bounds, then the i-th row remains unchanged
    else {
      if constexpr (ValuesViewRank == 1) {
        ASSERT_TRUE(v_before_sort_h(i) == v_after_sort_h(i));
      } else {
        for (size_t j = 0; j < v.extent(1); ++j) {
          ASSERT_TRUE(v_before_sort_h(i, j) == v_after_sort_h(i, j));
        }
      }
    }
  }
}

template <class ExecutionSpace, class KeyType, class ValueType>
void run_for_rank1() {
  constexpr int rank = 1;

  // trivial case
  test_on_view_with_stride<ExecutionSpace, KeyType, ValueType, rank>(1, 0, 1);

  // nontrivial cases
  for (std::size_t N : {311, 710017}) {
    // various cases for bounds
    test_on_view_with_stride<ExecutionSpace, KeyType, ValueType, rank>(N, 0, N);
    test_on_view_with_stride<ExecutionSpace, KeyType, ValueType, rank>(N, 3, N);
    test_on_view_with_stride<ExecutionSpace, KeyType, ValueType, rank>(N, 0,
                                                                       N - 4);
    test_on_view_with_stride<ExecutionSpace, KeyType, ValueType, rank>(N, 4,
                                                                       N - 3);
  }
}

template <class ExecutionSpace, class KeyType, class ValueType>
void run_for_rank2() {
  constexpr int rank = 2;

  // trivial case
  test_on_view_with_stride<ExecutionSpace, KeyType, ValueType, rank>(1, 0, 1,
                                                                     1);

  // nontrivial cases
  for (std::size_t Nr : {11, 1157, 710017}) {
    for (std::size_t Nc : {3, 51}) {
      // various cases for bounds
      test_on_view_with_stride<ExecutionSpace, KeyType, ValueType, rank>(
          Nr, 0, Nr, Nc);
      test_on_view_with_stride<ExecutionSpace, KeyType, ValueType, rank>(
          Nr, 3, Nr, Nc);
      test_on_view_with_stride<ExecutionSpace, KeyType, ValueType, rank>(
          Nr, 0, Nr - 4, Nc);
      test_on_view_with_stride<ExecutionSpace, KeyType, ValueType, rank>(
          Nr, 4, Nr - 3, Nc);
    }
  }
}

}  // namespace BinSortSetB

TEST(TEST_CATEGORY, BinSortUnsignedKeyLayoutStrideValues) {
  using ExeSpace = TEST_EXECSPACE;
  using key_type = unsigned;
  BinSortSetB::run_for_rank1<ExeSpace, key_type, int>();
  BinSortSetB::run_for_rank1<ExeSpace, key_type, double>();

  BinSortSetB::run_for_rank2<ExeSpace, key_type, int>();
  BinSortSetB::run_for_rank2<ExeSpace, key_type, double>();
}

}  // namespace Test
#endif
