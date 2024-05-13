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

#ifndef TEST_CONCURRENTBITSET_HPP
#define TEST_CONCURRENTBITSET_HPP

#include <gtest/gtest.h>

#include <sstream>
#include <iostream>

#include <impl/Kokkos_ConcurrentBitset.hpp>

namespace Test {

template <class DeviceType>
struct ConcurrentBitset {
  using view_unsigned_type = Kokkos::View<uint32_t*, DeviceType>;
  using view_int_type      = Kokkos::View<int*, DeviceType>;

  view_unsigned_type bitset;
  view_int_type acquired;
  uint32_t bitset_count_lg2;
  uint32_t bitset_count_mask;

  ConcurrentBitset(const uint32_t arg_bitset_count_lg2,
                   const view_unsigned_type& arg_bitset,
                   const view_int_type& arg_acquired)
      : bitset(arg_bitset),
        acquired(arg_acquired),
        bitset_count_lg2(arg_bitset_count_lg2),
        bitset_count_mask(uint32_t(1u << arg_bitset_count_lg2) - 1) {}

  struct TagAcquire {};
  struct TagRelease {};
  struct TagReacquire {};

  KOKKOS_INLINE_FUNCTION
  void operator()(TagAcquire, int i, long& update) const {
    unsigned hint = Kokkos::Impl::clock_tic() & bitset_count_mask;

    Kokkos::pair<int, int> result =
        Kokkos::Impl::concurrent_bitset::acquire_bounded_lg2(
            bitset.data(), bitset_count_lg2, hint);

    acquired(i) = result.first;

    if (0 <= result.first) ++update;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(TagRelease, int i, long& update) const {
    if (0 == (i % 3) && 0 <= acquired(i)) {
      Kokkos::Impl::concurrent_bitset::release(bitset.data(), acquired(i));
      acquired(i) = -1;
      ++update;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(TagReacquire, int i, long& update) const {
    if (acquired(i) < 0) {
      unsigned hint = Kokkos::Impl::clock_tic() & bitset_count_mask;

      Kokkos::pair<int, int> result =
          Kokkos::Impl::concurrent_bitset::acquire_bounded_lg2(
              bitset.data(), bitset_count_lg2, hint);

      acquired(i) = result.first;

      if (0 <= result.first) ++update;
    }
  }
};

template <class DeviceType>
void test_concurrent_bitset(int bit_count) {
  using Functor            = ConcurrentBitset<DeviceType>;
  using view_unsigned_type = typename Functor::view_unsigned_type;
  using view_int_type      = typename Functor::view_int_type;

  int bit_count_lg2 = 1;

  while ((1 << bit_count_lg2) < bit_count) ++bit_count_lg2;

  bit_count = 1 << bit_count_lg2;

  const int buffer_length =
      Kokkos::Impl::concurrent_bitset::buffer_bound_lg2(bit_count_lg2);

  view_unsigned_type bitset("bitset", buffer_length);

  // Try to acquire more than available:

  const size_t n = (bit_count * 3) / 2;

  view_int_type acquired("acquired", n);

  typename view_unsigned_type::HostMirror bitset_host =
      Kokkos::create_mirror_view(bitset);

  Kokkos::deep_copy(bitset, 0u);

  long total           = 0;
  long total_release   = 0;
  long total_reacquire = 0;

  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<DeviceType, typename Functor::TagAcquire>(0, n),
      Functor(bit_count_lg2, bitset, acquired), total);

  ASSERT_EQ(bit_count, total);

  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<DeviceType, typename Functor::TagRelease>(0, n),
      Functor(bit_count_lg2, bitset, acquired), total_release);

  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<DeviceType, typename Functor::TagReacquire>(0, n),
      Functor(bit_count_lg2, bitset, acquired), total_reacquire);

  ASSERT_EQ(total_release, total_reacquire);
}

}  // namespace Test

#endif /* #ifndef TEST_CONCURRENTBITSET_HPP */
