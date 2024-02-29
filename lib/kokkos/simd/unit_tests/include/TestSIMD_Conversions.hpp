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

#ifndef KOKKOS_TEST_SIMD_CONVERSIONS_HPP
#define KOKKOS_TEST_SIMD_CONVERSIONS_HPP

#include <Kokkos_SIMD.hpp>
#include <SIMDTesting_Utilities.hpp>

template <typename Abi>
inline void host_check_conversions() {
  {
    auto a = Kokkos::Experimental::simd<std::uint64_t, Abi>(1);
    auto b = Kokkos::Experimental::simd<std::int64_t, Abi>(a);
    EXPECT_TRUE(all_of(b == decltype(b)(1)));
  }
  {
    auto a = Kokkos::Experimental::simd<std::int32_t, Abi>(1);
    auto b = Kokkos::Experimental::simd<std::uint64_t, Abi>(a);
    EXPECT_TRUE(all_of(b == decltype(b)(1)));
  }
  {
    auto a = Kokkos::Experimental::simd<std::uint64_t, Abi>(1);
    auto b = Kokkos::Experimental::simd<std::int32_t, Abi>(a);
    EXPECT_TRUE(all_of(b == decltype(b)(1)));
  }
  {
    auto a = Kokkos::Experimental::simd_mask<double, Abi>(true);
    auto b = Kokkos::Experimental::simd_mask<std::int32_t, Abi>(a);
    EXPECT_TRUE(b == decltype(b)(true));
  }
  {
    auto a = Kokkos::Experimental::simd_mask<std::int32_t, Abi>(true);
    auto b = Kokkos::Experimental::simd_mask<std::uint64_t, Abi>(a);
    EXPECT_TRUE(b == decltype(b)(true));
  }
  {
    auto a = Kokkos::Experimental::simd_mask<std::int32_t, Abi>(true);
    auto b = Kokkos::Experimental::simd_mask<std::int64_t, Abi>(a);
    EXPECT_TRUE(b == decltype(b)(true));
  }
  {
    auto a = Kokkos::Experimental::simd_mask<std::int32_t, Abi>(true);
    auto b = Kokkos::Experimental::simd_mask<double, Abi>(a);
    EXPECT_TRUE(b == decltype(b)(true));
  }
}

template <typename... Abis>
inline void host_check_conversions_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  (host_check_conversions<Abis>(), ...);
}

template <typename Abi>
KOKKOS_INLINE_FUNCTION void device_check_conversions() {
  kokkos_checker checker;
  {
    auto a = Kokkos::Experimental::simd<std::uint64_t, Abi>(1);
    auto b = Kokkos::Experimental::simd<std::int64_t, Abi>(a);
    checker.truth(all_of(b == decltype(b)(1)));
  }
  {
    auto a = Kokkos::Experimental::simd<std::int32_t, Abi>(1);
    auto b = Kokkos::Experimental::simd<std::uint64_t, Abi>(a);
    checker.truth(all_of(b == decltype(b)(1)));
  }
  {
    auto a = Kokkos::Experimental::simd<std::uint64_t, Abi>(1);
    auto b = Kokkos::Experimental::simd<std::int32_t, Abi>(a);
    checker.truth(all_of(b == decltype(b)(1)));
  }
  {
    auto a = Kokkos::Experimental::simd_mask<double, Abi>(true);
    auto b = Kokkos::Experimental::simd_mask<std::int32_t, Abi>(a);
    checker.truth(b == decltype(b)(true));
  }
  {
    auto a = Kokkos::Experimental::simd_mask<std::int32_t, Abi>(true);
    auto b = Kokkos::Experimental::simd_mask<std::uint64_t, Abi>(a);
    checker.truth(b == decltype(b)(true));
  }
  {
    auto a = Kokkos::Experimental::simd_mask<std::int32_t, Abi>(true);
    auto b = Kokkos::Experimental::simd_mask<std::int64_t, Abi>(a);
    checker.truth(b == decltype(b)(true));
  }
  {
    auto a = Kokkos::Experimental::simd_mask<std::int32_t, Abi>(true);
    auto b = Kokkos::Experimental::simd_mask<double, Abi>(a);
    checker.truth(b == decltype(b)(true));
  }
}

template <typename... Abis>
KOKKOS_INLINE_FUNCTION void device_check_conversions_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  (device_check_conversions<Abis>(), ...);
}

class simd_device_conversions_functor {
 public:
  KOKKOS_INLINE_FUNCTION void operator()(int) const {
    device_check_conversions_all_abis(
        Kokkos::Experimental::Impl::device_abi_set());
  }
};

TEST(simd, host_conversions) {
  host_check_conversions_all_abis(Kokkos::Experimental::Impl::host_abi_set());
}

TEST(simd, device_conversions) {
  Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::IndexType<int>>(0, 1),
                       simd_device_conversions_functor());
}

#endif
