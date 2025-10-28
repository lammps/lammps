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

using Kokkos::Experimental::all_of;

template <typename DataType>
struct Gen_increment {
  KOKKOS_FUNCTION DataType operator()(std::size_t i) { return (DataType)(i); }
};

struct Gen_bool {
  KOKKOS_FUNCTION bool operator()(std::size_t i) { return 0 == (i & 1); }
};

template <typename Abi>
inline void host_check_conversions() {
  if constexpr (is_simd_avail_v<uint64_t, Abi>) {
    {
      auto a = Kokkos::Experimental::basic_simd<std::uint64_t, Abi>(1);
      auto b = Kokkos::Experimental::basic_simd<std::int64_t, Abi>(a);
      EXPECT_TRUE(all_of(b == decltype(b)(1)));
    }
    {
      auto a = Kokkos::Experimental::basic_simd<std::int32_t, Abi>(1);
      auto b = Kokkos::Experimental::basic_simd<std::uint64_t, Abi>(a);
      EXPECT_TRUE(all_of(b == decltype(b)(1)));
    }
    {
      auto a = Kokkos::Experimental::basic_simd<std::uint64_t, Abi>(1);
      auto b = Kokkos::Experimental::basic_simd<std::int32_t, Abi>(a);
      EXPECT_TRUE(all_of(b == decltype(b)(1)));
    }
    {
      auto a = Kokkos::Experimental::basic_simd<std::uint64_t, Abi>(
          Gen_increment<std::uint64_t>());
      auto b = Kokkos::Experimental::basic_simd<std::int64_t, Abi>(a);
      EXPECT_TRUE(all_of(b == decltype(b)(Gen_increment<std::int64_t>())));
    }
    {
      auto a = Kokkos::Experimental::basic_simd<std::int32_t, Abi>(
          Gen_increment<std::int32_t>());
      auto b = Kokkos::Experimental::basic_simd<std::uint64_t, Abi>(a);
      EXPECT_TRUE(all_of(b == decltype(b)(Gen_increment<std::uint64_t>())));
    }
    {
      auto a = Kokkos::Experimental::basic_simd<std::uint64_t, Abi>(
          Gen_increment<std::uint64_t>());
      auto b = Kokkos::Experimental::basic_simd<std::int32_t, Abi>(a);
      EXPECT_TRUE(all_of(b == decltype(b)(Gen_increment<std::int32_t>())));
    }
  }

  if constexpr (is_type_v<
                    Kokkos::Experimental::basic_simd_mask<int32_t, Abi>> &&
                is_type_v<
                    Kokkos::Experimental::basic_simd_mask<int64_t, Abi>>) {
    {
      auto a = Kokkos::Experimental::basic_simd_mask<double, Abi>(true);
      auto b = Kokkos::Experimental::basic_simd_mask<std::int32_t, Abi>(a);
      EXPECT_TRUE(all_of(b == decltype(b)(true)));
    }
    {
      auto a = Kokkos::Experimental::basic_simd_mask<std::int32_t, Abi>(true);
      auto b = Kokkos::Experimental::basic_simd_mask<std::uint64_t, Abi>(a);
      EXPECT_TRUE(all_of(b == decltype(b)(true)));
    }
    {
      auto a = Kokkos::Experimental::basic_simd_mask<std::int32_t, Abi>(true);
      auto b = Kokkos::Experimental::basic_simd_mask<std::int64_t, Abi>(a);
      EXPECT_TRUE(all_of(b == decltype(b)(true)));
    }
    {
      auto a = Kokkos::Experimental::basic_simd_mask<std::int32_t, Abi>(true);
      auto b = Kokkos::Experimental::basic_simd_mask<double, Abi>(a);
      EXPECT_TRUE(all_of(b == decltype(b)(true)));
    }
    {
      auto a = Kokkos::Experimental::basic_simd_mask<double, Abi>(Gen_bool());
      auto b = Kokkos::Experimental::basic_simd_mask<std::int32_t, Abi>(a);
      EXPECT_TRUE(all_of(b == decltype(b)(Gen_bool())));
    }
    {
      auto a =
          Kokkos::Experimental::basic_simd_mask<std::int32_t, Abi>(Gen_bool());
      auto b = Kokkos::Experimental::basic_simd_mask<std::uint64_t, Abi>(a);
      EXPECT_TRUE(all_of(b == decltype(b)(Gen_bool())));
    }
    {
      auto a =
          Kokkos::Experimental::basic_simd_mask<std::int32_t, Abi>(Gen_bool());
      auto b = Kokkos::Experimental::basic_simd_mask<std::int64_t, Abi>(a);
      EXPECT_TRUE(all_of(b == decltype(b)(Gen_bool())));
    }
    {
      auto a =
          Kokkos::Experimental::basic_simd_mask<std::int32_t, Abi>(Gen_bool());
      auto b = Kokkos::Experimental::basic_simd_mask<double, Abi>(a);
      EXPECT_TRUE(all_of(b == decltype(b)(Gen_bool())));
    }
  }
}

template <typename... Abis>
inline void host_check_conversions_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  (host_check_conversions<Abis>(), ...);
}

template <typename Abi>
KOKKOS_INLINE_FUNCTION void device_check_conversions() {
  if constexpr (is_type_v<Kokkos::Experimental::basic_simd<uint64_t, Abi>>) {
    kokkos_checker checker;
    {
      auto a = Kokkos::Experimental::basic_simd<std::uint64_t, Abi>(1);
      auto b = Kokkos::Experimental::basic_simd<std::int64_t, Abi>(a);
      checker.truth(all_of(b == decltype(b)(1)));
    }
    {
      auto a = Kokkos::Experimental::basic_simd<std::int32_t, Abi>(1);
      auto b = Kokkos::Experimental::basic_simd<std::uint64_t, Abi>(a);
      checker.truth(all_of(b == decltype(b)(1)));
    }
    {
      auto a = Kokkos::Experimental::basic_simd<std::uint64_t, Abi>(1);
      auto b = Kokkos::Experimental::basic_simd<std::int32_t, Abi>(a);
      checker.truth(all_of(b == decltype(b)(1)));
    }
    {
      auto a = Kokkos::Experimental::basic_simd<std::uint64_t, Abi>(
          Gen_increment<std::uint64_t>());
      auto b = Kokkos::Experimental::basic_simd<std::int64_t, Abi>(a);
      checker.truth(all_of(b == decltype(b)(Gen_increment<std::int64_t>())));
    }
    {
      auto a = Kokkos::Experimental::basic_simd<std::int32_t, Abi>(
          Gen_increment<std::int32_t>());
      auto b = Kokkos::Experimental::basic_simd<std::uint64_t, Abi>(a);
      checker.truth(all_of(b == decltype(b)(Gen_increment<std::uint64_t>())));
    }
    {
      auto a = Kokkos::Experimental::basic_simd<std::uint64_t, Abi>(
          Gen_increment<std::uint64_t>());
      auto b = Kokkos::Experimental::basic_simd<std::int32_t, Abi>(a);
      checker.truth(all_of(b == decltype(b)(Gen_increment<std::int32_t>())));
    }
    {
      auto a = Kokkos::Experimental::basic_simd_mask<double, Abi>(true);
      auto b = Kokkos::Experimental::basic_simd_mask<std::int32_t, Abi>(a);
      checker.truth(all_of(b == decltype(b)(true)));
    }
    {
      auto a = Kokkos::Experimental::basic_simd_mask<std::int32_t, Abi>(true);
      auto b = Kokkos::Experimental::basic_simd_mask<std::uint64_t, Abi>(a);
      checker.truth(all_of(b == decltype(b)(true)));
    }
    {
      auto a = Kokkos::Experimental::basic_simd_mask<std::int32_t, Abi>(true);
      auto b = Kokkos::Experimental::basic_simd_mask<std::int64_t, Abi>(a);
      checker.truth(all_of(b == decltype(b)(true)));
    }
    {
      auto a = Kokkos::Experimental::basic_simd_mask<std::int32_t, Abi>(true);
      auto b = Kokkos::Experimental::basic_simd_mask<double, Abi>(a);
      checker.truth(all_of(b == decltype(b)(true)));
    }
    {
      auto a = Kokkos::Experimental::basic_simd_mask<double, Abi>(Gen_bool());
      auto b = Kokkos::Experimental::basic_simd_mask<std::int32_t, Abi>(a);
      checker.truth(all_of(b == decltype(b)(Gen_bool())));
    }
    {
      auto a =
          Kokkos::Experimental::basic_simd_mask<std::int32_t, Abi>(Gen_bool());
      auto b = Kokkos::Experimental::basic_simd_mask<std::uint64_t, Abi>(a);
      checker.truth(all_of(b == decltype(b)(Gen_bool())));
    }
    {
      auto a =
          Kokkos::Experimental::basic_simd_mask<std::int32_t, Abi>(Gen_bool());
      auto b = Kokkos::Experimental::basic_simd_mask<std::int64_t, Abi>(a);
      checker.truth(all_of(b == decltype(b)(Gen_bool())));
    }
    {
      auto a =
          Kokkos::Experimental::basic_simd_mask<std::int32_t, Abi>(Gen_bool());
      auto b = Kokkos::Experimental::basic_simd_mask<double, Abi>(a);
      checker.truth(all_of(b == decltype(b)(Gen_bool())));
    }
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
