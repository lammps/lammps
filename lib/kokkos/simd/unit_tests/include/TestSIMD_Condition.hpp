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

#ifndef KOKKOS_TEST_SIMD_CONDITION_HPP
#define KOKKOS_TEST_SIMD_CONDITION_HPP

#include <Kokkos_SIMD.hpp>
#include <SIMDTesting_Utilities.hpp>

template <typename Abi, typename DataType>
inline void host_check_condition() {
  if constexpr (is_type_v<Kokkos::Experimental::simd<DataType, Abi>>) {
    using simd_type = typename Kokkos::Experimental::simd<DataType, Abi>;
    using mask_type = typename simd_type::mask_type;

    auto condition_op = [](mask_type const& mask, simd_type const& a,
                           simd_type const& b) {
      return Kokkos::Experimental::condition(mask, a, b);
    };

    simd_type value_a(16);
    simd_type value_b(20);

    auto condition_result = condition_op(mask_type(false), value_a, value_b);
    EXPECT_TRUE(all_of(condition_result == value_b));
    condition_result = condition_op(mask_type(true), value_a, value_b);
    EXPECT_TRUE(all_of(condition_result == value_a));
  }
}

template <typename Abi, typename... DataTypes>
inline void host_check_condition_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (host_check_condition<Abi, DataTypes>(), ...);
}

template <typename... Abis>
inline void host_check_condition_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (host_check_condition_all_types<Abis>(DataTypes()), ...);
}

template <typename Abi, typename DataType>
KOKKOS_INLINE_FUNCTION void device_check_condition() {
  if constexpr (is_type_v<Kokkos::Experimental::simd<DataType, Abi>>) {
    using simd_type = typename Kokkos::Experimental::simd<DataType, Abi>;
    using mask_type = typename simd_type::mask_type;
    kokkos_checker checker;

    auto condition_op = [](mask_type const& mask, simd_type const& a,
                           simd_type const& b) {
      return Kokkos::Experimental::condition(mask, a, b);
    };

    simd_type value_a(16);
    simd_type value_b(20);

    auto condition_result = condition_op(mask_type(false), value_a, value_b);
    checker.truth(all_of(condition_result == value_b));
    condition_result = condition_op(mask_type(true), value_a, value_b);
    checker.truth(all_of(condition_result == value_a));
  }
}

template <typename Abi, typename... DataTypes>
KOKKOS_INLINE_FUNCTION void device_check_condition_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (device_check_condition<Abi, DataTypes>(), ...);
}

template <typename... Abis>
KOKKOS_INLINE_FUNCTION void device_check_condition_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (device_check_condition_all_types<Abis>(DataTypes()), ...);
}

class simd_device_condition_functor {
 public:
  KOKKOS_INLINE_FUNCTION void operator()(int) const {
    device_check_condition_all_abis(
        Kokkos::Experimental::Impl::device_abi_set());
  }
};

TEST(simd, host_condition) {
  host_check_condition_all_abis(Kokkos::Experimental::Impl::host_abi_set());
}

TEST(simd, device_condition) {
  Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::IndexType<int>>(0, 1),
                       simd_device_condition_functor());
}

#endif
