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

#ifndef KOKKOS_TEST_SIMD_CONSTRUCTION_HPP
#define KOKKOS_TEST_SIMD_CONSTRUCTION_HPP

#include <Kokkos_SIMD.hpp>
#include <SIMDTesting_Utilities.hpp>

template <typename Abi, typename DataType>
inline void host_test_simd_traits() {
  using simd_type = Kokkos::Experimental::simd<DataType, Abi>;

  static_assert(std::is_nothrow_default_constructible_v<simd_type>);
  static_assert(std::is_nothrow_copy_assignable_v<simd_type>);
  static_assert(std::is_nothrow_copy_constructible_v<simd_type>);
  static_assert(std::is_nothrow_move_assignable_v<simd_type>);
  static_assert(std::is_nothrow_move_constructible_v<simd_type>);

  simd_type default_simd, result;
  simd_type test_simd(KOKKOS_LAMBDA(std::size_t i) { return (i % 2 == 0); });
  simd_type copy_simd(test_simd);
  simd_type move_simd(std::move(copy_simd));
  default_simd = std::move(move_simd);
  result       = default_simd;
  EXPECT_TRUE(all_of(test_simd == result));
}

template <typename Abi, typename DataType>
inline void host_test_mask_traits() {
  using mask_type = Kokkos::Experimental::simd_mask<DataType, Abi>;

  static_assert(std::is_nothrow_default_constructible_v<mask_type>);
  static_assert(std::is_nothrow_copy_assignable_v<mask_type>);
  static_assert(std::is_nothrow_copy_constructible_v<mask_type>);
  static_assert(std::is_nothrow_move_assignable_v<mask_type>);
  static_assert(std::is_nothrow_move_constructible_v<mask_type>);

  mask_type default_mask, result;
  mask_type test_mask(KOKKOS_LAMBDA(std::size_t i) { return (i % 2 == 0); });
  mask_type copy_mask(test_mask);
  mask_type move_mask(std::move(copy_mask));
  default_mask = std::move(move_mask);
  result       = default_mask;
  EXPECT_EQ(test_mask, result);
}

template <typename Abi, typename DataType>
inline void host_check_construction() {
  if constexpr (is_type_v<Kokkos::Experimental::simd<DataType, Abi>>) {
    host_test_simd_traits<Abi, DataType>();
    host_test_mask_traits<Abi, DataType>();
  }
}

template <typename Abi, typename... DataTypes>
inline void host_check_construction_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (host_check_construction<Abi, DataTypes>(), ...);
}

template <typename... Abis>
inline void host_check_construction_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (host_check_construction_all_types<Abis>(DataTypes()), ...);
}

template <typename Abi, typename DataType>
KOKKOS_INLINE_FUNCTION void device_test_simd_traits() {
  using simd_type = Kokkos::Experimental::simd<DataType, Abi>;

  simd_type default_simd, result;
  simd_type test_simd(KOKKOS_LAMBDA(std::size_t i) { return (i % 2 == 0); });
  simd_type copy_simd(test_simd);
  simd_type move_simd(std::move(copy_simd));
  default_simd = std::move(move_simd);
  result       = default_simd;

  kokkos_checker checker;
  checker.truth(all_of(test_simd == result));
}

template <typename Abi, typename DataType>
KOKKOS_INLINE_FUNCTION void device_test_mask_traits() {
  using mask_type = Kokkos::Experimental::simd_mask<DataType, Abi>;

  mask_type default_mask, result;
  mask_type test_mask(KOKKOS_LAMBDA(std::size_t i) { return (i % 2 == 0); });
  mask_type copy_mask(test_mask);
  mask_type move_mask(std::move(copy_mask));
  default_mask = std::move(move_mask);
  result       = default_mask;

  kokkos_checker checker;
  checker.truth(test_mask == result);
}

template <typename Abi, typename DataType>
KOKKOS_INLINE_FUNCTION void device_check_construction() {
  if constexpr (is_type_v<Kokkos::Experimental::simd<DataType, Abi>>) {
    device_test_simd_traits<Abi, DataType>();
    device_test_mask_traits<Abi, DataType>();
  }
}

template <typename Abi, typename... DataTypes>
KOKKOS_INLINE_FUNCTION void device_check_construction_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (device_check_construction<Abi, DataTypes>(), ...);
}

template <typename... Abis>
KOKKOS_INLINE_FUNCTION void device_check_construction_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (device_check_construction_all_types<Abis>(DataTypes()), ...);
}

class simd_device_construction_functor {
 public:
  KOKKOS_INLINE_FUNCTION void operator()(int) const {
    device_check_construction_all_abis(
        Kokkos::Experimental::Impl::device_abi_set());
  }
};

TEST(simd, host_construction) {
  host_check_construction_all_abis(Kokkos::Experimental::Impl::host_abi_set());
}

TEST(simd, device_construction) {
  Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::IndexType<int>>(0, 1),
                       simd_device_construction_functor());
}

#endif
