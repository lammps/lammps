// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_SIMD_MASK_OPS_HPP
#define KOKKOS_TEST_SIMD_MASK_OPS_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.simd;
import kokkos.simd_impl;
#else
#include <Kokkos_SIMD.hpp>
#endif
#include <SIMDTesting_Utilities.hpp>

template <typename Abi, typename DataType>
inline void host_check_mask_ops() {
  if constexpr (is_simd_avail_v<DataType, Abi>) {
    using mask_type = Kokkos::Experimental::basic_simd_mask<DataType, Abi>;

    EXPECT_FALSE(none_of(mask_type(true)));
    EXPECT_TRUE(none_of(mask_type(false)));
    EXPECT_TRUE(all_of(mask_type(true)));
    EXPECT_FALSE(all_of(mask_type(false)));
    EXPECT_TRUE(any_of(mask_type(true)));
    EXPECT_FALSE(any_of(mask_type(false)));

    for (std::size_t i = 0; i < mask_type::size(); ++i) {
      mask_type test_mask(KOKKOS_LAMBDA(std::size_t j) { return i == j; });

      EXPECT_TRUE(any_of(test_mask));
      EXPECT_FALSE(none_of(test_mask));

      if constexpr (mask_type::size() > 1) {
        EXPECT_FALSE(all_of(test_mask));
      } else {
        EXPECT_TRUE(all_of(test_mask));
      }
    }
  }
}

template <typename Abi, typename... DataTypes>
inline void host_check_mask_ops_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (host_check_mask_ops<Abi, DataTypes>(), ...);
}

template <typename... Abis>
inline void host_check_mask_ops_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (host_check_mask_ops_all_types<Abis>(DataTypes()), ...);
}

template <typename Abi, typename DataType>
KOKKOS_INLINE_FUNCTION void device_check_mask_ops() {
  if constexpr (is_type_v<Kokkos::Experimental::basic_simd<DataType, Abi>>) {
    using mask_type = Kokkos::Experimental::basic_simd_mask<DataType, Abi>;
    kokkos_checker checker;
    checker.truth(!none_of(mask_type(true)));
    checker.truth(none_of(mask_type(false)));
    checker.truth(all_of(mask_type(true)));
    checker.truth(!all_of(mask_type(false)));
    checker.truth(any_of(mask_type(true)));
    checker.truth(!any_of(mask_type(false)));

    for (std::size_t i = 0; i < mask_type::size(); ++i) {
      mask_type test_mask(KOKKOS_LAMBDA(std::size_t j) { return i == j; });

      checker.truth(any_of(test_mask));
      checker.truth(!none_of(test_mask));

      if constexpr (mask_type::size() > 1) {
        checker.truth(!all_of(test_mask));
      } else {
        checker.truth(all_of(test_mask));
      }
    }
  }
}

template <typename Abi, typename... DataTypes>
KOKKOS_INLINE_FUNCTION void device_check_mask_ops_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (device_check_mask_ops<Abi, DataTypes>(), ...);
}

template <typename... Abis>
KOKKOS_INLINE_FUNCTION void device_check_mask_ops_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (device_check_mask_ops_all_types<Abis>(DataTypes()), ...);
}

class simd_device_mask_ops_functor {
 public:
  KOKKOS_INLINE_FUNCTION void operator()(int) const {
    device_check_mask_ops_all_abis(
        Kokkos::Experimental::Impl::device_abi_set());
  }
};

TEST(simd, host_mask_ops) {
  host_check_mask_ops_all_abis(Kokkos::Experimental::Impl::host_abi_set());
}

TEST(simd, device_mask_ops) {
  Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::IndexType<int>>(0, 1),
                       simd_device_mask_ops_functor());
}

#endif
