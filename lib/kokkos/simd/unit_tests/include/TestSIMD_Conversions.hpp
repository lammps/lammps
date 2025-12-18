// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_SIMD_CONVERSIONS_HPP
#define KOKKOS_TEST_SIMD_CONVERSIONS_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.simd;
import kokkos.simd_impl;
#else
#include <Kokkos_SIMD.hpp>
#endif
#include <SIMDTesting_Utilities.hpp>

using Kokkos::Experimental::all_of;

template <typename Abi, typename DataTypeA, typename DataTypeB>
inline void host_check_conversions() {
  if constexpr (is_simd_avail_v<DataTypeA, Abi> &&
                is_simd_avail_v<DataTypeB, Abi>) {
    DataTypeA test_val =
        (std::is_signed_v<DataTypeA> && std::is_signed_v<DataTypeB>) ? -213
                                                                     : 213;
    bool test_mask_val = true;
    {
      auto from = Kokkos::Experimental::basic_simd<DataTypeA, Abi>(test_val);
      auto to   = Kokkos::Experimental::basic_simd<DataTypeB, Abi>(from);
      auto expected =
          Kokkos::Experimental::basic_simd<DataTypeB, Abi>(test_val);
      EXPECT_EQ(from.size(), to.size());
      host_check_equality(to, decltype(to)(test_val), to.size());
      host_check_equality(to, expected, to.size());
    }
    {
      auto from =
          Kokkos::Experimental::basic_simd_mask<DataTypeA, Abi>(test_mask_val);
      auto to = Kokkos::Experimental::basic_simd_mask<DataTypeB, Abi>(from);
      auto expected =
          Kokkos::Experimental::basic_simd_mask<DataTypeB, Abi>(test_mask_val);
      EXPECT_EQ(from.size(), to.size());
      EXPECT_TRUE(all_of(to == decltype(to)(test_mask_val)));
      EXPECT_TRUE(all_of(to == expected));
    }
  }
}

template <typename Abi, typename DataTypeA, typename... DataTypesB>
inline void host_check_conversions_all_types_to(
    Kokkos::Experimental::Impl::data_types<DataTypesB...>) {
  (host_check_conversions<Abi, DataTypeA, DataTypesB>(), ...);
}

template <typename Abi, typename... DataTypesA>
inline void host_check_conversions_all_types_from(
    Kokkos::Experimental::Impl::data_types<DataTypesA...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (host_check_conversions_all_types_to<Abi, DataTypesA>(DataTypes()), ...);
}

template <typename... Abis>
inline void host_check_conversions_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (host_check_conversions_all_types_from<Abis>(DataTypes()), ...);
}

template <typename Abi, typename DataTypeA, typename DataTypeB>
KOKKOS_INLINE_FUNCTION void device_check_conversions() {
  if constexpr (is_simd_avail_v<DataTypeA, Abi> &&
                is_simd_avail_v<DataTypeB, Abi>) {
    DataTypeA test_val =
        (std::is_signed_v<DataTypeA> && std::is_signed_v<DataTypeB>) ? -213
                                                                     : 213;
    bool test_mask_val = true;
    kokkos_checker checker;
    {
      auto from = Kokkos::Experimental::basic_simd<DataTypeA, Abi>(test_val);
      auto to   = Kokkos::Experimental::basic_simd<DataTypeB, Abi>(from);
      auto expected =
          Kokkos::Experimental::basic_simd<DataTypeB, Abi>(test_val);
      checker.truth(from.size() == to.size());
      device_check_equality(to, decltype(to)(test_val), to.size());
      device_check_equality(to, expected, to.size());
    }
    {
      auto from =
          Kokkos::Experimental::basic_simd_mask<DataTypeA, Abi>(test_mask_val);
      auto to = Kokkos::Experimental::basic_simd_mask<DataTypeB, Abi>(from);
      auto expected =
          Kokkos::Experimental::basic_simd_mask<DataTypeB, Abi>(test_mask_val);
      checker.truth(from.size() == to.size());
      checker.truth(all_of(to == decltype(to)(test_mask_val)));
      checker.truth(all_of(to == expected));
    }
  }
}

template <typename Abi, typename DataTypeA, typename... DataTypesB>
KOKKOS_INLINE_FUNCTION void device_check_conversions_all_types_to(
    Kokkos::Experimental::Impl::data_types<DataTypesB...>) {
  (device_check_conversions<Abi, DataTypeA, DataTypesB>(), ...);
}

template <typename Abi, typename... DataTypesA>
KOKKOS_INLINE_FUNCTION void device_check_conversions_all_types_from(
    Kokkos::Experimental::Impl::data_types<DataTypesA...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (device_check_conversions_all_types_to<Abi, DataTypesA>(DataTypes()), ...);
}

template <typename... Abis>
KOKKOS_INLINE_FUNCTION void device_check_conversions_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (device_check_conversions_all_types_from<Abis>(DataTypes()), ...);
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
