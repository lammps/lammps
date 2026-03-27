// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_SIMD_LOADSTORE_HPP
#define KOKKOS_TEST_SIMD_LOADSTORE_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.simd;
import kokkos.simd_impl;
#else
#include <Kokkos_SIMD.hpp>
#endif
#include <SIMDTesting_Utilities.hpp>

template <typename SimdType, typename... Args>
inline void host_test_simd_load(SimdType const& init, SimdType const& expected,
                                Args... args) {
  using data_type = typename SimdType::value_type;
  using abi_type  = typename SimdType::abi_type;

  constexpr size_t alignment =
      SimdType::size() * sizeof(typename SimdType::value_type);

  alignas(alignment) typename SimdType::value_type arr[SimdType::size()] = {0};
  simd_unchecked_store(init, arr, Kokkos::Experimental::simd_flag_default);

  SimdType result;

  if constexpr (sizeof...(args) > 1) {
    result = simd_partial_load(arr, args...);
  } else {
    if constexpr (std::is_same_v<Kokkos::Experimental::simd_abi::Impl::
                                     host_fixed_native<data_type>,
                                 abi_type>) {
      result = Kokkos::Experimental::simd_unchecked_load(arr, args...);
    } else {
      result =
          Kokkos::Experimental::simd_unchecked_load<SimdType>(arr, args...);
    }
  }

  auto mask = (result == expected);
  for (size_t i = 0; i < SimdType::size(); ++i) {
    EXPECT_TRUE(mask[i]);
  }
}

template <typename SimdType, typename... Args>
inline void host_test_simd_store(SimdType const& init, SimdType const& expected,
                                 Args... args) {
  constexpr size_t alignment =
      SimdType::size() * sizeof(typename SimdType::value_type);

  alignas(alignment) typename SimdType::value_type arr[SimdType::size()] = {0};

  if constexpr (sizeof...(args) > 1) {
    simd_partial_store(init, arr, args...);
  } else {
    simd_unchecked_store(init, arr, args...);
  }

  for (size_t i = 0; i < SimdType::size(); ++i) {
    EXPECT_EQ(arr[i], expected[i]);
  }
}

template <typename Abi, typename DataType>
inline void host_test_simd_loadstore() {
  using simd_type = Kokkos::Experimental::basic_simd<DataType, Abi>;
  using mask_type = typename simd_type::mask_type;

  mask_type mask(KOKKOS_LAMBDA(std::size_t i) { return i % 2 == 0; });
  simd_type expected(KOKKOS_LAMBDA(std::size_t i) { return (i + 1) * 12; });
  simd_type expected_masked(KOKKOS_LAMBDA(std::size_t i) {
    return (i % 2 == 0) ? (i + 1) * 12 : DataType();
  });

  host_test_simd_store(expected, expected,
                       Kokkos::Experimental::simd_flag_default);
  host_test_simd_store(expected, expected,
                       Kokkos::Experimental::simd_flag_aligned);
  host_test_simd_store(expected, expected_masked, mask,
                       Kokkos::Experimental::simd_flag_default);
  host_test_simd_store(expected, expected_masked, mask,
                       Kokkos::Experimental::simd_flag_aligned);

  host_test_simd_load(expected, expected,
                      Kokkos::Experimental::simd_flag_default);
  host_test_simd_load(expected, expected,
                      Kokkos::Experimental::simd_flag_aligned);
  host_test_simd_load(expected, expected_masked, mask,
                      Kokkos::Experimental::simd_flag_default);
  host_test_simd_load(expected, expected_masked, mask,
                      Kokkos::Experimental::simd_flag_aligned);
}

template <typename Abi, typename DataType>
inline void host_check_loadstore() {
  if constexpr (is_simd_avail_v<DataType, Abi>) {
    host_test_simd_loadstore<Abi, DataType>();
  }
}

template <typename Abi, typename... DataTypes>
inline void host_check_loadstore_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (host_check_loadstore<Abi, DataTypes>(), ...);
}

template <typename... Abis>
inline void host_check_loadstore_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (host_check_loadstore_all_types<Abis>(DataTypes()), ...);
}

template <typename SimdType, typename... Args>
KOKKOS_INLINE_FUNCTION void device_test_simd_load(SimdType const& init,
                                                  SimdType const& expected,
                                                  Args... args) {
  using data_type = typename SimdType::value_type;
  using abi_type  = typename SimdType::abi_type;

  constexpr size_t alignment =
      SimdType::size() * sizeof(typename SimdType::value_type);

  alignas(alignment) typename SimdType::value_type arr[SimdType::size()] = {0};
  simd_unchecked_store(init, arr, Kokkos::Experimental::simd_flag_default);

  SimdType result;

  if constexpr (sizeof...(args) > 1) {
    result = simd_partial_load(arr, args...);
  } else {
    if constexpr (std::is_same_v<Kokkos::Experimental::simd_abi::Impl::
                                     host_fixed_native<data_type>,
                                 abi_type>) {
      result = Kokkos::Experimental::simd_unchecked_load(arr, args...);
    } else {
      result =
          Kokkos::Experimental::simd_unchecked_load<SimdType>(arr, args...);
    }
  }

  device_check_equality(result, expected, SimdType::size());
}

template <typename SimdType, typename... Args>
KOKKOS_INLINE_FUNCTION void device_test_simd_store(SimdType const& init,
                                                   SimdType const& expected,
                                                   Args... args) {
  constexpr size_t alignment =
      SimdType::size() * sizeof(typename SimdType::value_type);

  alignas(alignment) typename SimdType::value_type arr[SimdType::size()] = {0};

  if constexpr (sizeof...(args) > 1) {
    simd_partial_store(init, arr, args...);
  } else {
    simd_unchecked_store(init, arr, args...);
  }

  kokkos_checker checker;
  for (size_t i = 0; i < SimdType::size(); ++i) {
    checker.equality(arr[i], expected[i]);
  }
}

template <typename Abi, typename DataType>
KOKKOS_INLINE_FUNCTION void device_test_simd_loadstore() {
  using simd_type = Kokkos::Experimental::basic_simd<DataType, Abi>;
  using mask_type = typename simd_type::mask_type;

  mask_type mask([=](std::size_t i) { return i % 2 == 0; });
  simd_type expected([=](std::size_t i) { return (i + 1) * 12; });
  simd_type expected_masked(
      [=](std::size_t i) { return (mask[i]) ? (i + 1) * 12 : DataType(); });

  device_test_simd_store(expected, expected,
                         Kokkos::Experimental::simd_flag_default);
  device_test_simd_store(expected, expected,
                         Kokkos::Experimental::simd_flag_aligned);
  device_test_simd_store(expected, expected_masked, mask,
                         Kokkos::Experimental::simd_flag_default);
  device_test_simd_store(expected, expected_masked, mask,
                         Kokkos::Experimental::simd_flag_aligned);

  device_test_simd_load(expected, expected,
                        Kokkos::Experimental::simd_flag_default);
  device_test_simd_load(expected, expected,
                        Kokkos::Experimental::simd_flag_aligned);
  device_test_simd_load(expected, expected_masked, mask,
                        Kokkos::Experimental::simd_flag_default);
  device_test_simd_load(expected, expected_masked, mask,
                        Kokkos::Experimental::simd_flag_aligned);
}

template <typename Abi, typename DataType>
KOKKOS_INLINE_FUNCTION void device_check_loadstore() {
  if constexpr (is_simd_avail_v<DataType, Abi>) {
    device_test_simd_loadstore<Abi, DataType>();
  }
}

template <typename Abi, typename... DataTypes>
KOKKOS_INLINE_FUNCTION void device_check_loadstore_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (device_check_loadstore<Abi, DataTypes>(), ...);
}

template <typename... Abis>
KOKKOS_INLINE_FUNCTION void device_check_loadstore_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (device_check_loadstore_all_types<Abis>(DataTypes()), ...);
}

class simd_device_loadstore_functor {
 public:
  KOKKOS_INLINE_FUNCTION void operator()(int) const {
    device_check_loadstore_all_abis(
        Kokkos::Experimental::Impl::device_abi_set());
  }
};

TEST(simd, host_loadstore) {
  host_check_loadstore_all_abis(Kokkos::Experimental::Impl::host_abi_set());
}

TEST(simd, device_loadstore) {
  Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::IndexType<int>>(0, 1),
                       simd_device_loadstore_functor());
}

#endif
