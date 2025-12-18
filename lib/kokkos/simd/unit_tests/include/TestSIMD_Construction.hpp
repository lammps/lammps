// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_SIMD_CONSTRUCTION_HPP
#define KOKKOS_TEST_SIMD_CONSTRUCTION_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.simd;
import kokkos.simd_impl;
#else
#include <Kokkos_SIMD.hpp>
#endif
#include <SIMDTesting_Utilities.hpp>

#include <climits>

using Kokkos::Experimental::all_of;

template <typename Abi, typename DataType>
inline void host_test_simd_traits() {
  using simd_type = Kokkos::Experimental::basic_simd<DataType, Abi>;

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
  using mask_type = Kokkos::Experimental::basic_simd_mask<DataType, Abi>;

  static_assert(std::is_nothrow_default_constructible_v<mask_type>);
  static_assert(std::is_nothrow_copy_assignable_v<mask_type>);
  static_assert(std::is_nothrow_copy_constructible_v<mask_type>);
  static_assert(std::is_nothrow_move_assignable_v<mask_type>);
  static_assert(std::is_nothrow_move_constructible_v<mask_type>);

  mask_type default_mask(false);
  mask_type result(false);
  mask_type test_mask(KOKKOS_LAMBDA(std::size_t i) { return (i % 2 == 0); });
  mask_type copy_mask(test_mask);
  mask_type move_mask(std::move(copy_mask));
  default_mask = std::move(move_mask);
  result       = default_mask;
  EXPECT_TRUE(all_of(test_mask == result));
}

template <typename Abi, typename DataType>
inline void host_test_simd_alias() {
  using basic_simd_type = Kokkos::Experimental::basic_simd<DataType, Abi>;
  using native_fixed_abi =
      Kokkos::Experimental::simd_abi::Impl::native_fixed_abi<DataType>;
  using native_abi =
      Kokkos::Experimental::simd_abi::Impl::native_abi<DataType,
                                                       basic_simd_type::size()>;

  if constexpr (std::is_same_v<Abi, native_fixed_abi>) {
    using simd_type =
        Kokkos::Experimental::simd<DataType, basic_simd_type::size()>;
    using simd_mask_type =
        Kokkos::Experimental::simd_mask<DataType, basic_simd_type::size()>;
    static_assert(std::is_same_v<basic_simd_type, simd_type>);
    static_assert(
        std::is_same_v<typename basic_simd_type::mask_type, simd_mask_type>);
  }
  if constexpr (std::is_same_v<Abi, native_abi>) {
    using simd_type =
        Kokkos::Experimental::simd<DataType, basic_simd_type::size()>;
    using simd_mask_type =
        Kokkos::Experimental::simd_mask<DataType, basic_simd_type::size()>;
    static_assert(std::is_same_v<basic_simd_type, simd_type>);
    static_assert(
        std::is_same_v<typename basic_simd_type::mask_type, simd_mask_type>);
  }
}

template <typename /*Abi*/, typename DataType>
inline void host_test_simd_default_abi() {
#if defined(KOKKOS_ENABLE_HPX) || defined(KOKKOS_ENABLE_OPENMPTARGET) || \
    defined(KOKKOS_ENABLE_OPENACC) || defined(KOKKOS_ENABLE_CUDA) ||     \
    defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_SYCL)
  constexpr int expected_size = 1;
#elif defined(KOKKOS_ARCH_AVX512XEON)
  constexpr int expected_size = 512 / (CHAR_BIT * sizeof(DataType));
#elif defined(KOKKOS_ARCH_AVX2)
  constexpr int expected_size = 256 / (CHAR_BIT * sizeof(DataType));
#elif defined(KOKKOS_ARCH_ARM_SVE)
  constexpr int expected_size =
      __ARM_FEATURE_SVE_BITS / (CHAR_BIT * sizeof(DataType));
#elif defined(KOKKOS_ARCH_ARM_NEON)
  constexpr int expected_size = 128 / (CHAR_BIT * sizeof(DataType));
#else
  constexpr int expected_size = 1;
#endif

  using simd_type      = Kokkos::Experimental::simd<DataType>;
  using simd_mask_type = typename simd_type::mask_type;

  static_assert(simd_type::size() == expected_size);
  static_assert(simd_mask_type::size() == expected_size);
}

template <typename Abi, typename DataType>
inline void host_check_construction() {
  if constexpr (is_simd_avail_v<DataType, Abi>) {
    host_test_simd_traits<Abi, DataType>();
    host_test_mask_traits<Abi, DataType>();
    host_test_simd_alias<Abi, DataType>();
    host_test_simd_default_abi<Abi, DataType>();
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
  using simd_type = Kokkos::Experimental::basic_simd<DataType, Abi>;

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
  using mask_type = Kokkos::Experimental::basic_simd_mask<DataType, Abi>;

  mask_type default_mask(false);
  mask_type result(false);
  mask_type test_mask(KOKKOS_LAMBDA(std::size_t i) { return (i % 2 == 0); });
  mask_type copy_mask(test_mask);
  mask_type move_mask(std::move(copy_mask));
  default_mask = std::move(move_mask);
  result       = default_mask;

  kokkos_checker checker;
  checker.truth(all_of(test_mask == result));
}

template <typename Abi, typename DataType>
KOKKOS_INLINE_FUNCTION void device_check_construction() {
  if constexpr (is_type_v<Kokkos::Experimental::basic_simd<DataType, Abi>>) {
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
