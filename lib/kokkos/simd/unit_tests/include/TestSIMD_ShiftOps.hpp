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

#ifndef KOKKOS_TEST_SIMD_SHIFT_OPS_HPP
#define KOKKOS_TEST_SIMD_SHIFT_OPS_HPP

#include <Kokkos_SIMD.hpp>
#include <SIMDTesting_Utilities.hpp>

template <typename Abi, typename Loader, typename ShiftOp, typename DataType>
inline void host_check_shift_on_one_loader(ShiftOp shift_op,
                                           DataType test_vals[],
                                           DataType shift_by[], std::size_t n) {
  using simd_type             = Kokkos::Experimental::simd<DataType, Abi>;
  constexpr std::size_t width = simd_type::size();
  Loader loader;

  for (std::size_t i = 0; i < n; ++i) {
    simd_type simd_vals;
    bool const loaded_arg = loader.host_load(test_vals, width, simd_vals);
    if (!loaded_arg) {
      continue;
    }

    simd_type expected_result;

    for (std::size_t lane = 0; lane < width; ++lane) {
      DataType value = simd_vals[lane];
      expected_result[lane] =
          shift_op.on_host(value, static_cast<int>(shift_by[i]));
      EXPECT_EQ(value, value);
    }

    simd_type const computed_result =
        shift_op.on_host(simd_vals, static_cast<int>(shift_by[i]));
    host_check_equality(expected_result, computed_result, width);
  }
}

template <typename Abi, typename Loader, typename ShiftOp, typename DataType>
inline void host_check_shift_by_lanes_on_one_loader(
    ShiftOp shift_op, DataType test_vals[],
    Kokkos::Experimental::simd<DataType, Abi>& shift_by) {
  using simd_type             = Kokkos::Experimental::simd<DataType, Abi>;
  constexpr std::size_t width = simd_type::size();
  Loader loader;

  simd_type simd_vals;
  bool const loaded_arg = loader.host_load(test_vals, width, simd_vals);
  ASSERT_TRUE(loaded_arg);

  simd_type expected_result;

  for (std::size_t lane = 0; lane < width; ++lane) {
    DataType value = simd_vals[lane];
    expected_result[lane] =
        shift_op.on_host(value, static_cast<int>(shift_by[lane]));
    EXPECT_EQ(value, value);
  }
  simd_type const computed_result = shift_op.on_host(simd_vals, shift_by);
  host_check_equality(expected_result, computed_result, width);
}

template <typename Abi, typename ShiftOp, typename DataType>
inline void host_check_shift_op_all_loaders(ShiftOp shift_op,
                                            DataType test_vals[],
                                            DataType shift_by[],
                                            std::size_t n) {
  host_check_shift_on_one_loader<Abi, load_element_aligned>(shift_op, test_vals,
                                                            shift_by, n);
  host_check_shift_on_one_loader<Abi, load_masked>(shift_op, test_vals,
                                                   shift_by, n);
  host_check_shift_on_one_loader<Abi, load_as_scalars>(shift_op, test_vals,
                                                       shift_by, n);
  host_check_shift_on_one_loader<Abi, load_vector_aligned>(shift_op, test_vals,
                                                           shift_by, n);

  Kokkos::Experimental::simd<DataType, Abi> shift_by_lanes;
  shift_by_lanes.copy_from(shift_by, Kokkos::Experimental::simd_flag_default);

  host_check_shift_by_lanes_on_one_loader<Abi, load_element_aligned>(
      shift_op, test_vals, shift_by_lanes);
  host_check_shift_by_lanes_on_one_loader<Abi, load_masked>(shift_op, test_vals,
                                                            shift_by_lanes);
  host_check_shift_by_lanes_on_one_loader<Abi, load_as_scalars>(
      shift_op, test_vals, shift_by_lanes);
  host_check_shift_by_lanes_on_one_loader<Abi, load_vector_aligned>(
      shift_op, test_vals, shift_by_lanes);
}

template <typename Abi, typename DataType>
inline void host_check_shift_ops() {
  if constexpr (std::is_integral_v<DataType>) {
    using simd_type                 = Kokkos::Experimental::simd<DataType, Abi>;
    constexpr std::size_t width     = simd_type::size();
    constexpr std::size_t num_cases = 8;
    constexpr size_t alignment =
        Kokkos::Experimental::simd<DataType, Abi>::size() * sizeof(DataType);

    DataType max = std::numeric_limits<DataType>::max();

    alignas(alignment) DataType shift_by[num_cases] = {
        0, 1, 3, width / 2, width / 2 + 1, width - 1, width, width + 1};
    alignas(alignment) DataType test_vals[width];
    for (std::size_t i = 0; i < width; ++i) {
      DataType inc = max / width;
      test_vals[i] = i * inc + 1;
    }

    host_check_shift_op_all_loaders<Abi>(shift_right(), test_vals, shift_by,
                                         num_cases);
    host_check_shift_op_all_loaders<Abi>(shift_left(), test_vals, shift_by,
                                         num_cases);

    if constexpr (std::is_signed_v<DataType>) {
      for (std::size_t i = 0; i < width; ++i) test_vals[i] *= -1;
      host_check_shift_op_all_loaders<Abi>(shift_right(), test_vals, shift_by,
                                           num_cases);
      host_check_shift_op_all_loaders<Abi>(shift_left(), test_vals, shift_by,
                                           num_cases);
    }
  }
}

template <typename Abi, typename... DataTypes>
inline void host_check_shift_ops_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (host_check_shift_ops<Abi, DataTypes>(), ...);
}

template <typename... Abis>
inline void host_check_shift_ops_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (host_check_shift_ops_all_types<Abis>(DataTypes()), ...);
}

template <typename Abi, typename Loader, typename ShiftOp, typename DataType>
KOKKOS_INLINE_FUNCTION void device_check_shift_on_one_loader(
    ShiftOp shift_op, DataType test_vals[], DataType shift_by[],
    std::size_t n) {
  using simd_type             = Kokkos::Experimental::simd<DataType, Abi>;
  constexpr std::size_t width = simd_type::size();
  Loader loader;

  for (std::size_t i = 0; i < n; ++i) {
    simd_type simd_vals;
    bool const loaded_arg = loader.device_load(test_vals, width, simd_vals);
    if (!loaded_arg) {
      continue;
    }

    simd_type expected_result;

    for (std::size_t lane = 0; lane < width; ++lane) {
      expected_result[lane] = shift_op.on_device(DataType(simd_vals[lane]),
                                                 static_cast<int>(shift_by[i]));
    }

    simd_type const computed_result =
        shift_op.on_device(simd_vals, static_cast<int>(shift_by[i]));
    device_check_equality(expected_result, computed_result, width);
  }
}

template <typename Abi, typename Loader, typename ShiftOp, typename DataType>
KOKKOS_INLINE_FUNCTION void device_check_shift_by_lanes_on_one_loader(
    ShiftOp shift_op, DataType test_vals[],
    Kokkos::Experimental::simd<DataType, Abi>& shift_by) {
  using simd_type             = Kokkos::Experimental::simd<DataType, Abi>;
  constexpr std::size_t width = simd_type::size();
  Loader loader;
  simd_type simd_vals;
  loader.device_load(test_vals, width, simd_vals);

  simd_type expected_result;

  for (std::size_t lane = 0; lane < width; ++lane) {
    expected_result[lane] = shift_op.on_device(
        DataType(simd_vals[lane]), static_cast<int>(shift_by[lane]));
  }
  simd_type const computed_result = shift_op.on_device(simd_vals, shift_by);
  device_check_equality(expected_result, computed_result, width);
}

template <typename Abi, typename ShiftOp, typename DataType>
KOKKOS_INLINE_FUNCTION void device_check_shift_op_all_loaders(
    ShiftOp shift_op, DataType test_vals[], DataType shift_by[],
    std::size_t n) {
  device_check_shift_on_one_loader<Abi, load_element_aligned>(
      shift_op, test_vals, shift_by, n);
  device_check_shift_on_one_loader<Abi, load_masked>(shift_op, test_vals,
                                                     shift_by, n);
  device_check_shift_on_one_loader<Abi, load_as_scalars>(shift_op, test_vals,
                                                         shift_by, n);
  device_check_shift_on_one_loader<Abi, load_vector_aligned>(
      shift_op, test_vals, shift_by, n);

  Kokkos::Experimental::simd<DataType, Abi> shift_by_lanes;
  shift_by_lanes.copy_from(shift_by, Kokkos::Experimental::simd_flag_default);

  device_check_shift_by_lanes_on_one_loader<Abi, load_element_aligned>(
      shift_op, test_vals, shift_by_lanes);
  device_check_shift_by_lanes_on_one_loader<Abi, load_masked>(
      shift_op, test_vals, shift_by_lanes);
  device_check_shift_by_lanes_on_one_loader<Abi, load_as_scalars>(
      shift_op, test_vals, shift_by_lanes);
  device_check_shift_by_lanes_on_one_loader<Abi, load_vector_aligned>(
      shift_op, test_vals, shift_by_lanes);
}

template <typename Abi, typename DataType>
KOKKOS_INLINE_FUNCTION void device_check_shift_ops() {
  if constexpr (std::is_integral_v<DataType>) {
    using simd_type                 = Kokkos::Experimental::simd<DataType, Abi>;
    constexpr std::size_t width     = simd_type::size();
    constexpr std::size_t num_cases = 8;

    DataType max = Kokkos::reduction_identity<DataType>::max();

    DataType shift_by[num_cases] = {
        0, 1, 3, width / 2, width / 2 + 1, width - 1, width, width + 1};
    DataType test_vals[width];

    for (std::size_t i = 0; i < width; ++i) {
      DataType inc = max / width;
      test_vals[i] = i * inc + 1;
    }

    device_check_shift_op_all_loaders<Abi>(shift_right(), test_vals, shift_by,
                                           num_cases);
    device_check_shift_op_all_loaders<Abi>(shift_left(), test_vals, shift_by,
                                           num_cases);

    if constexpr (std::is_signed_v<DataType>) {
      for (std::size_t i = 0; i < width; ++i) test_vals[i] *= -1;
      device_check_shift_op_all_loaders<Abi>(shift_right(), test_vals, shift_by,
                                             num_cases);
      device_check_shift_op_all_loaders<Abi>(shift_left(), test_vals, shift_by,
                                             num_cases);
    }
  }
}

template <typename Abi, typename... DataTypes>
KOKKOS_INLINE_FUNCTION void device_check_shift_ops_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (device_check_shift_ops<Abi, DataTypes>(), ...);
}

template <typename... Abis>
KOKKOS_INLINE_FUNCTION void device_check_shift_ops_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (device_check_shift_ops_all_types<Abis>(DataTypes()), ...);
}

class simd_device_shift_ops_functor {
 public:
  KOKKOS_INLINE_FUNCTION void operator()(int) const {
    device_check_shift_ops_all_abis(
        Kokkos::Experimental::Impl::device_abi_set());
  }
};

TEST(simd, host_shift_ops) {
  host_check_shift_ops_all_abis(Kokkos::Experimental::Impl::host_abi_set());
}

TEST(simd, device_shift_ops) {
  Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::IndexType<int>>(0, 1),
                       simd_device_shift_ops_functor());
}

#endif
