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

#ifndef KOKKOS_TEST_SIMD_GENERATOR_CTORS_HPP
#define KOKKOS_TEST_SIMD_GENERATOR_CTORS_HPP

#include <Kokkos_SIMD.hpp>
#include <SIMDTesting_Utilities.hpp>

template <typename Abi, typename DataType>
inline void host_check_gen_ctor() {
  using simd_type             = Kokkos::Experimental::simd<DataType, Abi>;
  using mask_type             = typename simd_type::mask_type;
  constexpr std::size_t lanes = simd_type::size();

  DataType init[lanes];
  DataType expected[lanes];
  mask_type init_mask(false);

  for (std::size_t i = 0; i < lanes; ++i) {
    if (i % 3 == 0) init_mask[i] = true;
    init[i]     = 7;
    expected[i] = (init_mask[i]) ? init[i] * 9 : init[i];
  }

  simd_type rhs;
  rhs.copy_from(init, Kokkos::Experimental::simd_flag_default);

  simd_type blend;
  blend.copy_from(expected, Kokkos::Experimental::simd_flag_default);

#if !(defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_COMPILER_MSVC))
  if constexpr (std::is_same_v<Abi, Kokkos::Experimental::simd_abi::scalar>) {
    simd_type basic(KOKKOS_LAMBDA(std::size_t i) { return init[i]; });
    host_check_equality(basic, rhs, lanes);

    simd_type lhs(KOKKOS_LAMBDA(std::size_t i) { return init[i] * 9; });
    mask_type mask(KOKKOS_LAMBDA(std::size_t i) { return init_mask[i]; });
    simd_type result(
        KOKKOS_LAMBDA(std::size_t i) { return (mask[i]) ? lhs[i] : rhs[i]; });

    host_check_equality(blend, result, lanes);
  } else {
    simd_type basic([=](std::size_t i) { return init[i]; });
    host_check_equality(basic, rhs, lanes);

    simd_type lhs([=](std::size_t i) { return init[i] * 9; });
    mask_type mask([=](std::size_t i) { return init_mask[i]; });
    simd_type result(
        [=](std::size_t i) { return (mask[i]) ? lhs[i] : rhs[i]; });

    host_check_equality(blend, result, lanes);
  }
#endif
}

template <typename Abi, typename... DataTypes>
inline void host_check_gen_ctors_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (host_check_gen_ctor<Abi, DataTypes>(), ...);
}

template <typename... Abis>
inline void host_check_gen_ctors_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (host_check_gen_ctors_all_types<Abis>(DataTypes()), ...);
}

template <typename Abi, typename DataType>
KOKKOS_INLINE_FUNCTION void device_check_gen_ctor() {
  using simd_type             = Kokkos::Experimental::simd<DataType, Abi>;
  using mask_type             = typename simd_type::mask_type;
  constexpr std::size_t lanes = simd_type::size();

  DataType init[lanes];
  DataType expected[lanes];
  mask_type mask(false);

  for (std::size_t i = 0; i < lanes; ++i) {
    if (i % 3 == 0) mask[i] = true;
    init[i]     = 7;
    expected[i] = (mask[i]) ? init[i] * 9 : init[i];
  }

  simd_type basic(KOKKOS_LAMBDA(std::size_t i) { return init[i]; });
  simd_type rhs;
  rhs.copy_from(init, Kokkos::Experimental::simd_flag_default);
  device_check_equality(basic, rhs, lanes);

  simd_type lhs(KOKKOS_LAMBDA(std::size_t i) { return init[i] * 9; });
  simd_type result(
      KOKKOS_LAMBDA(std::size_t i) { return (mask[i]) ? lhs[i] : rhs[i]; });

  simd_type blend;
  blend.copy_from(expected, Kokkos::Experimental::simd_flag_default);
  device_check_equality(result, blend, lanes);
}

template <typename Abi, typename... DataTypes>
KOKKOS_INLINE_FUNCTION void device_check_gen_ctors_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (device_check_gen_ctor<Abi, DataTypes>(), ...);
}

template <typename... Abis>
KOKKOS_INLINE_FUNCTION void device_check_gen_ctors_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (device_check_gen_ctors_all_types<Abis>(DataTypes()), ...);
}

class simd_device_gen_ctor_functor {
 public:
  KOKKOS_INLINE_FUNCTION void operator()(int) const {
    device_check_gen_ctors_all_abis(
        Kokkos::Experimental::Impl::device_abi_set());
  }
};

TEST(simd, host_gen_ctors) {
  host_check_gen_ctors_all_abis(Kokkos::Experimental::Impl::host_abi_set());
}

TEST(simd, device_gen_ctors) {
  Kokkos::parallel_for(1, simd_device_gen_ctor_functor());
}

#endif
