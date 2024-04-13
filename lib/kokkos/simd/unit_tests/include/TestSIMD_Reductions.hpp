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

#ifndef KOKKOS_TEST_SIMD_REDUCTIONS_HPP
#define KOKKOS_TEST_SIMD_REDUCTIONS_HPP

#include <Kokkos_SIMD.hpp>
#include <SIMDTesting_Utilities.hpp>

template <typename Abi, typename Loader, typename ReductionOp, typename T>
inline void host_check_reduction_one_loader(ReductionOp reduce_op,
                                            std::size_t n, T const* args) {
  Loader loader;
  using simd_type = Kokkos::Experimental::simd<T, Abi>;
  using mask_type = typename Kokkos::Experimental::simd<T, Abi>::mask_type;
  constexpr std::size_t width = simd_type::size();

  for (std::size_t i = 0; i < n; i += width) {
    std::size_t const nremaining = n - i;
    std::size_t const nlanes     = Kokkos::min(nremaining, width);
    simd_type arg;
    bool const loaded_arg = loader.host_load(args + i, nlanes, arg);
    if (!loaded_arg) continue;

    mask_type mask(false);
    for (std::size_t j = 0; j < n; ++j) {
      mask[j] = true;
    }
    auto value    = where(mask, arg);
    auto expected = reduce_op.on_host_serial(value);
    auto computed = reduce_op.on_host(value);

    gtest_checker().equality(expected, computed);
  }
}

template <typename Abi, typename ReductionOp, typename T>
inline void host_check_reduction_all_loaders(ReductionOp reduce_op,
                                             std::size_t n, T const* args) {
  host_check_reduction_one_loader<Abi, load_element_aligned>(reduce_op, n,
                                                             args);
  host_check_reduction_one_loader<Abi, load_masked>(reduce_op, n, args);
  host_check_reduction_one_loader<Abi, load_as_scalars>(reduce_op, n, args);
}

template <typename Abi, typename DataType, size_t n>
inline void host_check_all_reductions(const DataType (&args)[n]) {
  host_check_reduction_all_loaders<Abi>(hmin(), n, args);
  host_check_reduction_all_loaders<Abi>(hmax(), n, args);
  host_check_reduction_all_loaders<Abi>(reduce(), n, args);
}

template <typename Abi, typename DataType>
inline void host_check_reductions() {
  constexpr size_t n = 11;

  if constexpr (std::is_signed_v<DataType>) {
    DataType const args[n] = {1, 2, -1, 10, 0, 1, -2, 10, 0, 1, -2};
    host_check_all_reductions<Abi>(args);
  } else {
    DataType const args[n] = {1, 2, 1, 10, 0, 1, 2, 10, 0, 1, 2};
    host_check_all_reductions<Abi>(args);
  }
}

template <typename Abi, typename... DataTypes>
inline void host_check_reductions_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (host_check_reductions<Abi, DataTypes>(), ...);
}

template <typename... Abis>
inline void host_check_reductions_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (host_check_reductions_all_types<Abis>(DataTypes()), ...);
}

template <typename Abi, typename Loader, typename ReductionOp, typename T>
KOKKOS_INLINE_FUNCTION void device_check_reduction_one_loader(
    ReductionOp reduce_op, std::size_t n, T const* args) {
  Loader loader;
  using simd_type = Kokkos::Experimental::simd<T, Abi>;
  using mask_type = typename Kokkos::Experimental::simd<T, Abi>::mask_type;
  constexpr std::size_t width = simd_type::size();

  for (std::size_t i = 0; i < n; i += width) {
    std::size_t const nremaining = n - i;
    std::size_t const nlanes     = Kokkos::min(nremaining, width);
    simd_type arg;
    bool const loaded_arg = loader.device_load(args + i, nlanes, arg);
    if (!loaded_arg) continue;

    mask_type mask(false);
    for (std::size_t j = 0; j < n; ++j) {
      mask[j] = true;
    }
    auto value    = where(mask, arg);
    auto expected = reduce_op.on_device_serial(value);
    auto computed = reduce_op.on_device(value);

    kokkos_checker().equality(expected, computed);
  }
}

template <typename Abi, typename ReductionOp, typename T>
KOKKOS_INLINE_FUNCTION void device_check_reduction_all_loaders(
    ReductionOp reduce_op, std::size_t n, T const* args) {
  device_check_reduction_one_loader<Abi, load_element_aligned>(reduce_op, n,
                                                               args);
  device_check_reduction_one_loader<Abi, load_masked>(reduce_op, n, args);
  device_check_reduction_one_loader<Abi, load_as_scalars>(reduce_op, n, args);
}

template <typename Abi, typename DataType, size_t n>
KOKKOS_INLINE_FUNCTION void device_check_all_reductions(
    const DataType (&args)[n]) {
  device_check_reduction_all_loaders<Abi>(hmin(), n, args);
  device_check_reduction_all_loaders<Abi>(hmax(), n, args);
  device_check_reduction_all_loaders<Abi>(reduce(), n, args);
}

template <typename Abi, typename DataType>
KOKKOS_INLINE_FUNCTION void device_check_reductions() {
  constexpr size_t n = 11;

  if constexpr (std::is_signed_v<DataType>) {
    DataType const args[n] = {1, 2, -1, 10, 0, 1, -2, 10, 0, 1, -2};
    device_check_all_reductions<Abi>(args);
  } else {
    DataType const args[n] = {1, 2, 1, 10, 0, 1, 2, 10, 0, 1, 2};
    device_check_all_reductions<Abi>(args);
  }
}

template <typename Abi, typename... DataTypes>
KOKKOS_INLINE_FUNCTION void device_check_reductions_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (device_check_reductions<Abi, DataTypes>(), ...);
}

template <typename... Abis>
KOKKOS_INLINE_FUNCTION void device_check_reductions_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (device_check_reductions_all_types<Abis>(DataTypes()), ...);
}

class simd_device_reduction_functor {
 public:
  KOKKOS_INLINE_FUNCTION void operator()(int) const {
    device_check_reductions_all_abis(
        Kokkos::Experimental::Impl::device_abi_set());
  }
};

TEST(simd, host_reductions) {
  host_check_reductions_all_abis(Kokkos::Experimental::Impl::host_abi_set());
}

TEST(simd, device_reductions) {
#ifdef KOKKOS_ENABLE_OPENMPTARGET  // FIXME_OPENMPTARGET
  GTEST_SKIP()
      << "skipping because of a non-deterministic failure reporting: "
         "Failure to synchronize stream (nil): Error in "
         "cuStreamSynchronize: an illegal memory access was encountered";
#endif
  Kokkos::parallel_for(1, simd_device_reduction_functor());
}

#endif
