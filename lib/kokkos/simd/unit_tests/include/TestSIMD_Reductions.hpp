// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_SIMD_REDUCTIONS_HPP
#define KOKKOS_TEST_SIMD_REDUCTIONS_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.simd;
import kokkos.simd_impl;
#else
#include <Kokkos_SIMD.hpp>
#endif
#include <SIMDTesting_Utilities.hpp>

template <typename T, typename ReductionOp>
struct get_identity {
  KOKKOS_INLINE_FUNCTION T operator()() { return T(); }
};

template <typename T, typename BinaryOp>
struct get_identity<T, masked_reduce<BinaryOp>> {
  KOKKOS_INLINE_FUNCTION T operator()() {
    return Kokkos::Experimental::Impl::Identity<T, BinaryOp>();
  }
};

template <typename Abi, typename Loader, typename ReductionOp, typename T>
inline void host_check_reduction_one_loader(ReductionOp reduce_op,
                                            std::size_t n, T const* args) {
  Loader loader;
  using simd_type = Kokkos::Experimental::basic_simd<T, Abi>;
  using mask_type =
      typename Kokkos::Experimental::basic_simd<T, Abi>::mask_type;
  constexpr std::size_t width = simd_type::size();

  for (std::size_t i = 0; i < n; i += width) {
    std::size_t const nremaining = n - i;
    std::size_t const nlanes     = Kokkos::min(nremaining, width);
    simd_type arg;
    bool const loaded_arg = loader.host_load(args + i, nlanes, arg);
    if (!loaded_arg) continue;

    T true_identity = get_identity<T, ReductionOp>{}();
    T test_identity = 12;
    if constexpr (std::is_same_v<Abi, Kokkos::Experimental::simd_abi::scalar>) {
      mask_type mask_false(false);

      auto expected = reduce_op.on_host_serial(arg, test_identity, mask_false);
      auto computed = reduce_op.on_host(arg, test_identity, mask_false);
      gtest_checker().equality(expected, computed);

      mask_type mask_true(true);
      expected = reduce_op.on_host_serial(arg, true_identity, mask_true);
      computed = reduce_op.on_host(arg, true_identity, mask_true);

      gtest_checker().equality(expected, computed);
    } else {
      mask_type mask_false(false);
      auto expected = reduce_op.on_host_serial(arg, test_identity, mask_false);
      auto computed = reduce_op.on_host(arg, test_identity, mask_false);
      gtest_checker().equality(expected, computed);

      for (std::size_t j = 0; j < mask_type::size(); ++j) {
        mask_type mask([=](std::size_t idx) { return idx >= j; });
        expected = reduce_op.on_host_serial(arg, true_identity, mask);
        computed = reduce_op.on_host(arg, true_identity, mask);
        gtest_checker().equality(expected, computed);
      }
    }
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
  host_check_reduction_all_loaders<Abi>(reduce_min(), n, args);
  host_check_reduction_all_loaders<Abi>(reduce_max(), n, args);
  host_check_reduction_all_loaders<Abi>(reduce<std::plus<>>(), n, args);
  host_check_reduction_all_loaders<Abi>(reduce<std::multiplies<>>(), n, args);

  host_check_reduction_all_loaders<Abi>(masked_reduce_min(), n, args);
  host_check_reduction_all_loaders<Abi>(masked_reduce_max(), n, args);
  host_check_reduction_all_loaders<Abi>(masked_reduce<std::plus<>>(), n, args);
  host_check_reduction_all_loaders<Abi>(masked_reduce<std::multiplies<>>(), n,
                                        args);
}

template <typename Abi, typename DataType>
inline void host_check_reductions() {
  if constexpr (is_simd_avail_v<DataType, Abi>) {
    constexpr size_t n = 16;

    if constexpr (std::is_signed_v<DataType>) {
      DataType const args[n] = {1, 2, -1, 10,  0, 1,  -2,  10,
                                0, 1, -2, -15, 5, 17, -22, 20};
      host_check_all_reductions<Abi>(args);
    } else {
      DataType const args[n] = {1, 2, 1, 10, 0, 1,  2,  10,
                                0, 1, 2, 15, 5, 17, 22, 20};
      host_check_all_reductions<Abi>(args);
    }
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
  using simd_type = Kokkos::Experimental::basic_simd<T, Abi>;
  using mask_type =
      typename Kokkos::Experimental::basic_simd<T, Abi>::mask_type;
  constexpr std::size_t width = simd_type::size();

  T true_identity = get_identity<T, ReductionOp>{}();
  T test_identity = 12;
  for (std::size_t i = 0; i < n; i += width) {
    std::size_t const nremaining = n - i;
    std::size_t const nlanes     = Kokkos::min(nremaining, width);
    simd_type arg;
    bool const loaded_arg = loader.device_load(args + i, nlanes, arg);
    if (!loaded_arg) continue;

    mask_type mask_false(false);
    auto expected = reduce_op.on_device_serial(arg, test_identity, mask_false);
    auto computed = reduce_op.on_device(arg, test_identity, mask_false);
    kokkos_checker().equality(expected, computed);

    for (std::size_t j = 0; j < mask_type::size(); ++j) {
      mask_type mask(KOKKOS_LAMBDA(std::size_t idx) { return idx >= j; });
      expected = reduce_op.on_device_serial(arg, true_identity, mask);
      computed = reduce_op.on_device(arg, true_identity, mask);
      kokkos_checker().equality(expected, computed);
    }
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
  device_check_reduction_all_loaders<Abi>(reduce_min(), n, args);
  device_check_reduction_all_loaders<Abi>(reduce_max(), n, args);
  device_check_reduction_all_loaders<Abi>(reduce<std::plus<>>(), n, args);
  device_check_reduction_all_loaders<Abi>(reduce<std::multiplies<>>(), n, args);

  device_check_reduction_all_loaders<Abi>(masked_reduce_min(), n, args);
  device_check_reduction_all_loaders<Abi>(masked_reduce_max(), n, args);
  device_check_reduction_all_loaders<Abi>(masked_reduce<std::plus<>>(), n,
                                          args);
  device_check_reduction_all_loaders<Abi>(masked_reduce<std::multiplies<>>(), n,
                                          args);
}

template <typename Abi, typename DataType>
KOKKOS_INLINE_FUNCTION void device_check_reductions() {
  if constexpr (is_type_v<Kokkos::Experimental::basic_simd<DataType, Abi>>) {
    constexpr size_t n = 16;

    if constexpr (std::is_signed_v<DataType>) {
      DataType const args[n] = {1, 2, -1, 10,  0, 1,  -2,  10,
                                0, 1, -2, -15, 5, 17, -22, 20};
      device_check_all_reductions<Abi>(args);
    } else {
      DataType const args[n] = {1, 2, 1, 10, 0, 1,  2,  10,
                                0, 1, 2, 15, 5, 17, 22, 20};
      device_check_all_reductions<Abi>(args);
    }
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
#if defined(KOKKOS_ENABLE_OPENACC) && \
    defined(KOKKOS_COMPILER_CLANG)  // FIXME_CLACC
  GTEST_SKIP()
      << "skipping because of a non-deterministic failure reporting: "
         "Failure to synchronize stream (nil): Error in "
         "cuStreamSynchronize: an illegal memory access was encountered";
#endif
  Kokkos::parallel_for(1, simd_device_reduction_functor());
}

#endif
