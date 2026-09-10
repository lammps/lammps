// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_SIMD_COMPARISON_OPS_HPP
#define KOKKOS_TEST_SIMD_COMPARISON_OPS_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.simd;
#else
#include <Kokkos_SIMD.hpp>
#endif
#include <SIMDTesting_Utilities.hpp>

template <class Abi, class Loader, class BinaryOp, class T>
void host_check_comparison_op_one_loader(
    BinaryOp binary_op, Kokkos::Experimental::Impl::simd_size_t n,
    T const* first_args, T const* second_args) {
  Loader loader;
  using simd_type = Kokkos::Experimental::basic_simd<T, Abi>;
  using mask_type = typename simd_type::mask_type;
  using size_type = Kokkos::Experimental::Impl::simd_size_t;

  constexpr size_type width = simd_type::size();
  for (size_type i = 0; i < n; i += width) {
    size_type const nremaining = n - i;
    size_type const nlanes     = Kokkos::min(nremaining, width);

    simd_type first_arg;
    bool const loaded_first_arg =
        loader.host_load(first_args + i, nlanes, first_arg);
    simd_type second_arg;
    bool const loaded_second_arg =
        loader.host_load(second_args + i, nlanes, second_arg);
    if (!(loaded_first_arg && loaded_second_arg)) continue;

    bool expected_val[width];
    for (size_type lane = 0; lane < width; ++lane) {
      expected_val[lane] =
          binary_op.on_host(T(first_arg[lane]), T(second_arg[lane]));
    }

    mask_type const expected_result(
        KOKKOS_LAMBDA(size_type idx) { return expected_val[idx]; });
    mask_type const computed_result = binary_op.on_host(first_arg, second_arg);
    host_check_mask_equality(expected_result, computed_result, nlanes);
  }
}

template <class Abi, class Op, class T>
inline void host_check_comparison_op_all_loaders(
    Op op, Kokkos::Experimental::Impl::simd_size_t n, T const* first_args,
    T const* second_args) {
  host_check_comparison_op_one_loader<Abi, load_element_aligned>(
      op, n, first_args, second_args);
  host_check_comparison_op_one_loader<Abi, load_masked>(op, n, first_args,
                                                        second_args);
  host_check_comparison_op_one_loader<Abi, load_as_scalars>(op, n, first_args,
                                                            second_args);
  host_check_comparison_op_one_loader<Abi, load_vector_aligned>(
      op, n, first_args, second_args);
}

template <typename Abi, typename DataType,
          Kokkos::Experimental::Impl::simd_size_t n>
inline void host_check_all_comparison_ops(const DataType (&first_args)[n],
                                          const DataType (&second_args)[n]) {
  host_check_comparison_op_all_loaders<Abi>(equal(), n, first_args,
                                            second_args);
  host_check_comparison_op_all_loaders<Abi>(not_equal(), n, first_args,
                                            second_args);
  host_check_comparison_op_all_loaders<Abi>(less_than(), n, first_args,
                                            second_args);
  host_check_comparison_op_all_loaders<Abi>(less_equal(), n, first_args,
                                            second_args);
  host_check_comparison_op_all_loaders<Abi>(greater_than(), n, first_args,
                                            second_args);
  host_check_comparison_op_all_loaders<Abi>(greater_equal(), n, first_args,
                                            second_args);
}

template <typename Abi, typename DataType>
inline void host_check_comparison_ops() {
  if constexpr (is_simd_avail_v<DataType, Abi>) {
    constexpr size_t alignment =
        Kokkos::Experimental::basic_simd<DataType, Abi>::size() *
        sizeof(DataType);

    if constexpr (!std::is_integral_v<DataType>) {
      alignas(alignment) DataType const first_args[] = {
          0.1, 0.4, 0.5,  0.7, 1.0, 1.5,  -2.0, 10.0,
          0.0, 1.2, -2.8, 3.0, 4.0, -0.1, 5.0,  -0.2};
      alignas(alignment) DataType const second_args[] = {
          1.0,  0.2,  1.1,  1.8, -0.1,  -3.0, -2.4, 1.0,
          13.0, -3.2, -2.1, 3.0, -15.0, -0.5, -0.2, -0.2};
      host_check_all_comparison_ops<Abi>(first_args, second_args);
    } else {
      if constexpr (std::is_signed_v<DataType>) {
        alignas(alignment) DataType const first_args[] = {
            1, 2, -1, 10, 0, 1, -2, 10, 0, 1, -2, -3, 7, 4, -9, -15};
        alignas(alignment) DataType const second_args[] = {
            1, 2, 1, 1, 1, -3, -2, 1, 13, -3, -2, 10, -15, 7, 2, -10};
        host_check_all_comparison_ops<Abi>(first_args, second_args);
      } else {
        alignas(alignment) DataType const first_args[] = {
            1, 2, 1, 10, 0, 1, 2, 10, 0, 1, 2, 11, 5, 8, 2, 14};
        alignas(alignment) DataType const second_args[] = {
            1, 2, 1, 1, 1, 3, 2, 1, 13, 3, 2, 3, 6, 20, 5, 14};
        host_check_all_comparison_ops<Abi>(first_args, second_args);
      }
    }
  }
}

template <typename Abi, typename... DataTypes>
inline void host_check_comparison_ops_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (host_check_comparison_ops<Abi, DataTypes>(), ...);
}

template <typename... Abis>
inline void host_check_comparison_ops_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (host_check_comparison_ops_all_types<Abis>(DataTypes()), ...);
}

template <class Abi, class Loader, class BinaryOp, class T>
KOKKOS_INLINE_FUNCTION void device_check_comparison_op_one_loader(
    BinaryOp binary_op, Kokkos::Experimental::Impl::simd_size_t n,
    T const* first_args, T const* second_args) {
  Loader loader;
  using simd_type = Kokkos::Experimental::basic_simd<T, Abi>;
  using mask_type = typename simd_type::mask_type;
  using size_type = Kokkos::Experimental::Impl::simd_size_t;

  constexpr size_type width = simd_type::size();
  for (size_type i = 0; i < n; i += width) {
    size_type const nremaining = n - i;
    size_type const nlanes     = Kokkos::min(nremaining, width);

    simd_type first_arg;
    bool const loaded_first_arg =
        loader.device_load(first_args + i, nlanes, first_arg);
    simd_type second_arg;
    bool const loaded_second_arg =
        loader.device_load(second_args + i, nlanes, second_arg);
    if (!(loaded_first_arg && loaded_second_arg)) continue;

    T expected_val[width];
    for (size_type lane = 0; lane < width; ++lane) {
      expected_val[lane] =
          binary_op.on_device(T(first_arg[lane]), T(second_arg[lane]));
    }

    mask_type const expected_result(
        KOKKOS_LAMBDA(size_type idx) { return expected_val[idx]; });
    mask_type const computed_result =
        binary_op.on_device(first_arg, second_arg);
    device_check_mask_equality(expected_result, computed_result, nlanes);
  }
}

template <class Abi, class Op, class T>
KOKKOS_INLINE_FUNCTION void device_check_comparison_op_all_loaders(
    Op op, Kokkos::Experimental::Impl::simd_size_t n, T const* first_args,
    T const* second_args) {
  device_check_comparison_op_one_loader<Abi, load_element_aligned>(
      op, n, first_args, second_args);
  device_check_comparison_op_one_loader<Abi, load_masked>(op, n, first_args,
                                                          second_args);
  device_check_comparison_op_one_loader<Abi, load_as_scalars>(op, n, first_args,
                                                              second_args);
  device_check_comparison_op_one_loader<Abi, load_vector_aligned>(
      op, n, first_args, second_args);
}

template <typename Abi, typename DataType,
          Kokkos::Experimental::Impl::simd_size_t n>
KOKKOS_INLINE_FUNCTION void device_check_all_comparison_ops(
    const DataType (&first_args)[n], const DataType (&second_args)[n]) {
  device_check_comparison_op_all_loaders<Abi>(equal(), n, first_args,
                                              second_args);
  device_check_comparison_op_all_loaders<Abi>(not_equal(), n, first_args,
                                              second_args);
  device_check_comparison_op_all_loaders<Abi>(less_than(), n, first_args,
                                              second_args);
  device_check_comparison_op_all_loaders<Abi>(less_equal(), n, first_args,
                                              second_args);
  device_check_comparison_op_all_loaders<Abi>(greater_than(), n, first_args,
                                              second_args);
  device_check_comparison_op_all_loaders<Abi>(greater_equal(), n, first_args,
                                              second_args);
}

template <typename Abi, typename DataType>
KOKKOS_INLINE_FUNCTION void device_check_comparison_ops() {
  if constexpr (is_simd_avail_v<DataType, Abi>) {
    constexpr size_t alignment =
        Kokkos::Experimental::basic_simd<DataType, Abi>::size() *
        sizeof(DataType);

    if constexpr (!std::is_integral_v<DataType>) {
      alignas(alignment) DataType const first_args[] = {
          0.1, 0.4, 0.5,  0.7, 1.0, 1.5,  -2.0, 10.0,
          0.0, 1.2, -2.8, 3.0, 4.0, -0.1, 5.0,  -0.2};
      alignas(alignment) DataType const second_args[] = {
          1.0,  0.2,  1.1,  1.8, -0.1,  -3.0, -2.4, 1.0,
          13.0, -3.2, -2.1, 3.0, -15.0, -0.5, -0.2, -0.2};
      device_check_all_comparison_ops<Abi>(first_args, second_args);
    } else {
      if constexpr (std::is_signed_v<DataType>) {
        alignas(alignment) DataType const first_args[] = {
            1, 2, -1, 10, 0, 1, -2, 10, 0, 1, -2, -3, 7, 4, -9, -15};
        alignas(alignment) DataType const second_args[] = {
            1, 2, 1, 1, 1, -3, -2, 1, 13, -3, -2, 10, -15, 7, 2, -10};
        device_check_all_comparison_ops<Abi>(first_args, second_args);
      } else {
        alignas(alignment) DataType const first_args[] = {
            1, 2, 1, 10, 0, 1, 2, 10, 0, 1, 2, 11, 5, 8, 2, 14};
        alignas(alignment) DataType const second_args[] = {
            1, 2, 1, 1, 1, 3, 2, 1, 13, 3, 2, 3, 6, 20, 5, 14};
        device_check_all_comparison_ops<Abi>(first_args, second_args);
      }
    }
  }
}

template <typename Abi, typename... DataTypes>
KOKKOS_INLINE_FUNCTION void device_check_comparison_ops_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (device_check_comparison_ops<Abi, DataTypes>(), ...);
}

template <typename... Abis>
KOKKOS_INLINE_FUNCTION void device_check_comparison_ops_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (device_check_comparison_ops_all_types<Abis>(DataTypes()), ...);
}

class simd_device_comparison_ops_functor {
 public:
  KOKKOS_INLINE_FUNCTION void operator()(int) const {
    device_check_comparison_ops_all_abis(
        Kokkos::Experimental::Impl::device_abi_set());
  }
};

template <typename Abi, typename DataType>
inline void host_check_mask_comparison_ops() {
  if constexpr (is_simd_avail_v<DataType, Abi>) {
    using mask_type = Kokkos::Experimental::basic_simd_mask<DataType, Abi>;
    using size_type = Kokkos::Experimental::Impl::simd_size_t;

    using Kokkos::Experimental::all_of;
    using Kokkos::Experimental::none_of;

    constexpr size_type width = mask_type::size();

    mask_type const all(true);
    mask_type const none(false);

    gtest_checker checker;

    checker.truth(all_of(all == all));
    checker.truth(none_of(all != all));
    checker.truth(all_of(none == none));
    checker.truth(none_of(none != none));

    if constexpr (!std::is_same_v<Abi,
                                  Kokkos::Experimental::simd_abi::scalar>) {
      mask_type const hi([=](std::size_t const i) { return i >= width / 2; });
      mask_type const lo([=](std::size_t const i) { return i < width / 2; });
      mask_type const even([=](std::size_t const i) { return i % 2 == 0; });
      mask_type const odd([=](std::size_t const i) { return i % 2 == 1; });

      checker.truth(all_of(hi == hi));
      checker.truth(none_of(hi != hi));
      checker.truth(all_of(lo == lo));
      checker.truth(none_of(lo != lo));
      checker.truth(none_of(hi == lo));
      checker.truth(all_of(hi != lo));

      checker.truth(all_of(even == even));
      checker.truth(none_of(even != even));
      checker.truth(all_of(odd == odd));
      checker.truth(none_of(odd != odd));
      checker.truth(none_of(even == odd));
      checker.truth(all_of(even != odd));
    }
  }
}

template <typename Abi, typename... DataTypes>
inline void host_check_mask_comparison_ops_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (host_check_mask_comparison_ops<Abi, DataTypes>(), ...);
}

template <typename... Abis>
inline void host_check_mask_comparison_ops_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (host_check_mask_comparison_ops_all_types<Abis>(DataTypes()), ...);
}

template <typename Abi, typename DataType>
KOKKOS_INLINE_FUNCTION void device_check_mask_comparison_ops() {
  if constexpr (is_simd_avail_v<DataType, Abi>) {
    using mask_type = Kokkos::Experimental::basic_simd_mask<DataType, Abi>;

    using Kokkos::Experimental::all_of;
    using Kokkos::Experimental::none_of;

    mask_type const all(true);
    mask_type const none(false);

    kokkos_checker checker;

    checker.truth(all_of(all == all));
    checker.truth(none_of(all != all));
    checker.truth(all_of(none == none));
    checker.truth(none_of(none != none));
  }
}

template <typename Abi, typename... DataTypes>
KOKKOS_INLINE_FUNCTION void device_check_mask_comparison_ops_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (device_check_mask_comparison_ops<Abi, DataTypes>(), ...);
}

template <typename... Abis>
KOKKOS_INLINE_FUNCTION void device_check_mask_comparison_ops_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (device_check_mask_comparison_ops_all_types<Abis>(DataTypes()), ...);
}

class simd_device_mask_comparison_ops_functor {
 public:
  KOKKOS_INLINE_FUNCTION void operator()(int) const {
    device_check_mask_comparison_ops_all_abis(
        Kokkos::Experimental::Impl::device_abi_set());
  }
};

TEST(simd, host_comparison_ops) {
  host_check_comparison_ops_all_abis(
      Kokkos::Experimental::Impl::host_abi_set());
}

TEST(simd, host_mask_comparison_ops) {
  host_check_mask_comparison_ops_all_abis(
      Kokkos::Experimental::Impl::host_abi_set());
}

TEST(simd, device_comparison_ops) {
  Kokkos::parallel_for(1, simd_device_comparison_ops_functor());
  Kokkos::fence();
}

TEST(simd, device_mask_comparison_ops) {
  Kokkos::parallel_for(1, simd_device_mask_comparison_ops_functor());
  Kokkos::fence();
}

#endif
