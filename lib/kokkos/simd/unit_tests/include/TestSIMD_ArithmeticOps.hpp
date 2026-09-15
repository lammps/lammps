// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_SIMD_ARITHMETIC_OPS_HPP
#define KOKKOS_TEST_SIMD_ARITHMETIC_OPS_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.simd;
import kokkos.simd_impl;
#else
#include <Kokkos_SIMD.hpp>
#endif
#include <SIMDTesting_Utilities.hpp>

template <class Abi, class BinaryOp, class T>
void host_check_arithmetic_op_vector_vector(
    BinaryOp binary_op, Kokkos::Experimental::basic_simd<T, Abi> first_arg,
    Kokkos::Experimental::basic_simd<T, Abi> const& second_arg,
    Kokkos::Experimental::Impl::simd_size_t nlanes) {
  using simd_type = Kokkos::Experimental::basic_simd<T, Abi>;
  using size_type = Kokkos::Experimental::Impl::simd_size_t;

  constexpr size_type width = simd_type::size();
  T expected_val[width];
  for (size_type lane = 0; lane < width; ++lane) {
    expected_val[lane] = binary_op.on_host(first_arg[lane], second_arg[lane]);
  }

  simd_type expected_result =
      Kokkos::Experimental::simd_unchecked_load<simd_type>(
          expected_val, Kokkos::Experimental::simd_flag_default);
  simd_type const computed_result = binary_op.on_host(first_arg, second_arg);
  host_check_equality(expected_result, computed_result, nlanes);
}

template <class Abi, class BinaryOp, class T>
void host_check_arithmetic_op_vector_scalar(
    BinaryOp binary_op, Kokkos::Experimental::basic_simd<T, Abi> first_arg,
    T second_arg, Kokkos::Experimental::Impl::simd_size_t nlanes) {
  using simd_type = Kokkos::Experimental::basic_simd<T, Abi>;
  using size_type = Kokkos::Experimental::Impl::simd_size_t;

  constexpr size_type width = simd_type::size();
  T expected_val[width];
  for (size_type lane = 0; lane < width; ++lane) {
    expected_val[lane] = binary_op.on_host(first_arg[lane], second_arg);
  }

  simd_type expected_result =
      Kokkos::Experimental::simd_unchecked_load<simd_type>(
          expected_val, Kokkos::Experimental::simd_flag_default);
  simd_type const computed_result = binary_op.on_host(first_arg, second_arg);
  host_check_equality(expected_result, computed_result, nlanes);
}

template <class Abi, class BinaryOp, class T>
void host_check_arithmetic_op_scalar_vector(
    BinaryOp binary_op, T first_arg,
    Kokkos::Experimental::basic_simd<T, Abi> const& second_arg,
    Kokkos::Experimental::Impl::simd_size_t nlanes) {
  using simd_type = Kokkos::Experimental::basic_simd<T, Abi>;
  using size_type = Kokkos::Experimental::Impl::simd_size_t;

  constexpr size_type width = simd_type::size();
  T expected_val[width];
  for (size_type lane = 0; lane < width; ++lane) {
    expected_val[lane] = binary_op.on_host(first_arg, second_arg[lane]);
  }

  simd_type expected_result =
      Kokkos::Experimental::simd_unchecked_load<simd_type>(
          expected_val, Kokkos::Experimental::simd_flag_default);
  simd_type const computed_result = binary_op.on_host(first_arg, second_arg);
  host_check_equality(expected_result, computed_result, nlanes);
}

template <class Abi, class Loader, class BinaryOp, class T>
void host_check_arithmetic_op_one_loader(
    BinaryOp binary_op, Kokkos::Experimental::Impl::simd_size_t n,
    T const* first_args, T const* second_args) {
  Loader loader;
  using simd_type = Kokkos::Experimental::basic_simd<T, Abi>;
  using size_type = Kokkos::Experimental::Impl::simd_size_t;

  constexpr size_type width = simd_type::size();
  for (size_type i = 0; i < n; i += width) {
    const size_type nremaining = n - i;
    const size_type nlanes     = Kokkos::min(nremaining, width);
    if ((std::is_same_v<BinaryOp, divides> ||
         std::is_same_v<BinaryOp, divides_eq>)&&nremaining < width)
      continue;

    simd_type first_arg;
    bool const loaded_first_arg =
        loader.host_load(first_args + i, nlanes, first_arg);

    simd_type second_arg;
    bool const loaded_second_arg =
        loader.host_load(second_args + i, nlanes, second_arg);

    if (!(loaded_first_arg && loaded_second_arg)) continue;

    T first_scalar_arg  = first_args[i];
    T second_scalar_arg = second_args[i];
    host_check_arithmetic_op_vector_vector(binary_op, first_arg, second_arg,
                                           nlanes);
    host_check_arithmetic_op_vector_scalar(binary_op, first_arg,
                                           second_scalar_arg, nlanes);
    if constexpr (!(std::is_same_v<BinaryOp, plus_eq> ||
                    std::is_same_v<BinaryOp, minus_eq> ||
                    std::is_same_v<BinaryOp, multiplies_eq> ||
                    std::is_same_v<BinaryOp, divides_eq>)) {
      host_check_arithmetic_op_scalar_vector(binary_op, first_scalar_arg,
                                             second_arg, nlanes);
    }
  }
}

template <class Abi, class Op, class... T>
inline void host_check_arithmetic_op_all_loaders(
    Op op, Kokkos::Experimental::Impl::simd_size_t n, T const*... args) {
  host_check_arithmetic_op_one_loader<Abi, load_element_aligned>(op, n,
                                                                 args...);
  host_check_arithmetic_op_one_loader<Abi, load_masked>(op, n, args...);
  host_check_arithmetic_op_one_loader<Abi, load_as_scalars>(op, n, args...);
  host_check_arithmetic_op_one_loader<Abi, load_vector_aligned>(op, n, args...);
}

template <typename Abi, typename DataType, size_t n>
inline void host_check_all_arithmetic_ops(const DataType (&first_args)[n],
                                          const DataType (&second_args)[n]) {
  host_check_arithmetic_op_all_loaders<Abi>(plus(), n, first_args, second_args);
  host_check_arithmetic_op_all_loaders<Abi>(plus_eq(), n, first_args,
                                            second_args);
  host_check_arithmetic_op_all_loaders<Abi>(minus(), n, first_args,
                                            second_args);
  host_check_arithmetic_op_all_loaders<Abi>(minus_eq(), n, first_args,
                                            second_args);
  host_check_arithmetic_op_all_loaders<Abi>(multiplies(), n, first_args,
                                            second_args);
  host_check_arithmetic_op_all_loaders<Abi>(multiplies_eq(), n, first_args,
                                            second_args);
  host_check_arithmetic_op_all_loaders<Abi>(divides(), n, first_args,
                                            second_args);
  host_check_arithmetic_op_all_loaders<Abi>(divides_eq(), n, first_args,
                                            second_args);
}

template <typename Abi, typename DataType>
inline void host_check_arithmetic_ops() {
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
      host_check_all_arithmetic_ops<Abi>(first_args, second_args);
    } else {
      if constexpr (std::is_signed_v<DataType>) {
        alignas(alignment) DataType const first_args[] = {
            1, 2, -1, 10, 0, 1, -2, 10, 0, 1, -2, -3, 7, 4, -9, -15};
        alignas(alignment) DataType const second_args[] = {
            1, 2, 1, 1, 1, -3, -2, 1, 13, -3, -2, 10, -15, 7, 2, -10};
        host_check_all_arithmetic_ops<Abi>(first_args, second_args);
      } else {
        alignas(alignment) DataType const first_args[] = {
            1, 2, 1, 10, 0, 1, 2, 10, 0, 1, 2, 11, 5, 8, 2, 14};
        alignas(alignment) DataType const second_args[] = {
            1, 2, 1, 1, 1, 3, 2, 1, 13, 3, 2, 3, 6, 20, 5, 14};
        host_check_all_arithmetic_ops<Abi>(first_args, second_args);
      }
    }
  }
}

template <typename Abi, typename... DataTypes>
inline void host_check_arithmetic_ops_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (host_check_arithmetic_ops<Abi, DataTypes>(), ...);
}

template <typename... Abis>
inline void host_check_arithmetic_ops_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (host_check_arithmetic_ops_all_types<Abis>(DataTypes()), ...);
}

template <class Abi, class BinaryOp, class T>
KOKKOS_INLINE_FUNCTION void device_check_arithmetic_op_vector_vector(
    BinaryOp binary_op, Kokkos::Experimental::basic_simd<T, Abi> first_arg,
    Kokkos::Experimental::basic_simd<T, Abi> const& second_arg,
    Kokkos::Experimental::Impl::simd_size_t nlanes) {
  using simd_type = Kokkos::Experimental::basic_simd<T, Abi>;
  using size_type = Kokkos::Experimental::Impl::simd_size_t;

  constexpr size_type width = simd_type::size();
  T expected_val[width];
  for (size_type lane = 0; lane < width; ++lane) {
    expected_val[lane] = binary_op.on_device(first_arg[lane], second_arg[lane]);
  }

  simd_type expected_result =
      Kokkos::Experimental::simd_unchecked_load<simd_type>(
          expected_val, Kokkos::Experimental::simd_flag_default);
  simd_type const computed_result = binary_op.on_device(first_arg, second_arg);
  device_check_equality(expected_result, computed_result, nlanes);
}

template <class Abi, class BinaryOp, class T>
KOKKOS_INLINE_FUNCTION void device_check_arithmetic_op_vector_scalar(
    BinaryOp binary_op, Kokkos::Experimental::basic_simd<T, Abi> first_arg,
    T second_arg, Kokkos::Experimental::Impl::simd_size_t nlanes) {
  using simd_type = Kokkos::Experimental::basic_simd<T, Abi>;
  using size_type = Kokkos::Experimental::Impl::simd_size_t;

  constexpr size_type width = simd_type::size();
  T expected_val[width];
  for (size_type lane = 0; lane < width; ++lane) {
    expected_val[lane] = binary_op.on_device(first_arg[lane], second_arg);
  }

  simd_type expected_result =
      Kokkos::Experimental::simd_unchecked_load<simd_type>(
          expected_val, Kokkos::Experimental::simd_flag_default);
  simd_type const computed_result = binary_op.on_device(first_arg, second_arg);
  device_check_equality(expected_result, computed_result, nlanes);
}

template <class Abi, class BinaryOp, class T>
KOKKOS_INLINE_FUNCTION void device_check_arithmetic_op_scalar_vector(
    BinaryOp binary_op, T first_arg,
    Kokkos::Experimental::basic_simd<T, Abi> const& second_arg,
    Kokkos::Experimental::Impl::simd_size_t nlanes) {
  using simd_type = Kokkos::Experimental::basic_simd<T, Abi>;
  using size_type = Kokkos::Experimental::Impl::simd_size_t;

  constexpr size_type width = simd_type::size();
  T expected_val[width];
  for (size_type lane = 0; lane < width; ++lane) {
    expected_val[lane] = binary_op.on_device(first_arg, second_arg[lane]);
  }

  simd_type expected_result =
      Kokkos::Experimental::simd_unchecked_load<simd_type>(
          expected_val, Kokkos::Experimental::simd_flag_default);
  simd_type const computed_result = binary_op.on_device(first_arg, second_arg);
  device_check_equality(expected_result, computed_result, nlanes);
}

template <class Abi, class Loader, class BinaryOp, class T>
KOKKOS_INLINE_FUNCTION void device_check_arithmetic_op_one_loader(
    BinaryOp binary_op, Kokkos::Experimental::Impl::simd_size_t n,
    T const* first_args, T const* second_args) {
  Loader loader;
  using simd_type = Kokkos::Experimental::basic_simd<T, Abi>;
  using size_type = Kokkos::Experimental::Impl::simd_size_t;

  constexpr size_type width = simd_type::size();
  for (size_type i = 0; i < n; i += width) {
    const size_type nremaining = n - i;
    const size_type nlanes     = Kokkos::min(nremaining, width);
    if ((std::is_same_v<BinaryOp, divides> ||
         std::is_same_v<BinaryOp, divides_eq>)&&nremaining < width)
      continue;

    simd_type first_arg;
    bool const loaded_first_arg =
        loader.device_load(first_args + i, nlanes, first_arg);

    simd_type second_arg;
    bool const loaded_second_arg =
        loader.device_load(second_args + i, nlanes, second_arg);

    if (!(loaded_first_arg && loaded_second_arg)) continue;

    T first_scalar_arg  = first_args[i];
    T second_scalar_arg = second_args[i];
    device_check_arithmetic_op_vector_vector(binary_op, first_arg, second_arg,
                                             nlanes);
    device_check_arithmetic_op_vector_scalar(binary_op, first_arg,
                                             second_scalar_arg, nlanes);
    if constexpr (!(std::is_same_v<BinaryOp, plus_eq> ||
                    std::is_same_v<BinaryOp, minus_eq> ||
                    std::is_same_v<BinaryOp, multiplies_eq> ||
                    std::is_same_v<BinaryOp, divides_eq>)) {
      device_check_arithmetic_op_scalar_vector(binary_op, first_scalar_arg,
                                               second_arg, nlanes);
    }
  }
}

template <typename Abi, typename Op, typename... T>
KOKKOS_INLINE_FUNCTION void device_check_arithmetic_op_all_loaders(
    Op op, Kokkos::Experimental::Impl::simd_size_t n, T const*... args) {
  device_check_arithmetic_op_one_loader<Abi, load_element_aligned>(op, n,
                                                                   args...);
  device_check_arithmetic_op_one_loader<Abi, load_masked>(op, n, args...);
  device_check_arithmetic_op_one_loader<Abi, load_as_scalars>(op, n, args...);
  device_check_arithmetic_op_one_loader<Abi, load_vector_aligned>(op, n,
                                                                  args...);
}

template <typename Abi, typename DataType, size_t n>
KOKKOS_INLINE_FUNCTION void device_check_all_arithmetic_ops(
    const DataType (&first_args)[n], const DataType (&second_args)[n]) {
  device_check_arithmetic_op_all_loaders<Abi>(plus(), n, first_args,
                                              second_args);
  device_check_arithmetic_op_all_loaders<Abi>(plus_eq(), n, first_args,
                                              second_args);
  device_check_arithmetic_op_all_loaders<Abi>(minus(), n, first_args,
                                              second_args);
  device_check_arithmetic_op_all_loaders<Abi>(minus_eq(), n, first_args,
                                              second_args);
  device_check_arithmetic_op_all_loaders<Abi>(multiplies(), n, first_args,
                                              second_args);
  device_check_arithmetic_op_all_loaders<Abi>(multiplies_eq(), n, first_args,
                                              second_args);
  device_check_arithmetic_op_all_loaders<Abi>(divides(), n, first_args,
                                              second_args);
  device_check_arithmetic_op_all_loaders<Abi>(divides_eq(), n, first_args,
                                              second_args);
}

template <typename Abi, typename DataType>
KOKKOS_INLINE_FUNCTION void device_check_arithmetic_ops() {
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
      device_check_all_arithmetic_ops<Abi>(first_args, second_args);
    } else {
      if constexpr (std::is_signed_v<DataType>) {
        alignas(alignment) DataType const first_args[] = {
            1, 2, -1, 10, 0, 1, -2, 10, 0, 1, -2, -3, 7, 4, -9, -15};
        alignas(alignment) DataType const second_args[] = {
            1, 2, 1, 1, 1, -3, -2, 1, 13, -3, -2, 10, -15, 7, 2, -10};
        device_check_all_arithmetic_ops<Abi>(first_args, second_args);
      } else {
        alignas(alignment) DataType const first_args[] = {
            1, 2, 1, 10, 0, 1, 2, 10, 0, 1, 2, 11, 5, 8, 2, 14};
        alignas(alignment) DataType const second_args[] = {
            1, 2, 1, 1, 1, 3, 2, 1, 13, 3, 2, 3, 6, 20, 5, 14};
        device_check_all_arithmetic_ops<Abi>(first_args, second_args);
      }
    }
  }
}

template <typename Abi, typename... DataTypes>
KOKKOS_INLINE_FUNCTION void device_check_arithmetic_ops_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (device_check_arithmetic_ops<Abi, DataTypes>(), ...);
}

template <typename... Abis>
KOKKOS_INLINE_FUNCTION void device_check_arithmetic_ops_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (device_check_arithmetic_ops_all_types<Abis>(DataTypes()), ...);
}

class simd_device_arithmetic_ops_functor {
 public:
  KOKKOS_INLINE_FUNCTION void operator()(int) const {
    device_check_arithmetic_ops_all_abis(
        Kokkos::Experimental::Impl::device_abi_set());
  }
};

TEST(simd, host_arithmetic_ops) {
  host_check_arithmetic_ops_all_abis(
      Kokkos::Experimental::Impl::host_abi_set());
}

TEST(simd, device_arithmetic_ops) {
#if defined(KOKKOS_IMPL_32BIT)
  // FIXME_GCC: Floating point comparisons fail in a 32-bit build
  // See https://github.com/kokkos/kokkos/pull/7912#issuecomment-2747713935
  GTEST_SKIP()
      << "skipping due to a GCC bug associated with the computation of "
         "floating-point values in a 32-bit build.";
#endif
#if defined(KOKKOS_ENABLE_OPENACC) && \
    defined(KOKKOS_COMPILER_CLANG)  // FIXME_CLACC
  GTEST_SKIP()
      << "skipping because of a non-deterministic failure reporting: "
         "Failure to synchronize stream (nil): Error in "
         "cuStreamSynchronize: an illegal memory access was encountered";
#endif
  Kokkos::parallel_for(1, simd_device_arithmetic_ops_functor());
  Kokkos::fence();
}

#endif
