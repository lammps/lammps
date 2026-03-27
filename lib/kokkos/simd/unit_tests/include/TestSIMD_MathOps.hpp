// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_SIMD_MATH_OPS_HPP
#define KOKKOS_TEST_SIMD_MATH_OPS_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.simd;
import kokkos.simd_impl;
#else
#include <Kokkos_SIMD.hpp>
#endif
#include <SIMDTesting_Utilities.hpp>

template <class Abi, class Loader, class TernaryOp, class T>
void host_check_math_op_one_loader(TernaryOp ternary_op, std::size_t n,
                                   T const* first_args, T const* second_args,
                                   T const* third_args) {
  Loader loader;
  using simd_type             = Kokkos::Experimental::basic_simd<T, Abi>;
  constexpr std::size_t width = simd_type::size();
  for (std::size_t i = 0; i < n; i += width) {
    std::size_t const nremaining = n - i;
    std::size_t const nlanes     = Kokkos::min(nremaining, width);
    simd_type first_arg;
    bool const loaded_first_arg =
        loader.host_load(first_args + i, nlanes, first_arg);
    simd_type second_arg;
    bool const loaded_second_arg =
        loader.host_load(second_args + i, nlanes, second_arg);
    simd_type third_arg;
    bool const loaded_third_arg =
        loader.host_load(third_args + i, nlanes, third_arg);
    if (!(loaded_first_arg && loaded_second_arg && loaded_third_arg)) continue;

    T expected_val[width];
    for (std::size_t lane = 0; lane < width; ++lane) {
      expected_val[lane] = ternary_op.on_host(
          T(first_arg[lane]), T(second_arg[lane]), T(third_arg[lane]));
    }

    simd_type expected_result =
        Kokkos::Experimental::simd_unchecked_load<simd_type>(
            expected_val, Kokkos::Experimental::simd_flag_default);

    simd_type const computed_result =
        ternary_op.on_host(first_arg, second_arg, third_arg);
    host_check_equality(expected_result, computed_result, nlanes);
  }
}

template <class Abi, class Loader, class BinaryOp, class T>
void host_check_math_op_one_loader(BinaryOp binary_op, std::size_t n,
                                   T const* first_args, T const* second_args) {
  Loader loader;
  using simd_type             = Kokkos::Experimental::basic_simd<T, Abi>;
  constexpr std::size_t width = simd_type::size();
  for (std::size_t i = 0; i < n; i += width) {
    std::size_t const nremaining = n - i;
    std::size_t const nlanes     = Kokkos::min(nremaining, width);
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

    // The second argument of pow being negative and/or non-integer may provoke
    // a domain error
    if constexpr (std::is_same_v<BinaryOp, pow_op>) {
      second_arg = Kokkos::round(Kokkos::abs(second_arg));
    }

    T expected_val[width];
    for (std::size_t lane = 0; lane < width; ++lane) {
      expected_val[lane] =
          binary_op.on_host(T(first_arg[lane]), T(second_arg[lane]));
    }

    simd_type expected_result =
        Kokkos::Experimental::simd_unchecked_load<simd_type>(
            expected_val, Kokkos::Experimental::simd_flag_default);
    simd_type const computed_result = binary_op.on_host(first_arg, second_arg);
    host_check_equality(expected_result, computed_result, nlanes);
  }
}

template <class Abi, class Loader, class UnaryOp, class T>
void host_check_math_op_one_loader(UnaryOp unary_op, std::size_t n,
                                   T const* args) {
  Loader loader;
  using simd_type = Kokkos::Experimental::basic_simd<T, Abi>;

  constexpr std::size_t width = simd_type::size();
  for (std::size_t i = 0; i < n; i += width) {
    std::size_t const nremaining = n - i;
    std::size_t const nlanes     = Kokkos::min(nremaining, width);
    simd_type arg;
    bool const loaded_arg = loader.host_load(args + i, nlanes, arg);
    if (!loaded_arg) continue;

    if constexpr (std::is_same_v<UnaryOp, sqrt_op>) {
      arg = Kokkos::abs(arg);
    }

    if constexpr (std::is_same_v<UnaryOp, log_op> ||
                  std::is_same_v<UnaryOp, log10_op> ||
                  std::is_same_v<UnaryOp, log2_op> ||
                  std::is_same_v<UnaryOp, tgamma_op> ||
                  std::is_same_v<UnaryOp, lgamma_op>) {
      arg = Kokkos::abs(arg) + simd_type(0.1);
    }

    // These functions are defined for -1 < x < 1
    if constexpr (std::is_same_v<UnaryOp, asin_op> ||
                  std::is_same_v<UnaryOp, acos_op> ||
                  std::is_same_v<UnaryOp, atanh_op>) {
      arg /= simd_type(10.1);
    }

    // acosh is defined for x >= 1
    if constexpr (std::is_same_v<UnaryOp, acosh_op>) {
      arg = Kokkos::abs(arg) + simd_type(1.0);
    }

    auto unary_op_result   = unary_op.on_host(arg);
    using result_simd_type = decltype(unary_op_result);

    typename result_simd_type::value_type expected_val[width];
    for (std::size_t lane = 0; lane < width; ++lane) {
      expected_val[lane] = unary_op.on_host_serial(T(arg[lane]));
    }

    result_simd_type expected_result =
        Kokkos::Experimental::simd_unchecked_load<result_simd_type>(
            expected_val, Kokkos::Experimental::simd_flag_default);
    auto computed_result = unary_op.on_host(arg);
    host_check_equality(expected_result, computed_result, nlanes);
  }
}

template <class Abi, class Op, class... T>
inline void host_check_math_op_all_loaders(Op op, std::size_t n,
                                           T const*... args) {
  host_check_math_op_one_loader<Abi, load_element_aligned>(op, n, args...);
  host_check_math_op_one_loader<Abi, load_masked>(op, n, args...);
  host_check_math_op_one_loader<Abi, load_as_scalars>(op, n, args...);
  host_check_math_op_one_loader<Abi, load_vector_aligned>(op, n, args...);
}

template <typename Abi, typename DataType, size_t n>
inline void host_check_all_math_ops(const DataType (&first_args)[n],
                                    const DataType (&second_args)[n],
                                    const DataType (&third_args)[n]) {
  host_check_math_op_all_loaders<Abi>(plus(), n, first_args, second_args);
  host_check_math_op_all_loaders<Abi>(plus_eq(), n, first_args, second_args);
  host_check_math_op_all_loaders<Abi>(minus(), n, first_args, second_args);
  host_check_math_op_all_loaders<Abi>(minus_eq(), n, first_args, second_args);
  host_check_math_op_all_loaders<Abi>(multiplies(), n, first_args, second_args);
  host_check_math_op_all_loaders<Abi>(multiplies_eq(), n, first_args,
                                      second_args);
  host_check_math_op_all_loaders<Abi>(divides(), n, first_args, second_args);
  host_check_math_op_all_loaders<Abi>(divides_eq(), n, first_args, second_args);
  host_check_math_op_all_loaders<Abi>(absolutes(), n, first_args);

  host_check_math_op_all_loaders<Abi>(floors(), n, first_args);
  host_check_math_op_all_loaders<Abi>(ceils(), n, first_args);
  host_check_math_op_all_loaders<Abi>(rounds(), n, first_args);
  host_check_math_op_all_loaders<Abi>(truncates(), n, first_args);

  host_check_math_op_all_loaders<Abi>(minimum(), n, first_args, second_args);
  host_check_math_op_all_loaders<Abi>(maximum(), n, first_args, second_args);

  // TODO: Place fallback implementations for all simd integer types
  if constexpr (std::is_floating_point_v<DataType>) {
    host_check_math_op_all_loaders<Abi>(abs_op(), n, first_args);
    host_check_math_op_all_loaders<Abi>(exp_op(), n, first_args);
    host_check_math_op_all_loaders<Abi>(exp2_op(), n, first_args);
    host_check_math_op_all_loaders<Abi>(log_op(), n, first_args);
    host_check_math_op_all_loaders<Abi>(log10_op(), n, first_args);
    host_check_math_op_all_loaders<Abi>(log2_op(), n, first_args);
    host_check_math_op_all_loaders<Abi>(sqrt_op(), n, first_args);
    host_check_math_op_all_loaders<Abi>(cbrt_op(), n, first_args);
    host_check_math_op_all_loaders<Abi>(sin_op(), n, first_args);
    host_check_math_op_all_loaders<Abi>(cos_op(), n, first_args);
    host_check_math_op_all_loaders<Abi>(tan_op(), n, first_args);
    host_check_math_op_all_loaders<Abi>(asin_op(), n, first_args);
    host_check_math_op_all_loaders<Abi>(acos_op(), n, first_args);
    host_check_math_op_all_loaders<Abi>(atan_op(), n, first_args);
    host_check_math_op_all_loaders<Abi>(sinh_op(), n, first_args);
    host_check_math_op_all_loaders<Abi>(cosh_op(), n, first_args);
    host_check_math_op_all_loaders<Abi>(tanh_op(), n, first_args);
    host_check_math_op_all_loaders<Abi>(asinh_op(), n, first_args);
    host_check_math_op_all_loaders<Abi>(acosh_op(), n, first_args);
    host_check_math_op_all_loaders<Abi>(atanh_op(), n, first_args);
    host_check_math_op_all_loaders<Abi>(erf_op(), n, first_args);
    host_check_math_op_all_loaders<Abi>(erfc_op(), n, first_args);
    host_check_math_op_all_loaders<Abi>(tgamma_op(), n, first_args);
    host_check_math_op_all_loaders<Abi>(lgamma_op(), n, first_args);

    host_check_math_op_all_loaders<Abi>(pow_op(), n, first_args, second_args);
    host_check_math_op_all_loaders<Abi>(hypot_op(), n, first_args, second_args);
    host_check_math_op_all_loaders<Abi>(atan2_op(), n, first_args, second_args);
    host_check_math_op_all_loaders<Abi>(copysign_op(), n, first_args,
                                        second_args);

    host_check_math_op_all_loaders<Abi>(fma_op(), n, first_args, second_args,
                                        third_args);
    host_check_math_op_all_loaders<Abi>(ternary_hypot_op(), n, first_args,
                                        second_args, third_args);
  }
}

template <typename Abi, typename DataType>
inline void host_check_abi_size() {
  using simd_type = Kokkos::Experimental::basic_simd<DataType, Abi>;
  using mask_type = typename simd_type::mask_type;
  static_assert(simd_type::size() == mask_type::size());
}

template <typename Abi, typename DataType>
inline void host_check_math_ops() {
  if constexpr (is_simd_avail_v<DataType, Abi>) {
    constexpr size_t alignment =
        Kokkos::Experimental::basic_simd<DataType, Abi>::size() *
        sizeof(DataType);

    host_check_abi_size<Abi, DataType>();

    if constexpr (!std::is_integral_v<DataType>) {
      alignas(alignment) DataType const first_args[] = {
          0.1, 0.4, 0.5,  0.7, 1.0, 1.5,  -2.0, 10.0,
          0.0, 1.2, -2.8, 3.0, 4.0, -0.1, 5.0,  -0.2};
      alignas(alignment) DataType const second_args[] = {
          1.0,  0.2,  1.1,  1.8, -0.1,  -3.0, -2.4, 1.0,
          13.0, -3.2, -2.1, 3.0, -15.0, -0.5, -0.2, -0.2};
      alignas(alignment) DataType const third_args[] = {
          3.7,  2.6,  9.8, 9.3,  9.9, -5.3, 8.5,  1.6,
          -3.8, -0.4, 6.1, -5.3, 6.1, -8.9, -2.5, -5.2};
      host_check_all_math_ops<Abi>(first_args, second_args, third_args);
    } else {
      if constexpr (std::is_signed_v<DataType>) {
        alignas(alignment) DataType const first_args[] = {
            1, 2, -1, 10, 0, 1, -2, 10, 0, 1, -2, -3, 7, 4, -9, -15};
        alignas(alignment) DataType const second_args[] = {
            1, 2, 1, 1, 1, -3, -2, 1, 13, -3, -2, 10, -15, 7, 2, -10};
        alignas(alignment) DataType const third_args[] = {
            20, 3, -18, 3, 19, 11, 4, 20, 8, -8, 13, -18, -2, -5, -1, 11};
        host_check_all_math_ops<Abi>(first_args, second_args, third_args);
      } else {
        alignas(alignment) DataType const first_args[] = {
            1, 2, 1, 10, 0, 1, 2, 10, 0, 1, 2, 11, 5, 8, 2, 14};
        alignas(alignment) DataType const second_args[] = {
            1, 2, 1, 1, 1, 3, 2, 1, 13, 3, 2, 3, 6, 20, 5, 14};
        alignas(alignment) DataType const third_args[] = {
            10, 10, 6, 6, 16, 14, 18, 9, 19, 7, 0, 6, 2, 15, 10, 16};
        host_check_all_math_ops<Abi>(first_args, second_args, third_args);
      }
    }
  }
}

template <typename Abi, typename... DataTypes>
inline void host_check_math_ops_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (host_check_math_ops<Abi, DataTypes>(), ...);
}

template <typename... Abis>
inline void host_check_math_ops_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (host_check_math_ops_all_types<Abis>(DataTypes()), ...);
}

template <typename Abi, typename Loader, typename TernaryOp, typename T>
KOKKOS_INLINE_FUNCTION void device_check_math_op_one_loader(
    TernaryOp ternary_op, std::size_t n, T const* first_args,
    T const* second_args, T const* third_args) {
  Loader loader;
  using simd_type             = Kokkos::Experimental::basic_simd<T, Abi>;
  constexpr std::size_t width = simd_type::size();
  for (std::size_t i = 0; i < n; i += width) {
    std::size_t const nremaining = n - i;
    std::size_t const nlanes     = Kokkos::min(nremaining, width);
    simd_type first_arg;
    bool const loaded_first_arg =
        loader.device_load(first_args + i, nlanes, first_arg);
    simd_type second_arg;
    bool const loaded_second_arg =
        loader.device_load(second_args + i, nlanes, second_arg);
    simd_type third_arg;
    bool const loaded_third_arg =
        loader.device_load(third_args + i, nlanes, third_arg);
    if (!(loaded_first_arg && loaded_second_arg && loaded_third_arg)) continue;

    simd_type expected_result(KOKKOS_LAMBDA(std::size_t lane) {
      return ternary_op.on_device(first_arg[lane], second_arg[lane],
                                  third_arg[lane]);
    });

    simd_type const computed_result =
        ternary_op.on_device(first_arg, second_arg, third_arg);
    device_check_equality(expected_result, computed_result, nlanes);
  }
}

template <typename Abi, typename Loader, typename BinaryOp, typename T>
KOKKOS_INLINE_FUNCTION void device_check_math_op_one_loader(
    BinaryOp binary_op, std::size_t n, T const* first_args,
    T const* second_args) {
  Loader loader;
  using simd_type             = Kokkos::Experimental::basic_simd<T, Abi>;
  constexpr std::size_t width = simd_type::size();
  for (std::size_t i = 0; i < n; i += width) {
    std::size_t const nremaining = n - i;
    std::size_t const nlanes     = Kokkos::min(nremaining, width);
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

    // The second argument of pow being negative and/or non-integer may provoke
    // a domain error
    if constexpr (std::is_same_v<BinaryOp, pow_op>) {
      second_arg = Kokkos::round(Kokkos::abs(second_arg));
    }

    simd_type expected_result(KOKKOS_LAMBDA(std::size_t lane) {
      return binary_op.on_device(first_arg[lane], second_arg[lane]);
    });

    simd_type const computed_result =
        binary_op.on_device(first_arg, second_arg);
    device_check_equality(expected_result, computed_result, nlanes);
  }
}

template <typename Abi, typename Loader, typename UnaryOp, typename T>
KOKKOS_INLINE_FUNCTION void device_check_math_op_one_loader(UnaryOp unary_op,
                                                            std::size_t n,
                                                            T const* args) {
  Loader loader;
  using simd_type             = Kokkos::Experimental::basic_simd<T, Abi>;
  constexpr std::size_t width = simd_type::size();
  for (std::size_t i = 0; i < n; i += width) {
    std::size_t const nremaining = n - i;
    std::size_t const nlanes     = Kokkos::min(nremaining, width);
    simd_type arg;
    bool const loaded_arg = loader.device_load(args + i, nlanes, arg);
    if (!loaded_arg) continue;

    if constexpr (std::is_same_v<UnaryOp, sqrt_op>) {
      arg = Kokkos::abs(arg);
    }

    if constexpr (std::is_same_v<UnaryOp, log_op> ||
                  std::is_same_v<UnaryOp, log10_op> ||
                  std::is_same_v<UnaryOp, log2_op> ||
                  std::is_same_v<UnaryOp, tgamma_op> ||
                  std::is_same_v<UnaryOp, lgamma_op>) {
      arg = Kokkos::abs(arg) + simd_type(0.1);
    }

    // These functions are defined for -1 < x < 1
    if constexpr (std::is_same_v<UnaryOp, asin_op> ||
                  std::is_same_v<UnaryOp, acos_op> ||
                  std::is_same_v<UnaryOp, atanh_op>) {
      arg /= simd_type(10.1);
    }

    // acosh is defined for x >= 1
    if constexpr (std::is_same_v<UnaryOp, acosh_op>) {
      arg = Kokkos::abs(arg) + simd_type(1.0);
    }

    auto computed_result = unary_op.on_device(arg);

    decltype(computed_result) expected_result(KOKKOS_LAMBDA(std::size_t lane) {
      return unary_op.on_device_serial(arg[lane]);
    });

    device_check_equality(expected_result, computed_result, nlanes);
  }
}

template <typename Abi, typename Op, typename... T>
KOKKOS_INLINE_FUNCTION void device_check_math_op_all_loaders(Op op,
                                                             std::size_t n,
                                                             T const*... args) {
  device_check_math_op_one_loader<Abi, load_element_aligned>(op, n, args...);
  device_check_math_op_one_loader<Abi, load_masked>(op, n, args...);
  device_check_math_op_one_loader<Abi, load_as_scalars>(op, n, args...);
  device_check_math_op_one_loader<Abi, load_vector_aligned>(op, n, args...);
}

template <typename Abi, typename DataType, size_t n>
KOKKOS_INLINE_FUNCTION void device_check_all_math_ops(
    const DataType (&first_args)[n], const DataType (&second_args)[n],
    const DataType (&third_args)[n]) {
  device_check_math_op_all_loaders<Abi>(plus(), n, first_args, second_args);
  device_check_math_op_all_loaders<Abi>(plus_eq(), n, first_args, second_args);
  device_check_math_op_all_loaders<Abi>(minus(), n, first_args, second_args);
  device_check_math_op_all_loaders<Abi>(minus_eq(), n, first_args, second_args);
  device_check_math_op_all_loaders<Abi>(multiplies(), n, first_args,
                                        second_args);
  device_check_math_op_all_loaders<Abi>(multiplies_eq(), n, first_args,
                                        second_args);
  device_check_math_op_all_loaders<Abi>(divides(), n, first_args, second_args);
  device_check_math_op_all_loaders<Abi>(divides_eq(), n, first_args,
                                        second_args);
  device_check_math_op_all_loaders<Abi>(absolutes(), n, first_args);

  device_check_math_op_all_loaders<Abi>(floors(), n, first_args);
  device_check_math_op_all_loaders<Abi>(ceils(), n, first_args);
  device_check_math_op_all_loaders<Abi>(rounds(), n, first_args);
  device_check_math_op_all_loaders<Abi>(truncates(), n, first_args);

  device_check_math_op_all_loaders<Abi>(minimum(), n, first_args, second_args);
  device_check_math_op_all_loaders<Abi>(maximum(), n, first_args, second_args);

  if constexpr (std::is_floating_point_v<DataType>) {
    device_check_math_op_all_loaders<Abi>(abs_op(), n, first_args);
    device_check_math_op_all_loaders<Abi>(exp_op(), n, first_args);
    device_check_math_op_all_loaders<Abi>(exp2_op(), n, first_args);
    device_check_math_op_all_loaders<Abi>(log_op(), n, first_args);
    device_check_math_op_all_loaders<Abi>(log10_op(), n, first_args);
    device_check_math_op_all_loaders<Abi>(log2_op(), n, first_args);
    device_check_math_op_all_loaders<Abi>(sqrt_op(), n, first_args);
    device_check_math_op_all_loaders<Abi>(cbrt_op(), n, first_args);
    device_check_math_op_all_loaders<Abi>(sin_op(), n, first_args);
    device_check_math_op_all_loaders<Abi>(cos_op(), n, first_args);
    device_check_math_op_all_loaders<Abi>(tan_op(), n, first_args);
    device_check_math_op_all_loaders<Abi>(asin_op(), n, first_args);
    device_check_math_op_all_loaders<Abi>(acos_op(), n, first_args);
    device_check_math_op_all_loaders<Abi>(atan_op(), n, first_args);
    device_check_math_op_all_loaders<Abi>(sinh_op(), n, first_args);
    device_check_math_op_all_loaders<Abi>(cosh_op(), n, first_args);
    device_check_math_op_all_loaders<Abi>(tanh_op(), n, first_args);
    device_check_math_op_all_loaders<Abi>(asinh_op(), n, first_args);
    device_check_math_op_all_loaders<Abi>(acosh_op(), n, first_args);
    device_check_math_op_all_loaders<Abi>(atanh_op(), n, first_args);
    device_check_math_op_all_loaders<Abi>(erf_op(), n, first_args);
    device_check_math_op_all_loaders<Abi>(erfc_op(), n, first_args);
    device_check_math_op_all_loaders<Abi>(tgamma_op(), n, first_args);
    device_check_math_op_all_loaders<Abi>(lgamma_op(), n, first_args);

    device_check_math_op_all_loaders<Abi>(pow_op(), n, first_args, second_args);
    device_check_math_op_all_loaders<Abi>(hypot_op(), n, first_args,
                                          second_args);
    device_check_math_op_all_loaders<Abi>(atan2_op(), n, first_args,
                                          second_args);
    device_check_math_op_all_loaders<Abi>(copysign_op(), n, first_args,
                                          second_args);

    device_check_math_op_all_loaders<Abi>(fma_op(), n, first_args, second_args,
                                          third_args);
    device_check_math_op_all_loaders<Abi>(ternary_hypot_op(), n, first_args,
                                          second_args, third_args);
  }
}

template <typename Abi, typename DataType>
KOKKOS_INLINE_FUNCTION void device_check_abi_size() {
  using simd_type = Kokkos::Experimental::basic_simd<DataType, Abi>;
  using mask_type = typename simd_type::mask_type;
  static_assert(simd_type::size() == mask_type::size());
}

template <typename Abi, typename DataType>
KOKKOS_INLINE_FUNCTION void device_check_math_ops() {
  if constexpr (is_type_v<Kokkos::Experimental::basic_simd<DataType, Abi>>) {
    device_check_abi_size<Abi, DataType>();

    if constexpr (!std::is_integral_v<DataType>) {
      DataType const first_args[]  = {0.1,  0.4,  0.5, 0.7, 1.0,  1.5,
                                      -2.0, 10.0, 0.0, 1.2, -2.8, 3.0,
                                      4.0,  -0.1, 5.0, -0.2};
      DataType const second_args[] = {1.0,   0.2,  1.1,  1.8,  -0.1, -3.0,
                                      -2.4,  1.0,  13.0, -3.2, -2.1, 3.0,
                                      -15.0, -0.5, -0.2, -0.2};
      DataType const third_args[]  = {3.7, 2.6,  9.8,  9.3,  9.9, -5.3,
                                      8.5, 1.6,  -3.8, -0.4, 6.1, -5.3,
                                      6.1, -8.9, -2.5, -5.2};
      device_check_all_math_ops<Abi>(first_args, second_args, third_args);
    } else {
      if constexpr (std::is_signed_v<DataType>) {
        DataType const first_args[]  = {1, 2, -1, 10, 0, 1, -2, 10,
                                        0, 1, -2, -3, 7, 4, -9, -15};
        DataType const second_args[] = {1,  2,  1,  1,  1,   -3, -2, 1,
                                        13, -3, -2, 10, -15, 7,  2,  -10};
        DataType const third_args[]  = {20, 3,  -18, 3,   19, 11, 4,  20,
                                        8,  -8, 13,  -18, -2, -5, -1, 11};
        device_check_all_math_ops<Abi>(first_args, second_args, third_args);
      } else {
        DataType const first_args[]  = {1, 2, 1, 10, 0, 1, 2, 10,
                                        0, 1, 2, 11, 5, 8, 2, 14};
        DataType const second_args[] = {1,  2, 1, 1, 1, 3,  2, 1,
                                        13, 3, 2, 3, 6, 20, 5, 14};
        DataType const third_args[]  = {10, 10, 6, 6, 16, 14, 18, 9,
                                        19, 7,  0, 6, 2,  15, 10, 16};
        device_check_all_math_ops<Abi>(first_args, second_args, third_args);
      }
    }
  }
}

template <typename Abi, typename... DataTypes>
KOKKOS_INLINE_FUNCTION void device_check_math_ops_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (device_check_math_ops<Abi, DataTypes>(), ...);
}

template <typename... Abis>
KOKKOS_INLINE_FUNCTION void device_check_math_ops_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (device_check_math_ops_all_types<Abis>(DataTypes()), ...);
}

class simd_device_math_ops_functor {
 public:
  KOKKOS_INLINE_FUNCTION void operator()(int) const {
    device_check_math_ops_all_abis(
        Kokkos::Experimental::Impl::device_abi_set());
  }
};

TEST(simd, host_math_ops) {
  host_check_math_ops_all_abis(Kokkos::Experimental::Impl::host_abi_set());
}

TEST(simd, device_math_ops) {
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
  Kokkos::parallel_for(1, simd_device_math_ops_functor());
}

#endif
