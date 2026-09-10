// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_SIMD_PERFTEST_HOST_HPP
#define KOKKOS_SIMD_PERFTEST_HOST_HPP

#include <benchmark/benchmark.h>
#include <Kokkos_SIMD.hpp>
#include <cstdlib>

#include "Common.hpp"
#include "PerfTest_Operators.hpp"

template <class Abi, class UnaryOp, class T, class ExecSpace>
void host_bench_unary_op(benchmark::State& state,
                         const ArgsWrapper<T, ExecSpace>& wrapper) {
  constexpr bool force_serial = std::is_same_v<Abi, simd_abi_force_serial>;
  using RealAbi =
      std::conditional_t<force_serial, Kokkos::Experimental::simd_abi::scalar,
                         Abi>;
  using simd_type             = Kokkos::Experimental::basic_simd<T, RealAbi>;
  constexpr std::size_t width = simd_type::size();

  UnaryOp op;

  View<T*, ExecSpace> res("res", BENCH_SIZE);

  for (auto _ : state) {
    simd_type a, x;
    for (std::size_t i = 0; i < BENCH_SIZE; i += width) {
      a = simd_unchecked_load<simd_type>(
          wrapper.args->arg1.data() + i,
          Kokkos::Experimental::simd_flag_aligned);
      x = op.on_host(a);
      simd_unchecked_store(x, res.data() + i,
                           Kokkos::Experimental::simd_flag_aligned);
      if constexpr (force_serial) {
        benchmark::DoNotOptimize(x);
      }
    }
  }
}

template <class Abi, class BinaryOp, class T, class ExecSpace>
void host_bench_binary_op(benchmark::State& state,
                          const ArgsWrapper<T, ExecSpace>& wrapper) {
  constexpr bool force_serial = std::is_same_v<Abi, simd_abi_force_serial>;
  using RealAbi =
      std::conditional_t<force_serial, Kokkos::Experimental::simd_abi::scalar,
                         Abi>;
  using simd_type             = Kokkos::Experimental::basic_simd<T, RealAbi>;
  constexpr std::size_t width = simd_type::size();

  BinaryOp op;

  View<T*, ExecSpace> res("res", BENCH_SIZE);

  for (auto _ : state) {
    simd_type a, b, x;
    for (std::size_t i = 0; i < BENCH_SIZE; i += width) {
      a = simd_unchecked_load<simd_type>(
          wrapper.args->arg1.data() + i,
          Kokkos::Experimental::simd_flag_aligned);
      b = simd_unchecked_load<simd_type>(
          wrapper.args->arg2.data() + i,
          Kokkos::Experimental::simd_flag_aligned);
      x = op.on_host(a, b);
      simd_unchecked_store(x, res.data() + i,
                           Kokkos::Experimental::simd_flag_aligned);
      if constexpr (force_serial) {
        benchmark::DoNotOptimize(x);
      }
    }
  }
}

template <class Abi, class TernaryOp, class T, class ExecSpace>
void host_bench_ternary_op(benchmark::State& state,
                           const ArgsWrapper<T, ExecSpace>& wrapper) {
  constexpr bool force_serial = std::is_same_v<Abi, simd_abi_force_serial>;
  using RealAbi =
      std::conditional_t<force_serial, Kokkos::Experimental::simd_abi::scalar,
                         Abi>;
  using simd_type             = Kokkos::Experimental::basic_simd<T, RealAbi>;
  constexpr std::size_t width = simd_type::size();

  TernaryOp op;

  View<T*, ExecSpace> res("res", BENCH_SIZE);

  for (auto _ : state) {
    simd_type a, b, c, x;
    for (std::size_t i = 0; i < BENCH_SIZE; i += width) {
      a = simd_unchecked_load<simd_type>(
          wrapper.args->arg1.data() + i,
          Kokkos::Experimental::simd_flag_aligned);
      b = simd_unchecked_load<simd_type>(
          wrapper.args->arg2.data() + i,
          Kokkos::Experimental::simd_flag_aligned);
      c = simd_unchecked_load<simd_type>(
          wrapper.args->arg3.data() + i,
          Kokkos::Experimental::simd_flag_aligned);
      x = op.on_host(a, b, c);
      simd_unchecked_store(x, res.data() + i,
                           Kokkos::Experimental::simd_flag_aligned);
      if constexpr (force_serial) {
        benchmark::DoNotOptimize(x);
      }
    }
  }
}

template <class Abi, class ReductionOp, class T, class ExecSpace>
void host_bench_reduction_op(benchmark::State& state,
                             const ArgsWrapper<T, ExecSpace>& wrapper) {
  constexpr bool force_serial = std::is_same_v<Abi, simd_abi_force_serial>;
  using RealAbi =
      std::conditional_t<force_serial, Kokkos::Experimental::simd_abi::scalar,
                         Abi>;
  using simd_type             = Kokkos::Experimental::basic_simd<T, RealAbi>;
  using mask_type             = typename simd_type::mask_type;
  constexpr std::size_t width = simd_type::size();

  ReductionOp op;

  std::array<mask_type, BENCH_SIZE / width> masks;
  std::srand(58051);

  for (std::size_t i = 0; i < masks.size(); i++) {
    if constexpr (std::is_same_v<RealAbi,
                                 Kokkos::Experimental::simd_abi::scalar>) {
      masks[i] = mask_type(std::rand() % 2 == 0);
    } else {
      masks[i] = mask_type([=](std::size_t) { return std::rand() % 2 == 0; });
    }
  }

  View<typename simd_type::value_type*, ExecSpace> res("res",
                                                       BENCH_SIZE / width);

  for (auto _ : state) {
    simd_type a;
    for (std::size_t i = 0; i < BENCH_SIZE; i += width) {
      a = simd_unchecked_load<simd_type>(
          wrapper.args->arg1.data() + i,
          Kokkos::Experimental::simd_flag_aligned);
      res(i / width) = op.on_host(a, masks[i / width]);
      if constexpr (force_serial) {
        benchmark::DoNotOptimize(res(i / width));
      }
    }
  }
}

#define KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(prefix, name, op) \
  benchmark::RegisterBenchmark(                                      \
      benchmark_name<Abi, DataType>("host " #prefix, #name).data(),  \
      host_bench_unary_op<Abi, op, DataType, ExecSpace>, wrapper)    \
      ->Iterations(3)
#define KOKKOS_IMPL_SIMD_PERFTEST_HOST_BINARY_BENCH(prefix, name, op) \
  benchmark::RegisterBenchmark(                                       \
      benchmark_name<Abi, DataType>("host " #prefix, #name).data(),   \
      host_bench_binary_op<Abi, op, DataType, ExecSpace>, wrapper)    \
      ->Iterations(3)
#define KOKKOS_IMPL_SIMD_PERFTEST_HOST_TERNARY_BENCH(prefix, name, op) \
  benchmark::RegisterBenchmark(                                        \
      benchmark_name<Abi, DataType>("host " #prefix, #name).data(),    \
      host_bench_ternary_op<Abi, op, DataType, ExecSpace>, wrapper)    \
      ->Iterations(3)
#define KOKKOS_IMPL_SIMD_PERFTEST_HOST_REDUCTION_BENCH(prefix, name, op) \
  benchmark::RegisterBenchmark(                                          \
      benchmark_name<Abi, DataType>("host " #prefix, #name).data(),      \
      host_bench_reduction_op<Abi, op, DataType, ExecSpace>, wrapper)    \
      ->Iterations(3)

template <typename Abi, typename DataType>
inline void host_register_common_benchmarks() {
  using ExecSpace = Kokkos::DefaultHostExecutionSpace;

  const ArgsWrapper<DataType, ExecSpace> wrapper(BENCH_SIZE);

  if constexpr (is_simd_type_v<DataType, Abi>) {
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_BINARY_BENCH(common, add, plus);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_BINARY_BENCH(common, sub, minus);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_BINARY_BENCH(common, multiply, multiplies);

    if constexpr (std::is_floating_point_v<DataType>) {
      KOKKOS_IMPL_SIMD_PERFTEST_HOST_BINARY_BENCH(common, divide, divides);
      KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(common, floor, floors);
      KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(common, ceil, ceils);
      KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(common, round, rounds);
      KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(common, truncate, truncates);
    }

    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(common, abs, absolutes);

    if constexpr (std::is_integral_v<DataType>) {
      KOKKOS_IMPL_SIMD_PERFTEST_HOST_BINARY_BENCH(common, shift_right,
                                                  shift_right);
      KOKKOS_IMPL_SIMD_PERFTEST_HOST_BINARY_BENCH(common, shift_left,
                                                  shift_left);
    }

    KOKKOS_IMPL_SIMD_PERFTEST_HOST_BINARY_BENCH(common, min, minimum);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_BINARY_BENCH(common, max, maximum);

    KOKKOS_IMPL_SIMD_PERFTEST_HOST_REDUCTION_BENCH(common, reduce, reduce);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_REDUCTION_BENCH(common, reduce_min,
                                                   reduce_min);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_REDUCTION_BENCH(common, reduce_max,
                                                   reduce_max);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_REDUCTION_BENCH(common, masked_reduce,
                                                   masked_reduce);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_REDUCTION_BENCH(common, masked_reduce_min,
                                                   masked_reduce_min);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_REDUCTION_BENCH(common, masked_reduce_max,
                                                   masked_reduce_max);
  }
}

template <typename Abi, typename DataType>
inline void host_register_math_benchmarks() {
  using ExecSpace = Kokkos::DefaultHostExecutionSpace;

  const ArgsWrapper<DataType, ExecSpace> wrapper(BENCH_SIZE);

  if constexpr (std::is_floating_point_v<DataType> &&
                is_simd_type_v<DataType, Abi>) {
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(math, exp, exp_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(math, exp2, exp2_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(math, log, log_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(math, log10, log10_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(math, log2, log2_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(math, sqrt, sqrt_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(math, cbrt, cbrt_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(math, sin, sin_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(math, cos, cos_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(math, tan, tan_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(math, asin, asin_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(math, acos, acos_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(math, atan, atan_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(math, sinh, sinh_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(math, cosh, cosh_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(math, tanh, tanh_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(math, asinh, asinh_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(math, acosh, acosh_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(math, atanh, atanh_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(math, erf, erf_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(math, erfc, erfc_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(math, tgamma, tgamma_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_UNARY_BENCH(math, lgamma, lgamma_op);

    KOKKOS_IMPL_SIMD_PERFTEST_HOST_BINARY_BENCH(math, pow, pow_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_BINARY_BENCH(math, hypot, hypot_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_BINARY_BENCH(math, atan2, atan2_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_BINARY_BENCH(math, copysign, copysign_op);

    KOKKOS_IMPL_SIMD_PERFTEST_HOST_TERNARY_BENCH(math, fma, fma_op);
    KOKKOS_IMPL_SIMD_PERFTEST_HOST_TERNARY_BENCH(math, ternary_hypot,
                                                 ternary_hypot_op);
  }
}

template <typename Abi, typename... DataTypes>
inline void host_register_benchmarks_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (host_register_common_benchmarks<Abi, DataTypes>(), ...);
  (host_register_math_benchmarks<Abi, DataTypes>(), ...);
}

template <typename... Abis>
inline void host_register_benchmarks_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
#if defined(KOKKOS_ENABLE_BENCHMARKS_HEAVY)
  host_register_benchmarks_all_types<simd_abi_force_serial>(DataTypes());
#endif
  (host_register_benchmarks_all_types<Abis>(DataTypes()), ...);
}

inline void register_host_benchmarks() {
  host_register_benchmarks_all_abis(Kokkos::Experimental::Impl::host_abi_set());
}

#endif
