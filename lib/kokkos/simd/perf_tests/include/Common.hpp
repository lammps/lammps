// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_SIMD_PERFTEST_COMMON_HPP
#define KOKKOS_SIMD_PERFTEST_COMMON_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_SIMD.hpp>
#include <Kokkos_Random.hpp>
#include <cstdlib>
#include <benchmark/benchmark.h>

template <typename T, typename ExecSpace>
using View =
    Kokkos::View<T, ExecSpace,
                 Kokkos::MemoryTraits<Kokkos::Aligned | Kokkos::Restrict>>;

// Class wrapping the Arguments we use in the benchmarks.
// This avoids recreating the Views between each benches.
template <class T, class ExecSpace>
struct Args {
  using view_type = View<T*, ExecSpace>;

  view_type arg1;
  view_type arg2;
  view_type arg3;

  explicit Args(std::size_t size)
      : arg1("arg1", size), arg2("arg2", size), arg3("arg3", size) {
    Kokkos::Random_XorShift64_Pool<ExecSpace> random_pool(58051);

    Kokkos::fill_random(arg1, random_pool, 1, 100);
    Kokkos::fill_random(arg2, random_pool, 1, 15);
    Kokkos::fill_random(arg3, random_pool, 1, 100);
  }
};

// Class used to wrap the Args class, allowing us to destroy the Kokkos::Views
// inside it before Kokkos::finalize is called. Otherwise, Google Benchmark
// keeps a handle on Args and only destroys it after the call to
// Kokkos::finalize, causing a Kokkos::abort during main cleanup.
template <class T, class ExecSpace>
struct ArgsWrapper {
  Args<T, ExecSpace>* args;

  explicit ArgsWrapper(std::size_t size) {
    args = new Args<T, ExecSpace>(size);
    Kokkos::push_finalize_hook([ptr = args] { delete ptr; });
  }
};

// Additional abi used to bench the scalar abi without auto-vectorization.
class simd_abi_force_serial;

template <class Abi>
constexpr const char* abi_name() {
  using namespace Kokkos::Experimental::simd_abi;
  if constexpr (std::is_same_v<Abi, scalar>) {
    return "scalar";
  } else if constexpr (std::is_same_v<Abi, simd_abi_force_serial>) {
    return "scalar no-vec";
#if defined(KOKKOS_ARCH_AVX2)
  } else if constexpr (std::is_same_v<Abi, avx2_fixed_size<4>>) {
    return "avx2 x4";
  } else if constexpr (std::is_same_v<Abi, avx2_fixed_size<8>>) {
    return "avx2 x8";
#elif defined(KOKKOS_ARCH_AVX512XEON)
  } else if constexpr (std::is_same_v<Abi, avx512_fixed_size<8>>) {
    return "avx512 x8";
  } else if constexpr (std::is_same_v<Abi, avx512_fixed_size<16>>) {
    return "avx512 x16";
#elif defined(KOKKOS_ARCH_ARM_NEON)
  } else if constexpr (std::is_same_v<Abi, neon_fixed_size<2>>) {
    return "neon x2";
  } else if constexpr (std::is_same_v<Abi, neon_fixed_size<4>>) {
    return "neon x4";
#elif defined(KOKKOS_ARCH_ARM_SVE)
  } else if constexpr (std::is_same_v<Abi,
                                      sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>) {
    return "sve x" KOKKOS_IMPL_TOSTRING(SVE_DOUBLES_IN_VECTOR);
  } else if constexpr (std::is_same_v<Abi,
                                      sve_fixed_size<SVE_WORDS_IN_VECTOR>>) {
    return "sve x" KOKKOS_IMPL_TOSTRING(SVE_WORDS_IN_VECTOR);
#endif
  }
}

template <class DataType>
constexpr const char* datatype_name() {
  if constexpr (std::is_same_v<DataType, std::int32_t>) {
    return "i32";
  } else if constexpr (std::is_same_v<DataType, std::uint32_t>) {
    return "u32";
  } else if constexpr (std::is_same_v<DataType, std::int64_t>) {
    return "i64";
  } else if constexpr (std::is_same_v<DataType, std::uint64_t>) {
    return "u64";
  } else if constexpr (std::is_same_v<DataType, float>) {
    return "float";
  } else if constexpr (std::is_same_v<DataType, double>) {
    return "double";
  }
}

template <typename Abi, typename DataType>
std::string benchmark_name(const char* prefix, const char* name) {
  const std::string sp(" ");
  return prefix + sp + abi_name<Abi>() + sp + datatype_name<DataType>() + sp +
         name;
}

// This should be a multiple of the max vector size we can encounter, which is
// currently 16 (avx512 + 32bit datatype)
constexpr std::size_t BENCH_SIZE = 50'000;

// TODO: is_type_v and is_simd_avail_v are duplicated from
// SIMDTesting_Utilities.hpp, we should probably find a cleaner solution in the
// future
template <typename T, typename = void>
constexpr bool is_type_v = false;

template <typename T>
constexpr bool is_type_v<T, decltype(void(sizeof(T)))> = true;

template <typename DataType, typename Abi>
constexpr bool is_simd_avail_v =
    is_type_v<Kokkos::Experimental::basic_simd<DataType, Abi>> &&
    is_type_v<Kokkos::Experimental::basic_simd_mask<DataType, Abi>>;

template <typename T, typename Abi>
constexpr bool is_simd_type_v = is_simd_avail_v<
    T, std::conditional_t<std::is_same_v<Abi, simd_abi_force_serial>,
                          Kokkos::Experimental::simd_abi::scalar, Abi>>;

#endif
