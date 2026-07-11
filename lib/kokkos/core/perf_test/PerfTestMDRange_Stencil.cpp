// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <functional>
#include <iostream>
#include <numeric>
#include <limits>

#include <benchmark/benchmark.h>

#include "PerfTest_Category.hpp"
#include <Kokkos_Core.hpp>

namespace Test {

template <typename Layout>
struct LayoutToIterationPattern {};

template <>
struct LayoutToIterationPattern<Kokkos::LayoutRight> {
  static constexpr Kokkos::Iterate pattern = Kokkos::Iterate::Right;
};

template <>
struct LayoutToIterationPattern<Kokkos::LayoutLeft> {
  static constexpr Kokkos::Iterate pattern = Kokkos::Iterate::Left;
};

template <typename ScalarType, typename ViewType>
void check_computation(const ViewType& A, const ViewType& B) {
  int numErrors = 0;
  auto Ahost    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
  auto Bhost    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), B);

  // On KNL, this may vectorize - add print statement to prevent
  // Also, compare against epsilon, as vectorization can change bitwise
  // answer
  ScalarType epsilon = std::numeric_limits<ScalarType>::epsilon() * 100;
  if constexpr (ViewType::rank == 2) {
    const int n0 = Ahost.extent_int(0) - 2, n1 = Ahost.extent_int(1) - 2;
    for (int i0 = 1; i0 < n0 + 1; ++i0) {
      for (int i1 = 1; i1 < n1 + 1; ++i1) {
        ScalarType check =
            0.25 *
            (ScalarType)(Bhost(i0 + 1, i1) + Bhost(i0 - 1, i1) +
                         Bhost(i0, i1 + 1) + Bhost(i0, i1 - 1) + Bhost(i0, i1));
        if (Kokkos::abs(Ahost(i0, i1) - check) > epsilon) {
          ++numErrors;
          std::cerr << "Correctness error at index: " << i0 << "," << i1
                    << ", got " << Ahost(i0, i1) << ", expected " << check
                    << "\n";
        }
      }
    }
  } else if constexpr (ViewType::rank == 3) {
    const int n0 = Ahost.extent_int(0) - 2, n1 = Ahost.extent_int(1) - 2,
              n2 = Ahost.extent_int(2) - 2;
    for (int i0 = 1; i0 < n0 + 1; ++i0) {
      for (int i1 = 1; i1 < n1 + 1; ++i1) {
        for (int i2 = 1; i2 < n2 + 1; ++i2) {
          ScalarType check =
              0.25 *
              (ScalarType)(Bhost(i0 + 1, i1, i2) + Bhost(i0 - 1, i1, i2) +
                           Bhost(i0, i1 + 1, i2) + Bhost(i0, i1 - 1, i2) +
                           Bhost(i0, i1, i2 + 1) + Bhost(i0, i1, i2 - 1) +
                           Bhost(i0, i1, i2));
          if (Kokkos::abs(Ahost(i0, i1, i2) - check) > epsilon) {
            ++numErrors;
            std::cerr << "Correctness error at index: " << i0 << "," << i1
                      << "," << i2 << ", got " << Ahost(i0, i1, i2)
                      << ", expected " << check << "\n";
          }
        }
      }
    }
  } else if constexpr (ViewType::rank == 4) {
    const int n0 = Ahost.extent_int(0) - 2, n1 = Ahost.extent_int(1) - 2,
              n2 = Ahost.extent_int(2) - 2, n3 = Ahost.extent_int(3) - 2;
    for (int i0 = 1; i0 < n0 + 1; ++i0) {
      for (int i1 = 1; i1 < n1 + 1; ++i1) {
        for (int i2 = 1; i2 < n2 + 1; ++i2) {
          for (int i3 = 1; i3 < n3 + 1; ++i3) {
            ScalarType check = 0.25 * (ScalarType)(Bhost(i0 + 1, i1, i2, i3) +
                                                   Bhost(i0 - 1, i1, i2, i3) +
                                                   Bhost(i0, i1 + 1, i2, i3) +
                                                   Bhost(i0, i1 - 1, i2, i3) +
                                                   Bhost(i0, i1, i2 + 1, i3) +
                                                   Bhost(i0, i1, i2 - 1, i3) +
                                                   Bhost(i0, i1, i2, i3 + 1) +
                                                   Bhost(i0, i1, i2, i3 - 1) +
                                                   Bhost(i0, i1, i2, i3));
            if (Kokkos::abs(Ahost(i0, i1, i2, i3) - check) > epsilon) {
              ++numErrors;
              std::cerr << "Correctness error at index: " << i0 << "," << i1
                        << "," << i2 << "," << i3 << ", got "
                        << Ahost(i0, i1, i2, i3) << ", expected " << check
                        << "\n";
            }
          }
        }
      }
    }
    if (numErrors != 0) {
      std::cerr << "Detected some errors for a run with dimensions "
                << Ahost.extent(0);
      for (std::size_t i = 1; i < Ahost.rank(); i++) {
        std::cerr << "x" << Ahost.extent(i);
      }
      std::cerr << std::endl;
    }
  }
}

template <typename FunctorType, std::size_t... Idx>
void bench_mdrange(benchmark::State& state, std::index_sequence<Idx...>) {
  using execution_space = typename FunctorType::execution_space;
  using view_type       = typename FunctorType::view_type;

  Kokkos::Array<int, FunctorType::dimension> dims, tiles;
  for (std::size_t i = 0; i < dims.size(); i++) {
    dims[i]  = state.range(0);
    tiles[i] = state.range(1);
  }

  const auto policy = FunctorType::get_policy(dims, tiles);

  bool using_default_tiling = false;
  for (std::size_t i = 0; i < tiles.size(); i++) {
    state.counters[std::string("tile_") + std::to_string(i)] = tiles[i];
    using_default_tiling |= tiles[i] != state.range(1);
  }
  state.counters["default_tiling"] = using_default_tiling;

  view_type Atest("Atest", (dims[Idx] + 2)...);
  view_type Btest("Btest", (dims[Idx] + 2)...);

  Kokkos::deep_copy(Atest, 1.0);
  execution_space().fence();
  Kokkos::deep_copy(Btest, 1.0);
  execution_space().fence();

  for (auto _ : state) {
    Kokkos::Timer timer;
    Kokkos::parallel_for(policy, FunctorType(Atest, Btest, dims));
    execution_space().fence();
    const double dt = timer.seconds();
    state.SetIterationTime(dt);
  }
  // Correctness check
  check_computation<typename FunctorType::scalar_type>(Atest, Btest);
}

template <typename FunctorType>
void bench_mdrange(benchmark::State& state) {
  bench_mdrange<FunctorType>(
      state, std::make_index_sequence<FunctorType::dimension>());
}

template <typename T, std::size_t Rank>
struct add_pointer_n {
  using type = typename add_pointer_n<T*, Rank - 1>::type;
};

template <typename T>
struct add_pointer_n<T, 0> {
  using type = T;
};

template <typename T, std::size_t Rank>
using add_pointer_n_t = typename add_pointer_n<T, Rank>::type;

template <class DeviceType, int Dimension,
          typename TestLayout = Kokkos::LayoutRight,
          typename ScalarType = double>
struct MDRange {
  using execution_space = DeviceType;
  using scalar_type     = ScalarType;
  using size_type       = typename execution_space::size_type;
  using view_type       = Kokkos::View<add_pointer_n_t<ScalarType, Dimension>,
                                 TestLayout, DeviceType>;

  static constexpr int dimension = Dimension;

  view_type A;
  view_type B;
  const Kokkos::Array<int, dimension> ranges;

  template <typename... Dims>
  MDRange(const view_type& A_, const view_type& B_,
          const Kokkos::Array<int, dimension>& dims)
      : A(A_), B(B_), ranges(dims) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int i0, int i1) const
    requires(dimension == 2)
  {
    i0++;
    i1++;
    A(i0, i1) = 0.25 * (ScalarType)(B(i0 + 1, i1) + B(i0 - 1, i1) +
                                    B(i0, i1 + 1) + B(i0, i1 - 1) + B(i0, i1));
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(int i0, int i1, int i2) const
    requires(dimension == 3)
  {
    i0++;
    i1++;
    i2++;
    A(i0, i1, i2) = 0.25 * (ScalarType)(B(i0 + 1, i1, i2) + B(i0 - 1, i1, i2) +
                                        B(i0, i1 + 1, i2) + B(i0, i1 - 1, i2) +
                                        B(i0, i1, i2 + 1) + B(i0, i1, i2 - 1) +
                                        B(i0, i1, i2));
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(int i0, int i1, int i2, int i3) const
    requires(dimension == 4)
  {
    i0++;
    i1++;
    i2++;
    i3++;
    A(i0, i1, i2, i3) =
        0.25 * (ScalarType)(B(i0 + 1, i1, i2, i3) + B(i0 - 1, i1, i2, i3) +
                            B(i0, i1 + 1, i2, i3) + B(i0, i1 - 1, i2, i3) +
                            B(i0, i1, i2 + 1, i3) + B(i0, i1, i2 - 1, i3) +
                            B(i0, i1, i2, i3 + 1) + B(i0, i1, i2, i3 - 1) +
                            B(i0, i1, i2, i3));
  }

  static auto get_policy(const Kokkos::Array<int, dimension>& end,
                         Kokkos::Array<int, dimension>& tile) {
    constexpr Kokkos::Iterate iteration_pattern =
        LayoutToIterationPattern<TestLayout>::pattern;
    const Kokkos::MDRangePolicy<
        Kokkos::Rank<dimension, iteration_pattern, iteration_pattern>,
        execution_space>
        policy(Kokkos::Array<int, dimension>{}, end, tile);

    for (int i = 0; i < dimension; i++) {
      tile[i] = policy.m_tile[i];
    }

    return policy;
  }
};

template <class DeviceType, int Dimension,
          typename TestLayout = Kokkos::LayoutRight,
          typename ScalarType = double>
struct CollapseTwo {
  // RangePolicy for ND range, but will collapse only 2 dims; unroll 2 dims in
  // one-dim

  using execution_space = DeviceType;
  using scalar_type     = ScalarType;
  using size_type       = typename execution_space::size_type;
  using view_type       = Kokkos::View<add_pointer_n_t<ScalarType, Dimension>,
                                 TestLayout, DeviceType>;

  static constexpr int dimension = Dimension;

  view_type A;
  view_type B;
  const Kokkos::Array<int, dimension> ranges;

  CollapseTwo(view_type& A_, const view_type& B_,
              const Kokkos::Array<int, dimension>& dims)
      : A(A_), B(B_), ranges(dims) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int r) const
    requires(dimension == 3)
  {
    if constexpr (std::is_same_v<TestLayout, Kokkos::LayoutLeft>) {
      const int i0 = r % ranges[0] + 1, i1 = r / ranges[0] + 1;
      for (int i2 = 1; i2 < ranges[2] + 1; ++i2) {
        A(i0, i1, i2) =
            0.25 *
            (ScalarType)(B(i0 + 1, i1, i2) + B(i0 - 1, i1, i2) +
                         B(i0, i1 + 1, i2) + B(i0, i1 - 1, i2) +
                         B(i0, i1, i2 + 1) + B(i0, i1, i2 - 1) + B(i0, i1, i2));
      }
    } else {
      const int i2 = r % ranges[2] + 1, i1 = r / ranges[2] + 1;
      for (int i0 = 1; i0 < ranges[0] + 1; ++i0) {
        A(i0, i1, i2) =
            0.25 *
            (ScalarType)(B(i0 + 1, i1, i2) + B(i0 - 1, i1, i2) +
                         B(i0, i1 + 1, i2) + B(i0, i1 - 1, i2) +
                         B(i0, i1, i2 + 1) + B(i0, i1, i2 - 1) + B(i0, i1, i2));
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int r) const
    requires(dimension == 4)
  {
    if constexpr (std::is_same_v<TestLayout, Kokkos::LayoutLeft>) {
      const int i0 = r % ranges[0] + 1, i12 = r / ranges[0];
      const int i1 = i12 % ranges[1] + 1, i2 = i12 / ranges[1] + 1;
      for (int i3 = 1; i3 < ranges[3] + 1; ++i3) {
        A(i0, i1, i2, i3) =
            0.25 * (ScalarType)(B(i0 + 1, i1, i2, i3) + B(i0 - 1, i1, i2, i3) +
                                B(i0, i1 + 1, i2, i3) + B(i0, i1 - 1, i2, i3) +
                                B(i0, i1, i2 + 1, i3) + B(i0, i1, i2 - 1, i3) +
                                B(i0, i1, i2, i3 + 1) + B(i0, i1, i2, i3 - 1) +
                                B(i0, i1, i2, i3));
      }
    } else {
      const int i3 = r % ranges[3] + 1, i21 = r / ranges[3];
      const int i2 = i21 % ranges[2] + 1, i1 = i21 / ranges[2] + 1;
      for (int i0 = 1; i0 < ranges[0] + 1; ++i0) {
        A(i0, i1, i2, i3) =
            0.25 * (ScalarType)(B(i0 + 1, i1, i2, i3) + B(i0 - 1, i1, i2, i3) +
                                B(i0, i1 + 1, i2, i3) + B(i0, i1 - 1, i2, i3) +
                                B(i0, i1, i2 + 1, i3) + B(i0, i1, i2 - 1, i3) +
                                B(i0, i1, i2, i3 + 1) + B(i0, i1, i2, i3 - 1) +
                                B(i0, i1, i2, i3));
      }
    }
  }

  static auto get_policy(const Kokkos::Array<int, dimension>& dims,
                         const Kokkos::Array<int, dimension>&) {
    int collapse_index_rangeA = 0;
    if constexpr (std::is_same_v<TestLayout, Kokkos::LayoutLeft>) {
      collapse_index_rangeA = std::reduce(Kokkos::begin(dims),
                                          Kokkos::begin(dims) + (dimension - 1),
                                          1, std::multiplies<int>{});
    } else if constexpr (std::is_same_v<TestLayout, Kokkos::LayoutRight>) {
      collapse_index_rangeA =
          std::reduce(Kokkos::begin(dims) + 1, Kokkos::end(dims), 1,
                      std::multiplies<int>{});
    } else {
      static_assert(!(std::is_same_v<TestLayout, Kokkos::LayoutRight> ||
                      std::is_same_v<TestLayout, Kokkos::LayoutLeft>),
                    "LayoutRight or LayoutLeft required");
    }

    return Kokkos::RangePolicy<execution_space>(0, collapse_index_rangeA);
  }
};

template <class DeviceType, int Dimension,
          typename TestLayout = Kokkos::LayoutRight,
          typename ScalarType = double>
struct CollapseAll {
  // RangePolicy for ND range, but will collapse all dims

  using execution_space = DeviceType;
  using scalar_type     = ScalarType;
  using size_type       = typename execution_space::size_type;
  using view_type       = Kokkos::View<add_pointer_n_t<ScalarType, Dimension>,
                                 TestLayout, DeviceType>;

  static constexpr int dimension = Dimension;

  view_type A;
  view_type B;
  const Kokkos::Array<int, dimension> ranges;

  template <typename... Dims>
  CollapseAll(view_type& A_, const view_type& B_,
              const Kokkos::Array<int, dimension>& dims)
      : A(A_), B(B_), ranges(dims) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int r) const
    requires(dimension == 2)
  {
    if constexpr (std::is_same_v<TestLayout, Kokkos::LayoutLeft>) {
      const int i0 = r % ranges[0] + 1, i1 = r / ranges[0] + 1;
      A(i0, i1) =
          0.25 * (ScalarType)(B(i0 + 1, i1) + B(i0 - 1, i1) + B(i0, i1 + 1) +
                              B(i0, i1 - 1) + B(i0, i1));
    } else {
      const int i1 = r % ranges[1] + 1, i0 = r / ranges[1] + 1;
      A(i0, i1) =
          0.25 * (ScalarType)(B(i0 + 1, i1) + B(i0 - 1, i1) + B(i0, i1 + 1) +
                              B(i0, i1 - 1) + B(i0, i1));
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int r) const
    requires(dimension == 3)
  {
    if constexpr (std::is_same_v<TestLayout, Kokkos::LayoutLeft>) {
      const int i0 = r % ranges[0] + 1, i12 = r / ranges[0];
      const int i1 = i12 % ranges[1] + 1, i2 = i12 / ranges[1] + 1;
      A(i0, i1, i2) =
          0.25 *
          (ScalarType)(B(i0 + 1, i1, i2) + B(i0 - 1, i1, i2) +
                       B(i0, i1 + 1, i2) + B(i0, i1 - 1, i2) +
                       B(i0, i1, i2 + 1) + B(i0, i1, i2 - 1) + B(i0, i1, i2));
    } else {
      const int i2 = r % ranges[2] + 1, i10 = r / ranges[2];
      const int i1 = i10 % ranges[1] + 1, i0 = i10 / ranges[1] + 1;
      A(i0, i1, i2) =
          0.25 *
          (ScalarType)(B(i0 + 1, i1, i2) + B(i0 - 1, i1, i2) +
                       B(i0, i1 + 1, i2) + B(i0, i1 - 1, i2) +
                       B(i0, i1, i2 + 1) + B(i0, i1, i2 - 1) + B(i0, i1, i2));
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int r) const
    requires(dimension == 4)
  {
    if constexpr (std::is_same_v<TestLayout, Kokkos::LayoutLeft>) {
      const int i0 = r % ranges[0] + 1, i123 = r / ranges[0];
      const int i1 = i123 % ranges[1] + 1, i23 = i123 / ranges[1];
      const int i2 = i23 % ranges[2] + 1, i3 = i23 / ranges[2] + 1;
      A(i0, i1, i2, i3) =
          0.25 * (ScalarType)(B(i0 + 1, i1, i2, i3) + B(i0 - 1, i1, i2, i3) +
                              B(i0, i1 + 1, i2, i3) + B(i0, i1 - 1, i2, i3) +
                              B(i0, i1, i2 + 1, i3) + B(i0, i1, i2 - 1, i3) +
                              B(i0, i1, i2, i3 + 1) + B(i0, i1, i2, i3 - 1) +
                              B(i0, i1, i2, i3));
    } else {
      const int i3 = r % ranges[3] + 1, i210 = r / ranges[3];
      const int i2 = i210 % ranges[2] + 1, i10 = i210 / ranges[2];
      const int i1 = i10 % ranges[1] + 1, i0 = i10 / ranges[1] + 1;
      A(i0, i1, i2, i3) =
          0.25 * (ScalarType)(B(i0 + 1, i1, i2, i3) + B(i0 - 1, i1, i2, i3) +
                              B(i0, i1 + 1, i2, i3) + B(i0, i1 - 1, i2, i3) +
                              B(i0, i1, i2 + 1, i3) + B(i0, i1, i2 - 1, i3) +
                              B(i0, i1, i2, i3 + 1) + B(i0, i1, i2, i3 - 1) +
                              B(i0, i1, i2, i3));
    }
  }

  static auto get_policy(const Kokkos::Array<int, dimension>& dims,
                         const Kokkos::Array<int, dimension>&) {
    const int flat_index_range = std::reduce(
        Kokkos::begin(dims), Kokkos::end(dims), 1, std::multiplies<int>{});
    return Kokkos::RangePolicy<execution_space>(0, flat_index_range);
  }
};

#if !defined(KOKKOS_ENABLE_BENCHMARKS_HEAVY)
#define MDRANGE_STENCIL_BENCHMARK(functor, dim, layout, sizes, ...)      \
  BENCHMARK(bench_mdrange<functor<TEST_EXECSPACE, dim, Kokkos::layout>>) \
      ->UseManualTime()                                                  \
      ->Unit(benchmark::kMillisecond)                                    \
      ->Name("MDRangeStencil_" #dim "D_" #functor "_" #layout)           \
      ->ArgNames({"size", "tile_size"})                                  \
      ->ArgsProduct({sizes, __VA_ARGS__})                                \
      ->Iterations(1);

MDRANGE_STENCIL_BENCHMARK(MDRange, 2, LayoutRight, {512}, {0})
MDRANGE_STENCIL_BENCHMARK(MDRange, 3, LayoutLeft, {128}, {0})
MDRANGE_STENCIL_BENCHMARK(MDRange, 4, LayoutRight, {32}, {0})
#else
#define MDRANGE_STENCIL_BENCHMARK(functor, dim, layout, sizes, ...)      \
  BENCHMARK(bench_mdrange<functor<TEST_EXECSPACE, dim, Kokkos::layout>>) \
      ->UseManualTime()                                                  \
      ->Unit(benchmark::kMillisecond)                                    \
      ->Name("MDRangeStencil_" #dim "D_" #functor "_" #layout)           \
      ->ArgNames({"size", "tile_size"})                                  \
      ->ArgsProduct({sizes, __VA_ARGS__});

#define SIZES_2D \
  { 512, 1024, 2048, 4096, 8192 }
#define SIZES_3D \
  { 128, 192, 256, 512 }
#define SIZES_4D \
  { 32, 64, 96 }
MDRANGE_STENCIL_BENCHMARK(MDRange, 2, LayoutRight, SIZES_2D, {0, 1})
MDRANGE_STENCIL_BENCHMARK(MDRange, 3, LayoutRight, SIZES_3D, {0, 1})
MDRANGE_STENCIL_BENCHMARK(MDRange, 4, LayoutRight, SIZES_4D, {0, 1})
MDRANGE_STENCIL_BENCHMARK(MDRange, 2, LayoutLeft, SIZES_2D, {0, 1})
MDRANGE_STENCIL_BENCHMARK(MDRange, 3, LayoutLeft, SIZES_3D, {0, 1})
MDRANGE_STENCIL_BENCHMARK(MDRange, 4, LayoutLeft, SIZES_4D, {0, 1})

MDRANGE_STENCIL_BENCHMARK(CollapseTwo, 3, LayoutRight, SIZES_3D, {-1})
MDRANGE_STENCIL_BENCHMARK(CollapseTwo, 3, LayoutLeft, SIZES_3D, {-1})
MDRANGE_STENCIL_BENCHMARK(CollapseTwo, 4, LayoutRight, SIZES_4D, {-1})
MDRANGE_STENCIL_BENCHMARK(CollapseTwo, 4, LayoutLeft, SIZES_4D, {-1})

MDRANGE_STENCIL_BENCHMARK(CollapseAll, 2, LayoutRight, SIZES_2D, {-1})
MDRANGE_STENCIL_BENCHMARK(CollapseAll, 2, LayoutLeft, SIZES_2D, {-1})
MDRANGE_STENCIL_BENCHMARK(CollapseAll, 3, LayoutRight, SIZES_3D, {-1})
MDRANGE_STENCIL_BENCHMARK(CollapseAll, 3, LayoutLeft, SIZES_3D, {-1})
MDRANGE_STENCIL_BENCHMARK(CollapseAll, 4, LayoutRight, SIZES_4D, {-1})
MDRANGE_STENCIL_BENCHMARK(CollapseAll, 4, LayoutLeft, SIZES_4D, {-1})
#undef SIZES_2D
#undef SIZES_3D
#undef SIZES_4D
#endif

#undef MDRANGE_STENCIL_BENCHMARK

}  // end namespace Test
