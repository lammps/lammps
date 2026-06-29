// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>

#include <benchmark/benchmark.h>

#include "Benchmark_Context.hpp"
#include "PerfTest_Category.hpp"

namespace Benchmark {

struct Tag_Set {};
struct Tag_Copy {};
struct Tag_Scale {};
struct Tag_Add {};
struct Tag_Triad {};

template <typename Tag>
int tag_to_data_ratio() {
  if constexpr (std::is_same_v<Tag, Tag_Set>) return 1;
  if constexpr (std::is_same_v<Tag, Tag_Copy>) return 2;
  if constexpr (std::is_same_v<Tag, Tag_Scale>) return 2;
  if constexpr (std::is_same_v<Tag, Tag_Add>) return 3;
  if constexpr (std::is_same_v<Tag, Tag_Triad>) return 3;
  return 0;
}

template <int Rank, typename ScalarType, typename Layout, typename MemorySpace>
struct ViewTypeRank {
  using type = void;
};

template <typename ScalarType, typename Layout, typename MemorySpace>
struct ViewTypeRank<1, ScalarType, Layout, MemorySpace> {
  using type = Kokkos::View<ScalarType *, Layout, MemorySpace>;
};

template <typename ScalarType, typename Layout, typename MemorySpace>
struct ViewTypeRank<2, ScalarType, Layout, MemorySpace> {
  using type = Kokkos::View<ScalarType **, Layout, MemorySpace>;
};

template <typename ScalarType, typename Layout, typename MemorySpace>
struct ViewTypeRank<3, ScalarType, Layout, MemorySpace> {
  using type = Kokkos::View<ScalarType ***, Layout, MemorySpace>;
};

template <typename ScalarType, typename Layout, typename MemorySpace>
struct ViewTypeRank<4, ScalarType, Layout, MemorySpace> {
  using type = Kokkos::View<ScalarType ****, Layout, MemorySpace>;
};

template <typename ScalarType, typename Layout, typename MemorySpace>
struct ViewTypeRank<5, ScalarType, Layout, MemorySpace> {
  using type = Kokkos::View<ScalarType *****, Layout, MemorySpace>;
};

template <typename ScalarType, typename Layout, typename MemorySpace>
struct ViewTypeRank<6, ScalarType, Layout, MemorySpace> {
  using type = Kokkos::View<ScalarType ******, Layout, MemorySpace>;
};

// Select between MDRangePolicy and RangePolicy
template <int Rank, typename ExecutionSpace,
          typename IndexType = Kokkos::IndexType<uint32_t>, typename Tag = void>
struct policy_selector {
  using preferred_layout = typename ExecutionSpace::array_layout;
  static const Kokkos::Iterate outer_iter =
      Kokkos::Impl::layout_iterate_type_selector<
          preferred_layout>::outer_iteration_pattern;
  static const Kokkos::Iterate inner_iter =
      Kokkos::Impl::layout_iterate_type_selector<
          preferred_layout>::inner_iteration_pattern;
  using type = Kokkos::MDRangePolicy<ExecutionSpace,
                                     Kokkos::Rank<Rank, outer_iter, inner_iter>,
                                     IndexType, Tag>;
};

// Specialization for 1D RangePolicy
template <typename ExecutionSpace, typename IndexType, typename Tag>
struct policy_selector<1, ExecutionSpace, IndexType, Tag> {
  using type = Kokkos::RangePolicy<ExecutionSpace, IndexType, Tag>;
};

// Choose between 1d and Nd bound types
template <typename PolicyType>
struct bound_type_selector {
  using type = typename PolicyType::index_type;
};

// Specialization for MDRangePolicy
template <typename... Args>
struct bound_type_selector<Kokkos::MDRangePolicy<Args...>> {
  using type = typename Kokkos::MDRangePolicy<Args...>::point_type;
};

// Functor for stream test (set, copy, scale, add, triad).
// The problem size is N^6, meaning that each view will have the size of N^6
// whatever the rank is.
template <typename ExecutionSpace, int Rank, typename ScalarType = double,
          typename IndexType = Kokkos::IndexType<uint32_t>>
struct MDStreamTest {
  using scalar_type      = ScalarType;
  using execution_space  = ExecutionSpace;
  using memory_space     = typename ExecutionSpace::memory_space;
  using preferred_layout = typename ExecutionSpace::array_layout;
  using view_type = typename ViewTypeRank<Rank, scalar_type, preferred_layout,
                                          memory_space>::type;
  using policy_init_type =
      typename policy_selector<Rank, ExecutionSpace, IndexType, void>::type;
  using bound_type = typename bound_type_selector<policy_init_type>::type;

  view_type m_view_A;
  view_type m_view_B;
  view_type m_view_C;
  ScalarType m_scalar;

  bound_type m_lower_bounds;
  bound_type m_upper_bounds;

  // Functor for initialization
  struct Init {
    view_type m_tensor;
    scalar_type m_value;

    Init(const view_type &tensor, const scalar_type &value)
        : m_tensor(tensor), m_value(value) {}

    template <typename... Indices>
    KOKKOS_INLINE_FUNCTION void operator()(Indices... indices) const {
      m_tensor(indices...) = m_value;
    }
  };

  MDStreamTest(const int N) {
    static_assert(Rank >= 1 && Rank <= 6,
                  "MDStreamTest: Only ranks 1 to 6 supported");

    m_view_A = create_test_view("MDStreamTest::view_A", N);
    m_view_B = create_test_view("MDStreamTest::view_B", N);
    m_view_C = create_test_view("MDStreamTest::view_C", N);
    m_scalar = static_cast<scalar_type>(2.718281828);

    if constexpr (Rank == 1) {
      m_lower_bounds = 0;
      m_upper_bounds = m_view_A.extent(0);
    } else {
      for (int i = 0; i < Rank; ++i) {
        m_lower_bounds[i] = 0;
        m_upper_bounds[i] = m_view_A.extent(i);
      }
    }

    policy_init_type init_policy(m_lower_bounds, m_upper_bounds);
    Kokkos::parallel_for(init_policy,
                         Init(m_view_A, static_cast<ScalarType>(1.0)));
    Kokkos::parallel_for(init_policy,
                         Init(m_view_B, static_cast<ScalarType>(2.0)));
    Kokkos::parallel_for(init_policy,
                         Init(m_view_C, static_cast<ScalarType>(3.0)));
    execution_space().fence();
  }

  // Tagged operator()
  template <typename... Args>
  KOKKOS_INLINE_FUNCTION void operator()(Tag_Set, Args... args) const {
    m_view_A(args...) = static_cast<ScalarType>(m_scalar);
  }

  template <typename... Args>
  KOKKOS_INLINE_FUNCTION void operator()(Tag_Copy, Args... args) const {
    m_view_B(args...) = m_view_A(args...);
  }

  template <typename... Args>
  KOKKOS_INLINE_FUNCTION void operator()(Tag_Scale, Args... args) const {
    m_view_B(args...) = m_scalar * m_view_A(args...);
  }

  template <typename... Args>
  KOKKOS_INLINE_FUNCTION void operator()(Tag_Add, Args... args) const {
    m_view_C(args...) = m_view_A(args...) + m_view_B(args...);
  }

  template <typename... Args>
  KOKKOS_INLINE_FUNCTION void operator()(Tag_Triad, Args... args) const {
    m_view_C(args...) = m_view_A(args...) + m_scalar * m_view_B(args...);
  }

  // Create test views of size N^6
  view_type create_test_view(const char *name, int dim) {
    long N1 = dim;
    long N2 = N1 * N1;
    long N3 = N2 * N1;
    long N6 = N3 * N3;
    std::string view_name(name);
    if constexpr (Rank == 1) {
      return view_type(view_name, N6);
    } else if constexpr (Rank == 2) {
      return view_type(view_name, N3, N3);
    } else if constexpr (Rank == 3) {
      return view_type(view_name, N2, N2, N2);
    } else if constexpr (Rank == 4) {
      return view_type(view_name, N2, N1, N1, N2);
    } else if constexpr (Rank == 5) {
      return view_type(view_name, N1, N1, N2, N1, N1);
    } else if constexpr (Rank == 6) {
      return view_type(view_name, N1, N1, N1, N1, N1, N1);
    }
  }

  template <typename Tag>
  void run_test(benchmark::State &state) {
    using policy_test_type =
        typename policy_selector<Rank, ExecutionSpace, IndexType, Tag>::type;

    policy_test_type compute_policy(m_lower_bounds, m_upper_bounds);

    const int data_ratio = tag_to_data_ratio<Tag>();

    for (auto _ : state) {
      Kokkos::Timer timer;
      Kokkos::parallel_for(compute_policy, *this);
      execution_space().fence();
      KokkosBenchmark::report_results(state, m_view_A, data_ratio,
                                      timer.seconds());
    }
  }

  void test_set(benchmark::State &state) { run_test<Tag_Set>(state); }

  void test_copy(benchmark::State &state) { run_test<Tag_Copy>(state); }

  void test_scale(benchmark::State &state) { run_test<Tag_Scale>(state); }

  void test_add(benchmark::State &state) { run_test<Tag_Add>(state); }

  void test_triad(benchmark::State &state) { run_test<Tag_Triad>(state); }
};

}  // namespace Benchmark
