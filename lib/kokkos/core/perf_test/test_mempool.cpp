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

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <limits>

#include <benchmark/benchmark.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Timer.hpp>

#include "Benchmark_Context.hpp"
#include "PerfTest_Category.hpp"

using ExecSpace   = Kokkos::DefaultExecutionSpace;
using MemorySpace = Kokkos::DefaultExecutionSpace::memory_space;

using MemoryPool = Kokkos::MemoryPool<ExecSpace>;

struct TestFunctor {
  using ptrs_type = Kokkos::View<uintptr_t*, ExecSpace>;

  enum : unsigned { chunk = 32 };

  MemoryPool pool;
  ptrs_type ptrs;
  unsigned chunk_span;
  unsigned fill_stride;
  unsigned range_iter;
  unsigned repeat_inner;

  TestFunctor(size_t total_alloc_size, unsigned min_superblock_size,
              unsigned number_alloc, unsigned arg_stride_alloc,
              unsigned arg_chunk_span, unsigned arg_repeat)
      : pool(), ptrs(), chunk_span(0), fill_stride(0), repeat_inner(0) {
    MemorySpace m;

    const unsigned min_block_size = chunk;
    const unsigned max_block_size = chunk * arg_chunk_span;
    pool = MemoryPool(m, total_alloc_size, min_block_size, max_block_size,
                      min_superblock_size);

    ptrs         = ptrs_type(Kokkos::view_alloc(m, "ptrs"), number_alloc);
    fill_stride  = arg_stride_alloc;
    chunk_span   = arg_chunk_span;
    range_iter   = fill_stride * number_alloc;
    repeat_inner = arg_repeat;
  }

  //----------------------------------------

  using value_type = long;

  //----------------------------------------

  struct TagFill {};

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFill, int i, value_type& update) const noexcept {
    if (0 == i % fill_stride) {
      const int j = i / fill_stride;

      const unsigned size_alloc = chunk * (1 + (j % chunk_span));

      ptrs(j) = reinterpret_cast<uintptr_t>(pool.allocate(size_alloc));

      if (ptrs(j)) ++update;
    }
  }

  bool test_fill() {
    using policy = Kokkos::RangePolicy<ExecSpace, TagFill>;

    long result = 0;

    Kokkos::parallel_reduce(policy(0, range_iter), *this, result);

    if (result == long(ptrs.extent(0))) return true;
    pool.print_state(std::cerr);
    return false;
  }

  //----------------------------------------

  struct TagDel {};

  KOKKOS_INLINE_FUNCTION
  void operator()(TagDel, int i) const noexcept {
    if (0 == i % fill_stride) {
      const int j = i / fill_stride;

      const unsigned size_alloc = chunk * (1 + (j % chunk_span));

      pool.deallocate(reinterpret_cast<void*>(ptrs(j)), size_alloc);
    }
  }

  void test_del() {
    using policy = Kokkos::RangePolicy<ExecSpace, TagDel>;

    Kokkos::parallel_for(policy(0, range_iter), *this);
    Kokkos::fence();
  }

  //----------------------------------------

  struct TagAllocDealloc {};

  KOKKOS_INLINE_FUNCTION
  void operator()(TagAllocDealloc, int i, long& update) const noexcept {
    if (0 == i % fill_stride) {
      const int j = i / fill_stride;

      if (0 == j % 3) {
        for (unsigned k = 0; k < repeat_inner; ++k) {
          const unsigned size_alloc = chunk * (1 + (j % chunk_span));

          pool.deallocate(reinterpret_cast<void*>(ptrs(j)), size_alloc);

          ptrs(j) = reinterpret_cast<uintptr_t>(pool.allocate(size_alloc));

          if (0 == ptrs(j)) update++;
        }
      }
    }
  }

  bool test_alloc_dealloc() {
    using policy = Kokkos::RangePolicy<ExecSpace, TagAllocDealloc>;

    long error_count = 0;

    Kokkos::parallel_reduce(policy(0, range_iter), *this, error_count);

    return 0 == error_count;
  }
};

int get_number_alloc(int chunk_span, int min_superblock_size,
                     long total_alloc_size, int fill_level) {
  int chunk_span_bytes = 0;
  for (int i = 0; i < chunk_span; ++i) {
    auto chunk_bytes = TestFunctor::chunk * (1 + i);
    if (chunk_bytes < 64) chunk_bytes = 64;
    auto block_bytes_lg2 =
        Kokkos::Impl::integral_power_of_two_that_contains(chunk_bytes);
    auto block_bytes = (1 << block_bytes_lg2);
    chunk_span_bytes += block_bytes;
  }
  auto actual_superblock_bytes_lg2 =
      Kokkos::Impl::integral_power_of_two_that_contains(min_superblock_size);
  auto actual_superblock_bytes = (1 << actual_superblock_bytes_lg2);
  auto superblock_mask         = actual_superblock_bytes - 1;
  auto nsuperblocks =
      (total_alloc_size + superblock_mask) >> actual_superblock_bytes_lg2;
  auto actual_total_bytes = nsuperblocks * actual_superblock_bytes;
  auto bytes_wanted       = (actual_total_bytes * fill_level) / 100;
  auto chunk_spans        = bytes_wanted / chunk_span_bytes;
  auto number_alloc       = int(chunk_spans * chunk_span);
  return number_alloc;
}

template <class T>
T get_parameter(const char flag[], T default_value) {
  auto argc  = Test::command_line_num_args();
  auto value = default_value;

  for (int i = 1; i < argc; i++) {
    const char* const a = Test::command_line_arg(i);

    if (!strncmp(a, flag, strlen(flag))) value = std::stoi(a + strlen(flag));
  }

  return value;
}

static void Mempool_Fill(benchmark::State& state) {
  long total_alloc_size =
      get_parameter("--alloc_size=", static_cast<long>(state.range(0)));
  int min_superblock_size = get_parameter("--super_size=", state.range(1));
  int chunk_span          = get_parameter("--chunk_span=", state.range(2));
  int fill_stride         = get_parameter("--fill_stride=", state.range(3));
  int fill_level          = get_parameter("--fill_level=", state.range(4));
  int repeat_inner        = get_parameter("--repeat_inner=", state.range(5));
  int number_alloc        = get_number_alloc(chunk_span, min_superblock_size,
                                      total_alloc_size, fill_level);

  for (auto _ : state) {
    TestFunctor functor(total_alloc_size, min_superblock_size, number_alloc,
                        fill_stride, chunk_span, repeat_inner);
    Kokkos::Timer timer;

    if (!functor.test_fill()) {
      Kokkos::abort("fill ");
    }

    state.SetIterationTime(timer.seconds());
    state.counters[KokkosBenchmark::benchmark_fom("fill ops per second")] =
        benchmark::Counter(number_alloc,
                           benchmark::Counter::kIsIterationInvariantRate);
  }
}

static void Mempool_Alloc_Dealloc(benchmark::State& state) {
  long total_alloc_size =
      get_parameter("--alloc_size=", static_cast<long>(state.range(0)));
  int min_superblock_size = get_parameter("--super_size=", state.range(1));
  int chunk_span          = get_parameter("--chunk_span=", state.range(2));
  int fill_stride         = get_parameter("--fill_stride=", state.range(3));
  int fill_level          = get_parameter("--fill_level=", state.range(4));
  int repeat_inner        = get_parameter("--repeat_inner=", state.range(5));
  int number_alloc        = get_number_alloc(chunk_span, min_superblock_size,
                                      total_alloc_size, fill_level);

  for (auto _ : state) {
    TestFunctor functor(total_alloc_size, min_superblock_size, number_alloc,
                        fill_stride, chunk_span, repeat_inner);
    Kokkos::Timer timer;

    if (!functor.test_alloc_dealloc()) {
      Kokkos::abort("alloc/dealloc ");
    }

    state.SetIterationTime(timer.seconds());
    state.counters[KokkosBenchmark::benchmark_fom("cycle ops per second")] =
        benchmark::Counter(2 * number_alloc * repeat_inner,
                           benchmark::Counter::kIsIterationInvariantRate);
  }
}

const std::vector<std::string> ARG_NAMES = {
    "total_alloc_size", "min_superblock_size", "chunk_span",
    "fill_stride",      "fill_level",          "repeat_inner"};
const std::vector<int64_t> ARGS = {1'000'000, 10'000, 5, 1, 70, 1};

BENCHMARK(Mempool_Fill)->ArgNames(ARG_NAMES)->Args(ARGS)->UseManualTime();

BENCHMARK(Mempool_Alloc_Dealloc)
    ->ArgNames(ARG_NAMES)
    ->Args(ARGS)
    ->UseManualTime();
