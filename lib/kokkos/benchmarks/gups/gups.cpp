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

/*! \brief file gups.cpp

    An implementation of something like HPCC RandomAccess.
*/

#include "Kokkos_Core.hpp"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <numeric>
#include <algorithm>
#include <random>

#define HLINE "-------------------------------------------------------------\n"

using Index = int;
using Datum = int64_t;

using IndexView = Kokkos::View<Index*>;
using DataView  = Kokkos::View<Datum*>;

using Clock    = std::chrono::steady_clock;
using Duration = std::chrono::duration<double>;

using RandomDevice = std::random_device;
using RNG          = std::mt19937;

IndexView randomized_indices(const Index indicesCount, const Index dataCount,
                             RNG& rng) {
  // generate random indices 0..dataCount
  std::uniform_int_distribution<Index> uid(0, dataCount);
  std::vector<Index> indices(indicesCount);
  std::generate(indices.begin(), indices.end(), [&]() { return uid(rng); });

  // Copy to the default space and return
  Kokkos::View<Index*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
      unmanaged_indices(indices.data(), indices.size());
  IndexView dev_indices("dev_indices", indicesCount);
  Kokkos::deep_copy(dev_indices, unmanaged_indices);
  return dev_indices;
}

IndexView permuted_indices(const Index indicesCount, const Index dataCount,
                           RNG& rng) {
  // create a permutation array of offsets into the data
  std::vector<Index> perm(dataCount);
  std::iota(perm.begin(), perm.end(), 0);
  std::shuffle(perm.begin(), perm.end(), rng);

  // indices is repeated copies of the permutation array
  // (or the first entries of the permutation array if there
  // are fewer indices than data elements)
  IndexView dev_indices("dev_indices", indicesCount);
  auto indices = Kokkos::create_mirror_view(dev_indices);
  for (Index i = 0; i < Index(indices.extent(0)); ++i) {
    indices(i) = perm[i % perm.size()];
  }

  // Copy to the default space and return

  Kokkos::deep_copy(dev_indices, indices);
  return dev_indices;
}

void run_gups(IndexView& indices, DataView& data, const Datum datum,
              const bool performAtomics) {
  if (performAtomics) {
    Kokkos::parallel_for(
        "bench-gups-atomic", indices.extent(0), KOKKOS_LAMBDA(const Index i) {
          Kokkos::atomic_fetch_xor(&data[indices[i]], datum);
        });
  } else {
    Kokkos::parallel_for(
        "bench-gups-non-atomic", indices.extent(0),
        KOKKOS_LAMBDA(const Index i) { data[indices[i]] ^= datum; });
  }

  Kokkos::fence();
}

enum class AccessPattern { random, permutation };

int run_benchmark(const Index indicesCount, const Index dataCount,
                  const int repeats, const bool useAtomics,
                  const AccessPattern pattern) {
  constexpr auto arbitrary_seed = 20230913;
  RNG rng(arbitrary_seed);

  printf("Reports fastest timing per kernel\n");
  printf("Creating Views...\n");

  printf("Memory Sizes:\n");
  printf("- Elements:      %15" PRIu64 " (%12.4f MB)\n",
         static_cast<uint64_t>(dataCount),
         1.0e-6 * ((double)dataCount * (double)sizeof(Datum)));
  printf("- Indices:       %15" PRIu64 " (%12.4f MB)\n",
         static_cast<uint64_t>(indicesCount),
         1.0e-6 * ((double)indicesCount * (double)sizeof(Index)));
  printf(" - Atomics:      %15s\n", (useAtomics ? "Yes" : "No"));
  printf("Benchmark kernels will be performed for %d iterations.\n", repeats);

  printf(HLINE);

  printf("Initializing Data...\n");
  DataView data("data", dataCount);
  Kokkos::parallel_for(
      "init-data",
      Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, dataCount),
      KOKKOS_LAMBDA(const int i) { data[i] = 10101010101; });

  printf("Starting benchmarking...\n");
  double gupsTime       = 0.0;
  constexpr Datum datum = -1;
  for (Index k = 0; k < repeats; ++k) {
    IndexView indices;
    switch (pattern) {
      case AccessPattern::random: {
        indices = randomized_indices(indicesCount, dataCount, rng);
        break;
      }
      case AccessPattern::permutation: {
        indices = permuted_indices(indicesCount, dataCount, rng);
        break;
      }
      default: {
        throw std::runtime_error("unexpected mode");
      }
    }

    auto start = Clock::now();
    run_gups(indices, data, datum, useAtomics);
    gupsTime += Duration(Clock::now() - start).count();
  }

  printf(HLINE);
  printf("GUP/s Random:      %18.6f\n",
         (1.0e-9 * ((double)repeats) * (double)indicesCount) / gupsTime);
  printf(HLINE);

  return 0;
}

int main(int argc, char* argv[]) {
  printf(HLINE);
  printf("Kokkos GUPS Benchmark\n");
  printf(HLINE);

  Kokkos::initialize(argc, argv);

  int64_t indices       = 8192;
  int64_t data          = 33554432;
  int64_t repeats       = 10;
  bool useAtomics       = false;
  AccessPattern pattern = AccessPattern::random;

  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "--indices") == 0) {
      indices = std::atoll(argv[i + 1]);
      ++i;
    } else if (strcmp(argv[i], "--data") == 0) {
      data = std::atoll(argv[i + 1]);
      ++i;
    } else if (strcmp(argv[i], "--repeats") == 0) {
      repeats = std::atoll(argv[i + 1]);
      ++i;
    } else if (strcmp(argv[i], "--atomics") == 0) {
      useAtomics = true;
    } else if (strcmp(argv[i], "--pattern-permutation") == 0) {
      pattern = AccessPattern::permutation;
    }
  }

  const int rc = run_benchmark(indices, data, repeats, useAtomics, pattern);

  Kokkos::finalize();

  return rc;
}
