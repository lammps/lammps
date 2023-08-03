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

#include "Kokkos_Core.hpp"
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <sys/time.h>

#define HLINE "-------------------------------------------------------------\n"

#if defined(KOKKOS_ENABLE_CUDA)
using GUPSHostArray   = Kokkos::View<int64_t*, Kokkos::CudaSpace>::HostMirror;
using GUPSDeviceArray = Kokkos::View<int64_t*, Kokkos::CudaSpace>;
#else
using GUPSHostArray   = Kokkos::View<int64_t*, Kokkos::HostSpace>::HostMirror;
using GUPSDeviceArray = Kokkos::View<int64_t*, Kokkos::HostSpace>;
#endif

using GUPSIndex = int;

double now() {
  struct timeval now;
  gettimeofday(&now, nullptr);

  return (double)now.tv_sec + ((double)now.tv_usec * 1.0e-6);
}

void randomize_indices(GUPSHostArray& indices, GUPSDeviceArray& dev_indices,
                       const int64_t dataCount) {
  for (GUPSIndex i = 0; i < indices.extent(0); ++i) {
    indices[i] = lrand48() % dataCount;
  }

  Kokkos::deep_copy(dev_indices, indices);
}

void run_gups(GUPSDeviceArray& indices, GUPSDeviceArray& data,
              const int64_t datum, const bool performAtomics) {
  if (performAtomics) {
    Kokkos::parallel_for(
        "bench-gups-atomic", indices.extent(0),
        KOKKOS_LAMBDA(const GUPSIndex i) {
          Kokkos::atomic_fetch_xor(&data[indices[i]], datum);
        });
  } else {
    Kokkos::parallel_for(
        "bench-gups-non-atomic", indices.extent(0),
        KOKKOS_LAMBDA(const GUPSIndex i) { data[indices[i]] ^= datum; });
  }

  Kokkos::fence();
}

int run_benchmark(const GUPSIndex indicesCount, const GUPSIndex dataCount,
                  const int repeats, const bool useAtomics) {
  printf("Reports fastest timing per kernel\n");
  printf("Creating Views...\n");

  printf("Memory Sizes:\n");
  printf("- Elements:      %15" PRIu64 " (%12.4f MB)\n",
         static_cast<uint64_t>(dataCount),
         1.0e-6 * ((double)dataCount * (double)sizeof(int64_t)));
  printf("- Indices:       %15" PRIu64 " (%12.4f MB)\n",
         static_cast<uint64_t>(indicesCount),
         1.0e-6 * ((double)indicesCount * (double)sizeof(int64_t)));
  printf(" - Atomics:      %15s\n", (useAtomics ? "Yes" : "No"));
  printf("Benchmark kernels will be performed for %d iterations.\n", repeats);

  printf(HLINE);

  GUPSDeviceArray dev_indices("indices", indicesCount);
  GUPSDeviceArray dev_data("data", dataCount);
  int64_t datum = -1;

  GUPSHostArray indices = Kokkos::create_mirror_view(dev_indices);
  GUPSHostArray data    = Kokkos::create_mirror_view(dev_data);

  double gupsTime = 0.0;

  printf("Initializing Views...\n");

#if defined(KOKKOS_HAVE_OPENMP)
  Kokkos::parallel_for(
      "init-data", Kokkos::RangePolicy<Kokkos::OpenMP>(0, dataCount),
#else
  Kokkos::parallel_for(
      "init-data", Kokkos::RangePolicy<Kokkos::Serial>(0, dataCount),
#endif
      KOKKOS_LAMBDA(const int i) { data[i] = 10101010101; });

#if defined(KOKKOS_HAVE_OPENMP)
  Kokkos::parallel_for(
      "init-indices", Kokkos::RangePolicy<Kokkos::OpenMP>(0, indicesCount),
#else
  Kokkos::parallel_for(
      "init-indices", Kokkos::RangePolicy<Kokkos::Serial>(0, indicesCount),
#endif
      KOKKOS_LAMBDA(const int i) { indices[i] = 0; });

  Kokkos::deep_copy(dev_data, data);
  Kokkos::deep_copy(dev_indices, indices);
  double start;

  printf("Starting benchmarking...\n");

  for (GUPSIndex k = 0; k < repeats; ++k) {
    randomize_indices(indices, dev_indices, data.extent(0));

    start = now();
    run_gups(dev_indices, dev_data, datum, useAtomics);
    gupsTime += now() - start;
  }

  Kokkos::deep_copy(indices, dev_indices);
  Kokkos::deep_copy(data, dev_data);

  printf(HLINE);
  printf(
      "GUP/s Random:      %18.6f\n",
      (1.0e-9 * ((double)repeats) * (double)dev_indices.extent(0)) / gupsTime);
  printf(HLINE);

  return 0;
}

int main(int argc, char* argv[]) {
  printf(HLINE);
  printf("Kokkos GUPS Benchmark\n");
  printf(HLINE);

  srand48(1010101);

  Kokkos::initialize(argc, argv);

  int64_t indices = 8192;
  int64_t data    = 33554432;
  int64_t repeats = 10;
  bool useAtomics = false;

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
    }
  }

  const int rc = run_benchmark(indices, data, repeats, useAtomics);

  Kokkos::finalize();

  return rc;
}
