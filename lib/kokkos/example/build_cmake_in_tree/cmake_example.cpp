// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>

#include <cstdio>
#include <iostream>

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  Kokkos::DefaultExecutionSpace{}.print_configuration(std::cout);

  if (argc < 2) {
    fprintf(stderr, "Usage: %s [<kokkos_options>] <size>\n", argv[0]);
    Kokkos::finalize();
    exit(1);
  }

  const long n = strtol(argv[1], nullptr, 10);

  printf("Number of even integers from 0 to %ld\n", n - 1);

  Kokkos::Timer timer;
  timer.reset();

  // Compute the number of even integers from 0 to n-1, in parallel.
  long count = 0;
  Kokkos::parallel_reduce(
      n, KOKKOS_LAMBDA(const long i, long& lcount) { lcount += (i % 2) == 0; },
      count);

  double count_time = timer.seconds();
  printf("  Parallel: %ld    %10.6f\n", count, count_time);

  timer.reset();

  // Compare to a sequential loop.
  long seq_count = 0;
  for (long i = 0; i < n; ++i) {
    seq_count += (i % 2) == 0;
  }

  count_time = timer.seconds();
  printf("Sequential: %ld    %10.6f\n", seq_count, count_time);

  Kokkos::finalize();

  return (count == seq_count) ? 0 : -1;
}
