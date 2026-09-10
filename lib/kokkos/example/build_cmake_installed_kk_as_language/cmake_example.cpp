// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>

#include <cstdio>
#include <iostream>

extern "C" void print_fortran_();
void print_cxx();

struct CountEvenIntegers {
  KOKKOS_FUNCTION void operator()(const long i, long& lcount) const {
    lcount += (i % 2) == 0;
  }
};

int main(int argc, char* argv[]) {
  Kokkos::ScopeGuard guard(argc, argv);
  Kokkos::DefaultExecutionSpace().print_configuration(std::cout);

  const long n = argc > 1 ? atoi(argv[1]) : 10;

  printf("Number of even integers from 0 to %ld\n", n - 1);

  Kokkos::Timer timer;
  timer.reset();

  // Compute the number of even integers from 0 to n-1, in parallel.
  long count = 0;
  CountEvenIntegers functor;
  Kokkos::parallel_reduce(n, functor, count);

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

  print_fortran_();

  return (count == seq_count) ? 0 : -1;
}
