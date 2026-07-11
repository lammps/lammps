// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

// Replaces the usual #include <Kokkos_Core.hpp>
import kokkos.core;
// We don't get any transitive includes or macro definitions when importing
// C++20 modules. Since pretty much every Kokkos-based code needs configuration
// macros, it's recommended to always include Kokkos_Macros.hpp. In this
// example, we need it for KOKKOS_LAMBDA in particular.
#include <Kokkos_Macros.hpp>

#include <iostream>

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  Kokkos::DefaultExecutionSpace{}.print_configuration(std::cout);

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " [<kokkos_options>] <size>\n";
    Kokkos::finalize();
    exit(1);
  }

  const long n = strtol(argv[1], nullptr, 10);

  std::cout << "Number of even integers from 0 to " << n - 1 << '\n';

  Kokkos::Timer timer;

  // Compute the number of even integers from 0 to n-1, in parallel.
  long count = 0;
  Kokkos::parallel_reduce(
      n, KOKKOS_LAMBDA(const long i, long& lcount) { lcount += (i % 2) == 0; },
      count);

  double count_time = timer.seconds();
  std::cout << "  Parallel: " << count << "    " << count_time << '\n';

  timer.reset();

  // Compare to a sequential loop.
  long seq_count = 0;
  for (long i = 0; i < n; ++i) {
    seq_count += (i % 2) == 0;
  }

  count_time = timer.seconds();
  std::cout << "Sequential: " << seq_count << "    " << count_time << '\n';

  Kokkos::finalize();

  return (count == seq_count) ? 0 : -1;
}
