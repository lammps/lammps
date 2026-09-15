// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>

//
// First reduction (parallel_reduce) example:
//   1. Start up Kokkos
//   2. Execute a parallel_reduce loop in the default execution space,
//      using a functor to define the loop body
//   3. Shut down Kokkos
//
struct collision {
  // Reduction functor
  // For each i, we generate 10 hashes, look for and count collisions
  // We use parallel_reduce to count the total collisions
  // Note that we're just counting collisions within the 10 generated
  // one i.
  // The body keeps many scalars active (illustrative workload only; how that
  // maps to occupancy or register use is backend- and toolchain-specific).
  using value_type = int;

  KOKKOS_INLINE_FUNCTION
  int hash(int q) const {
    // A simple hash by Justin Sobel
    // Thanks to Arash Partow (partow.net)
    char* fourchars =
        (char*)&q;  // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
    int hash = 1315423911;
    for (int i = 0; i < 4; fourchars++, i++) {
      hash ^= ((hash << 5) + *fourchars + (hash >> 2));
    }
    return hash;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, int& lsum) const {
    // This is a silly function which generates 10 hashes
    // then checks for collisions
    int a = hash(i) % 64;
    int b = hash(i * 3) % 64;
    int c = hash(i * 5) % 64;
    int d = hash(i * 7) % 64;
    int e = hash(i * 11) % 64;
    int f = hash(i * 17) % 64;
    int g = hash(i * 23) % 64;
    int h = hash(i * 29) % 64;
    int j = hash(i * 31) % 64;
    int k = hash(i * 37) % 64;

    if (a == b) lsum++;
    if (a == c) lsum++;
    if (a == d) lsum++;
    if (a == e) lsum++;
    if (a == f) lsum++;
    if (a == g) lsum++;
    if (a == h) lsum++;
    if (a == j) lsum++;
    if (a == k) lsum++;
    if (b == c) lsum++;
    if (b == d) lsum++;
    if (b == e) lsum++;
    if (b == f) lsum++;
    if (b == g) lsum++;
    if (b == h) lsum++;
    if (b == j) lsum++;
    if (b == k) lsum++;
    if (c == d) lsum++;
    if (c == e) lsum++;
    if (c == f) lsum++;
    if (c == g) lsum++;
    if (c == h) lsum++;
    if (c == j) lsum++;
    if (c == k) lsum++;
    if (d == e) lsum++;
    if (d == f) lsum++;
    if (d == g) lsum++;
    if (d == h) lsum++;
    if (d == j) lsum++;
    if (d == k) lsum++;
    if (e == f) lsum++;
    if (e == g) lsum++;
    if (e == h) lsum++;
    if (e == j) lsum++;
    if (e == k) lsum++;
    if (f == g) lsum++;
    if (f == h) lsum++;
    if (f == j) lsum++;
    if (f == k) lsum++;
    if (g == h) lsum++;
    if (g == j) lsum++;
    if (g == k) lsum++;
    if (h == j) lsum++;
    if (h == k) lsum++;
    if (j == k) lsum++;
  }
};

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  const int n = 10000;

  // Compute and count hash collisions in parallel, using Kokkos.
  // This is not really a useful algorithm, but it demonstrates the
  // LaunchBounds functionality
  int sum_default = 0;
  int sum_bounds  = 0;

  // Default reduction (no LaunchBounds trait on the range policy).
  Kokkos::parallel_reduce(n, collision(), sum_default);

  // Same reduction with LaunchBounds<512, min_blocks> carried on the policy.
  // That is only a hint; resource usage and correctness should be verified on
  // your target GPU and Kokkos build.
  const int min_blocks =
#if defined(KOKKOS_ARCH_TURING75) || defined(KOKKOS_ARCH_ADA89) || \
    defined(KOKKOS_ARCH_AMPERE86)
      // We only use 2 blocks to stay below 1024 threads per SM.
      2;
#else
      // We have a maximum of 2048 resident threads per SM.
      4;
#endif

  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<Kokkos::LaunchBounds<512, min_blocks>>(0, n),
      collision(), sum_bounds);

  Kokkos::printf(
      "Number of collisions, "
      "computed in parallel, is %i\n",
      sum_default);

  if (sum_default != sum_bounds) {
    Kokkos::printf("Uh-oh! Results do not match\n");
  }

  Kokkos::finalize();
  return 0;
}
