// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <concepts>

// consteval specifier
consteval int sqr(int n) { return n * n; }
static_assert(sqr(100) == 10000);

// conditional explicit
struct S {
  explicit(sizeof(int) > 0) S(int) {}
};

// concepts library
constexpr std::floating_point auto x2(std::floating_point auto x) {
  return x + x;
}
constexpr std::integral auto x2(std::integral auto x) { return x << 1; }

int main() { return 0; }
