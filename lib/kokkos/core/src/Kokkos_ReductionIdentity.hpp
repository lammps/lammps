// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_REDUCTION_IDENTITY_HPP
#define KOKKOS_REDUCTION_IDENTITY_HPP

#include <Kokkos_Macros.hpp>
#include <concepts>
#include <limits>

namespace Kokkos {

template <class T>
struct reduction_identity;

template <typename Integral>
  requires(std::integral<Integral>)
struct reduction_identity<Integral> {
 private:
  static constexpr auto max_ = std::numeric_limits<Integral>::max();
  static constexpr auto min_ = std::numeric_limits<Integral>::min();

 public:
  KOKKOS_FUNCTION constexpr static Integral sum() noexcept { return 0; }
  KOKKOS_FUNCTION constexpr static Integral prod() noexcept { return 1; }
  KOKKOS_FUNCTION constexpr static Integral max() noexcept { return min_; }
  KOKKOS_FUNCTION constexpr static Integral min() noexcept { return max_; }
  KOKKOS_FUNCTION constexpr static Integral bor() noexcept { return 0x0; }
  KOKKOS_FUNCTION constexpr static Integral band() noexcept { return 0x0; }
  KOKKOS_FUNCTION constexpr static Integral lor() noexcept { return 0; }
  KOKKOS_FUNCTION constexpr static Integral land() noexcept { return 1; }
};

template <typename Floating>
  requires(std::floating_point<Floating> && sizeof(Floating) <= sizeof(double))
struct reduction_identity<Floating> {
 private:
  static constexpr auto inf =
#if __FINITE_MATH_ONLY__
      std::numeric_limits<Floating>::max();
#else
      std::numeric_limits<Floating>::infinity();
#endif
 public:
  KOKKOS_FUNCTION constexpr static Floating sum() noexcept { return 0; }
  KOKKOS_FUNCTION constexpr static Floating prod() noexcept { return 1; }
  KOKKOS_FUNCTION constexpr static Floating max() noexcept { return -inf; }
  KOKKOS_FUNCTION constexpr static Floating min() noexcept { return +inf; }
};

// No __host__ __device__ annotation because long double treated as double in
// device code.  May be revisited later if that is not true any more.
template <typename Floating>
  requires(std::floating_point<Floating> && sizeof(Floating) > sizeof(double))
struct reduction_identity<Floating> {
 private:
  static constexpr auto inf =
#if __FINITE_MATH_ONLY__
      std::numeric_limits<Floating>::max();
#else
      std::numeric_limits<Floating>::infinity();
#endif

 public:
  constexpr static Floating sum() noexcept { return 0; }
  constexpr static Floating prod() noexcept { return 1; }
  constexpr static Floating max() noexcept { return -inf; }
  constexpr static Floating min() noexcept { return +inf; }
};

}  // namespace Kokkos

#endif
