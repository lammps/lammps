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

#include <Kokkos_Core.hpp>
#include "Kokkos_Core_fwd.hpp"

namespace {

struct TestRangePolicyCTAD {
  struct SomeExecutionSpace {
    using execution_space = SomeExecutionSpace;
    using size_type       = size_t;

    [[maybe_unused]] static int concurrency() { return 0; }
  };
  static_assert(Kokkos::is_execution_space_v<SomeExecutionSpace>);

  struct ImplicitlyConvertibleToDefaultExecutionSpace {
    [[maybe_unused]] operator Kokkos::DefaultExecutionSpace() const {
      return Kokkos::DefaultExecutionSpace();
    }
  };
  static_assert(!Kokkos::is_execution_space_v<
                ImplicitlyConvertibleToDefaultExecutionSpace>);

  [[maybe_unused]] static inline auto i64 = int64_t();
  [[maybe_unused]] static inline auto i32 = int32_t();
  [[maybe_unused]] static inline auto cs  = Kokkos::ChunkSize(0);
  [[maybe_unused]] static inline auto des = Kokkos::DefaultExecutionSpace();
  [[maybe_unused]] static inline auto nes =
      ImplicitlyConvertibleToDefaultExecutionSpace();
  [[maybe_unused]] static inline auto ses = SomeExecutionSpace();

  // RangePolicy()

  [[maybe_unused]] static inline auto rp = Kokkos::RangePolicy{};
  static_assert(std::is_same_v<Kokkos::RangePolicy<>, decltype(rp)>);

  // RangePolicy(index_type, index_type)

  [[maybe_unused]] static inline auto rpi64i64 = Kokkos::RangePolicy(i64, i64);
  static_assert(std::is_same_v<Kokkos::RangePolicy<>, decltype(rpi64i64)>);

  [[maybe_unused]] static inline auto rpi64i32 = Kokkos::RangePolicy(i64, i32);
  static_assert(std::is_same_v<Kokkos::RangePolicy<>, decltype(rpi64i32)>);

  [[maybe_unused]] static inline auto rpi32i64 = Kokkos::RangePolicy(i32, i64);
  static_assert(std::is_same_v<Kokkos::RangePolicy<>, decltype(rpi32i64)>);

  [[maybe_unused]] static inline auto rpi32i32 = Kokkos::RangePolicy(i32, i32);
  static_assert(std::is_same_v<Kokkos::RangePolicy<>, decltype(rpi32i32)>);

  // RangePolicy(index_type, index_type, ChunkSize)

  [[maybe_unused]] static inline auto rpi64i64cs =
      Kokkos::RangePolicy(i64, i64, cs);
  static_assert(std::is_same_v<Kokkos::RangePolicy<>, decltype(rpi64i64cs)>);

  [[maybe_unused]] static inline auto rpi64i32cs =
      Kokkos::RangePolicy(i64, i32, cs);
  static_assert(std::is_same_v<Kokkos::RangePolicy<>, decltype(rpi64i32cs)>);

  [[maybe_unused]] static inline auto rpi32i64cs =
      Kokkos::RangePolicy(i32, i64, cs);
  static_assert(std::is_same_v<Kokkos::RangePolicy<>, decltype(rpi32i64cs)>);

  [[maybe_unused]] static inline auto rpi32i32cs =
      Kokkos::RangePolicy(i32, i32, cs);
  static_assert(std::is_same_v<Kokkos::RangePolicy<>, decltype(rpi32i32cs)>);

  // RangePolicy(execution_space, index_type, index_type)

  [[maybe_unused]] static inline auto rpdesi64i64 =
      Kokkos::RangePolicy(des, i64, i64);
  static_assert(std::is_same_v<Kokkos::RangePolicy<>, decltype(rpdesi64i64)>);

  [[maybe_unused]] static inline auto rpdesi32i32 =
      Kokkos::RangePolicy(des, i32, i32);
  static_assert(std::is_same_v<Kokkos::RangePolicy<>, decltype(rpdesi32i32)>);

  [[maybe_unused]] static inline auto rpnesi64i64 =
      Kokkos::RangePolicy(nes, i64, i64);
  static_assert(std::is_same_v<Kokkos::RangePolicy<>, decltype(rpnesi64i64)>);

  [[maybe_unused]] static inline auto rpnesi32i32 =
      Kokkos::RangePolicy(nes, i32, i32);
  static_assert(std::is_same_v<Kokkos::RangePolicy<>, decltype(rpnesi32i32)>);

  [[maybe_unused]] static inline auto rpsesi64i64 =
      Kokkos::RangePolicy(ses, i64, i64);
  static_assert(std::is_same_v<Kokkos::RangePolicy<SomeExecutionSpace>,
                               decltype(rpsesi64i64)>);

  [[maybe_unused]] static inline auto rpsesi32i32 =
      Kokkos::RangePolicy(ses, i32, i32);
  static_assert(std::is_same_v<Kokkos::RangePolicy<SomeExecutionSpace>,
                               decltype(rpsesi32i32)>);

  // RangePolicy(execution_space, index_type, index_type, ChunkSize)

  [[maybe_unused]] static inline auto rpdesi64i64cs =
      Kokkos::RangePolicy(des, i64, i64, cs);
  static_assert(std::is_same_v<Kokkos::RangePolicy<>, decltype(rpdesi64i64cs)>);

  [[maybe_unused]] static inline auto rpdesi32i32cs =
      Kokkos::RangePolicy(des, i32, i32, cs);
  static_assert(std::is_same_v<Kokkos::RangePolicy<>, decltype(rpdesi32i32cs)>);

  [[maybe_unused]] static inline auto rpnesi64i64cs =
      Kokkos::RangePolicy(nes, i64, i64, cs);
  static_assert(std::is_same_v<Kokkos::RangePolicy<>, decltype(rpnesi64i64cs)>);

  [[maybe_unused]] static inline auto rpnesi32i32cs =
      Kokkos::RangePolicy(nes, i32, i32, cs);
  static_assert(std::is_same_v<Kokkos::RangePolicy<>, decltype(rpnesi32i32cs)>);

  [[maybe_unused]] static inline auto rpsesi64i64cs =
      Kokkos::RangePolicy(ses, i64, i64, cs);
  static_assert(std::is_same_v<Kokkos::RangePolicy<SomeExecutionSpace>,
                               decltype(rpsesi64i64cs)>);

  [[maybe_unused]] static inline auto rpsesi32i32cs =
      Kokkos::RangePolicy(ses, i32, i32, cs);
  static_assert(std::is_same_v<Kokkos::RangePolicy<SomeExecutionSpace>,
                               decltype(rpsesi32i32cs)>);

};  // TestRangePolicyCTAD struct

// To eliminate maybe_unused warning on some compilers

[[maybe_unused]] const Kokkos::DefaultExecutionSpace nestodes =
    TestRangePolicyCTAD::ImplicitlyConvertibleToDefaultExecutionSpace();

[[maybe_unused]] const auto sesconcurrency =
    TestRangePolicyCTAD::ses.concurrency();

}  // namespace
