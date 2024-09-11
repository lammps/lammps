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

namespace {

struct TestTeamPolicyCTAD {
  template <typename... Ts>
  static void maybe_unused(Ts&&...) {}

  struct SomeExecutionSpace {
    using execution_space = SomeExecutionSpace;
    using size_type       = size_t;
  };
  static_assert(Kokkos::is_execution_space_v<SomeExecutionSpace>);

  struct ImplicitlyConvertibleToDefaultExecutionSpace {
    [[maybe_unused]] operator Kokkos::DefaultExecutionSpace() const {
      return Kokkos::DefaultExecutionSpace();
    }
  };
  static_assert(!Kokkos::is_execution_space_v<
                ImplicitlyConvertibleToDefaultExecutionSpace>);

  [[maybe_unused]] static inline Kokkos::DefaultExecutionSpace des;
  [[maybe_unused]] static inline ImplicitlyConvertibleToDefaultExecutionSpace
      notEs;
  [[maybe_unused]] static inline SomeExecutionSpace ses;

  [[maybe_unused]] static inline int i;

  // Workaround for nvc++ (CUDA-11.7-NVHPC) ignoring [[maybe_unused]] on
  // ImplicitlyConvertibleToDefaultExecutionSpace::operator
  // Kokkos::DefaultExecutionSpace() const
  [[maybe_unused]] static inline Kokkos::DefaultExecutionSpace notEsToDes =
      notEs;

  // Workaround for HIP-ROCm-5.2 warning about was declared but never referenced
  TestTeamPolicyCTAD() { maybe_unused(des, notEs, ses, i, notEsToDes); }

  // Default construction deduces to TeamPolicy<>
  static_assert(
      std::is_same_v<Kokkos::TeamPolicy<>, decltype(Kokkos::TeamPolicy{})>);

  // Execution space not provided deduces to TeamPolicy<>

  static_assert(
      std::is_same_v<Kokkos::TeamPolicy<>, decltype(Kokkos::TeamPolicy(i, i))>);
  static_assert(std::is_same_v<Kokkos::TeamPolicy<>,
                               decltype(Kokkos::TeamPolicy(i, i, i))>);
  static_assert(std::is_same_v<Kokkos::TeamPolicy<>,
                               decltype(Kokkos::TeamPolicy(i, Kokkos::AUTO))>);
  static_assert(
      std::is_same_v<Kokkos::TeamPolicy<>,
                     decltype(Kokkos::TeamPolicy(i, Kokkos::AUTO, i))>);
  static_assert(std::is_same_v<Kokkos::TeamPolicy<>,
                               decltype(Kokkos::TeamPolicy(i, Kokkos::AUTO,
                                                           Kokkos::AUTO))>);
  static_assert(
      std::is_same_v<Kokkos::TeamPolicy<>,
                     decltype(Kokkos::TeamPolicy(i, i, Kokkos::AUTO))>);

  // DefaultExecutionSpace deduces to TeamPolicy<>

  static_assert(std::is_same_v<Kokkos::TeamPolicy<>,
                               decltype(Kokkos::TeamPolicy(des, i, i))>);
  static_assert(std::is_same_v<Kokkos::TeamPolicy<>,
                               decltype(Kokkos::TeamPolicy(des, i, i, i))>);
  static_assert(
      std::is_same_v<Kokkos::TeamPolicy<>,
                     decltype(Kokkos::TeamPolicy(des, i, Kokkos::AUTO))>);
  static_assert(
      std::is_same_v<Kokkos::TeamPolicy<>,
                     decltype(Kokkos::TeamPolicy(des, i, Kokkos::AUTO, i))>);
  static_assert(std::is_same_v<Kokkos::TeamPolicy<>,
                               decltype(Kokkos::TeamPolicy(des, i, Kokkos::AUTO,
                                                           Kokkos::AUTO))>);
  static_assert(
      std::is_same_v<Kokkos::TeamPolicy<>,
                     decltype(Kokkos::TeamPolicy(des, i, i, Kokkos::AUTO))>);

  // Convertible to DefaultExecutionSpace deduces to TeamPolicy<>

  static_assert(std::is_same_v<Kokkos::TeamPolicy<>,
                               decltype(Kokkos::TeamPolicy(notEs, i, i))>);
  static_assert(std::is_same_v<Kokkos::TeamPolicy<>,
                               decltype(Kokkos::TeamPolicy(notEs, i, i, i))>);
  static_assert(
      std::is_same_v<Kokkos::TeamPolicy<>,
                     decltype(Kokkos::TeamPolicy(notEs, i, Kokkos::AUTO))>);
  static_assert(
      std::is_same_v<Kokkos::TeamPolicy<>,
                     decltype(Kokkos::TeamPolicy(notEs, i, Kokkos::AUTO, i))>);
  static_assert(std::is_same_v<Kokkos::TeamPolicy<>,
                               decltype(Kokkos::TeamPolicy(
                                   notEs, i, Kokkos::AUTO, Kokkos::AUTO))>);
  static_assert(
      std::is_same_v<Kokkos::TeamPolicy<>,
                     decltype(Kokkos::TeamPolicy(notEs, i, i, Kokkos::AUTO))>);

  // SES != DefaultExecutionSpace deduces to TeamPolicy<SES>

  static_assert(std::is_same_v<Kokkos::TeamPolicy<SomeExecutionSpace>,
                               decltype(Kokkos::TeamPolicy(ses, i, i))>);
  static_assert(std::is_same_v<Kokkos::TeamPolicy<SomeExecutionSpace>,
                               decltype(Kokkos::TeamPolicy(ses, i, i, i))>);
  static_assert(
      std::is_same_v<Kokkos::TeamPolicy<SomeExecutionSpace>,
                     decltype(Kokkos::TeamPolicy(ses, i, Kokkos::AUTO))>);
  static_assert(
      std::is_same_v<Kokkos::TeamPolicy<SomeExecutionSpace>,
                     decltype(Kokkos::TeamPolicy(ses, i, Kokkos::AUTO, i))>);
  static_assert(std::is_same_v<Kokkos::TeamPolicy<SomeExecutionSpace>,
                               decltype(Kokkos::TeamPolicy(ses, i, Kokkos::AUTO,
                                                           Kokkos::AUTO))>);
  static_assert(
      std::is_same_v<Kokkos::TeamPolicy<SomeExecutionSpace>,
                     decltype(Kokkos::TeamPolicy(ses, i, i, Kokkos::AUTO))>);
};

}  // namespace
