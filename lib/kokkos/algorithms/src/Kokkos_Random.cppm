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

module;

#include <Kokkos_Random.hpp>

export module kokkos.random;

export {
  namespace Kokkos {
  using ::Kokkos::fill_random;
  using ::Kokkos::rand;
  using ::Kokkos::Random_XorShift1024_Pool;
  using ::Kokkos::Random_XorShift64_Pool;
  }  // namespace Kokkos
}
