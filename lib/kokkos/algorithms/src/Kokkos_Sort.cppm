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

#include <Kokkos_NestedSort.hpp>
#include <Kokkos_Sort.hpp>

export module kokkos.sort;

export {
  namespace Kokkos {
  using ::Kokkos::BinOp1D;
  using ::Kokkos::BinOp3D;
  using ::Kokkos::BinSort;

  using ::Kokkos::sort;

  namespace Experimental {
  using ::Kokkos::Experimental::sort_by_key;
  using ::Kokkos::Experimental::sort_by_key_team;
  using ::Kokkos::Experimental::sort_by_key_thread;
  using ::Kokkos::Experimental::sort_team;
  using ::Kokkos::Experimental::sort_thread;
  }  // namespace Experimental
  }  // namespace Kokkos
}
