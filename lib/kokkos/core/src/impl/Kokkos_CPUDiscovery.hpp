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
namespace Kokkos {
namespace Impl {

int mpi_ranks_per_node();
int mpi_local_rank_on_node();
// returns true if MPI execution environment is detected, false otherwise.
bool mpi_detected();

}  // namespace Impl
}  // namespace Kokkos
