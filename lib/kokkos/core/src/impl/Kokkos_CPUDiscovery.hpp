// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
namespace Kokkos {
namespace Impl {

int mpi_ranks_per_node();
int mpi_local_rank_on_node();
// returns true if MPI execution environment is detected, false otherwise.
bool mpi_detected();

}  // namespace Impl
}  // namespace Kokkos
