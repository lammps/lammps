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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <impl/Kokkos_CPUDiscovery.hpp>

#include <cstdlib>  // getenv
#include <string>

int Kokkos::Impl::mpi_ranks_per_node() {
  for (char const* env_var : {
           "OMPI_COMM_WORLD_LOCAL_SIZE",  // OpenMPI
           "MV2_COMM_WORLD_LOCAL_SIZE",   // MVAPICH2
           "MPI_LOCALNRANKS",             // MPICH
                                          // SLURM???
           "PMI_LOCAL_SIZE",              // PMI
       }) {
    char const* str = std::getenv(env_var);
    if (str) {
      return std::stoi(str);
    }
  }
  return -1;
}

int Kokkos::Impl::mpi_local_rank_on_node() {
  for (char const* env_var : {
           "OMPI_COMM_WORLD_LOCAL_RANK",  // OpenMPI
           "MV2_COMM_WORLD_LOCAL_RANK",   // MVAPICH2
           "MPI_LOCALRANKID",             // MPICH
           "SLURM_LOCALID",               // SLURM
           "PMI_LOCAL_RANK",              // PMI
       }) {
    char const* str = std::getenv(env_var);
    if (str) {
      return std::stoi(str);
    }
  }
  return -1;
}

bool Kokkos::Impl::mpi_detected() { return mpi_local_rank_on_node() != -1; }
