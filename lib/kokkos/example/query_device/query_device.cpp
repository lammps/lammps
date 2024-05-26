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

#include <iostream>
#include <sstream>

#include <Kokkos_Macros.hpp>

//#define USE_MPI
#if defined(USE_MPI)
#include <mpi.h>
#endif

#include <Kokkos_Core.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

int main(int argc, char** argv) {
  std::ostringstream msg;

  (void)argc;
  (void)argv;
#if defined(USE_MPI)

  MPI_Init(&argc, &argv);

  int mpi_rank = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  msg << "MPI rank(" << mpi_rank << ") ";

#endif
  Kokkos::initialize(argc, argv);
  msg << "{" << std::endl;

  if (Kokkos::hwloc::available()) {
    msg << "hwloc( NUMA[" << Kokkos::hwloc::get_available_numa_count()
        << "] x CORE[" << Kokkos::hwloc::get_available_cores_per_numa()
        << "] x HT[" << Kokkos::hwloc::get_available_threads_per_core() << "] )"
        << std::endl;
  }

  Kokkos::print_configuration(msg);

  msg << "}" << std::endl;

  std::cout << msg.str();
  Kokkos::finalize();
#if defined(USE_MPI)

  MPI_Finalize();

#endif

  return 0;
}
