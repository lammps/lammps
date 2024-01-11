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
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif
#ifndef KOKKOS_HWLOC_HPP
#define KOKKOS_HWLOC_HPP

#include <Kokkos_Macros.hpp>

#include <utility>

namespace Kokkos {

/** \brief  Minimal subset of logical 'hwloc' functionality available
 *          from http://www.open-mpi.org/projects/hwloc/.
 *
 *  The calls are NOT thread safe in order to avoid mutexes,
 *  memory allocations, or other actions which could give the
 *  runtime system an opportunity to migrate the threads or
 *  touch allocated memory during the function calls.
 *
 *  All calls to these functions should be performed by a thread
 *  when it has guaranteed exclusive access; e.g., for OpenMP
 *  within a 'critical' region.
 */
namespace hwloc {

/** \brief  Query if hwloc is available */
bool available();

/** \brief  Query number of available NUMA regions.
 *          This will be less than the hardware capacity
 *          if the MPI process is pinned to a NUMA region.
 */
unsigned get_available_numa_count();

/** \brief  Query number of available cores per NUMA regions.
 *          This will be less than the hardware capacity
 *          if the MPI process is pinned to a set of cores.
 */
unsigned get_available_cores_per_numa();

/** \brief  Query number of available "hard" threads per core; i.e.,
 * hyperthreads */
unsigned get_available_threads_per_core();

} /* namespace hwloc */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Internal functions for binding persistent spawned threads.

namespace Kokkos {
namespace hwloc {

/** \brief  Recommend mapping of threads onto cores.
 *
 * If thread_count == 0 then choose and set a value.
 * If use_numa_count == 0 then choose and set a value.
 * If use_cores_per_numa == 0 then choose and set a value.
 *
 * Return 0 if asynchronous,
 * Return 1 if synchronous and threads_coord[0] is process core
 */
unsigned thread_mapping(const char* const label, const bool allow_async,
                        unsigned& thread_count, unsigned& use_numa_count,
                        unsigned& use_cores_per_numa,
                        std::pair<unsigned, unsigned> threads_coord[]);

/** \brief  Query core-coordinate of the current thread
 *          with respect to the core_topology.
 *
 *  As long as the thread is running within the
 *  process binding the following condition holds.
 *
 *  core_coordinate.first  < core_topology.first
 *  core_coordinate.second < core_topology.second
 */
std::pair<unsigned, unsigned> get_this_thread_coordinate();

/** \brief  Bind the current thread to a core. */
bool bind_this_thread(const std::pair<unsigned, unsigned>);

/** \brief Can hwloc bind threads? */
bool can_bind_threads();

/** \brief  Bind the current thread to one of the cores in the list.
 *          Set that entry to (~0,~0) and return the index.
 *          If binding fails return ~0.
 */
unsigned bind_this_thread(const unsigned coordinate_count,
                          std::pair<unsigned, unsigned> coordinate[]);

/** \brief  Unbind the current thread back to the original process binding */
bool unbind_this_thread();

} /* namespace hwloc */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #define KOKKOS_HWLOC_HPP */
