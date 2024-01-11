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

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_OPENMP) && defined(KOKKOS_ENABLE_TASKDAG)

#include <Kokkos_Core.hpp>

#include <impl/Kokkos_TaskQueue_impl.hpp>
#include <impl/Kokkos_HostThreadTeam.hpp>
#include <OpenMP/Kokkos_OpenMP_Task.hpp>
#include <cassert>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template class TaskQueue<Kokkos::OpenMP, typename Kokkos::OpenMP::memory_space>;

HostThreadTeamData& HostThreadTeamDataSingleton::singleton() {
  static HostThreadTeamDataSingleton s;
  return s;
}

HostThreadTeamDataSingleton::HostThreadTeamDataSingleton()
    : HostThreadTeamData() {
  Kokkos::OpenMP::memory_space space;
  const size_t num_pool_reduce_bytes  = 32;
  const size_t num_team_reduce_bytes  = 32;
  const size_t num_team_shared_bytes  = 1024;
  const size_t num_thread_local_bytes = 1024;
  const size_t alloc_bytes            = HostThreadTeamData::scratch_size(
      num_pool_reduce_bytes, num_team_reduce_bytes, num_team_shared_bytes,
      num_thread_local_bytes);

  void* ptr = nullptr;
  try {
    ptr = space.allocate(alloc_bytes);
  } catch (Kokkos::Experimental::RawMemoryAllocationFailure const& f) {
    // For now, just rethrow the error message with a note
    // Note that this could, in turn, trigger an out of memory exception,
    // but it's pretty unlikely, so we won't worry about it for now.
    // TODO reasonable error message when `std::string` causes OOM error
    Kokkos::Impl::throw_runtime_exception(
        std::string("Failure to allocate scratch memory:  ") +
        f.get_error_message());
  }

  HostThreadTeamData::scratch_assign(
      ptr, alloc_bytes, num_pool_reduce_bytes, num_team_reduce_bytes,
      num_team_shared_bytes, num_thread_local_bytes);
}

HostThreadTeamDataSingleton::~HostThreadTeamDataSingleton() {
  Kokkos::OpenMP::memory_space space;
  space.deallocate(HostThreadTeamData::scratch_buffer(),
                   static_cast<size_t>(HostThreadTeamData::scratch_bytes()));
}

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
#else
void KOKKOS_CORE_SRC_OPENMP_KOKKOS_OPENMP_TASK_PREVENT_LINK_ERROR() {}
#endif /* #if defined( KOKKOS_ENABLE_OPENMP ) && defined( \
          KOKKOS_ENABLE_TASKDAG ) */
