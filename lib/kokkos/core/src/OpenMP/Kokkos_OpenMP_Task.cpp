/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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
