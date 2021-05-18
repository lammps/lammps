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

#include <stdio.h>
#include <limits>
#include <iostream>
#include <vector>
#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Error.hpp>
#include <iostream>
#include <impl/Kokkos_CPUDiscovery.hpp>
#include <impl/Kokkos_Tools.hpp>

#ifdef KOKKOS_ENABLE_OPENMPTARGET

// FIXME_OPENMPTARGET currently unused
/*
namespace Kokkos {
namespace Impl {
namespace {

KOKKOS_INLINE_FUNCTION
int kokkos_omp_in_parallel();

KOKKOS_INLINE_FUNCTION
int kokkos_omp_in_parallel() { return omp_in_parallel(); }

bool s_using_hwloc = false;

}  // namespace
}  // namespace Impl
}  // namespace Kokkos
*/

namespace Kokkos {
namespace Impl {

void OpenMPTargetExec::verify_is_process(const char* const label) {
  if (omp_in_parallel()) {
    std::string msg(label);
    msg.append(" ERROR: in parallel");
    Kokkos::Impl::throw_runtime_exception(msg);
  }
}

void OpenMPTargetExec::verify_initialized(const char* const label) {
  if (0 == Kokkos::Experimental::OpenMPTarget().impl_is_initialized()) {
    std::string msg(label);
    msg.append(" ERROR: not initialized");
    Kokkos::Impl::throw_runtime_exception(msg);
  }
}

void* OpenMPTargetExec::m_scratch_ptr         = nullptr;
int64_t OpenMPTargetExec::m_scratch_size      = 0;
int* OpenMPTargetExec::m_lock_array           = nullptr;
int64_t OpenMPTargetExec::m_lock_size         = 0;
uint32_t* OpenMPTargetExec::m_uniquetoken_ptr = nullptr;

void OpenMPTargetExec::clear_scratch() {
  Kokkos::Experimental::OpenMPTargetSpace space;
  space.deallocate(m_scratch_ptr, m_scratch_size);
  m_scratch_ptr  = nullptr;
  m_scratch_size = 0;
}

void OpenMPTargetExec::clear_lock_array() {
  if (m_lock_array != nullptr) {
    Kokkos::Experimental::OpenMPTargetSpace space;
    space.deallocate(m_lock_array, m_lock_size);
    m_lock_array = nullptr;
    m_lock_size  = 0;
  }
}

void* OpenMPTargetExec::get_scratch_ptr() { return m_scratch_ptr; }

void OpenMPTargetExec::resize_scratch(int64_t team_size, int64_t shmem_size_L0,
                                      int64_t shmem_size_L1) {
  Kokkos::Experimental::OpenMPTargetSpace space;
  const int64_t shmem_size =
      shmem_size_L0 + shmem_size_L1;  // L0 + L1 scratch memory per team.
  const int64_t padding = shmem_size * 10 / 100;  // Padding per team.
  // Total amount of scratch memory allocated is depenedent
  // on the maximum number of in-flight teams possible.
  int64_t total_size =
      (shmem_size + OpenMPTargetExecTeamMember::TEAM_REDUCE_SIZE + padding) *
      (MAX_ACTIVE_THREADS / team_size);

  if (total_size > m_scratch_size) {
    space.deallocate(m_scratch_ptr, m_scratch_size);
    m_scratch_size = total_size;
    m_scratch_ptr  = space.allocate(total_size);
  }
}

int* OpenMPTargetExec::get_lock_array(int num_teams) {
  Kokkos::Experimental::OpenMPTargetSpace space;
  int max_active_league_size = MAX_ACTIVE_THREADS / 32;
  int lock_array_elem =
      (num_teams > max_active_league_size) ? num_teams : max_active_league_size;
  if (m_lock_size < (lock_array_elem * sizeof(int))) {
    space.deallocate(m_lock_array, m_lock_size);
    m_lock_size  = lock_array_elem * sizeof(int);
    m_lock_array = static_cast<int*>(space.allocate(m_lock_size));

    // FIXME_OPENMPTARGET - Creating a target region here to initialize the
    // lock_array with 0's fails. Hence creating an equivalent host array to
    // achieve the same. Value of host array are then copied to the lock_array.
    int* h_lock_array = static_cast<int*>(
        omp_target_alloc(m_lock_size, omp_get_initial_device()));

    for (int i = 0; i < lock_array_elem; ++i) h_lock_array[i] = 0;

    OMPT_SAFE_CALL(omp_target_memcpy(m_lock_array, h_lock_array, m_lock_size, 0,
                                     0, omp_get_default_device(),
                                     omp_get_initial_device()));

    omp_target_free(h_lock_array, omp_get_initial_device());
  }

  return m_lock_array;
}

}  // namespace Impl
}  // namespace Kokkos

#endif  // KOKKOS_ENABLE_OPENMPTARGET
