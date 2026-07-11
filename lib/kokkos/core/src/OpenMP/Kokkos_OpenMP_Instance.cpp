// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif

#include <OpenMP/Kokkos_OpenMP_Instance.hpp>
#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_CPUDiscovery.hpp>
#include <impl/Kokkos_Tools.hpp>
#include <impl/Kokkos_ExecSpaceManager.hpp>

#include <cstdlib>
#include <iostream>
#include <new>
#include <sstream>
#include <thread>

namespace Kokkos {
namespace Impl {

std::vector<OpenMPInternal *> OpenMPInternal::all_instances;
std::mutex OpenMPInternal::all_instances_mutex;
HostSharedPtr<OpenMPInternal> OpenMPInternal::default_instance;
int OpenMPInternal::hardware_max_threads;

int OpenMPInternal::max_hardware_threads() noexcept {
  return hardware_max_threads;
}

void OpenMPInternal::clear_thread_data() {
  const size_t member_bytes =
      sizeof(int64_t) *
      HostThreadTeamData::align_to_int64(sizeof(HostThreadTeamData));

  const int old_alloc_bytes =
      m_pool[0] ? (member_bytes + m_pool[0]->scratch_bytes()) : 0;

  OpenMP::memory_space space;

  for (int rank = 0; rank < m_pool_size; ++rank) {
    if (nullptr != m_pool[rank]) {
      m_pool[rank]->disband_pool();

      m_pool[rank]->~HostThreadTeamData();

      space.deallocate(m_pool[rank], old_alloc_bytes);

      m_pool[rank] = nullptr;
    }
  }
}

void OpenMPInternal::resize_thread_data(size_t pool_reduce_bytes,
                                        size_t team_reduce_bytes,
                                        size_t team_shared_bytes,
                                        size_t thread_local_bytes) {
  const size_t member_bytes =
      sizeof(int64_t) *
      HostThreadTeamData::align_to_int64(sizeof(HostThreadTeamData));

  HostThreadTeamData *root = m_pool[0];

  const size_t old_pool_reduce  = root ? root->pool_reduce_bytes() : 0;
  const size_t old_team_reduce  = root ? root->team_reduce_bytes() : 0;
  const size_t old_team_shared  = root ? root->team_shared_bytes() : 0;
  const size_t old_thread_local = root ? root->thread_local_bytes() : 0;
  const size_t old_alloc_bytes =
      root ? (member_bytes + root->scratch_bytes()) : 0;

  // Allocate if any of the old allocation is tool small:

  const bool allocate = (old_pool_reduce < pool_reduce_bytes) ||
                        (old_team_reduce < team_reduce_bytes) ||
                        (old_team_shared < team_shared_bytes) ||
                        (old_thread_local < thread_local_bytes);

  if (allocate) {
    if (pool_reduce_bytes < old_pool_reduce) {
      pool_reduce_bytes = old_pool_reduce;
    }
    if (team_reduce_bytes < old_team_reduce) {
      team_reduce_bytes = old_team_reduce;
    }
    if (team_shared_bytes < old_team_shared) {
      team_shared_bytes = old_team_shared;
    }
    if (thread_local_bytes < old_thread_local) {
      thread_local_bytes = old_thread_local;
    }

    const size_t alloc_bytes =
        member_bytes +
        HostThreadTeamData::scratch_size(pool_reduce_bytes, team_reduce_bytes,
                                         team_shared_bytes, thread_local_bytes);

    OpenMP::memory_space space;

    memory_fence();

    for (int rank = 0; rank < m_pool_size; ++rank) {
      if (nullptr != m_pool[rank]) {
        m_pool[rank]->disband_pool();

        // impl_deallocate to not fence here
        space.impl_deallocate("[unlabeled]", m_pool[rank], old_alloc_bytes);
      }

      void *ptr = space.allocate("Kokkos::OpenMP::scratch_mem", alloc_bytes);

      m_pool[rank] = new (ptr) HostThreadTeamData();

      m_pool[rank]->scratch_assign(static_cast<char *>(ptr) + member_bytes,
                                   alloc_bytes, pool_reduce_bytes,
                                   team_reduce_bytes, team_shared_bytes,
                                   thread_local_bytes);
    }

    HostThreadTeamData::organize_pool(m_pool, m_pool_size);
  }
}

int OpenMPInternal::get_current_max_threads() noexcept {
  // Using omp_get_max_threads(); is problematic in conjunction with
  // Hwloc on Intel (essentially an initial call to the OpenMP runtime
  // without a parallel region before will set a process mask for a single core
  // The runtime will than bind threads for a parallel region to other cores on
  // the entering the first parallel region and make the process mask the
  // aggregate of the thread masks. The intend seems to be to make serial code
  // run fast, if you compile with OpenMP enabled but don't actually use
  // parallel regions or so static int omp_max_threads = omp_get_max_threads();

  int count = 0;
#pragma omp parallel
  {
#pragma omp atomic
    ++count;
  }
  return count;
}

OpenMPInternal::OpenMPInternal(int arg_pool_size)
    : m_pool_size{arg_pool_size}, m_level{omp_get_level()}, m_pool() {
  // guard pushing to all_instances
  {
    if (omp_get_level() != 0) {
      constexpr char msg[] =
          "Kokkos::OpenMP instances can only be created outside OpenMP "
          "regions!";
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
      if (Kokkos::show_warnings()) {
        std::cerr << msg << '\n';
      }
#else
      Kokkos::abort(msg);
#endif
    }
    std::scoped_lock lock(all_instances_mutex);
    all_instances.push_back(this);
  }

  // New, unified host thread team data:
  {
    size_t pool_reduce_bytes  = static_cast<size_t>(32) * arg_pool_size;
    size_t team_reduce_bytes  = static_cast<size_t>(32) * arg_pool_size;
    size_t team_shared_bytes  = static_cast<size_t>(1024) * arg_pool_size;
    size_t thread_local_bytes = 1024;

    resize_thread_data(pool_reduce_bytes, team_reduce_bytes, team_shared_bytes,
                       thread_local_bytes);
  }
}

void OpenMPInternal::fence(const std::string &name) {
  Kokkos::Tools::Experimental::Impl::profile_fence_event<Kokkos::OpenMP>(
      name, Kokkos::Tools::Experimental::Impl::DirectFenceIDHandle{1},
      [this]() { std::lock_guard<std::mutex> lock(m_instance_mutex); });
}

OpenMPInternal::~OpenMPInternal() {
  if (omp_in_parallel()) {
    std::string msg("Kokkos::OpenMP::finalize ERROR : in parallel");
    Kokkos::Impl::throw_runtime_exception(msg);
  }

  fence("Kokkos::OpenMPInternal: fence on destruction");

  // guard erasing from all_instances
  {
    std::scoped_lock lock(all_instances_mutex);

    auto it = std::find(all_instances.begin(), all_instances.end(), this);
    if (it == all_instances.end())
      Kokkos::abort(
          "Execution space instance to be removed couldn't be found!");
    *it = all_instances.back();
    all_instances.pop_back();
  }

  clear_thread_data();
}

void OpenMPInternal::print_configuration(std::ostream &s) const {
  s << "Kokkos::OpenMP";

  const int numa_count      = 1;
  const int core_per_numa   = hardware_max_threads;
  const int thread_per_core = 1;

  s << " thread_pool_topology[ " << numa_count << " x " << core_per_numa
    << " x " << thread_per_core << " ]" << std::endl;
}

}  // namespace Impl
}  // namespace Kokkos
