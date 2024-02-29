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

#include <Kokkos_Core.hpp>

#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_CPUDiscovery.hpp>
#include <impl/Kokkos_Tools.hpp>
#include <impl/Kokkos_ExecSpaceManager.hpp>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <thread>

namespace Kokkos {
namespace Impl {

void OpenMPInternal::acquire_lock() {
  while (1 == desul::atomic_compare_exchange(&m_pool_mutex, 0, 1,
                                             desul::MemoryOrderAcquire(),
                                             desul::MemoryScopeDevice())) {
    // do nothing
  }
}

void OpenMPInternal::release_lock() {
  desul::atomic_store(&m_pool_mutex, 0, desul::MemoryOrderRelease(),
                      desul::MemoryScopeDevice());
}

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
void OpenMPInternal::validate_partition_impl(const int nthreads,
                                             int &num_partitions,
                                             int &partition_size) {
  if (nthreads == 1) {
    num_partitions = 1;
    partition_size = 1;
  } else if (num_partitions < 1 && partition_size < 1) {
    int idle = nthreads;
    for (int np = 2; np <= nthreads; ++np) {
      for (int ps = 1; ps <= nthreads / np; ++ps) {
        if (nthreads - np * ps < idle) {
          idle           = nthreads - np * ps;
          num_partitions = np;
          partition_size = ps;
        }
        if (idle == 0) {
          break;
        }
      }
    }
  } else if (num_partitions < 1 && partition_size > 0) {
    if (partition_size <= nthreads) {
      num_partitions = nthreads / partition_size;
    } else {
      num_partitions = 1;
      partition_size = nthreads;
    }
  } else if (num_partitions > 0 && partition_size < 1) {
    if (num_partitions <= nthreads) {
      partition_size = nthreads / num_partitions;
    } else {
      num_partitions = nthreads;
      partition_size = 1;
    }
  } else if (num_partitions * partition_size > nthreads) {
    int idle     = nthreads;
    const int NP = num_partitions;
    const int PS = partition_size;
    for (int np = NP; np > 0; --np) {
      for (int ps = PS; ps > 0; --ps) {
        if ((np * ps <= nthreads) && (nthreads - np * ps < idle)) {
          idle           = nthreads - np * ps;
          num_partitions = np;
          partition_size = ps;
        }
        if (idle == 0) {
          break;
        }
      }
    }
  }
}
#endif

void OpenMPInternal::clear_thread_data() {
  const size_t member_bytes =
      sizeof(int64_t) *
      HostThreadTeamData::align_to_int64(sizeof(HostThreadTeamData));

  const int old_alloc_bytes =
      m_pool[0] ? (member_bytes + m_pool[0]->scratch_bytes()) : 0;

  OpenMP::memory_space space;

#pragma omp parallel num_threads(m_pool_size)
  {
    const int rank = omp_get_thread_num();

    if (nullptr != m_pool[rank]) {
      m_pool[rank]->disband_pool();

      space.deallocate(m_pool[rank], old_alloc_bytes);

      m_pool[rank] = nullptr;
    }
  }
  /* END #pragma omp parallel */
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

        space.deallocate(m_pool[rank], old_alloc_bytes);
      }

      void *ptr = nullptr;
      try {
        ptr = space.allocate(alloc_bytes);
      } catch (
          Kokkos::Experimental::RawMemoryAllocationFailure const &failure) {
        // For now, just rethrow the error message the existing way
        Kokkos::Impl::throw_runtime_exception(failure.get_error_message());
      }

      m_pool[rank] = new (ptr) HostThreadTeamData();

      m_pool[rank]->scratch_assign(((char *)ptr) + member_bytes, alloc_bytes,
                                   pool_reduce_bytes, team_reduce_bytes,
                                   team_shared_bytes, thread_local_bytes);
    }

    HostThreadTeamData::organize_pool(m_pool, m_pool_size);
  }
}

OpenMPInternal &OpenMPInternal::singleton() {
  static OpenMPInternal *self = nullptr;
  if (self == nullptr) {
    self = new OpenMPInternal(get_current_max_threads());
  }

  return *self;
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

void OpenMPInternal::initialize(int thread_count) {
  if (m_initialized) {
    Kokkos::abort(
        "Calling OpenMP::initialize after OpenMP::finalize is illegal\n");
  }

  if (omp_in_parallel()) {
    std::string msg("Kokkos::OpenMP::initialize ERROR : in parallel");
    Kokkos::Impl::throw_runtime_exception(msg);
  }

  {
    if (Kokkos::show_warnings() && !std::getenv("OMP_PROC_BIND")) {
      std::cerr
          << R"WARNING(Kokkos::OpenMP::initialize WARNING: OMP_PROC_BIND environment variable not set
  In general, for best performance with OpenMP 4.0 or better set OMP_PROC_BIND=spread and OMP_PLACES=threads
  For best performance with OpenMP 3.1 set OMP_PROC_BIND=true
  For unit testing set OMP_PROC_BIND=false
)WARNING" << std::endl;

      if (mpi_detected()) {
        std::cerr
            << R"WARNING(MPI detected: For OpenMP binding to work as intended, MPI ranks must be bound to exclusive CPU sets.
)WARNING" << std::endl;
      }
    }

    // Before any other call to OMP query the maximum number of threads
    // and save the value for re-initialization unit testing.

    Impl::g_openmp_hardware_max_threads = get_current_max_threads();

    int process_num_threads = Impl::g_openmp_hardware_max_threads;

    if (Kokkos::hwloc::available()) {
      process_num_threads = Kokkos::hwloc::get_available_numa_count() *
                            Kokkos::hwloc::get_available_cores_per_numa() *
                            Kokkos::hwloc::get_available_threads_per_core();
    }

    // if thread_count  < 0, use g_openmp_hardware_max_threads;
    // if thread_count == 0, set g_openmp_hardware_max_threads to
    // process_num_threads if thread_count  > 0, set
    // g_openmp_hardware_max_threads to thread_count
    if (thread_count < 0) {
      thread_count = Impl::g_openmp_hardware_max_threads;
    } else if (thread_count == 0) {
      if (Impl::g_openmp_hardware_max_threads != process_num_threads) {
        Impl::g_openmp_hardware_max_threads = process_num_threads;
        omp_set_num_threads(Impl::g_openmp_hardware_max_threads);
      }
    } else {
      if (Kokkos::show_warnings() && thread_count > process_num_threads) {
        std::cerr << "Kokkos::OpenMP::initialize WARNING: You are likely "
                     "oversubscribing your CPU cores.\n"
                  << "  process threads available : " << std::setw(3)
                  << process_num_threads
                  << ",  requested thread : " << std::setw(3) << thread_count
                  << std::endl;
      }
      Impl::g_openmp_hardware_max_threads = thread_count;
      omp_set_num_threads(Impl::g_openmp_hardware_max_threads);
    }

// setup thread local
#pragma omp parallel num_threads(Impl::g_openmp_hardware_max_threads)
    { Impl::SharedAllocationRecord<void, void>::tracking_enable(); }

    auto &instance       = OpenMPInternal::singleton();
    instance.m_pool_size = Impl::g_openmp_hardware_max_threads;

    // New, unified host thread team data:
    {
      size_t pool_reduce_bytes  = 32 * thread_count;
      size_t team_reduce_bytes  = 32 * thread_count;
      size_t team_shared_bytes  = 1024 * thread_count;
      size_t thread_local_bytes = 1024;

      instance.resize_thread_data(pool_reduce_bytes, team_reduce_bytes,
                                  team_shared_bytes, thread_local_bytes);
    }
  }

  // Check for over-subscription
  auto const reported_ranks = mpi_ranks_per_node();
  auto const mpi_local_size = reported_ranks < 0 ? 1 : reported_ranks;
  int const procs_per_node  = std::thread::hardware_concurrency();
  if (Kokkos::show_warnings() &&
      (mpi_local_size * long(thread_count) > procs_per_node)) {
    std::cerr << "Kokkos::OpenMP::initialize WARNING: You are likely "
                 "oversubscribing your CPU cores."
              << std::endl;
    std::cerr << "                                    Detected: "
              << procs_per_node << " cores per node." << std::endl;
    std::cerr << "                                    Detected: "
              << mpi_local_size << " MPI_ranks per node." << std::endl;
    std::cerr << "                                    Requested: "
              << thread_count << " threads per process." << std::endl;
  }

  m_initialized = true;
}

void OpenMPInternal::finalize() {
  if (omp_in_parallel()) {
    std::string msg("Kokkos::OpenMP::finalize ERROR ");
    if (this != &singleton()) msg.append(": not initialized");
    if (omp_in_parallel()) msg.append(": in parallel");
    Kokkos::Impl::throw_runtime_exception(msg);
  }

  if (this == &singleton()) {
    auto const &instance = singleton();
    // Silence Cuda Warning
    const int nthreads =
        instance.m_pool_size <= Impl::g_openmp_hardware_max_threads
            ? Impl::g_openmp_hardware_max_threads
            : instance.m_pool_size;
    (void)nthreads;

#pragma omp parallel num_threads(nthreads)
    { Impl::SharedAllocationRecord<void, void>::tracking_disable(); }

    // allow main thread to track
    Impl::SharedAllocationRecord<void, void>::tracking_enable();

    Impl::g_openmp_hardware_max_threads = 1;
  }

  m_initialized = false;

  Kokkos::Profiling::finalize();
}

void OpenMPInternal::print_configuration(std::ostream &s) const {
  s << "Kokkos::OpenMP";

  if (m_initialized) {
    const int numa_count      = 1;
    const int core_per_numa   = Impl::g_openmp_hardware_max_threads;
    const int thread_per_core = 1;

    s << " thread_pool_topology[ " << numa_count << " x " << core_per_numa
      << " x " << thread_per_core << " ]" << std::endl;
  } else {
    s << " not initialized" << std::endl;
  }
}

bool OpenMPInternal::verify_is_initialized(const char *const label) const {
  if (!m_initialized) {
    std::cerr << "Kokkos::OpenMP " << label
              << " : ERROR OpenMP is not initialized" << std::endl;
  }
  return m_initialized;
}
}  // namespace Impl
}  // namespace Kokkos
