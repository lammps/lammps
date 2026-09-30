// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <iostream>

#include <OpenMP/Kokkos_OpenMP.hpp>
#include <OpenMP/Kokkos_OpenMP_Instance.hpp>

#include <impl/Kokkos_CheckUsage.hpp>
#include <impl/Kokkos_ExecSpaceManager.hpp>
#include <impl/Kokkos_CPUDiscovery.hpp>
#include <Kokkos_hwloc.hpp>
#include <thread>
#include <iomanip>

namespace Kokkos {

OpenMP::~OpenMP() {
  Impl::check_execution_space_destructor_precondition(name());
}

OpenMP::OpenMP()
    : m_space_instance(
          (Impl::check_execution_space_constructor_precondition(name()),
           Impl::OpenMPInternal::default_instance)) {}

OpenMP::OpenMP(int pool_size)
    : m_space_instance(
          (Impl::check_execution_space_constructor_precondition(name()),
           Impl::HostSharedPtr(new Impl::OpenMPInternal(pool_size)))) {}

int OpenMP::impl_get_current_max_threads() noexcept {
  return Impl::OpenMPInternal::get_current_max_threads();
}

void OpenMP::impl_initialize(InitializationSettings const &settings) {
  int thread_count =
      settings.has_num_threads() ? settings.get_num_threads() : -1;
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

      if (Impl::mpi_detected()) {
        std::cerr
            << R"WARNING(MPI detected: For OpenMP binding to work as intended, MPI ranks must be bound to exclusive CPU sets.
)WARNING" << std::endl;
      }
    }

    // Before any other call to OMP query the maximum number of threads
    // and save the value for re-initialization unit testing.

    Impl::OpenMPInternal::hardware_max_threads =
        Impl::OpenMPInternal::get_current_max_threads();

    int process_num_threads = Impl::OpenMPInternal::hardware_max_threads;

    if (Kokkos::hwloc::available()) {
      process_num_threads = Kokkos::hwloc::get_available_numa_count() *
                            Kokkos::hwloc::get_available_cores_per_numa() *
                            Kokkos::hwloc::get_available_threads_per_core();
    }

    // if thread_count  < 0, use hardware_max_threads;
    // if thread_count == 0, set hardware_max_threads to
    // process_num_threads if thread_count  > 0, set
    // hardware_max_threads to thread_count
    if (thread_count < 0) {
      thread_count = Impl::OpenMPInternal::hardware_max_threads;
    } else if (thread_count == 0) {
      if (Impl::OpenMPInternal::hardware_max_threads != process_num_threads) {
        Impl::OpenMPInternal::hardware_max_threads = process_num_threads;
        omp_set_num_threads(Impl::OpenMPInternal::hardware_max_threads);
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
      Impl::OpenMPInternal::hardware_max_threads = thread_count;
      omp_set_num_threads(Impl::OpenMPInternal::hardware_max_threads);
    }

// setup thread local
#pragma omp parallel num_threads(Impl::OpenMPInternal::hardware_max_threads)
    { Impl::SharedAllocationRecord<void, void>::tracking_enable(); }
  }

  // Create the default instance.
  Impl::OpenMPInternal::default_instance =
      Impl::HostSharedPtr<Impl::OpenMPInternal>(
          new Impl::OpenMPInternal(Impl::OpenMPInternal::hardware_max_threads));

  // Check for over-subscription
  auto const reported_ranks = Impl::mpi_ranks_per_node();
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
}

void OpenMP::impl_finalize() {
  // Destroy the default instance.
  Impl::OpenMPInternal::default_instance = nullptr;
}

void OpenMP::print_configuration(std::ostream &os, bool /*verbose*/) const {
  os << "Host Parallel Execution Space:\n";
  os << "  KOKKOS_ENABLE_OPENMP: yes\n";

  os << "\nOpenMP Runtime Configuration:\n";

  m_space_instance->print_configuration(os);
}

int OpenMP::concurrency() const { return impl_thread_pool_size(); }

void OpenMP::impl_static_fence(std::string const &name) {
  Kokkos::Tools::Experimental::Impl::profile_fence_event<Kokkos::OpenMP>(
      name,
      Kokkos::Tools::Experimental::SpecialSynchronizationCases::
          GlobalDeviceSynchronization,
      []() {
        std::lock_guard<std::mutex> lock_all_instances(
            Impl::OpenMPInternal::all_instances_mutex);
        for (auto *instance_ptr : Impl::OpenMPInternal::all_instances) {
          std::lock_guard<std::mutex> lock_instance(
              instance_ptr->m_instance_mutex);
        }
      });
}

void OpenMP::fence(const std::string &name) const {
  impl_internal_space_instance()->fence(name);
}

int OpenMP::impl_thread_pool_size() const noexcept {
  return (impl_internal_space_instance()->get_level() < omp_get_level())
             ? omp_get_num_threads()
             : impl_internal_space_instance()->m_pool_size;
}

int OpenMP::impl_max_hardware_threads() noexcept {
  return Impl::OpenMPInternal::max_hardware_threads();
}

namespace Impl {

int g_openmp_space_factory_initialized =
    initialize_space_factory<OpenMP>("050_OpenMP");

}  // namespace Impl

}  // namespace Kokkos
