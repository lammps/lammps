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

#include <OpenMP/Kokkos_OpenMP.hpp>
#include <OpenMP/Kokkos_OpenMP_Instance.hpp>

#include <impl/Kokkos_ExecSpaceManager.hpp>

namespace Kokkos {

OpenMP::OpenMP()
    : m_space_instance(&Impl::OpenMPInternal::singleton(),
                       [](Impl::OpenMPInternal *) {}) {
  Impl::OpenMPInternal::singleton().verify_is_initialized(
      "OpenMP instance constructor");
}

OpenMP::OpenMP(int pool_size)
    : m_space_instance(new Impl::OpenMPInternal(pool_size),
                       [](Impl::OpenMPInternal *ptr) {
                         ptr->finalize();
                         delete ptr;
                       }) {
  Impl::OpenMPInternal::singleton().verify_is_initialized(
      "OpenMP instance constructor");
}

int OpenMP::impl_get_current_max_threads() noexcept {
  return Impl::OpenMPInternal::get_current_max_threads();
}

void OpenMP::impl_initialize(InitializationSettings const &settings) {
  Impl::OpenMPInternal::singleton().initialize(
      settings.has_num_threads() ? settings.get_num_threads() : -1);
}

void OpenMP::impl_finalize() { Impl::OpenMPInternal::singleton().finalize(); }

void OpenMP::print_configuration(std::ostream &os, bool /*verbose*/) const {
  os << "Host Parallel Execution Space:\n";
  os << "  KOKKOS_ENABLE_OPENMP: yes\n";

  os << "\nOpenMP Runtime Configuration:\n";

  m_space_instance->print_configuration(os);
}

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
int OpenMP::concurrency(OpenMP const &instance) {
  return instance.impl_thread_pool_size();
}
#else
int OpenMP::concurrency() const { return impl_thread_pool_size(); }
#endif

void OpenMP::fence(const std::string &name) const {
  Kokkos::Tools::Experimental::Impl::profile_fence_event<Kokkos::OpenMP>(
      name, Kokkos::Tools::Experimental::Impl::DirectFenceIDHandle{1}, []() {});
}

bool OpenMP::impl_is_initialized() noexcept {
  return Impl::OpenMPInternal::singleton().is_initialized();
}

bool OpenMP::in_parallel(OpenMP const &exec_space) noexcept {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
  return (
      (exec_space.impl_internal_space_instance()->m_level < omp_get_level()) &&
      (!Impl::t_openmp_instance ||
       Impl::t_openmp_instance->m_level < omp_get_level()));
#else
  return exec_space.impl_internal_space_instance()->m_level < omp_get_level();
#endif
}

int OpenMP::impl_thread_pool_size() const noexcept {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
  return OpenMP::in_parallel(*this)
             ? omp_get_num_threads()
             : (Impl::t_openmp_instance
                    ? Impl::t_openmp_instance->m_pool_size
                    : impl_internal_space_instance()->m_pool_size);
#else
  return OpenMP::in_parallel(*this)
             ? omp_get_num_threads()
             : impl_internal_space_instance()->m_pool_size;
#endif
}

int OpenMP::impl_max_hardware_threads() noexcept {
  return Impl::g_openmp_hardware_max_threads;
}

namespace Impl {

int g_openmp_space_factory_initialized =
    initialize_space_factory<OpenMP>("050_OpenMP");

}  // namespace Impl

}  // namespace Kokkos
