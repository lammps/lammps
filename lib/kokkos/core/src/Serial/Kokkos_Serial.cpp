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

#include <Serial/Kokkos_Serial.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_ExecSpaceManager.hpp>
#include <impl/Kokkos_SharedAlloc.hpp>

#include <cstdlib>
#include <iostream>
#include <sstream>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

bool SerialInternal::is_initialized() { return m_is_initialized; }

void SerialInternal::initialize() {
  if (is_initialized()) return;

  Impl::SharedAllocationRecord<void, void>::tracking_enable();

  m_is_initialized = true;
}

void SerialInternal::finalize() {
  if (m_thread_team_data.scratch_buffer()) {
    m_thread_team_data.disband_team();
    m_thread_team_data.disband_pool();

    Kokkos::HostSpace space;

    space.deallocate(m_thread_team_data.scratch_buffer(),
                     m_thread_team_data.scratch_bytes());

    m_thread_team_data.scratch_assign(nullptr, 0, 0, 0, 0, 0);
  }

  m_is_initialized = false;
}

SerialInternal& SerialInternal::singleton() {
  static SerialInternal* self = nullptr;
  if (!self) {
    self = new SerialInternal();
  }
  return *self;
}

// Resize thread team data scratch memory
void SerialInternal::resize_thread_team_data(size_t pool_reduce_bytes,
                                             size_t team_reduce_bytes,
                                             size_t team_shared_bytes,
                                             size_t thread_local_bytes) {
  if (pool_reduce_bytes < 512) pool_reduce_bytes = 512;
  if (team_reduce_bytes < 512) team_reduce_bytes = 512;

  const size_t old_pool_reduce  = m_thread_team_data.pool_reduce_bytes();
  const size_t old_team_reduce  = m_thread_team_data.team_reduce_bytes();
  const size_t old_team_shared  = m_thread_team_data.team_shared_bytes();
  const size_t old_thread_local = m_thread_team_data.thread_local_bytes();
  const size_t old_alloc_bytes  = m_thread_team_data.scratch_bytes();

  // Allocate if any of the old allocation is tool small:

  const bool allocate = (old_pool_reduce < pool_reduce_bytes) ||
                        (old_team_reduce < team_reduce_bytes) ||
                        (old_team_shared < team_shared_bytes) ||
                        (old_thread_local < thread_local_bytes);

  if (allocate) {
    Kokkos::HostSpace space;

    if (old_alloc_bytes) {
      m_thread_team_data.disband_team();
      m_thread_team_data.disband_pool();

      space.deallocate("Kokkos::Serial::scratch_mem",
                       m_thread_team_data.scratch_buffer(),
                       m_thread_team_data.scratch_bytes());
    }

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
        HostThreadTeamData::scratch_size(pool_reduce_bytes, team_reduce_bytes,
                                         team_shared_bytes, thread_local_bytes);

    void* ptr = nullptr;
    try {
      ptr = space.allocate("Kokkos::Serial::scratch_mem", alloc_bytes);
    } catch (Kokkos::Experimental::RawMemoryAllocationFailure const& failure) {
      // For now, just rethrow the error message the existing way
      Kokkos::Impl::throw_runtime_exception(failure.get_error_message());
    }

    m_thread_team_data.scratch_assign(static_cast<char*>(ptr), alloc_bytes,
                                      pool_reduce_bytes, team_reduce_bytes,
                                      team_shared_bytes, thread_local_bytes);

    HostThreadTeamData* pool[1] = {&m_thread_team_data};

    m_thread_team_data.organize_pool(pool, 1);
    m_thread_team_data.organize_team(1);
  }
}
}  // namespace Impl

Serial::Serial()
    : m_space_instance(&Impl::SerialInternal::singleton(),
                       [](Impl::SerialInternal*) {}) {}

Serial::Serial(NewInstance)
    : m_space_instance(new Impl::SerialInternal, [](Impl::SerialInternal* ptr) {
        ptr->finalize();
        delete ptr;
      }) {}

void Serial::print_configuration(std::ostream& os, bool /*verbose*/) const {
  os << "Host Serial Execution Space:\n";
  os << "  KOKKOS_ENABLE_SERIAL: yes\n";

#ifdef KOKKOS_ENABLE_ATOMICS_BYPASS
  os << "Kokkos atomics disabled\n";
#endif

  os << "\nSerial Runtime Configuration:\n";
}

bool Serial::impl_is_initialized() {
  return Impl::SerialInternal::singleton().is_initialized();
}

void Serial::impl_initialize(InitializationSettings const&) {
  Impl::SerialInternal::singleton().initialize();
}

void Serial::impl_finalize() { Impl::SerialInternal::singleton().finalize(); }

const char* Serial::name() { return "Serial"; }

namespace Impl {

int g_serial_space_factory_initialized =
    initialize_space_factory<Serial>("100_Serial");

}  // namespace Impl

}  // namespace Kokkos
