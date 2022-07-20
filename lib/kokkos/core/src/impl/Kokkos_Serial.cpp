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

#include <Kokkos_Core.hpp>
#if defined(KOKKOS_ENABLE_SERIAL)

#include <cstdlib>
#include <sstream>
#include <Kokkos_Serial.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Error.hpp>

#include <impl/Kokkos_SharedAlloc.hpp>
#include <sstream>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

bool SerialInternal::is_initialized() { return m_is_initialized; }

void SerialInternal::initialize() {
  if (is_initialized()) return;

  Impl::SharedAllocationRecord<void, void>::tracking_enable();

  // Init the array of locks used for arbitrarily sized atomics
  Impl::init_lock_array_host_space();

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

  Kokkos::Profiling::finalize();

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
#ifdef KOKKOS_IMPL_WORKAROUND_ICE_IN_TRILINOS_WITH_OLD_INTEL_COMPILERS
    : m_space_instance(&Impl::SerialInternal::singleton()) {
}
#else
    : m_space_instance(&Impl::SerialInternal::singleton(),
                       [](Impl::SerialInternal*) {}) {
}
#endif

bool Serial::impl_is_initialized() {
  return Impl::SerialInternal::singleton().is_initialized();
}

void Serial::impl_initialize() {
  Impl::SerialInternal::singleton().initialize();
}

void Serial::impl_finalize() { Impl::SerialInternal::singleton().finalize(); }

const char* Serial::name() { return "Serial"; }

namespace Impl {

int g_serial_space_factory_initialized =
    initialize_space_factory<SerialSpaceInitializer>("100_Serial");

void SerialSpaceInitializer::initialize(const InitArguments& args) {
  // Prevent "unused variable" warning for 'args' input struct.  If
  // Serial::initialize() ever needs to take arguments from the input
  // struct, you may remove this line of code.
  (void)args;

  // Always initialize Serial if it is configure time enabled
  Kokkos::Serial::impl_initialize();
}

void SerialSpaceInitializer::finalize(const bool) {
  if (Kokkos::Serial::impl_is_initialized()) Kokkos::Serial::impl_finalize();
}

void SerialSpaceInitializer::fence() { Kokkos::Serial::impl_static_fence(); }
void SerialSpaceInitializer::fence(const std::string& name) {
  Kokkos::Serial::impl_static_fence(name);
}

void SerialSpaceInitializer::print_configuration(std::ostream& msg,
                                                 const bool detail) {
  msg << "Host Serial Execution Space:" << std::endl;
  msg << "  KOKKOS_ENABLE_SERIAL: ";
  msg << "yes" << std::endl;

  msg << "Serial Atomics:" << std::endl;
  msg << "  KOKKOS_ENABLE_SERIAL_ATOMICS: ";
#ifdef KOKKOS_ENABLE_SERIAL_ATOMICS
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif

  msg << "\nSerial Runtime Configuration:" << std::endl;
  Serial::print_configuration(msg, detail);
}

}  // namespace Impl

#ifdef KOKKOS_ENABLE_CXX14
namespace Tools {
namespace Experimental {
constexpr DeviceType DeviceTypeTraits<Serial>::id;
}
}  // namespace Tools
#endif

}  // namespace Kokkos

#else
void KOKKOS_CORE_SRC_IMPL_SERIAL_PREVENT_LINK_ERROR() {}
#endif  // defined( KOKKOS_ENABLE_SERIAL )
