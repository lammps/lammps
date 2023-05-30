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

#define KOKKOS_IMPL_PUBLIC_INCLUDE

#include <OpenACC/Kokkos_OpenACC_Instance.hpp>
#include <OpenACC/Kokkos_OpenACC.hpp>
#include <OpenACC/Kokkos_OpenACC_Traits.hpp>
#include <impl/Kokkos_Profiling.hpp>
#include <impl/Kokkos_DeviceManagement.hpp>

#include <openacc.h>

#include <iostream>

namespace Kokkos {
bool show_warnings() noexcept;
}

Kokkos::Experimental::Impl::OpenACCInternal*
Kokkos::Experimental::Impl::OpenACCInternal::singleton() {
  static OpenACCInternal self;
  return &self;
}

void Kokkos::Experimental::Impl::OpenACCInternal::initialize(
    InitializationSettings const& settings) {
  if (OpenACC_Traits::may_fallback_to_host &&
      acc_get_num_devices(OpenACC_Traits::dev_type) == 0 &&
      !settings.has_device_id()) {
    if (show_warnings()) {
      std::cerr << "Warning: No GPU available for execution, falling back to"
                   " using the host!"
                << std::endl;
    }
    acc_set_device_type(acc_device_host);
    // FIXME_OPENACC if multiple execution space instances are supported,
    // device id variable should be explicitly set to the value returned by
    // acc_get_device_num(acc_device_host).
  } else {
    using Kokkos::Impl::get_gpu;
    int const dev_num = get_gpu(settings);
    acc_set_device_num(dev_num, OpenACC_Traits::dev_type);
  }
  m_is_initialized = true;
}

void Kokkos::Experimental::Impl::OpenACCInternal::finalize() {
  m_is_initialized = false;
}

bool Kokkos::Experimental::Impl::OpenACCInternal::is_initialized() const {
  return m_is_initialized;
}

void Kokkos::Experimental::Impl::OpenACCInternal::print_configuration(
    std::ostream& os, bool /*verbose*/) const {
  os << "Using OpenACC\n";  // FIXME_OPENACC
}

void Kokkos::Experimental::Impl::OpenACCInternal::fence(
    std::string const& name) const {
  Kokkos::Tools::Experimental::Impl::profile_fence_event<
      Kokkos::Experimental::OpenACC>(
      name,
      Kokkos::Tools::Experimental::Impl::DirectFenceIDHandle{instance_id()},
      [&]() {
        //[DEBUG] disabled due to synchronous behaviors of the current
        // parallel construct implementations. acc_wait_all();
      });
}

uint32_t Kokkos::Experimental::Impl::OpenACCInternal::instance_id() const
    noexcept {
  return Kokkos::Tools::Experimental::Impl::idForInstance<
      Kokkos::Experimental::OpenACC>(reinterpret_cast<uintptr_t>(this));
}
