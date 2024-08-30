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

#define KOKKOS_IMPL_PUBLIC_INCLUDE

#include <OpenACC/Kokkos_OpenACC.hpp>
#include <OpenACC/Kokkos_OpenACC_Instance.hpp>
#include <impl/Kokkos_Profiling.hpp>
#include <impl/Kokkos_DeviceManagement.hpp>

#include <openacc.h>

#include <iostream>

// Arbitrary value to denote that we don't know yet what device to use.
int Kokkos::Experimental::Impl::OpenACCInternal::m_acc_device_num = -1;

Kokkos::Experimental::Impl::OpenACCInternal&
Kokkos::Experimental::Impl::OpenACCInternal::singleton() {
  static OpenACCInternal self;
  return self;
}

bool Kokkos::Experimental::Impl::OpenACCInternal::verify_is_initialized(
    const char* const label) const {
  if (!m_is_initialized) {
    Kokkos::abort((std::string("Kokkos::Experimental::OpenACC::") + label +
                   " : ERROR device not initialized\n")
                      .c_str());
  }
  return m_is_initialized;
}

void Kokkos::Experimental::Impl::OpenACCInternal::initialize(int async_arg) {
  if ((async_arg < 0) && (async_arg != acc_async_sync) &&
      (async_arg != acc_async_noval)) {
    Kokkos::abort((std::string("Kokkos::Experimental::OpenACC::initialize()") +
                   " : ERROR async_arg should be a non-negative integer" +
                   " unless being a special value defined in OpenACC\n")
                      .c_str());
  }
  m_async_arg      = async_arg;
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
      [&]() { acc_wait(m_async_arg); });
}

uint32_t Kokkos::Experimental::Impl::OpenACCInternal::instance_id() const
    noexcept {
  return Kokkos::Tools::Experimental::Impl::idForInstance<OpenACC>(
      reinterpret_cast<uintptr_t>(this));
}
