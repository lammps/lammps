// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <Kokkos_Core.hpp>
#include <mutex>
#include <ostream>
#include <cstdint>
#include <cstddef>

namespace Kokkos::Experimental::Impl {

Kokkos::Impl::HostSharedPtr<Kokkos::Experimental::Impl::NextSiliconInternal>
    Kokkos::Experimental::Impl::NextSiliconInternal::default_instance;

NextSiliconInternal::NextSiliconInternal() {}

NextSiliconInternal::~NextSiliconInternal() {
  fence(
      "Kokkos::Experimental::Impl::NextSiliconInternal: fence on destruction");
}

void NextSiliconInternal::print_configuration(std::ostream &os) const {
  os << "macro  KOKKOS_ENABLE_NEXTSILICON      : defined\n";
  // FIXME_NEXTSILICON print_configuration doesn't do anything useful, fix
  // once device properties nextapi is available
}

void NextSiliconInternal::fence(std::string const &name) const {
  // FIXME_NEXTSILICON fence doesn't do anything in the OpenMP-style interface

  Kokkos::Tools::Experimental::Impl::profile_fence_event<NextSilicon>(
      name,
      Kokkos::Tools::Experimental::Impl::DirectFenceIDHandle{instance_id()},
      [&]() {});
}

uint32_t NextSiliconInternal::instance_id() const noexcept {
  return Kokkos::Tools::Experimental::Impl::idForInstance<
      Kokkos::Experimental::NextSilicon>(reinterpret_cast<uintptr_t>(this));
}

std::byte *NextSiliconInternal::resize_functor_buffer(size_t requested) {
  constexpr static size_t MIN_FUNCTOR_BUFFER_SIZE = 4 * 1024 * 1024;  // 4 MB
  requested = std::max(requested, MIN_FUNCTOR_BUFFER_SIZE);

  return functorBuffer_.ensure("functor heap buffer", requested);
}

std::lock_guard<std::mutex>
Kokkos::Experimental::Impl::NextSiliconInternal::lock_device() {
  KOKKOS_IF_ON_DEVICE(
      (KOKKOS_ASSERT(false && "lock_device should never be called on device");))
  return std::lock_guard<std::mutex>(this->device_mutex_);
}

}  // namespace Kokkos::Experimental::Impl
