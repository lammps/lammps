// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <NextSilicon/Kokkos_NextSilicon_ThreadSpaceGuard.hpp>
#include <Kokkos_Assert.hpp>
#include <nextapi/intrinsics.h>

namespace Kokkos::Impl {

thread_local PageAlignedData<bool, PageLocation::Host>
    NextSiliconThreadSpaceGuard::thread_is_on_device = false;

NextSiliconThreadSpaceGuard::NextSiliconThreadSpaceGuard() noexcept {
  // Touching thread_local variables cannot be done on device.
  // We can always determine if we are actually handed off by calling
  // __next_is_in_handed_off_code.
  if (!__next_is_in_handed_off_code()) {
    KOKKOS_ASSERT(!thread_is_on_device);
    thread_is_on_device = true;
  }
}

NextSiliconThreadSpaceGuard::~NextSiliconThreadSpaceGuard() noexcept {
  // We assume the invariant that __next_is_in_handed_off_code's result
  // cannot change mid-scope, otherwise it will break the counter.
  if (!__next_is_in_handed_off_code()) {
    // Catch if the invariant doesn't hold
    KOKKOS_ASSERT(thread_is_on_device);
    thread_is_on_device = false;
  }
}

}  // namespace Kokkos::Impl
