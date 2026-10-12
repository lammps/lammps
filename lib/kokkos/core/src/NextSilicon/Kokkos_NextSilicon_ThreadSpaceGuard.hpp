// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_NEXTSILICON_THREAD_SPACE_GUARD_HPP
#define KOKKOS_NEXTSILICON_THREAD_SPACE_GUARD_HPP

#include <nextapi/intrinsics.h>
#include <NextSilicon/Kokkos_NextSilicon_PageAlignedData.hpp>

namespace Kokkos::Impl {
// used to implement KOKKOS_IF_ON_HOST / KOKKOS_IF_ON_DEVICE

// NextSiliconThreadSpaceGuard maintains a thread-local flag that tracks whether
// the current thread is executing "within" the NextSilicon device execution
// space. This is needed to correctly implement KOKKOS_IF_ON_DEVICE /
// KOKKOS_IF_ON_HOST semantics: unlike other backends, NextSilicon kernels may
// run on the host during training/telemetry collection while still being
// "within" the device execution space, so we cannot simply equate "on device"
// with "offloaded to device".

// The guard interface is modeled on the interface of std::lock_guard -- RAII,
// non-copyable, non-movable. Unlike std::lock_guard, because the underlying
// static flag is thread_local, no actual lock type is necessary.
//
// The class provides a static member function is_on_device() that returns true
// if the current thread should be considered "on device" for the purposes of
// KOKKOS_IF_ON_DEVICE / KOKKOS_IF_ON_HOST, and false otherwise. There are two
// distinct scenarios where a thread should be considered "on device":
//   1) The thread is actually executing on the device (i.e. after optimization
//      and handoff), as indicated by the __next_is_in_handed_off_code()
//      intrinsic returning true. This short circuit allows us to avoid touching
//      the thread_local variable at all from device code, thus avoiding
//      migration of the host thread's stack to device.
//   2) The thread is executing on the host, but has been "logically" handed off
//      to the device execution space. This can occur during training/telemetry
//      collection, where we want to execute the same code paths as we would on
//      the device, or if the optimizer has decided to schedule the kernel on
//      the host. Because the thread is still executing on the host,
//      __next_is_in_handed_off_code() will return false. In this case, the
//      thread_local flag being true indicates that we should still treat the
//      thread as logically "on device".
//
// Whenever a thread enters a scope where it should be considered "on device"
// (e.g. at the beginning of a parallel_function/microtask), a
// NextSiliconThreadSpaceGuard should be created to set the thread_local flag to
// true for the duration of that scope. Because the guard is non-copyable and
// non-movable, this ensures that the thread_local flag will be correctly
// reset to false when the scope is exited, even if an exception is thrown.
//
// If this thread is executing on the device, the guard must be a no-op to avoid
// touching thread_local flags from device code, which would cause the host
// thread's stack to be migrated to the device. This is achieved by
// short-circuiting on the __next_is_in_handed_off_code() intrinsic, which
// allows us to determine whether we are actually executing on the device
// without touching any thread_local flags.
class [[nodiscard]] NextSiliconThreadSpaceGuard {
 private:
  static thread_local PageAlignedData<bool, PageLocation::Host>
      thread_is_on_device;

  // Containment for the TLS read. If `thread_is_on_device > 0` were inlined
  // into is_on_device()'s callers, the underlying `llvm.threadlocal.address`
  // intrinsic (declared speculatable) could be hoisted past the `||`
  // short-circuit in is_on_device(), breaking semantics. Keeping the read
  // behind a non-inlined wrapper makes the call opaque at every use site,
  // so the speculation never reaches the load. `weak` is what prevents
  // inlining (its primary side effect is that the linker may interpose the
  // body, so the optimizer must treat it as opaque).
  //
  // FIXME_NEXTSILICON: __attribute__((noinline)) prevents scheduling to grid;
  // when ticket https://nextsilicon.atlassian.net/browse/SW-25677 is resolved,
  // __attribute__((weak)) can be switched to __attribute__((noinline)).
  static __attribute__((weak)) bool host_thread_is_on_device() {
    return thread_is_on_device;
  }

 public:
  explicit NextSiliconThreadSpaceGuard() noexcept;
  ~NextSiliconThreadSpaceGuard() noexcept;

  static __attribute__((weak)) bool is_on_device() noexcept;

  // Copy and move operations are deleted to emulate semantics of
  // std::lock_guard.
  NextSiliconThreadSpaceGuard(const NextSiliconThreadSpaceGuard&) = delete;
  NextSiliconThreadSpaceGuard(NextSiliconThreadSpaceGuard&&)      = delete;
  NextSiliconThreadSpaceGuard& operator=(const NextSiliconThreadSpaceGuard&) =
      delete;
  NextSiliconThreadSpaceGuard& operator=(NextSiliconThreadSpaceGuard&&) =
      delete;
};

inline __attribute__((weak)) bool
NextSiliconThreadSpaceGuard::is_on_device() noexcept {
  // We are either really handed off to the device, or should pretend to be
  // (for training mode).
  return __next_is_in_handed_off_code() || host_thread_is_on_device();
}

}  // namespace Kokkos::Impl

#endif
