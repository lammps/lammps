// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_KOKKOS_DEVICEHANDLE_HPP
#define KOKKOS_IMPL_KOKKOS_DEVICEHANDLE_HPP

#include "View/Kokkos_ViewCtor.hpp"

namespace Kokkos::Impl {
template <Kokkos::ExecutionSpace Exec>
struct DeviceHandle {
  // A device is implicitly contained in an execution space instance.
  static DeviceHandle from(Exec exec) {
    return DeviceHandle{.m_exec = std::move(exec)};
  }

  auto operator<=>(const DeviceHandle&) const = default;

  // For now, let's store the execution space instance.
  // It is the best portable way to ensure that the device handle has enough
  // information. It is an implementation detail, let's keep this reference
  // counted member extractible easily (no need to make it private yet).
  // In the future, it could be treated as the default execution queue for the
  // device.
  Exec m_exec;
};

template <typename>
struct is_device_handle : public std::false_type {};

template <typename T>
struct is_device_handle<DeviceHandle<T>> : public std::true_type {};

template <typename T>
constexpr bool is_device_handle_v = is_device_handle<T>::value;
}  // namespace Kokkos::Impl

namespace Kokkos::Experimental {
// The user shall treat the return type as an opaque type.
template <typename Exec>
  requires Kokkos::ExecutionSpace<std::remove_cvref_t<Exec>>
auto get_device_handle(Exec&& exec) {
  return Kokkos::Impl::DeviceHandle<std::remove_cvref_t<Exec>>::from(
      std::forward<Exec>(exec));
}
}  // namespace Kokkos::Experimental

#endif  // KOKKOS_IMPL_KOKKOS_DEVICEHANDLE_HPP
