// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_DEVICE_MANAGEMENT_HPP
#define KOKKOS_DEVICE_MANAGEMENT_HPP

#include <optional>
#include <vector>

namespace Kokkos {
class InitializationSettings;
namespace Impl {
std::optional<int> get_gpu(const Kokkos::InitializationSettings& settings);
// This declaration is provided for testing purposes only
int get_ctest_gpu(int local_rank);
std::vector<int> get_visible_devices(int device_count);  // test-only
std::vector<int> const& get_visible_devices();           // use this instead
}  // namespace Impl
}  // namespace Kokkos

#endif
