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
