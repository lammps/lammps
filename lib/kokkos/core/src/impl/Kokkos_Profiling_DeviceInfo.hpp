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

#ifndef KOKKOSP_DEVICE_INFO_HPP
#define KOKKOSP_DEVICE_INFO_HPP

#include <cstdint>
#include <impl/Kokkos_Profiling_C_Interface.h>
namespace Kokkos {
namespace Profiling {
using KokkosPDeviceInfo = Kokkos_Profiling_KokkosPDeviceInfo;
}  // namespace Profiling
}  // namespace Kokkos

#endif
