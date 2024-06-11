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

#ifndef KOKKOS_DECLARE_HIP_HPP
#define KOKKOS_DECLARE_HIP_HPP

#if defined(KOKKOS_ENABLE_HIP)
#include <HIP/Kokkos_HIP.hpp>
#include <HIP/Kokkos_HIP_Space.hpp>
#include <HIP/Kokkos_HIP_DeepCopy.hpp>
#include <HIP/Kokkos_HIP_Half_Impl_Type.hpp>
#include <HIP/Kokkos_HIP_Half_Conversion.hpp>
#include <HIP/Kokkos_HIP_Instance.hpp>
#include <HIP/Kokkos_HIP_MDRangePolicy.hpp>
#include <HIP/Kokkos_HIP_ParallelFor_Range.hpp>
#include <HIP/Kokkos_HIP_ParallelFor_MDRange.hpp>
#include <HIP/Kokkos_HIP_ParallelFor_Team.hpp>
#include <HIP/Kokkos_HIP_ParallelReduce_Range.hpp>
#include <HIP/Kokkos_HIP_ParallelReduce_MDRange.hpp>
#include <HIP/Kokkos_HIP_ParallelReduce_Team.hpp>
#include <HIP/Kokkos_HIP_ParallelScan_Range.hpp>
#include <HIP/Kokkos_HIP_SharedAllocationRecord.hpp>
#include <HIP/Kokkos_HIP_UniqueToken.hpp>
#include <HIP/Kokkos_HIP_ZeroMemset.hpp>

namespace Kokkos {
namespace Experimental {
using HIPSpace           = ::Kokkos::HIPSpace;
using HIPHostPinnedSpace = ::Kokkos::HIPHostPinnedSpace;
using HIPManagedSpace    = ::Kokkos::HIPManagedSpace;
using HIP                = ::Kokkos::HIP;
}  // namespace Experimental
}  // namespace Kokkos
#endif

#endif
