// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_DECLARE_SYCL_HPP
#define KOKKOS_DECLARE_SYCL_HPP

#if defined(KOKKOS_ENABLE_SYCL)
#include <SYCL/Kokkos_SYCL.hpp>
#ifdef KOKKOS_IMPL_SYCL_GRAPH_SUPPORT
#include <SYCL/Kokkos_SYCL_GraphNodeKernel.hpp>
#endif
#include <SYCL/Kokkos_SYCL_DeepCopy.hpp>
#include <SYCL/Kokkos_SYCL_MDRangePolicy.hpp>
#include <SYCL/Kokkos_SYCL_ParallelFor_Range.hpp>
#include <SYCL/Kokkos_SYCL_ParallelFor_MDRange.hpp>
#include <SYCL/Kokkos_SYCL_ParallelFor_Team.hpp>
#include <SYCL/Kokkos_SYCL_ParallelReduce_Range.hpp>
#include <SYCL/Kokkos_SYCL_ParallelReduce_MDRange.hpp>
#include <SYCL/Kokkos_SYCL_ParallelReduce_Team.hpp>
#include <SYCL/Kokkos_SYCL_ParallelScan_Range.hpp>
#include <SYCL/Kokkos_SYCL_UniqueToken.hpp>
#include <SYCL/Kokkos_SYCL_ZeroMemset.hpp>

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_5
namespace Kokkos {
namespace Experimental {
using SYCLDeviceUSMSpace KOKKOS_DEPRECATED_WITH_COMMENT(
    "Use Kokkos::SYCLDeviceUSMSpace instead!") = ::Kokkos::SYCLDeviceUSMSpace;
using SYCLHostUSMSpace KOKKOS_DEPRECATED_WITH_COMMENT(
    "Use Kokkos::SYCLHostUSMSpace instead!") = ::Kokkos::SYCLHostUSMSpace;
using SYCLSharedUSMSpace KOKKOS_DEPRECATED_WITH_COMMENT(
    "Use Kokkos::SYCLSharedUSMSpace instead!") = ::Kokkos::SYCLSharedUSMSpace;
using SYCL KOKKOS_DEPRECATED_WITH_COMMENT("Use Kokkos::SYCL instead!") =
    ::Kokkos::SYCL;
}  // namespace Experimental
}  // namespace Kokkos
#endif

#endif

#endif
