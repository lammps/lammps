// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_SYCL_FWD_HPP_
#define KOKKOS_SYCL_FWD_HPP_

#if defined(KOKKOS_ENABLE_SYCL)
namespace Kokkos {
class SYCLDeviceUSMSpace;  ///< Memory space on SYCL device, not accessible from
                           ///< the host
class SYCLSharedUSMSpace;  ///< Memory space accessible from both the SYCL
                           ///< device and the host
class SYCLHostUSMSpace;    ///< Memory space accessible from both the SYCL
                           ///< device and the host (host pinned)
class SYCL;                ///< Execution space for SYCL
}  // namespace Kokkos
#endif
#endif
