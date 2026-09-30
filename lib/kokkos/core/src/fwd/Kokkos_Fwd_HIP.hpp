// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_HIP_FWD_HPP_
#define KOKKOS_HIP_FWD_HPP_

#if defined(KOKKOS_ENABLE_HIP)
namespace Kokkos {
class HIPSpace;            ///< Memory space on HIP GPU
class HIPHostPinnedSpace;  ///< Memory space on Host accessible to HIP GPU
class HIPManagedSpace;     ///< Memory migratable between Host and HIP GPU
class HIP;                 ///< Execution space for HIP GPU
}  // namespace Kokkos
#endif
#endif
