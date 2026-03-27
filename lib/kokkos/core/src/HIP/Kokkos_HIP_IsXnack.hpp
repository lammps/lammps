// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

/*--------------------------------------------------------------------------*/

#ifndef KOKKOS_HIP_ISXNACK_HPP
#define KOKKOS_HIP_ISXNACK_HPP

#include <Kokkos_Macros.hpp>

namespace Kokkos::Impl {

/*Based on AMD's ROCm 6.3.1 documentation:
https://github.com/ROCm/HIP/blob/2c240cacff16c2bb18ce9e5b4c1b937ab17a0199/docs/how-to/hip_runtime_api/memory_management/unified_memory.rst?plain=1#L141-L146

    To ensure the proper functioning of system allocated unified memory on
    supported GPUs, it is essential to configure the environment variable
    ``XNACK=1`` (sic) and use a kernel that supports HMM. Without this
    configuration, the behavior will be similar to that of systems without HMM
    support.

This clearly states two things are required:
* HSA_XNACK=1 is set in the environment
* The kernel must support HMM

Across a couple Nvidia and AMD systems, we have observed that
CONFIG_HMM_MIRROR=y was set in /boot/config-$(uname -r). This test may need to
be modified if a better way is determined to check for HMM support in Linux.
Checking for CONFIG_HMM was considered, but it was not present on El Capitan,
so we infer its presence is not necessary.
*/

// Returns true iff we detect HSA_XNACK=1 in the environment.
bool xnack_environment_enabled();
// Returns true iff we detect CONFIG_HMM_MIROR=y in /boot/config-$(uname -r).
bool xnack_boot_config_has_hmm_mirror();
// Returns true iff the architecture of the gpu supports accessing system
// allocated memory
constexpr bool gpu_arch_can_access_system_allocations() {
#if defined(KOKKOS_ARCH_AMD_GFX908) || defined(KOKKOS_ARCH_AMD_GFX90A) || \
    defined(KOKKOS_ARCH_AMD_GFX942) || defined(KOKKOS_ARCH_AMD_GFX942_APU)
  return true;
#elif defined(KOKKOS_ARCH_AMD_GFX906) || defined(KOKKOS_ARCH_AMD_GFX1103) || \
    defined(KOKKOS_ARCH_AMD_GFX1100) || defined(KOKKOS_ARCH_AMD_GFX1030) ||  \
    defined(KOKKOS_ARCH_AMD_GFX1201)
  return false;
#endif
}
}  // namespace Kokkos::Impl

#endif  // KOKKOS_HIP_ISXNACK_HPP
