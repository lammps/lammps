// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_HIP_SHARED_ALLOCATION_RECORD_HPP
#define KOKKOS_HIP_SHARED_ALLOCATION_RECORD_HPP

#include <HIP/Kokkos_HIP_Space.hpp>
#include <impl/Kokkos_SharedAlloc.hpp>

#if defined(KOKKOS_IMPL_HIP_UNIFIED_MEMORY)
KOKKOS_IMPL_SHARED_ALLOCATION_SPECIALIZATION(Kokkos::HIPSpace);
#else
KOKKOS_IMPL_HOST_INACCESSIBLE_SHARED_ALLOCATION_SPECIALIZATION(
    Kokkos::HIPSpace);
#endif
KOKKOS_IMPL_SHARED_ALLOCATION_SPECIALIZATION(Kokkos::HIPHostPinnedSpace);
KOKKOS_IMPL_SHARED_ALLOCATION_SPECIALIZATION(Kokkos::HIPManagedSpace);

#endif
