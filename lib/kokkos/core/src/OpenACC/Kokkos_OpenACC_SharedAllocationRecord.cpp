// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#define KOKKOS_IMPL_PUBLIC_INCLUDE

#include <OpenACC/Kokkos_OpenACC.hpp>
#include <OpenACC/Kokkos_OpenACC_DeepCopy.hpp>
#include <OpenACC/Kokkos_OpenACC_SharedAllocationRecord.hpp>

#include <impl/Kokkos_SharedAlloc_timpl.hpp>

KOKKOS_IMPL_HOST_INACCESSIBLE_SHARED_ALLOCATION_RECORD_EXPLICIT_INSTANTIATION(
    Kokkos::Experimental::OpenACCSpace);
