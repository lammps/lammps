// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_OPENACC_SHARED_ALLOCATION_RECORD_HPP
#define KOKKOS_OPENACC_SHARED_ALLOCATION_RECORD_HPP

#include <OpenACC/Kokkos_OpenACCSpace.hpp>
#include <impl/Kokkos_SharedAlloc.hpp>

KOKKOS_IMPL_HOST_INACCESSIBLE_SHARED_ALLOCATION_SPECIALIZATION(
    Kokkos::Experimental::OpenACCSpace);

#endif
