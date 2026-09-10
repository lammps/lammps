// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_DECLARE_OPENACC_HPP
#define KOKKOS_DECLARE_OPENACC_HPP

#if defined(KOKKOS_ENABLE_OPENACC)
#include <OpenACC/Kokkos_OpenACC.hpp>
#include <OpenACC/Kokkos_OpenACC_Instance.hpp>
#include <OpenACC/Kokkos_OpenACCSpace.hpp>
#include <OpenACC/Kokkos_OpenACC_DeepCopy.hpp>
#include <OpenACC/Kokkos_OpenACC_SharedAllocationRecord.hpp>
#include <OpenACC/Kokkos_OpenACC_ParallelFor_Range.hpp>
#include <OpenACC/Kokkos_OpenACC_ParallelReduce_Range.hpp>
#include <OpenACC/Kokkos_OpenACC_ParallelScan_Range.hpp>
#include <OpenACC/Kokkos_OpenACC_MDRangePolicy.hpp>
#include <OpenACC/Kokkos_OpenACC_ParallelFor_MDRange.hpp>
#include <OpenACC/Kokkos_OpenACC_ParallelReduce_MDRange.hpp>
#include <OpenACC/Kokkos_OpenACC_ParallelFor_Team.hpp>
#include <OpenACC/Kokkos_OpenACC_ParallelReduce_Team.hpp>
#endif

#endif
