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

#ifndef KOKKOS_DECLARE_OPENACC_HPP
#define KOKKOS_DECLARE_OPENACC_HPP

#if defined(KOKKOS_ENABLE_OPENACC)
#include <OpenACC/Kokkos_OpenACC.hpp>
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
