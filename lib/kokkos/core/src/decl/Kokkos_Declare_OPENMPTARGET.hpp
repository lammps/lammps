// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_DECLARE_OPENMPTARGET_HPP
#define KOKKOS_DECLARE_OPENMPTARGET_HPP

#if defined(KOKKOS_ENABLE_OPENMPTARGET)
#include <OpenMPTarget/Kokkos_OpenMPTarget.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTargetSpace.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_Reducer.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_MDRangePolicy.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_UniqueToken.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_ParallelFor_Range.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_ParallelFor_Team.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_ParallelReduce_Range.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_ParallelReduce_Team.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_ParallelScan_Range.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_ParallelScan_Team.hpp>
#endif

#endif
