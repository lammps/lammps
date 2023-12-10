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

#ifndef KOKKOS_DECLARE_THREADS_HPP
#define KOKKOS_DECLARE_THREADS_HPP

#if defined(KOKKOS_ENABLE_THREADS)
#include <Threads/Kokkos_Threads.hpp>
#include <Threads/Kokkos_ThreadsExec.hpp>
#include <Threads/Kokkos_Threads_MDRangePolicy.hpp>
#include <Threads/Kokkos_Threads_ParallelFor_Range.hpp>
#include <Threads/Kokkos_Threads_ParallelFor_MDRange.hpp>
#include <Threads/Kokkos_Threads_ParallelFor_Team.hpp>
#include <Threads/Kokkos_Threads_ParallelReduce_Range.hpp>
#include <Threads/Kokkos_Threads_ParallelReduce_MDRange.hpp>
#include <Threads/Kokkos_Threads_ParallelReduce_Team.hpp>
#include <Threads/Kokkos_Threads_ParallelScan_Range.hpp>
#include <Threads/Kokkos_ThreadsTeam.hpp>
#include <Threads/Kokkos_Threads_UniqueToken.hpp>
#endif

#endif
