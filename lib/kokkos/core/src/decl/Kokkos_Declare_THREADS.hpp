// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_DECLARE_THREADS_HPP
#define KOKKOS_DECLARE_THREADS_HPP

#if defined(KOKKOS_ENABLE_THREADS)
#include <Threads/Kokkos_Threads.hpp>
#include <Threads/Kokkos_Threads_Instance.hpp>
#include <Threads/Kokkos_Threads_MDRangePolicy.hpp>
#include <Threads/Kokkos_Threads_ParallelFor_Range.hpp>
#include <Threads/Kokkos_Threads_ParallelFor_MDRange.hpp>
#include <Threads/Kokkos_Threads_ParallelFor_Team.hpp>
#include <Threads/Kokkos_Threads_ParallelReduce_Range.hpp>
#include <Threads/Kokkos_Threads_ParallelReduce_MDRange.hpp>
#include <Threads/Kokkos_Threads_ParallelReduce_Team.hpp>
#include <Threads/Kokkos_Threads_ParallelScan_Range.hpp>
#include <Threads/Kokkos_Threads_Team.hpp>
#include <Threads/Kokkos_Threads_UniqueToken.hpp>
#include <Threads/Kokkos_Threads_ZeroMemset.hpp>
#endif

#endif
