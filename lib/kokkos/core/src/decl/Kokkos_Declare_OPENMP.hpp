// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_DECLARE_OPENMP_HPP
#define KOKKOS_DECLARE_OPENMP_HPP

#if defined(KOKKOS_ENABLE_OPENMP)
#include <OpenMP/Kokkos_OpenMP.hpp>
#include <OpenMP/Kokkos_OpenMP_MDRangePolicy.hpp>
#include <OpenMP/Kokkos_OpenMP_UniqueToken.hpp>
#include <OpenMP/Kokkos_OpenMP_Parallel_For.hpp>
#include <OpenMP/Kokkos_OpenMP_Parallel_Reduce.hpp>
#include <OpenMP/Kokkos_OpenMP_Parallel_Scan.hpp>
#include <OpenMP/Kokkos_OpenMP_ZeroMemset.hpp>
#endif

#endif
