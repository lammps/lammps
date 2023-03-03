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

#ifndef KOKKOS_TEST_OPENMP_HPP
#define KOKKOS_TEST_OPENMP_HPP

#include <gtest/gtest.h>

#include <Kokkos_Macros.hpp>

#ifdef KOKKOS_LAMBDA
#undef KOKKOS_LAMBDA
#endif
#define KOKKOS_LAMBDA [=]

#include <Kokkos_Core.hpp>

#include <TestViewMapping.hpp>
#include <TestViewAPI.hpp>
#include <TestViewOfClass.hpp>
#include <TestViewSubview.hpp>
#include <TestAtomic.hpp>
#include <TestAtomicOperations.hpp>
#include <TestAtomicViews.hpp>
#include <TestRange.hpp>
#include <TestTeam.hpp>
#include <TestReduce.hpp>
#include <TestScan.hpp>
#include <TestAggregate.hpp>
#include <TestCompilerMacros.hpp>
#include <TestTaskScheduler.hpp>
#include <TestMemoryPool.hpp>
#include <TestCXX11.hpp>
#include <TestCXX11Deduction.hpp>
#include <TestTeamVector.hpp>
#include <TestPolicyConstruction.hpp>
#include <TestMDRange.hpp>
#include <TestConcurrentBitset.hpp>

namespace Test {

class openmp : public ::testing::Test {
 protected:
  static void SetUpTestCase() {
    int threads_count = 0;
#pragma omp parallel
    {
#pragma omp atomic
      ++threads_count;
    }

    if (threads_count > 3) {
      threads_count /= 2;
    }

    Kokkos::OpenMP::initialize(threads_count);
    Kokkos::print_configuration(std::cout, true);

    srand(10231);
  }

  static void TearDownTestCase() { Kokkos::OpenMP::finalize(); }
};

}  // namespace Test

#endif
