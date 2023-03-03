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

#ifndef KOKKOS_TEST_OPENMPTARGET_HPP
#define KOKKOS_TEST_OPENMPTARGET_HPP

#include <gtest/gtest.h>

#include <Kokkos_Macros.hpp>

#ifdef KOKKOS_LAMBDA
#undef KOKKOS_LAMBDA
#endif
#define KOKKOS_LAMBDA [=]

#include <Kokkos_Core.hpp>

//#include <TestViewAPI.hpp>
//#include <TestViewOfClass.hpp>
//#include <TestViewSubview.hpp>
//#include <TestAtomic.hpp>
//#include <TestAtomicOperations.hpp>
//#include <TestAtomicViews.hpp>
#include <TestRange.hpp>
#include <TestTeam.hpp>
//#include <TestReduce.hpp>
//#include <TestScan.hpp>
//#include <TestAggregate.hpp>
//#include <TestCompilerMacros.hpp>

// TODO enable task scheduler tests for openmptarget
//#include <TestTaskScheduler.hpp>

//#include <TestMemoryPool.hpp>
//#include <TestCXX11.hpp>
//#include <TestCXX11Deduction.hpp>
#include <TestTeamVector.hpp>
//#include <TestPolicyConstruction.hpp>
//#include <TestMDRange.hpp>

namespace Test {

class openmptarget : public ::testing::Test {
 protected:
  static void SetUpTestCase() {
    const unsigned numa_count = Kokkos::hwloc::get_available_numa_count();
    const unsigned cores_per_numa =
        Kokkos::hwloc::get_available_cores_per_numa();
    const unsigned openmptarget_per_core =
        Kokkos::hwloc::get_available_openmptarget_per_core();

    unsigned openmptarget_count = 0;

    openmptarget_count = std::max(1u, numa_count) *
                         std::max(2u, cores_per_numa * openmptarget_per_core);

    Kokkos::OpenMPTarget::initialize(openmptarget_count);
    Kokkos::print_configuration(std::cout, true /* detailed */);
  }

  static void TearDownTestCase() { Kokkos::OpenMPTarget::finalize(); }
};

}  // namespace Test

#endif
