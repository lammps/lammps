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

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

#include <cstdlib>
#include <type_traits>

#include "KokkosExecutionEnvironmentNeverInitializedFixture.hpp"

namespace {

using ExecutionEnvironmentNonInitializedOrFinalized_DeathTest =
    KokkosExecutionEnvironmentNeverInitialized;

struct NonTrivial {
  KOKKOS_FUNCTION NonTrivial() {}
};
static_assert(!std::is_trivial_v<NonTrivial>);

TEST_F(ExecutionEnvironmentNonInitializedOrFinalized_DeathTest,
       default_constructed_views) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  auto make_views = [] {
    Kokkos::View<int> v0;
    Kokkos::View<float*> v1;
    Kokkos::View<NonTrivial**> v2;
    return std::make_tuple(v0, v1, v2);
  };
  EXPECT_EXIT(
      {
        { auto views = make_views(); }
        std::exit(EXIT_SUCCESS);
      },
      ::testing::ExitedWithCode(EXIT_SUCCESS), "");
  EXPECT_EXIT(
      {
        {
          Kokkos::initialize();
          auto views =
              make_views();  // views outlive the Kokkos execution environment
          Kokkos::finalize();
        }
        std::exit(EXIT_SUCCESS);
      },
      ::testing::ExitedWithCode(EXIT_SUCCESS), "");
  EXPECT_EXIT(
      {
        {
          Kokkos::initialize();
          Kokkos::finalize();
          auto views = make_views();
        }
        std::exit(EXIT_SUCCESS);
      },
      ::testing::ExitedWithCode(EXIT_SUCCESS), "");
}

TEST_F(ExecutionEnvironmentNonInitializedOrFinalized_DeathTest, views) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  EXPECT_EXIT(
      {
        {
          Kokkos::View<int*> v;
          Kokkos::initialize();
          v = Kokkos::View<int*>("v", 10);
          v = Kokkos::View<int*>();
          Kokkos::finalize();
        }
        std::exit(EXIT_SUCCESS);
      },
      ::testing::ExitedWithCode(EXIT_SUCCESS), "");
  EXPECT_EXIT(
      {
        {
          Kokkos::initialize();
          Kokkos::View<int*> v("v", 10);
          v = {};  // assign default constructed view
          Kokkos::finalize();
        }
        std::exit(EXIT_SUCCESS);
      },
      ::testing::ExitedWithCode(EXIT_SUCCESS), "");
  EXPECT_DEATH(
      {
        Kokkos::initialize();
        Kokkos::View<int*> v("v", 0);
        Kokkos::finalize();
      },
      "Kokkos allocation \"v\" is being deallocated after Kokkos::finalize was "
      "called");
  // NOLINTBEGIN(bugprone-unused-local-non-trivial-variable)
  [[maybe_unused]] std::string error_constructing_view_with_unitialized_exec =
      "Constructing View and initializing data with uninitialized execution "
      "space";
  [[maybe_unused]] std::string error_constructing_exec_space_instance =
      std::string("Kokkos::") +
#ifdef KOKKOS_ENABLE_OPENACC
      "Experimental::" +
#endif
      Kokkos::DefaultExecutionSpace::name() +
      "::" + Kokkos::DefaultExecutionSpace::name() +
      " instance constructor : ERROR device not initialized";
  // NOLINTEND(bugprone-unused-local-non-trivial-variable)
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || \
    defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ENABLE_OPENACC)
  std::string matcher1 = error_constructing_exec_space_instance;
#else
  std::string matcher1 = error_constructing_view_with_unitialized_exec;
#endif
#if defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ENABLE_OPENACC)
  std::string matcher2 = error_constructing_exec_space_instance;
#else
  std::string matcher2 = error_constructing_view_with_unitialized_exec;
#endif
  EXPECT_DEATH({ Kokkos::View<int*> v("v", 0); }, matcher1);
  EXPECT_DEATH(
      {
        Kokkos::initialize();
        Kokkos::finalize();
        Kokkos::View<int*> v("v", 0);
      },
      matcher2);
}

}  // namespace
