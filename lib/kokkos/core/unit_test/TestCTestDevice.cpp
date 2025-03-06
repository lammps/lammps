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

#include <impl/Kokkos_DeviceManagement.hpp>  // get_ctest_gpu

#ifdef _WIN32
int setenv(const char *name, const char *value, int overwrite) {
  int errcode = 0;
  if (!overwrite) {
    size_t envsize = 0;
    errcode        = getenv_s(&envsize, NULL, 0, name);
    if (errcode || envsize) return errcode;
  }
  return _putenv_s(name, value);
}

int unsetenv(const char *name) { return _putenv_s(name, ""); }
#endif

class ctest_environment : public ::testing::Test {
 protected:
  void SetUp() override;
};

void ctest_environment::SetUp() {
  setenv("CTEST_KOKKOS_DEVICE_TYPE", "gpus", 1);
  setenv("CTEST_RESOURCE_GROUP_COUNT", "10", 1);
  unsetenv("CTEST_RESOURCE_GROUP_0");
  setenv("CTEST_RESOURCE_GROUP_1", "threads", 1);
  setenv("CTEST_RESOURCE_GROUP_2", "threads,cores", 1);

  setenv("CTEST_RESOURCE_GROUP_3", "gpus", 1);
  unsetenv("CTEST_RESOURCE_GROUP_3_GPUS");

  setenv("CTEST_RESOURCE_GROUP_4", "gpus", 1);
  setenv("CTEST_RESOURCE_GROUP_4_GPUS", "id:2", 1);

  setenv("CTEST_RESOURCE_GROUP_5", "gpus", 1);
  setenv("CTEST_RESOURCE_GROUP_5_GPUS", "slots:1,id:2", 1);

  setenv("CTEST_RESOURCE_GROUP_6", "gpus", 1);
  setenv("CTEST_RESOURCE_GROUP_6_GPUS", "id:2,slots:1", 1);

  setenv("CTEST_RESOURCE_GROUP_7", "threads,gpus", 1);
  setenv("CTEST_RESOURCE_GROUP_7_GPUS", "id:3,slots:1", 1);

  setenv("CTEST_RESOURCE_GROUP_8", "gpus,threads", 1);
  setenv("CTEST_RESOURCE_GROUP_8_GPUS", "id:1,slots:1", 1);

  setenv("CTEST_RESOURCE_GROUP_9", "cores,gpus,threads", 1);
  setenv("CTEST_RESOURCE_GROUP_9_GPUS", "id:4,slots:1", 1);
}

struct ctest_environment_DeathTest : public ctest_environment {};

TEST_F(ctest_environment, no_device_type) {
  unsetenv("CTEST_KOKKOS_DEVICE_TYPE");
  EXPECT_EQ(Kokkos::Impl::get_ctest_gpu(0), 0);
}

TEST_F(ctest_environment, no_process_count) {
  unsetenv("CTEST_RESOURCE_GROUP_COUNT");
  EXPECT_EQ(Kokkos::Impl::get_ctest_gpu(0), 0);
}

TEST_F(ctest_environment_DeathTest, invalid_rank) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH(Kokkos::Impl::get_ctest_gpu(10),
               "Error: local rank 10 is outside the bounds of resource groups "
               "provided by CTest.");
}

TEST_F(ctest_environment_DeathTest, no_type_str) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH(Kokkos::Impl::get_ctest_gpu(0),
               "Error: CTEST_RESOURCE_GROUP_0 is not specified. Raised by "
               "Kokkos::Impl::get_ctest_gpu\\(\\).");
}

TEST_F(ctest_environment_DeathTest, missing_type) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH(
      Kokkos::Impl::get_ctest_gpu(1),
      "Error: device type 'gpus' not included in CTEST_RESOURCE_GROUP_1. "
      "Raised by Kokkos::Impl::get_ctest_gpu\\(\\).");
  EXPECT_DEATH(
      Kokkos::Impl::get_ctest_gpu(2),
      "Error: device type 'gpus' not included in CTEST_RESOURCE_GROUP_2. "
      "Raised by Kokkos::Impl::get_ctest_gpu\\(\\).");
}

TEST_F(ctest_environment_DeathTest, no_id_str) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH(Kokkos::Impl::get_ctest_gpu(3),
               "Error: CTEST_RESOURCE_GROUP_3_GPUS is not specified. Raised by "
               "Kokkos::Impl::get_ctest_gpu\\(\\).");
}

TEST_F(ctest_environment_DeathTest, invalid_id_str) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH(Kokkos::Impl::get_ctest_gpu(4),
               "Error: invalid value of CTEST_RESOURCE_GROUP_4_GPUS: 'id:2'. "
               "Raised by Kokkos::Impl::get_ctest_gpu\\(\\).");
  EXPECT_DEATH(Kokkos::Impl::get_ctest_gpu(5),
               "Error: invalid value of CTEST_RESOURCE_GROUP_5_GPUS: "
               "'slots:1,id:2'. Raised by Kokkos::Impl::get_ctest_gpu\\(\\).");
}

TEST_F(ctest_environment, good) {
  EXPECT_EQ(Kokkos::Impl::get_ctest_gpu(6), 2);
  EXPECT_EQ(Kokkos::Impl::get_ctest_gpu(7), 3);
  EXPECT_EQ(Kokkos::Impl::get_ctest_gpu(8), 1);
  EXPECT_EQ(Kokkos::Impl::get_ctest_gpu(9), 4);
}
