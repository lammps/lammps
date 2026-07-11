// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_HIP_HPP
#define KOKKOS_TEST_HIP_HPP

#include <gtest/gtest.h>

#define TEST_CATEGORY hip
#define TEST_CATEGORY_NUMBER 6
#define TEST_CATEGORY_DEATH hip_DeathTest
#define TEST_EXECSPACE Kokkos::HIP
#define TEST_CATEGORY_FIXTURE(name) hip_##name

#endif
