// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_HIPUNIFIED_HPP
#define KOKKOS_TEST_HIPUNIFIED_HPP

#include <gtest/gtest.h>

#define TEST_CATEGORY hip_managed
#define TEST_CATEGORY_DEATH hip_managed_DeathTest
#define TEST_EXECSPACE Kokkos::HIPManagedSpace

#endif
