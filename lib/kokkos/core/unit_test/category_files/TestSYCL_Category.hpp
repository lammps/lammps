// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_SYCL_HPP
#define KOKKOS_TEST_SYCL_HPP

#include <gtest/gtest.h>

#define TEST_CATEGORY sycl
#define TEST_CATEGORY_NUMBER 7
#define TEST_CATEGORY_DEATH sycl_DeathTest
#define TEST_EXECSPACE Kokkos::SYCL
#define TEST_CATEGORY_FIXTURE(name) sycl_##name

#endif
