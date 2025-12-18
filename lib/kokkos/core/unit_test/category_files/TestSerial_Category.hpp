// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_SERIAL_HPP
#define KOKKOS_TEST_SERIAL_HPP

#include <gtest/gtest.h>

#define TEST_CATEGORY serial
#define TEST_CATEGORY_NUMBER 0
#define TEST_CATEGORY_DEATH serial_DeathTest
#define TEST_EXECSPACE Kokkos::Serial
#define TEST_CATEGORY_FIXTURE(name) serial_##name

#endif
