// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_THREADS_HPP
#define KOKKOS_TEST_THREADS_HPP

#include <gtest/gtest.h>

#define TEST_CATEGORY threads
#define TEST_CATEGORY_NUMBER 1
#define TEST_CATEGORY_DEATH threads_DeathTest
#define TEST_EXECSPACE Kokkos::Threads
#define TEST_CATEGORY_FIXTURE(name) threads_##name

#endif
