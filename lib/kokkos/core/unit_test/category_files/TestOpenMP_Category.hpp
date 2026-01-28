// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_OMP_HPP
#define KOKKOS_TEST_OMP_HPP

#include <gtest/gtest.h>

#define TEST_CATEGORY openmp
#define TEST_CATEGORY_NUMBER 2
#define TEST_CATEGORY_DEATH openmp_DeathTest
#define TEST_EXECSPACE Kokkos::OpenMP
#define TEST_CATEGORY_FIXTURE(name) openmp_##name

#endif
