// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_CUDA_HPP
#define KOKKOS_TEST_CUDA_HPP

#include <gtest/gtest.h>

#define TEST_CATEGORY cuda
#define TEST_CATEGORY_NUMBER 5
#define TEST_CATEGORY_DEATH cuda_DeathTest
#define TEST_EXECSPACE Kokkos::Cuda
#define TEST_CATEGORY_FIXTURE(name) cuda_##name

#endif
