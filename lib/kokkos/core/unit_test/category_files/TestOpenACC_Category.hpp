// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_OACC_HPP
#define KOKKOS_TEST_OACC_HPP

#include <gtest/gtest.h>

#define TEST_CATEGORY openacc
#define TEST_CATEGORY_NUMBER 8
#define TEST_CATEGORY_DEATH openacc_DeathTest
#define TEST_EXECSPACE Kokkos::Experimental::OpenACC
#define TEST_CATEGORY_FIXTURE(name) openacc_##name

#endif
