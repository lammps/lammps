// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_HPX_HPP
#define KOKKOS_TEST_HPX_HPP

#include <gtest/gtest.h>

#define TEST_CATEGORY hpx
#define TEST_CATEGORY_NUMBER 3
#define TEST_CATEGORY_DEATH hpx_DeathTest
#define TEST_EXECSPACE Kokkos::Experimental::HPX
#define TEST_CATEGORY_FIXTURE(name) hpx_##name

#endif
