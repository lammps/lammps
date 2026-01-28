// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_DEFAULTDEVICETYPE_HPP
#define KOKKOS_TEST_DEFAULTDEVICETYPE_HPP

#include <gtest/gtest.h>

#define TEST_CATEGORY defaultdevicetype
#define TEST_CATEGORY_DEATH defaultdevicetype_DeathTest
#define TEST_EXECSPACE Kokkos::DefaultExecutionSpace
#define TEST_CATEGORY_FIXTURE(name) default_##name

#endif
