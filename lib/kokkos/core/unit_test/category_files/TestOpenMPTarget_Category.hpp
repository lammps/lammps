// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_OMPTARGET_HPP
#define KOKKOS_TEST_OMPTARGET_HPP

#include <gtest/gtest.h>

#define TEST_CATEGORY openmptarget
#define TEST_CATEGORY_NUMBER 4
#define TEST_CATEGORY_DEATH openmptarget_DeathTest
#define TEST_EXECSPACE Kokkos::Experimental::OpenMPTarget
#define TEST_CATEGORY_FIXTURE(name) openmptarget_##name

#endif
