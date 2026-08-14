// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_HIPHOSTPINNED_HPP
#define KOKKOS_TEST_HIPHOSTPINNED_HPP

#include <gtest/gtest.h>

#define TEST_CATEGORY hip_hostpinned
#define TEST_CATEGORY_DEATH hip_hostpinned_DeathTest
#define TEST_EXECSPACE Kokkos::HIPHostPinnedSpace

#endif
