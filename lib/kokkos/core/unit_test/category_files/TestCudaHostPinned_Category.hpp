// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_CUDAHOSTPINNED_HPP
#define KOKKOS_TEST_CUDAHOSTPINNED_HPP

#include <gtest/gtest.h>

#define TEST_CATEGORY cuda_hostpinned
#define TEST_CATEGORY_DEATH cuda_hostpinned_DeathTest
// #define TEST_EXECSPACE
//  Kokkos::Device<Kokkos::Cuda,Kokkos::CudaHostPinnedSpace>
#define TEST_EXECSPACE Kokkos::CudaHostPinnedSpace

#endif
