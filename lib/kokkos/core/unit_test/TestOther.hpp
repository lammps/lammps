//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_TEST_OTHER_HPP
#define KOKKOS_TEST_OTHER_HPP
#include <TestAggregate.hpp>
#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC:
// NVC++-F-0000-Internal compiler error. Basic LLVM base data type required 23
// (/ascldap/users/crtrott/Kokkos/kokkos/build/core/unit_test/cuda/TestCuda_Other.cpp:
// 204) NVC++/x86-64 Linux 22.3-0: compilation aborted
#include <TestMemoryPool.hpp>
#endif
#include <TestCXX11.hpp>

#include <TestViewCtorPropEmbeddedDim.hpp>
// with VS 16.11.3 and CUDA 11.4.2 getting cudafe stackoverflow crash
#if !(defined(_WIN32) && defined(KOKKOS_ENABLE_CUDA))
#include <TestViewLayoutTiled.hpp>
#endif
#endif
