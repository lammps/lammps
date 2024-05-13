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

#include <Kokkos_Core.hpp>

#ifndef KOKKOS_VERSION
static_assert(false, "KOKKOS_VERSION macro is not defined!");
#endif

#ifndef KOKKOS_VERSION_MAJOR
static_assert(false, "KOKKOS_VERSION_MAJOR macro is not defined!");
#endif

#ifndef KOKKOS_VERSION_MINOR
static_assert(false, "KOKKOS_VERSION_MINOR macro is not defined!");
#endif

#ifndef KOKKOS_VERSION_PATCH
static_assert(false, "KOKKOS_VERSION_PATCH macro is not defined!");
#endif

static_assert(KOKKOS_VERSION == KOKKOS_VERSION_MAJOR * 10000 +
                                    KOKKOS_VERSION_MINOR * 100 +
                                    KOKKOS_VERSION_PATCH);

// clang-format off
static_assert(!KOKKOS_VERSION_LESS            (KOKKOS_VERSION_MAJOR    , KOKKOS_VERSION_MINOR, KOKKOS_VERSION_PATCH));
static_assert(!KOKKOS_VERSION_LESS            (KOKKOS_VERSION_MAJOR - 1, KOKKOS_VERSION_MINOR, KOKKOS_VERSION_PATCH));
static_assert( KOKKOS_VERSION_LESS            (KOKKOS_VERSION_MAJOR + 1, KOKKOS_VERSION_MINOR, KOKKOS_VERSION_PATCH));

static_assert( KOKKOS_VERSION_LESS_EQUAL      (KOKKOS_VERSION_MAJOR    , KOKKOS_VERSION_MINOR, KOKKOS_VERSION_PATCH));
static_assert(!KOKKOS_VERSION_LESS_EQUAL      (KOKKOS_VERSION_MAJOR - 1, KOKKOS_VERSION_MINOR, KOKKOS_VERSION_PATCH));
static_assert( KOKKOS_VERSION_LESS_EQUAL      (KOKKOS_VERSION_MAJOR + 1, KOKKOS_VERSION_MINOR, KOKKOS_VERSION_PATCH));

static_assert(!KOKKOS_VERSION_GREATER         (KOKKOS_VERSION_MAJOR    , KOKKOS_VERSION_MINOR, KOKKOS_VERSION_PATCH));
static_assert( KOKKOS_VERSION_GREATER         (KOKKOS_VERSION_MAJOR - 1, KOKKOS_VERSION_MINOR, KOKKOS_VERSION_PATCH));
static_assert(!KOKKOS_VERSION_GREATER         (KOKKOS_VERSION_MAJOR + 1, KOKKOS_VERSION_MINOR, KOKKOS_VERSION_PATCH));

static_assert( KOKKOS_VERSION_GREATER_EQUAL   (KOKKOS_VERSION_MAJOR    , KOKKOS_VERSION_MINOR, KOKKOS_VERSION_PATCH));
static_assert( KOKKOS_VERSION_GREATER_EQUAL   (KOKKOS_VERSION_MAJOR - 1, KOKKOS_VERSION_MINOR, KOKKOS_VERSION_PATCH));
static_assert(!KOKKOS_VERSION_GREATER_EQUAL   (KOKKOS_VERSION_MAJOR + 1, KOKKOS_VERSION_MINOR, KOKKOS_VERSION_PATCH));

static_assert( KOKKOS_VERSION_EQUAL           (KOKKOS_VERSION_MAJOR    , KOKKOS_VERSION_MINOR, KOKKOS_VERSION_PATCH));
static_assert(!KOKKOS_VERSION_EQUAL           (KOKKOS_VERSION_MAJOR - 1, KOKKOS_VERSION_MINOR, KOKKOS_VERSION_PATCH));
static_assert(!KOKKOS_VERSION_EQUAL           (KOKKOS_VERSION_MAJOR + 1, KOKKOS_VERSION_MINOR, KOKKOS_VERSION_PATCH));
// clang-format on
