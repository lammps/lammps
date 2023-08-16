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

#include <classes.hpp>

KOKKOS_FUNCTION
Foo::Foo() { val = 0; }

KOKKOS_FUNCTION
Foo_1::Foo_1() { val = 1; }

KOKKOS_FUNCTION
int Foo_1::value() { return val; }

KOKKOS_FUNCTION
Foo_2::Foo_2() { val = 2; }

KOKKOS_FUNCTION
int Foo_2::value() { return val; }
