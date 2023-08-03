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

#ifndef KOKKOS_EXAMPLE_VIRTUAL_FUNCTIONS_CLASSES_HPP
#define KOKKOS_EXAMPLE_VIRTUAL_FUNCTIONS_CLASSES_HPP

#include <Kokkos_Core.hpp>

class Foo {
 protected:
  int val;

 public:
  KOKKOS_FUNCTION
  Foo();

  KOKKOS_FUNCTION
  virtual int value() { return 0; };

  KOKKOS_FUNCTION
  virtual ~Foo() {}
};

class Foo_1 : public Foo {
 public:
  KOKKOS_FUNCTION
  Foo_1();

  KOKKOS_FUNCTION
  int value();
};

class Foo_2 : public Foo {
 public:
  KOKKOS_FUNCTION
  Foo_2();

  KOKKOS_FUNCTION
  int value();
};

#endif  // KOKKOS_EXAMPLE_VIRTUAL_FUNCTIONS_CLASSES_HPP
