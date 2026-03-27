// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_EXAMPLE_VIRTUAL_FUNCTIONS_CLASSES_HPP
#define KOKKOS_EXAMPLE_VIRTUAL_FUNCTIONS_CLASSES_HPP

#include <Kokkos_Core.hpp>

class Foo {
 protected:
  int val;

 public:
  KOKKOS_FUNCTION Foo();

  KOKKOS_FUNCTION virtual int value() { return 0; }

  KOKKOS_FUNCTION virtual ~Foo() {}
};

class Foo_1 : public Foo {
 public:
  KOKKOS_FUNCTION Foo_1();

  KOKKOS_FUNCTION int value() override;
};

class Foo_2 : public Foo {
 public:
  KOKKOS_FUNCTION Foo_2();

  KOKKOS_FUNCTION int value() override;
};

#endif  // KOKKOS_EXAMPLE_VIRTUAL_FUNCTIONS_CLASSES_HPP
