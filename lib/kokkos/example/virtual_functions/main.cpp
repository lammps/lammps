// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <classes.hpp>

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);

  {
    Foo* f_1 = static_cast<Foo*>(Kokkos::kokkos_malloc(sizeof(Foo_1)));
    Foo* f_2 = static_cast<Foo*>(Kokkos::kokkos_malloc(sizeof(Foo_2)));

    Kokkos::parallel_for(
        "CreateObjects", 1, KOKKOS_LAMBDA(const int&) {
          new (f_1) Foo_1();
          new (f_2) Foo_2();
        });

    int value_1, value_2;
    Kokkos::parallel_reduce(
        "CheckValues", 1,
        KOKKOS_LAMBDA(const int&, int& lsum) { lsum = f_1->value(); }, value_1);

    Kokkos::parallel_reduce(
        "CheckValues", 1,
        KOKKOS_LAMBDA(const int&, int& lsum) { lsum = f_2->value(); }, value_2);

    printf("Values: %i %i\n", value_1, value_2);

    Kokkos::parallel_for(
        "DestroyObjects", 1, KOKKOS_LAMBDA(const int&) {
          f_1->~Foo();
          f_2->~Foo();
        });

    Kokkos::kokkos_free(f_1);
    Kokkos::kokkos_free(f_2);
  }

  Kokkos::finalize();
}
