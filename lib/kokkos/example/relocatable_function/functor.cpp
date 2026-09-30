// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>

KOKKOS_RELOCATABLE_FUNCTION void count_even(const long i, long& lcount) {
  lcount += (i % 2) == 0;
}
