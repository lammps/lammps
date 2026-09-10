// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <omp.h>

int main(int, char**) {
  int thr = omp_get_num_threads();
  if (thr > 0)
    return thr;
  else
    return 0;
}
