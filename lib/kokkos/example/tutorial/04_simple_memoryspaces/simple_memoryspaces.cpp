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
#include <cstdio>

// The type of a two-dimensional N x 3 array of double.
// It lives in Kokkos' default memory space.
using view_type = Kokkos::View<double * [3]>;

// The "HostMirror" type corresponding to view_type above is also a
// two-dimensional N x 3 array of double.  However, it lives in the
// host memory space corresponding to view_type's memory space.  For
// example, if view_type lives in CUDA device memory, host_view_type
// lives in host (CPU) memory.  Furthermore, declaring host_view_type
// as the host mirror of view_type means that host_view_type has the
// same layout as view_type.  This makes it easier to copy between the
// two Views.
// Advanced issues: If a memory space is accessible from the host without
// performance penalties then it is its own host_mirror_space. This is
// the case for HostSpace, CudaUVMSpace and CudaHostPinnedSpace.

using host_view_type = view_type::HostMirror;

struct ReduceFunctor {
  view_type a;
  ReduceFunctor(view_type a_) : a(a_) {}
  using value_type = int;  // Specify type for reduction value, lsum

  KOKKOS_INLINE_FUNCTION
  void operator()(int i, int &lsum) const {
    lsum += a(i, 0) - a(i, 1) + a(i, 2);
  }
};

int main() {
  Kokkos::initialize();

  {
    view_type a("A", 10);
    // If view_type and host_mirror_type live in the same memory space,
    // a "mirror view" is just an alias, and deep_copy does nothing.
    // Otherwise, a mirror view of a device View lives in host memory,
    // and deep_copy does a deep copy.
    host_view_type h_a = Kokkos::create_mirror_view(a);

    // The View h_a lives in host (CPU) memory, so it's legal to fill
    // the view sequentially using ordinary code, like this.
    for (int i = 0; i < 10; i++) {
      for (int j = 0; j < 3; j++) {
        h_a(i, j) = i * 10 + j;
      }
    }
    Kokkos::deep_copy(a, h_a);  // Copy from host to device.

    int sum = 0;
    Kokkos::parallel_reduce(10, ReduceFunctor(a), sum);
    printf("Result is %i\n", sum);
  }

  Kokkos::finalize();
}
