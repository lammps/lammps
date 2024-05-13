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
#include <Kokkos_DualView.hpp>
#include <Kokkos_Timer.hpp>
#include <cstdio>
#include <cstdlib>

#ifdef KOKKOS_ENABLE_CUDA
using view_type = Kokkos::View<double*, Kokkos::CudaUVMSpace>;
using idx_type  = Kokkos::View<int**, Kokkos::CudaUVMSpace>;
#else
using view_type = Kokkos::View<double*, Kokkos::HostSpace>;
using idx_type  = Kokkos::View<int**, Kokkos::HostSpace>;
#endif

template <class Device>
struct localsum {
  // Define the execution space for the functor (overrides the
  // DefaultExecutionSpace)
  using execution_space = Device;

  // Get the view types on the particular device the functor is instantiated for
  idx_type::const_type idx;
  view_type dest;
  Kokkos::View<view_type::const_data_type, view_type::array_layout,
               view_type::device_type, Kokkos::MemoryRandomAccess>
      src;

  localsum(idx_type idx_, view_type dest_, view_type src_)
      : idx(idx_), dest(dest_), src(src_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const {
    double tmp = 0.0;
    for (int j = 0; j < int(idx.extent(1)); j++) {
      const double val = src(idx(i, j));
      tmp += val * val + 0.5 * (idx.extent(0) * val - idx.extent(1) * val);
    }
    dest(i) += tmp;
  }
};

int main(int narg, char* arg[]) {
  Kokkos::initialize(narg, arg);

  {
    int size = 1000000;

    // Create Views
    idx_type idx("Idx", size, 64);
    view_type dest("Dest", size);
    view_type src("Src", size);

    srand(134231);

    Kokkos::fence();

    // When using UVM Cuda views can be accessed on the Host directly
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < int(idx.extent(1)); j++)
        idx(i, j) = (size + i + (rand() % 500 - 250)) % size;
    }

    Kokkos::fence();
    // Run on the device
    // This will cause a sync of idx to the device since it was modified on the
    // host
    Kokkos::Timer timer;
    Kokkos::parallel_for(size,
                         localsum<view_type::execution_space>(idx, dest, src));
    Kokkos::fence();
    double sec1_dev = timer.seconds();

    // No data transfer will happen now, since nothing is accessed on the host
    timer.reset();
    Kokkos::parallel_for(size,
                         localsum<view_type::execution_space>(idx, dest, src));
    Kokkos::fence();
    double sec2_dev = timer.seconds();

    // Run on the host
    // This will cause a sync back to the host of dest which was changed on the
    // device Compare runtime here with the dual_view example: dest will be
    // copied back in 4k blocks when they are accessed the first time during the
    // parallel_for. Due to the latency of a memcpy this gives lower effective
    // bandwidth when doing a manual copy via dual views
    timer.reset();
    Kokkos::parallel_for(
        size, localsum<Kokkos::HostSpace::execution_space>(idx, dest, src));
    Kokkos::fence();
    double sec1_host = timer.seconds();

    // No data transfers will happen now
    timer.reset();
    Kokkos::parallel_for(
        size, localsum<Kokkos::HostSpace::execution_space>(idx, dest, src));
    Kokkos::fence();
    double sec2_host = timer.seconds();

    printf("Device Time with Sync: %e without Sync: %e \n", sec1_dev, sec2_dev);
    printf("Host   Time with Sync: %e without Sync: %e \n", sec1_host,
           sec2_host);
  }

  Kokkos::finalize();
}
