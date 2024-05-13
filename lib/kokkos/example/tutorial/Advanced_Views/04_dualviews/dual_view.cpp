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

// DualView helps you manage data and computations that take place on
// two different memory spaces.  Examples include CUDA device memory
// and (CPU) host memory (currently implemented), or Intel Knights
// Landing MCDRAM and DRAM (not yet implemented).  For example, if you
// have ported only some parts of you application to run in CUDA,
// DualView can help manage moving data between the parts of your
// application that work best with CUDA, and the parts that work
// better on the CPU.
//
// A DualView takes the same template parameters as a View, but
// contains two Views: One that lives in the DualView's memory space,
// and one that lives in that memory space's host mirror space.  If
// both memory spaces are the same, then the two Views just alias one
// another.  This means that you can use DualView all the time, even
// when not running in a memory space like CUDA.  DualView's
// operations to help you manage memory take almost no time in that
// case.  This makes your code even more performance portable.

using view_type = Kokkos::DualView<double*>;
using idx_type  = Kokkos::DualView<int**>;

template <class ExecutionSpace>
struct localsum {
  // If the functor has a public 'execution_space' alias, that defines
  // the functor's execution space (where it runs in parallel).  This
  // overrides Kokkos' default execution space.
  using execution_space = ExecutionSpace;

  using memory_space = typename Kokkos::Impl::if_c<
      std::is_same<ExecutionSpace, Kokkos::DefaultExecutionSpace>::value,
      idx_type::memory_space, idx_type::host_mirror_space>::type;

  // Get the view types on the particular device for which the functor
  // is instantiated.
  //
  // "const_data_type" is an alias in View (and DualView) which is
  // the const version of the first template parameter of the View.
  // For example, the const_data_type version of double** is const
  // double**.
  Kokkos::View<idx_type::const_data_type, idx_type::array_layout, memory_space>
      idx;
  // "scalar_array_type" is an alias in ViewTraits (and DualView) which is the
  // array version of the value(s) stored in the View.
  Kokkos::View<view_type::scalar_array_type, view_type::array_layout,
               memory_space>
      dest;
  Kokkos::View<view_type::const_data_type, view_type::array_layout,
               memory_space, Kokkos::MemoryRandomAccess>
      src;

  // Constructor takes DualViews, synchronizes them to the device,
  // then marks them as modified on the device.
  localsum(idx_type dv_idx, view_type dv_dest, view_type dv_src) {
    // Extract the view on the correct Device (i.e., the correct
    // memory space) from the DualView.  DualView has a template
    // method, view(), which is templated on the memory space.  If the
    // DualView has a View from that memory space, view() returns the
    // View in that space.
    idx  = dv_idx.view<memory_space>();
    dest = dv_dest.template view<memory_space>();
    src  = dv_src.template view<memory_space>();

    // Synchronize the DualView to the correct Device.
    //
    // DualView's sync() method is templated on a memory space, and
    // synchronizes the DualView in a one-way fashion to that memory
    // space.  "Synchronizing" means copying, from the other memory
    // space to the Device memory space.  sync() does _nothing_ if the
    // Views on the two memory spaces are in sync.  DualView
    // determines this by the user manually marking one side or the
    // other as modified; see the modify() call below.

    dv_idx.sync<memory_space>();
    dv_dest.template sync<memory_space>();
    dv_src.template sync<memory_space>();

    // Mark dest as modified on Device.
    dv_dest.template modify<memory_space>();
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    double tmp = 0.0;
    for (int j = 0; j < (int)idx.extent(1); ++j) {
      const double val = src(idx(i, j));
      tmp += val * val + 0.5 * (idx.extent(0) * val - idx.extent(1) * val);
    }
    dest(i) += tmp;
  }
};

class ParticleType {
 public:
  double q;
  double m;
  double q_over_m;
  KOKKOS_INLINE_FUNCTION
  ParticleType(double q_ = -1, double m_ = 1) : q(q_), m(m_), q_over_m(q / m) {}

 protected:
};

using ParticleTypes = Kokkos::DualView<ParticleType[10]>;
int main(int narg, char* arg[]) {
  Kokkos::initialize(narg, arg);

  // If View is non-trivial constructible type then add braces so it is out of
  // scope before Kokkos::finalize() call
  {
    ParticleTypes test("Test");
    Kokkos::fence();
    test.h_view(0) = ParticleType(-1e4, 1);
    Kokkos::fence();

    int size = 1000000;

    // Create DualViews. This will allocate on both the device and its
    // host_mirror_device.
    idx_type idx("Idx", size, 64);
    view_type dest("Dest", size);
    view_type src("Src", size);

    srand(134231);

    // Get a reference to the host view of idx directly (equivalent to
    // idx.view<idx_type::host_mirror_space>() )
    idx_type::t_host h_idx = idx.h_view;
    using size_type        = view_type::size_type;
    for (int i = 0; i < size; ++i) {
      for (size_type j = 0; j < static_cast<size_type>(h_idx.extent(1)); ++j) {
        h_idx(i, j) = (size + i + (rand() % 500 - 250)) % size;
      }
    }

    // Mark idx as modified on the host_mirror_space so that a
    // sync to the device will actually move data.  The sync happens in
    // the functor's constructor.
    idx.modify<idx_type::host_mirror_space>();

    // Run on the device.  This will cause a sync of idx to the device,
    // since it was marked as modified on the host.
    Kokkos::Timer timer;
    Kokkos::parallel_for(size,
                         localsum<view_type::execution_space>(idx, dest, src));
    Kokkos::fence();
    double sec1_dev = timer.seconds();

    timer.reset();
    Kokkos::parallel_for(size,
                         localsum<view_type::execution_space>(idx, dest, src));
    Kokkos::fence();
    double sec2_dev = timer.seconds();

    // Run on the host's default execution space (could be the same as device).
    // This will cause a sync back to the host of dest.  Note that if the Device
    // is CUDA, the data layout will not be optimal on host, so performance is
    // lower than what it would be for a pure host compilation.
    timer.reset();
    Kokkos::parallel_for(
        size, localsum<Kokkos::HostSpace::execution_space>(idx, dest, src));
    Kokkos::fence();
    double sec1_host = timer.seconds();

    timer.reset();
    Kokkos::parallel_for(
        size, localsum<Kokkos::HostSpace::execution_space>(idx, dest, src));
    Kokkos::fence();
    double sec2_host = timer.seconds();

    printf("Device Time with Sync: %f without Sync: %f \n", sec1_dev, sec2_dev);
    printf("Host   Time with Sync: %f without Sync: %f \n", sec1_host,
           sec2_host);
  }

  Kokkos::finalize();
}
