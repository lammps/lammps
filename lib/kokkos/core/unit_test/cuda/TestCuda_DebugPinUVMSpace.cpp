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
#include <TestCuda_Category.hpp>

namespace Test {

template <class View>
struct CopyFunctor {
  View a;
  View b;

  CopyFunctor(int N) : a(View("A", N)), b(View("B", N)) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const { a(i) = b(i); }

  double time_copy(int R) {
    Kokkos::parallel_for("CopyFunctor::time_copy", a.extent(0), *this);
    Kokkos::fence();

    Kokkos::Timer timer;
    for (int r = 0; r < R; r++)
      Kokkos::parallel_for("CopyFunctor::time_copy", a.extent(0), *this);
    Kokkos::fence();
    return timer.seconds();
  }
};

TEST(cuda, debug_pin_um_to_host) {
  double time_cuda_space;
  double time_cuda_host_pinned_space;
  double time_cuda_uvm_space_not_pinned_1;
  double time_cuda_uvm_space_pinned;
  double time_cuda_uvm_space_not_pinned_2;

  int N = 10000000;
  int R = 100;
  {
    CopyFunctor<Kokkos::View<int*, Kokkos::CudaSpace>> f(N);
    time_cuda_space = f.time_copy(R);
  }
  {
    CopyFunctor<Kokkos::View<int*, Kokkos::CudaHostPinnedSpace>> f(N);
    time_cuda_host_pinned_space = f.time_copy(R);
  }
  {
    CopyFunctor<Kokkos::View<int*, Kokkos::CudaUVMSpace>> f(N);
    time_cuda_uvm_space_not_pinned_1 = f.time_copy(R);
  }
  {
#ifdef KOKKOS_IMPL_DEBUG_CUDA_PIN_UVM_TO_HOST
    kokkos_impl_cuda_set_pin_uvm_to_host(true);
#endif
    CopyFunctor<Kokkos::View<int*, Kokkos::CudaUVMSpace>> f(N);
    time_cuda_uvm_space_pinned = f.time_copy(R);
#ifdef KOKKOS_IMPL_DEBUG_CUDA_PIN_UVM_TO_HOST
    kokkos_impl_cuda_set_pin_uvm_to_host(false);
#endif
  }
  {
    CopyFunctor<Kokkos::View<int*, Kokkos::CudaUVMSpace>> f(N);
    time_cuda_uvm_space_not_pinned_2 = f.time_copy(R);
  }
  bool uvm_approx_cuda_1 =
      time_cuda_uvm_space_not_pinned_1 < time_cuda_space * 2.0;
  bool uvm_approx_cuda_2 =
      time_cuda_uvm_space_not_pinned_2 < time_cuda_space * 2.0;
  bool pinned_slower_cuda = time_cuda_host_pinned_space > time_cuda_space * 2.0;
  bool uvm_pinned_slower_cuda =
      time_cuda_uvm_space_pinned > time_cuda_space * 2.0;

  bool passed = uvm_approx_cuda_1 && uvm_approx_cuda_2 && pinned_slower_cuda &&
#ifdef KOKKOS_IMPL_DEBUG_CUDA_PIN_UVM_TO_HOST
                uvm_pinned_slower_cuda;
#else
                !uvm_pinned_slower_cuda;
#endif
  if (!passed)
    printf(
        "Time CudaSpace: %lf CudaUVMSpace_1: %lf CudaUVMSpace_2: %lf "
        "CudaPinnedHostSpace: %lf CudaUVMSpace_Pinned: %lf\n",
        time_cuda_space, time_cuda_uvm_space_not_pinned_1,
        time_cuda_uvm_space_not_pinned_2, time_cuda_host_pinned_space,
        time_cuda_uvm_space_pinned);
  ASSERT_TRUE(passed);
}

}  // namespace Test
