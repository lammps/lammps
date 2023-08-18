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
#include <typeinfo>
#include <cmath>
#include <Kokkos_Timer.hpp>

struct FillDevice {
  double value;
  Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::CudaSpace> a;
  FillDevice(
      const double& val,
      const Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::CudaSpace>& d_a)
      : value(val), a(d_a) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i) const { a(i) = value; }
};

struct ComputeADevice {
  int iter;
  Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::CudaSpace> a;
  Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::CudaSpace> b;
  ComputeADevice(
      const int& iter_,
      const Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::CudaSpace>& d_a,
      const Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::CudaSpace>& d_b)
      : iter(iter_), a(d_a), b(d_b) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i) const {
    for (int j = 1; j < iter; j++) {
      a(i) += std::pow(b(i), 1.0 + 1.0 / iter);
    }
  }
};

struct ComputeAHost {
  Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::CudaHostPinnedSpace> a;
  Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::CudaHostPinnedSpace> b;
  ComputeAHost(const Kokkos::View<double*, Kokkos::LayoutLeft,
                                  Kokkos::CudaHostPinnedSpace>& d_a,
               const Kokkos::View<double*, Kokkos::LayoutLeft,
                                  Kokkos::CudaHostPinnedSpace>& d_b)
      : a(d_a), b(d_b) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i) const { a(i) += b(i); }
};

struct MergeDevice {
  Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::CudaSpace> a;
  Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::CudaSpace> b;
  MergeDevice(
      const Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::CudaSpace>& d_a,
      const Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::CudaSpace>& d_b)
      : a(d_a), b(d_b) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i) const { a(i) += b(i); }
};

int main(int argc, char* argv[]) {
  int size = 100000000;
  Kokkos::initialize();
  int synch = std::stoi(argv[1]);
  Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::CudaSpace> d_a("Device A",
                                                                   size);
  Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::CudaSpace> d_b("Device B",
                                                                   size);
  Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::CudaSpace> d_tmp(
      "Device tmp", size);
  Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::CudaHostPinnedSpace> h_a(
      "Host A", size);
  Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::CudaHostPinnedSpace> h_b(
      "Host B", size);

  Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::Cuda>(0, size),
                       FillDevice(0.0, d_a));
  Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::Cuda>(0, size),
                       FillDevice(1.3513, d_b));
  Kokkos::fence();
  Kokkos::Timer timer;
  Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::Cuda>(0, size),
                       ComputeADevice(20, d_a, d_b));

  if (synch == 1) Kokkos::deep_copy(Kokkos::OpenMP(), h_b, d_b);
  if (synch == 2) Kokkos::deep_copy(h_b, d_b);

  Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::OpenMP>(0, size),
                       [=](const int& i) { h_a(i) = 0.0; });
  Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::OpenMP>(0, size),
                       ComputeAHost(h_a, h_b));
  Kokkos::OpenMP().fence();
  if (synch == 1) Kokkos::deep_copy(Kokkos::OpenMP(), d_tmp, h_a);
  if (synch == 2) Kokkos::deep_copy(d_tmp, h_a);
  Kokkos::fence();

  std::cout << "Time " << timer.seconds() << std::endl;
  Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::Cuda>(0, size),
                       MergeDevice(d_a, d_tmp));

  Kokkos::deep_copy(h_a, d_a);
  std::cout << "h_a(0): " << h_a(0) << " ( Correct: 27.4154 )" << std::endl;
  Kokkos::finalize();
}
