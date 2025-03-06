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
#include <iostream>
#include <string>
#include <thread>

int get_num_devices() {
  int num_devices;
#if defined(KOKKOS_ENABLE_CUDA)
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaGetDeviceCount(&num_devices));
#elif defined(KOKKOS_ENABLE_HIP)
  KOKKOS_IMPL_HIP_SAFE_CALL(hipGetDeviceCount(&num_devices));
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
  num_devices = omp_get_num_devices();
#elif defined(KOKKOS_ENABLE_OPENACC)
  num_devices = acc_get_num_devices(acc_get_device_type());
#elif defined(KOKKOS_ENABLE_SYCL)
  num_devices = sycl::device::get_devices(sycl::info::device_type::gpu).size();
#else
  num_devices = -1;
#endif
  assert(num_devices == Kokkos::num_devices());
  return num_devices;
}

int get_device_id() {
  int device_id;
#if defined(KOKKOS_ENABLE_CUDA)
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaGetDevice(&device_id));
#elif defined(KOKKOS_ENABLE_HIP)
  KOKKOS_IMPL_HIP_SAFE_CALL(hipGetDevice(&device_id));
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
  device_id   = omp_get_default_device();
#elif defined(KOKKOS_ENABLE_OPENACC)
  device_id   = acc_get_device_num(acc_get_device_type());
#elif defined(KOKKOS_ENABLE_SYCL)
  // Not able to query the underlying runtime because there is no such thing as
  // device currently being used with SYCL.  We go through the Kokkos runtime
  // which makes the assert below pointless but it still let us check that
  // Kokkos selected the device we asked for from the Python tests.
  device_id = Kokkos::device_id();
#else
  device_id   = -1;
#endif
  assert(device_id == Kokkos::device_id());
  return device_id;
}

int get_max_threads() {
#if defined(KOKKOS_ENABLE_OPENMP)
  return omp_get_max_threads();
#elif defined(KOKKOS_ENABLE_THREADS)
  return std::thread::hardware_concurrency();
#else
  return 1;
#endif
}

int get_hwloc_enabled() {
#ifdef KOKKOS_ENABLE_HWLOC
  return 1;
#else
  return 0;
#endif
}

int get_num_threads() {
  int const num_threads = Kokkos::DefaultHostExecutionSpace().concurrency();
  assert(num_threads == Kokkos::num_threads());
  return num_threads;
}

int get_disable_warnings() { return !Kokkos::show_warnings(); }

int get_tune_internals() { return Kokkos::tune_internals(); }

int print_flag(std::string const& flag) {
  std::vector<std::string> valid_flags;
#define KOKKOS_TEST_PRINT_FLAG(NAME)   \
  if (flag == #NAME) {                 \
    std::cout << get_##NAME() << '\n'; \
    return EXIT_SUCCESS;               \
  }                                    \
  valid_flags.push_back(#NAME)

  KOKKOS_TEST_PRINT_FLAG(num_threads);
  KOKKOS_TEST_PRINT_FLAG(max_threads);
  KOKKOS_TEST_PRINT_FLAG(device_id);
  KOKKOS_TEST_PRINT_FLAG(num_devices);
  KOKKOS_TEST_PRINT_FLAG(disable_warnings);
  KOKKOS_TEST_PRINT_FLAG(tune_internals);
  KOKKOS_TEST_PRINT_FLAG(hwloc_enabled);

#undef KOKKOS_TEST_PRINT_FLAG

  std::cerr << "Invalid flag name " << flag << ".  Valid names are ";
  for (int i = 0; i < (int)valid_flags.size() - 1; ++i) {
    std::cerr << valid_flags[i] << ", ";
  }
  std::cerr << "and " << valid_flags.back() << ".\n";
  return EXIT_FAILURE;
}

int main(int argc, char* argv[]) {
  Kokkos::ScopeGuard guard(argc, argv);
  if (argc != 2) {
    std::cerr << "Usage: <executable> NAME_OF_FLAG\n";
    return EXIT_FAILURE;
  }
  int exit_code = print_flag(argv[1]);
  return exit_code;
}
