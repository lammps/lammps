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

int get_device_count() {
#if defined(KOKKOS_ENABLE_CUDA)
  int count;
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaGetDeviceCount(&count));
  return count;
#elif defined(KOKKOS_ENABLE_HIP)
  int count;
  KOKKOS_IMPL_HIP_SAFE_CALL(hipGetDevice(&count));
  return count;
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
  return omp_get_num_devices();
#elif defined(KOKKOS_ENABLE_OPENACC)
  return acc_get_num_devices(acc_get_device_type());
#else
  return 0;
#endif
}

int get_device_id() {
  int device_id;
#if defined(KOKKOS_ENABLE_CUDA)
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaGetDevice(&device_id));
#elif defined(KOKKOS_ENABLE_HIP)
  KOKKOS_IMPL_HIP_SAFE_CALL(hipGetDevice(&device_id));
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
  device_id = omp_get_device_num();
#elif defined(KOKKOS_ENABLE_OPENACC)
  device_id = acc_get_device_num(acc_get_device_type());
#elif defined(KOKKOS_ENABLE_SYCL)
  // FIXME_SYCL ?
  assert(false);
  return -2;
#else
  device_id = -1;
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
  KOKKOS_TEST_PRINT_FLAG(device_count);
  KOKKOS_TEST_PRINT_FLAG(disable_warnings);
  KOKKOS_TEST_PRINT_FLAG(tune_internals);

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
