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

#include <iostream>
#include <hip/hip_runtime_api.h>

int main() {
  hipDeviceProp_t hipProp;
  hipError_t error = hipGetDeviceProperties(&hipProp, 0);

  if (error != hipSuccess) {
    std::cout << hipGetErrorString(error) << '\n';
    return error;
  }

  if (hipProp.integrated == 1) {
    // We detected an APU
    std::cout << "ON";
  } else {
    // We detected a discrete GPU
    std::cout << "OFF";
  }

  return 0;
}
