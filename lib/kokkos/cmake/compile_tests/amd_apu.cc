// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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
