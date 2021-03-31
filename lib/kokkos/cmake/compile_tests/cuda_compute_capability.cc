/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <iostream>

int main() {
  cudaDeviceProp device_properties;
  const cudaError_t error = cudaGetDeviceProperties(&device_properties,
                                                    /*device*/ 0);
  if (error != cudaSuccess) {
    std::cout << "CUDA error: " << cudaGetErrorString(error) << '\n';
    return error;
  }
  unsigned int const compute_capability =
      device_properties.major * 10 + device_properties.minor;
#ifdef SM_ONLY
  std::cout << compute_capability;
#else
  switch (compute_capability) {
      // clang-format off
    case 30: std::cout << "Set -DKokkos_ARCH_KEPLER30=ON ." << std::endl; break;
    case 32: std::cout << "Set -DKokkos_ARCH_KEPLER32=ON ." << std::endl; break;
    case 35: std::cout << "Set -DKokkos_ARCH_KEPLER35=ON ." << std::endl; break;
    case 37: std::cout << "Set -DKokkos_ARCH_KEPLER37=ON ." << std::endl; break;
    case 50: std::cout << "Set -DKokkos_ARCH_MAXWELL50=ON ." << std::endl; break;
    case 52: std::cout << "Set -DKokkos_ARCH_MAXWELL52=ON ." << std::endl; break;
    case 53: std::cout << "Set -DKokkos_ARCH_MAXWELL53=ON ." << std::endl; break;
    case 60: std::cout << "Set -DKokkos_ARCH_PASCAL60=ON ." << std::endl; break;
    case 61: std::cout << "Set -DKokkos_ARCH_PASCAL61=ON ." << std::endl; break;
    case 70: std::cout << "Set -DKokkos_ARCH_VOLTA70=ON ." << std::endl; break;
    case 72: std::cout << "Set -DKokkos_ARCH_VOLTA72=ON ." << std::endl; break;
    case 75: std::cout << "Set -DKokkos_ARCH_TURING75=ON ." << std::endl; break;
    case 80: std::cout << "Set -DKokkos_ARCH_AMPERE80=ON ." << std::endl; break;
    default:
      std::cout << "Compute capability " << compute_capability
                << " is not supported" << std::endl;
      // clang-format on
  }
#endif
  return 0;
}
