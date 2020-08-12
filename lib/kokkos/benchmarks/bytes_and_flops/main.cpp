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

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>
#include <bench.hpp>
#include <cstdlib>

int main(int argc, char* argv[]) {
  Kokkos::initialize();

  if (argc < 10) {
    printf("Arguments: N K R D U F T S\n");
    printf("  P:   Precision (1==float, 2==double)\n");
    printf("  N,K: dimensions of the 2D array to allocate\n");
    printf("  R:   how often to loop through the K dimension with each team\n");
    printf("  D:   distance between loaded elements (stride)\n");
    printf("  U:   how many independent flops to do per load\n");
    printf(
        "  F:   how many times to repeat the U unrolled operations before "
        "reading next element\n");
    printf("  T:   team size\n");
    printf(
        "  S:   shared memory per team (used to control occupancy on GPUs)\n");
    printf("Example Input GPU:\n");
    printf("  Bandwidth Bound : 2 100000 1024 1 1 1 1 256 6000\n");
    printf("  Cache Bound     : 2 100000 1024 64 1 1 1 512 20000\n");
    printf("  Compute Bound   : 2 100000 1024 1 1 8 64 256 6000\n");
    printf("  Load Slots Used : 2 20000 256 32 16 1 1 256 6000\n");
    printf("  Inefficient Load: 2 20000 256 32 2 1 1 256 20000\n");
    Kokkos::finalize();
    return 0;
  }

  int P = atoi(argv[1]);
  int N = atoi(argv[2]);
  int K = atoi(argv[3]);
  int R = atoi(argv[4]);
  int D = atoi(argv[5]);
  int U = atoi(argv[6]);
  int F = atoi(argv[7]);
  int T = atoi(argv[8]);
  int S = atoi(argv[9]);

  if (U > 8) {
    printf("U must be 1-8\n");
    return 0;
  }
  if ((D != 1) && (D != 2) && (D != 4) && (D != 8) && (D != 16) && (D != 32)) {
    printf("D must be one of 1,2,4,8,16,32\n");
    return 0;
  }
  if ((P != 1) && (P != 2)) {
    printf("P must be one of 1,2\n");
    return 0;
  }

  if (P == 1) {
    run_stride_unroll<float>(N, K, R, D, U, F, T, S);
  }
  if (P == 2) {
    run_stride_unroll<double>(N, K, R, D, U, F, T, S);
  }

  Kokkos::finalize();
}
