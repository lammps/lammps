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
#include <cstdio>

extern "C" void print_fortran_();

struct CountFunctor {
  KOKKOS_FUNCTION void operator()(const long i, long& lcount) const {
    lcount += (i % 2) == 0;
  }
};

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  Kokkos::DefaultExecutionSpace().print_configuration(std::cout);

  if (argc < 2) {
    fprintf(stderr, "Usage: %s [<kokkos_options>] <size>\n", argv[0]);
    Kokkos::finalize();
    exit(1);
  }

  const long n = strtol(argv[1], nullptr, 10);

  printf("Number of even integers from 0 to %ld\n", n - 1);

  Kokkos::Timer timer;
  timer.reset();

  // Compute the number of even integers from 0 to n-1, in parallel.
  long count = 0;
  CountFunctor functor;
  Kokkos::parallel_reduce(n, functor, count);

  double count_time = timer.seconds();
  printf("  Parallel: %ld    %10.6f\n", count, count_time);

  timer.reset();

  // Compare to a sequential loop.
  long seq_count = 0;
  for (long i = 0; i < n; ++i) {
    seq_count += (i % 2) == 0;
  }

  count_time = timer.seconds();
  printf("Sequential: %ld    %10.6f\n", seq_count, count_time);

  print_fortran_();

  Kokkos::finalize();

  return (count == seq_count) ? 0 : -1;
}
