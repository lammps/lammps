/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>
#include <PerfTest_Category.hpp>
#include <Kokkos_Random.hpp>

#ifdef KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA
namespace Test {
template<class Scalar>
void custom_reduction_test(int N, int R, int num_trials) {
  Kokkos::Random_XorShift64_Pool<> rand_pool(183291);
  Kokkos::View<Scalar*> a("A",N);
  Kokkos::fill_random(a,rand_pool,1.0);

  Scalar max;

  // Warm up
  Kokkos::parallel_reduce(Kokkos::TeamPolicy<>(N/1024,32), KOKKOS_LAMBDA( const Kokkos::TeamPolicy<>::member_type& team, Scalar& lmax) {
    Scalar team_max = Scalar(0);
    for(int rr = 0; rr<R; rr++) {
    int i = team.league_rank();
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team,32), [&] (const int& j, Scalar& thread_max) {
      Scalar t_max = Scalar(0);
      Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team,32), [&] (const int& k, Scalar& max_) {
        const Scalar val =  a((i*32 + j)*32 + k);
        if(val>lmax) lmax = val;
        if((k == 11) && (j==17) && (i==2)) lmax = 11.5;
      },Kokkos::Experimental::Max<Scalar>(t_max));
      if(t_max>thread_max) thread_max = t_max;
    },Kokkos::Experimental::Max<Scalar>(team_max));
    }
    if(team_max>lmax) lmax = team_max;
  },Kokkos::Experimental::Max<Scalar>(max));

  // Timing
  Kokkos::Timer timer;
  for(int r = 0; r<num_trials; r++) {
    Kokkos::parallel_reduce(Kokkos::TeamPolicy<>(N/1024,32), KOKKOS_LAMBDA( const Kokkos::TeamPolicy<>::member_type& team, Scalar& lmax) {
      Scalar team_max = Scalar(0);
      for(int rr = 0; rr<R; rr++) {
      int i = team.league_rank();
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team,32), [&] (const int& j, Scalar& thread_max) {
        Scalar t_max = Scalar(0);
        Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team,32), [&] (const int& k, Scalar& max_) {
          const Scalar val =  a((i*32 + j)*32 + k);
          if(val>lmax) lmax = val;
          if((k == 11) && (j==17) && (i==2)) lmax = 11.5;
        },Kokkos::Experimental::Max<Scalar>(t_max));
        if(t_max>thread_max) thread_max = t_max;
      },Kokkos::Experimental::Max<Scalar>(team_max));
      }
      if(team_max>lmax) lmax = team_max;
    },Kokkos::Experimental::Max<Scalar>(max));
  }
  double time = timer.seconds();
  printf("%e %e %e\n",time,1.0*N*R*num_trials*sizeof(Scalar)/time/1024/1024/1024,max);
}

TEST_F( default_exec, custom_reduction ) {
  int N = 100000;
  int R = 1000;
  int num_trials = 1;

  if(command_line_num_args()>1)
    N = atoi(command_line_arg(1));
  if(command_line_num_args()>2)
    R = atoi(command_line_arg(2));
  if(command_line_num_args()>3)
    num_trials = atoi(command_line_arg(3));
  custom_reduction_test<double>(N,R,num_trials);
}
}
#endif
