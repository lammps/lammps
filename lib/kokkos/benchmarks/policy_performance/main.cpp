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
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <Kokkos_Core.hpp>
#include "policy_perf_test.hpp"

int main(int argc, char* argv[] ) {
  Kokkos::initialize(argc,argv);

  if(argc<10) {
    printf("  Ten arguments are needed to run this program:\n");
    printf("    (1)team_range, (2)thread_range, (3)vector_range, (4)outer_repeat, (5)thread_repeat, (6)vector_repeat, (7)team_size, (8)vector_size, (9)schedule, (10)test_type\n");
    printf("  team_range:     number of teams (league_size)\n");
    printf("  thread_range:   range for nested TeamThreadRange parallel_*\n");
    printf("  vector_range:   range for nested ThreadVectorRange parallel_*\n");
    printf("  outer_repeat:   number of repeats for outer parallel_* call\n");
    printf("  thread_repeat:  number of repeats for TeamThreadRange parallel_* call\n");
    printf("  vector_repeat:  number of repeats for ThreadVectorRange parallel_* call\n");
    printf("  team_size:      number of team members (team_size)\n");
    printf("  vector_size:    desired vectorization (if possible)\n");
    printf("  schedule:       1 == Static  2 == Dynamic\n");
    printf("  test_type:      3-digit code XYZ for testing (nested) parallel_*\n");
    printf("  code key:       XYZ    X in {1,2,3,4,5}, Y in {0,1,2}, Z in {0,1,2}\n");
    printf("                  TeamPolicy:\n");
    printf("                    X: 0 = none (never used, makes no sense); 1 = parallel_for; 2 = parallel_reduce\n");
    printf("                    Y: 0 = none; 1 = parallel_for; 2 = parallel_reduce\n");
    printf("                    Z: 0 = none; 1 = parallel_for; 2 = parallel_reduce\n");
    printf("                  RangePolicy:\n");
    printf("                    X: 3 = parallel_for; 4 = parallel_reduce; 5 = parallel_scan\n");
    printf("                    Y: 0 = none\n");
    printf("                    Z: 0 = none\n");
    printf("  Example Input:\n");
    printf("  100000 32 32 100 100 100 8 1 1 100\n"); 
    Kokkos::finalize();
    return 0;
  }

  int team_range = atoi(argv[1]);
  int thread_range = atoi(argv[2]);
  int vector_range = atoi(argv[3]);

  int outer_repeat = atoi(argv[4]);
  int thread_repeat = atoi(argv[5]);
  int vector_repeat = atoi(argv[6]);

  int team_size = atoi(argv[7]);
  int vector_size = atoi(argv[8]);
  int schedule = atoi(argv[9]);
  int test_type = atoi(argv[10]);

  int disable_verbose_output = 0; 
  if ( argc > 11 ) {
    disable_verbose_output = atoi(argv[11]);
  }

  if ( schedule != 1 && schedule != 2 ) {
    printf("schedule: %d\n", schedule);
    printf("Options for schedule are: 1 == Static  2 == Dynamic\n");
    Kokkos::finalize();
    return -1;
  }

  if ( test_type != 100 && test_type != 110 && test_type != 111 && test_type != 112 && test_type != 120  && test_type != 121  && test_type != 122
     && test_type != 200 && test_type != 210 && test_type != 211 && test_type != 212 && test_type != 220  && test_type != 221  && test_type != 222
     && test_type != 300 && test_type != 400 && test_type != 500
     )
  {
    printf("Incorrect test_type option\n");
    Kokkos::finalize();
    return -2;
  }

  double result = 0.0;

  Kokkos::parallel_reduce( "parallel_reduce warmup", Kokkos::TeamPolicy<>(10,1), 
    KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type team, double& lval) {
      lval += 1;
    }, result);

  typedef Kokkos::View<double*, Kokkos::LayoutRight>   view_type_1d;
  typedef Kokkos::View<double**, Kokkos::LayoutRight>  view_type_2d;
  typedef Kokkos::View<double***, Kokkos::LayoutRight> view_type_3d;

  // Allocate view without initializing
  // Call a 'warmup' test with 1 repeat - this will initialize the corresponding view appropriately for test and should obey first-touch etc
  // Second call to test is the one we actually care about and time
  view_type_1d v_1( Kokkos::ViewAllocateWithoutInitializing("v_1"), team_range*team_size);
  view_type_2d v_2( Kokkos::ViewAllocateWithoutInitializing("v_2"), team_range*team_size, thread_range);
  view_type_3d v_3( Kokkos::ViewAllocateWithoutInitializing("v_3"), team_range*team_size, thread_range, vector_range);

  double result_computed = 0.0;
  double result_expect = 0.0;
  double time = 0.0;

  if(schedule==1) {
    if ( test_type != 500 ) {
      // warmup - no repeat of loops
      test_policy<Kokkos::Schedule<Kokkos::Static>,int>(team_range,thread_range,vector_range,1,1,1,team_size,vector_size,test_type,v_1,v_2,v_3,result_computed,result_expect,time);
      test_policy<Kokkos::Schedule<Kokkos::Static>,int>(team_range,thread_range,vector_range,outer_repeat,thread_repeat,vector_repeat,team_size,vector_size,test_type,v_1,v_2,v_3,result_computed,result_expect,time);
    }
    else {
      // parallel_scan: initialize 1d view for parallel_scan
      test_policy<Kokkos::Schedule<Kokkos::Static>,int>(team_range,thread_range,vector_range,1,1,1,team_size,vector_size,100,v_1,v_2,v_3,result_computed,result_expect,time);
      test_policy<Kokkos::Schedule<Kokkos::Static>,int>(team_range,thread_range,vector_range,outer_repeat,thread_repeat,vector_repeat,team_size,vector_size,test_type,v_1,v_2,v_3,result_computed,result_expect,time);
    }
  }
  if(schedule==2) {
    if ( test_type != 500 ) {
      // warmup - no repeat of loops
      test_policy<Kokkos::Schedule<Kokkos::Dynamic>,int>(team_range,thread_range,vector_range,1,1,1,team_size,vector_size,test_type,v_1,v_2,v_3,result_computed,result_expect,time);
      test_policy<Kokkos::Schedule<Kokkos::Dynamic>,int>(team_range,thread_range,vector_range,outer_repeat,thread_repeat,vector_repeat,team_size,vector_size,test_type,v_1,v_2,v_3,result_computed,result_expect,time);
    }
    else {
      // parallel_scan: initialize 1d view for parallel_scan
      test_policy<Kokkos::Schedule<Kokkos::Static>,int>(team_range,thread_range,vector_range,1,1,1,team_size,vector_size,100,v_1,v_2,v_3,result_computed,result_expect,time);
      test_policy<Kokkos::Schedule<Kokkos::Static>,int>(team_range,thread_range,vector_range,outer_repeat,thread_repeat,vector_repeat,team_size,vector_size,test_type,v_1,v_2,v_3,result_computed,result_expect,time);
    }
  }

  if ( disable_verbose_output == 0 ) {
    printf("%7i %4i %2i %9i %4i %4i %4i %2i %1i %3i %e %e %lf\n",team_range,thread_range,vector_range,outer_repeat,thread_repeat,vector_repeat,team_size,vector_size,schedule,test_type,result_computed,result_expect,time);
  }
  else {
    printf("%lf\n",time);
  }

  Kokkos::finalize();

  return 0;
}
