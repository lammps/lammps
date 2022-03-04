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

#if !defined(KOKKOS_ENABLE_TASKDAG) || \
    defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_THREADS)

int main() { return 0; }

#else

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <limits>

#include <Kokkos_Timer.hpp>

using ExecSpace = Kokkos::DefaultExecutionSpace;

inline long eval_fib(long n) {
  constexpr long mask = 0x03;

  long fib[4] = {0, 1, 0, 0};

  for (long i = 2; i <= n; ++i) {
    fib[i & mask] = fib[(i - 1) & mask] + fib[(i - 2) & mask];
  }

  return fib[n & mask];
}

inline long fib_alloc_count(long n) {
  constexpr long mask = 0x03;

  long count[4] = {1, 1, 0, 0};

  for (long i = 2; i <= n; ++i) {
    count[i & mask] = 2  // this task plus the 'when_all' task
                      + count[(i - 1) & mask] + count[(i - 2) & mask];
  }

  return count[n & mask];
}

template <class Scheduler>
struct TestFib {
  using MemorySpace = typename Scheduler::memory_space;
  using MemberType  = typename Scheduler::member_type;
  using FutureType  = Kokkos::BasicFuture<long, Scheduler>;

  using value_type = long;

  FutureType dep[2];
  const value_type n;

  KOKKOS_INLINE_FUNCTION
  TestFib(const value_type arg_n) : dep{}, n(arg_n) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(MemberType& member, value_type& result) noexcept {
    auto& sched = member.scheduler();
    if (n < 2) {
      result = n;
    } else if (!dep[0].is_null() && !dep[1].is_null()) {
      result = dep[0].get() + dep[1].get();
    } else {
      // Spawn new children and respawn myself to sum their results.
      // Spawn lower value at higher priority as it has a shorter
      // path to completion.

      dep[1] = Kokkos::task_spawn(
          Kokkos::TaskSingle(sched, Kokkos::TaskPriority::High),
          TestFib(n - 2));

      dep[0] = Kokkos::task_spawn(Kokkos::TaskSingle(sched), TestFib(n - 1));

      auto fib_all = sched.when_all(dep, 2);

      if (!dep[0].is_null() && !dep[1].is_null() && !fib_all.is_null()) {
        // High priority to retire this branch.
        Kokkos::respawn(this, fib_all, Kokkos::TaskPriority::High);
      } else {
        Kokkos::abort("Failed nested task spawn (allocation)");
      }
    }
  }
};

int main(int argc, char* argv[]) {
  static const char help[]         = "--help";
  static const char alloc_size[]   = "--alloc_size=";
  static const char super_size[]   = "--super_size=";
  static const char repeat_outer[] = "--repeat_outer=";
  static const char input_value[]  = "--input=";

  long total_alloc_size   = 1000000;
  int min_superblock_size = 10000;
  int test_repeat_outer   = 1;
  int fib_input           = 4;

  int ask_help = 0;

  for (int i = 1; i < argc; i++) {
    const char* const a = argv[i];

    if (!strncmp(a, help, strlen(help))) ask_help = 1;

    if (!strncmp(a, alloc_size, strlen(alloc_size)))
      total_alloc_size = atol(a + strlen(alloc_size));

    if (!strncmp(a, super_size, strlen(super_size)))
      min_superblock_size = std::stoi(a + strlen(super_size));

    if (!strncmp(a, repeat_outer, strlen(repeat_outer)))
      test_repeat_outer = std::stoi(a + strlen(repeat_outer));

    if (!strncmp(a, input_value, strlen(input_value)))
      fib_input = std::stoi(a + strlen(input_value));
  }

  const long fib_output   = eval_fib(fib_input);
  const long number_alloc = fib_alloc_count(fib_input);

  const unsigned min_block_size = 32;
  const unsigned max_block_size = 128;

  long task_count_max   = 0;
  long task_count_accum = 0;
  long test_result      = 0;

  if (ask_help) {
    std::cout << "command line options:"
              << " " << help << " " << alloc_size << "##"
              << " " << super_size << "##"
              << " " << input_value << "##"
              << " " << repeat_outer << "##" << std::endl;
    return -1;
  }

  using Scheduler = Kokkos::TaskSchedulerMultiple<ExecSpace>;

  using Functor = TestFib<Scheduler>;

  Kokkos::initialize(argc, argv);

  {
    Scheduler sched(Functor::MemorySpace(), total_alloc_size, min_block_size,
                    max_block_size, min_superblock_size);

    Functor::FutureType f =
        Kokkos::host_spawn(Kokkos::TaskSingle(sched), Functor(fib_input));

    Kokkos::wait(sched);

    test_result = f.get();

    // task_count_max   = sched.allocated_task_count_max();
    // task_count_accum = sched.allocated_task_count_accum();

    // if ( number_alloc != task_count_accum ) {
    //  std::cout << " number_alloc( " << number_alloc << " )"
    //            << " != task_count_accum( " << task_count_accum << " )"
    //            << std::endl ;
    //}

    if (fib_output != test_result) {
      std::cout << " answer( " << fib_output << " )"
                << " != result( " << test_result << " )" << std::endl;
    }

    if (fib_output != test_result) {  // || number_alloc != task_count_accum ) {
      printf("  TEST FAILED\n");
      return -1;
    }

    double min_time = std::numeric_limits<double>::max();
    double time_sum = 0;

    for (int i = 0; i < test_repeat_outer; ++i) {
      Kokkos::Timer timer;

      Functor::FutureType ftmp =
          Kokkos::host_spawn(Kokkos::TaskSingle(sched), Functor(fib_input));

      Kokkos::wait(sched);
      auto this_time = timer.seconds();
      min_time       = std::min(min_time, this_time);
      time_sum += this_time;
    }

    auto avg_time = time_sum / test_repeat_outer;

    printf(
        "\"taskdag: alloc super repeat input output task-accum task-max\" %ld "
        "%d %d %d %ld %ld %ld\n",
        total_alloc_size, min_superblock_size, test_repeat_outer, fib_input,
        fib_output, task_count_accum, task_count_max);

    printf("\"taskdag: time (min, avg)\" %g %g\n", min_time, avg_time);
    printf("\"taskdag: tasks per second (max, avg)\" %g %g\n",
           number_alloc / min_time, number_alloc / avg_time);
  }  // end scope to destroy scheduler prior to finalize

  Kokkos::finalize();

  return 0;
}

#endif
