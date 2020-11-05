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
#include <Kokkos_Random.hpp>
#include <Kokkos_DualView.hpp>
#include <impl/Kokkos_Timer.hpp>
#include <cstdlib>

using DefaultHostType = Kokkos::HostSpace::execution_space;

// Kokkos provides two different random number generators with a 64 bit and a
// 1024 bit state. These generators are based on Vigna, Sebastiano (2014). "An
// experimental exploration of Marsaglia's xorshift generators, scrambled" See:
// http://arxiv.org/abs/1402.6246 The generators can be used fully independently
// on each thread and have been tested to produce good statistics for both inter
// and intra thread numbers. Note that within a kernel NO random number
// operations are (team) collective operations. Everything can be called within
// branches. This is a difference to the curand library where certain operations
// are required to be called by all threads in a block.
//
// In Kokkos you are required to create a pool of generator states, so that
// threads can grep their own. On CPU architectures the pool size is equal to
// the thread number, on CUDA about 128k states are generated (enough to give
// every potentially simultaneously running thread its own state). With a kernel
// a thread is required to acquire a state from the pool and later return it. On
// CPUs the Random number generator is deterministic if using the same number of
// threads. On GPUs (i.e. using the CUDA backend it is not deterministic because
// threads acquire states via atomics.

// A Functor for generating uint64_t random numbers templated on the
// GeneratorPool type
template <class GeneratorPool>
struct generate_random {
  // Output View for the random numbers
  Kokkos::View<uint64_t*> vals;

  // The GeneratorPool
  GeneratorPool rand_pool;

  int samples;

  // Initialize all members
  generate_random(Kokkos::View<uint64_t*> vals_, GeneratorPool rand_pool_,
                  int samples_)
      : vals(vals_), rand_pool(rand_pool_), samples(samples_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const {
    // Get a random number state from the pool for the active thread
    typename GeneratorPool::generator_type rand_gen = rand_pool.get_state();

    // Draw samples numbers from the pool as urand64 between 0 and
    // rand_pool.MAX_URAND64 Note there are function calls to get other type of
    // scalars, and also to specify Ranges or get a normal distributed float.
    for (int k = 0; k < samples; k++)
      vals(i * samples + k) = rand_gen.urand64();

    // Give the state back, which will allow another thread to acquire it
    rand_pool.free_state(rand_gen);
  }
};

int main(int argc, char* args[]) {
  if (argc != 3) {
    printf("Please pass two integers on the command line\n");
  } else {
    // Initialize Kokkos
    Kokkos::initialize(argc, args);
    int size    = std::stoi(args[1]);
    int samples = std::stoi(args[2]);

    // Create two random number generator pools one for 64bit states and one for
    // 1024 bit states Both take an 64 bit unsigned integer seed to initialize a
    // Random_XorShift64 generator which is used to fill the generators of the
    // pool.
    Kokkos::Random_XorShift64_Pool<> rand_pool64(5374857);
    Kokkos::Random_XorShift1024_Pool<> rand_pool1024(5374857);
    Kokkos::DualView<uint64_t*> vals("Vals", size * samples);

    // Run some performance comparisons
    Kokkos::Timer timer;
    Kokkos::parallel_for(size,
                         generate_random<Kokkos::Random_XorShift64_Pool<> >(
                             vals.d_view, rand_pool64, samples));
    Kokkos::fence();

    timer.reset();
    Kokkos::parallel_for(size,
                         generate_random<Kokkos::Random_XorShift64_Pool<> >(
                             vals.d_view, rand_pool64, samples));
    Kokkos::fence();
    double time_64 = timer.seconds();

    Kokkos::parallel_for(size,
                         generate_random<Kokkos::Random_XorShift1024_Pool<> >(
                             vals.d_view, rand_pool1024, samples));
    Kokkos::fence();

    timer.reset();
    Kokkos::parallel_for(size,
                         generate_random<Kokkos::Random_XorShift1024_Pool<> >(
                             vals.d_view, rand_pool1024, samples));
    Kokkos::fence();
    double time_1024 = timer.seconds();

    printf("#Time XorShift64*:   %e %e\n", time_64,
           1.0e-9 * samples * size / time_64);
    printf("#Time XorShift1024*: %e %e\n", time_1024,
           1.0e-9 * samples * size / time_1024);

    Kokkos::deep_copy(vals.h_view, vals.d_view);

    Kokkos::finalize();
  }
  return 0;
}
