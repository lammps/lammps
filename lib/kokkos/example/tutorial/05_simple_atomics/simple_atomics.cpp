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
#include <cstdlib>
#include <cmath>

// Type of a one-dimensional length-N array of int.
typedef Kokkos::View<int*> view_type;
typedef view_type::HostMirror host_view_type;
// This is a "zero-dimensional" View, that is, a View of a single
// value (an int, in this case).  Access the value using operator()
// with no arguments: e.g., 'count()'.
//
// Zero-dimensional Views are useful for reduction results that stay
// resident in device memory, as well as for irregularly updated
// shared state.  We use it for the latter in this example.
typedef Kokkos::View<int> count_type;
typedef count_type::HostMirror host_count_type;

// Functor for finding a list of primes in a given set of numbers.  If
// run in parallel, the order of results is nondeterministic, because
// hardware atomic updates do not guarantee an order of execution.
struct findprimes {
  view_type data;
  view_type result;
  count_type count;

  findprimes(view_type data_, view_type result_, count_type count_)
      : data(data_), result(result_), count(count_) {}

  // Test if data(i) is prime.  If it is, increment the count of
  // primes (stored in the zero-dimensional View 'count') and add the
  // value to the current list of primes 'result'.
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    const int number = data(i);  // the current number

    // Test all numbers from 3 to ceiling(sqrt(data(i))), to see if
    // they are factors of data(i).  It's not the most efficient prime
    // test, but it works.
    const int upper_bound = std::sqrt(1.0 * number) + 1;
    bool is_prime         = !(number % 2 == 0);
    int k                 = 3;
    while (k < upper_bound && is_prime) {
      is_prime = !(number % k == 0);
      k += 2;  // don't have to test even numbers
    }

    if (is_prime) {
      // Use an atomic update both to update the current count of
      // primes, and to find a place in the current list of primes for
      // the new result.
      //
      // atomic_fetch_add results the _current_ count, but increments
      // it (by 1 in this case).  The current count of primes indexes
      // into the first unoccupied position of the 'result' array.
      const int idx = Kokkos::atomic_fetch_add(&count(), 1);
      result(idx)   = number;
    }
  }
};

int main() {
  Kokkos::initialize();

  {
    srand(61391);  // Set the random seed

    int nnumbers = 100000;
    view_type data("RND", nnumbers);
    view_type result("Prime", nnumbers);
    count_type count("Count");

    host_view_type h_data   = Kokkos::create_mirror_view(data);
    host_view_type h_result = Kokkos::create_mirror_view(result);
    host_count_type h_count = Kokkos::create_mirror_view(count);

    typedef view_type::size_type size_type;
    // Fill the 'data' array on the host with random numbers.  We assume
    // that they come from some process which is only implemented on the
    // host, via some library.  (That's true in this case.)
    for (size_type i = 0; i < data.extent(0); ++i) {
      h_data(i) = rand() % nnumbers;
    }
    Kokkos::deep_copy(data, h_data);  // copy from host to device

    Kokkos::parallel_for(data.extent(0), findprimes(data, result, count));
    Kokkos::deep_copy(h_count, count);  // copy from device to host

    printf("Found %i prime numbers in %i random numbers\n", h_count(),
           nnumbers);
  }
  Kokkos::finalize();
}
