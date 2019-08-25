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
#include <cstdio>

// The type of a two-dimensional N x 3 array of double.
// It lives in Kokkos' default memory space.
typedef Kokkos::View<double*[3]> view_type;

// The "HostMirror" type corresponding to view_type above is also a
// two-dimensional N x 3 array of double.  However, it lives in the
// host memory space corresponding to view_type's memory space.  For
// example, if view_type lives in CUDA device memory, host_view_type
// lives in host (CPU) memory.  Furthermore, declaring host_view_type
// as the host mirror of view_type means that host_view_type has the
// same layout as view_type.  This makes it easier to copy between the
// two Views.
// Advanced issues: If a memory space is accessible from the host without
// performance penalties then it is its own host_mirror_space. This is
// the case for HostSpace, CudaUVMSpace and CudaHostPinnedSpace.

typedef view_type::HostMirror host_view_type;

struct ReduceFunctor {
  view_type a;
  ReduceFunctor (view_type a_) : a (a_) {}
  typedef int value_type; //Specify type for reduction value, lsum

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int &lsum) const {
    lsum += a(i,0)-a(i,1)+a(i,2);
  }
};

int main() {
  Kokkos::initialize();

  {
    view_type a ("A", 10);
    // If view_type and host_mirror_type live in the same memory space,
    // a "mirror view" is just an alias, and deep_copy does nothing.
    // Otherwise, a mirror view of a device View lives in host memory,
    // and deep_copy does a deep copy.
    host_view_type h_a = Kokkos::create_mirror_view (a);

    // The View h_a lives in host (CPU) memory, so it's legal to fill
    // the view sequentially using ordinary code, like this.
    for (int i = 0; i < 10; i++) {
      for (int j = 0; j < 3; j++) {
        h_a(i,j) = i*10 + j;
      }
    }
    Kokkos::deep_copy (a, h_a); // Copy from host to device.

    int sum = 0;
    Kokkos::parallel_reduce (10, ReduceFunctor (a), sum);
    printf ("Result is %i\n",sum);
  }

  Kokkos::finalize ();
}

