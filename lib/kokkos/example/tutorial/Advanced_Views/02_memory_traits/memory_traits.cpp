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
#include <impl/Kokkos_Timer.hpp>
#include <cstdio>
#include <cstdlib>

typedef Kokkos::View<double*> view_type;
// Kokkos::Views have an MemoryTraits template parameter which
// allows users to specify usage scenarios of a View.
// Some of those act simply as hints, which can be used to insert
// optimal load and store paths, others change the symantics of the
// access. The trait Kokkos::Atomic is one of the latter. A view with
// that MemoryTrait will perform any access atomicly (read, write, update).
//
// In this example we use a view with a usage hint for RandomAccess.
// Kokkos::RandomAccess means that we expect to use this view
// with indirect indexing.
//
// In CUDA, RandomAccess allows accesses through the texture
// cache.  This only works if the View is read-only, which we enforce
// through the first template parameter.
//
// Note that we are still talking about views of the data, its not a new allocation.
// For example you can have an atomic view of a default view. While you even
// could use both in the same kernel, this could lead to undefined behaviour because
// one of your access paths is not atomic. Think of it in the same way as you think of
// pointers to const data and pointers to non-const data (i.e. const double* and double*).
// While these pointers can point to the same data you should not use them together if that
// brakes the const guarantee of the first pointer.
typedef Kokkos::View<const double*, Kokkos::MemoryTraits<Kokkos::RandomAccess> > view_type_rnd;
typedef Kokkos::View<int**> idx_type;
typedef idx_type::HostMirror idx_type_host;

// We template this functor on the ViewTypes to show the effect of the RandomAccess trait.
template<class DestType, class SrcType>
struct localsum {
  idx_type::const_type idx;
  DestType dest;
  SrcType src;
  localsum (idx_type idx_, DestType dest_, SrcType src_) :
    idx (idx_), dest (dest_), src (src_)
  {}

  // Calculate a local sum of values
  KOKKOS_INLINE_FUNCTION
  void operator() (const int i) const {
    double tmp = 0.0;
    for (int j = 0; j < (int) idx.dimension_1 (); ++j) {
      // This is an indirect access on src
      const double val = src(idx(i,j));
      tmp += val*val + 0.5*(idx.dimension_0()*val -idx.dimension_1()*val);
    }
    dest(i) = tmp;
  }
};

int main(int narg, char* arg[]) {
  Kokkos::initialize (narg, arg);

  int size = 1000000;

  idx_type idx("Idx",size,64);
  idx_type_host h_idx = Kokkos::create_mirror_view (idx);

  view_type dest ("Dest", size);
  view_type src ("Src", size);

  srand(134231);

  for (int i = 0; i < size; i++) {
    for (view_type::size_type j = 0; j < h_idx.dimension_1 (); ++j) {
      h_idx(i,j) = (size + i + (rand () % 500 - 250)) % size;
    }
  }

  // Deep copy the initial data to the device
  Kokkos::deep_copy(idx,h_idx);
  // Run the first kernel to warmup caches
  Kokkos::parallel_for(size,localsum<view_type,view_type_rnd>(idx,dest,src));
  Kokkos::fence();

  // Run the localsum functor using the RandomAccess trait. On CPUs there should
  // not be any different in performance to not using the RandomAccess trait.
  // On GPUs where can be a dramatic difference
  Kokkos::Timer time1;
  Kokkos::parallel_for(size,localsum<view_type,view_type_rnd>(idx,dest,src));
  Kokkos::fence();
  double sec1 = time1.seconds();

  Kokkos::Timer time2;
  Kokkos::parallel_for(size,localsum<view_type,view_type>(idx,dest,src));
  Kokkos::fence();
  double sec2 = time2.seconds();

  printf("Time with Trait RandomAccess: %f with Plain: %f \n",sec1,sec2);

  Kokkos::finalize();
}

