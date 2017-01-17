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

#include<Kokkos_Core.hpp>
#include<Kokkos_Random.hpp>

template<class Scalar>
struct RunGather<Scalar,UNROLL> {
static void run(int N, int K, int D, int R, int F) {
  Kokkos::View<int**> connectivity("Connectivity",N,K);
  Kokkos::View<Scalar*> A_in("Input",N);
  Kokkos::View<Scalar*> B_in("Input",N);
  Kokkos::View<Scalar*> C("Output",N);

  Kokkos::Random_XorShift64_Pool<> rand_pool(12313);

  Kokkos::deep_copy(A_in,1.5);
  Kokkos::deep_copy(B_in,2.0);

  Kokkos::View<const Scalar*, Kokkos::MemoryTraits<Kokkos::RandomAccess> > A(A_in);
  Kokkos::View<const Scalar*, Kokkos::MemoryTraits<Kokkos::RandomAccess> > B(B_in);

  Kokkos::parallel_for("InitKernel",N,
      KOKKOS_LAMBDA (const int& i) {
    auto rand_gen = rand_pool.get_state();
    for( int jj=0; jj<K; jj++) {
      connectivity(i,jj) = (rand_gen.rand(D) + i - D/2 + N)%N;
    }
    rand_pool.free_state(rand_gen);
  });
  Kokkos::fence();


  Kokkos::Timer timer;
  for(int r = 0; r<R; r++) {
  Kokkos::parallel_for("BenchmarkKernel",N,
      KOKKOS_LAMBDA (const int& i) {
      Scalar c = Scalar(0.0);
      for( int jj=0; jj<K; jj++) {
        const int j = connectivity(i,jj);
        Scalar a1 = A(j);
        const Scalar b = B(j);
#if(UNROLL>1)
        Scalar a2 = a1*Scalar(1.3);
#endif
#if(UNROLL>2)
        Scalar a3 = a2*Scalar(1.1);
#endif
#if(UNROLL>3)
        Scalar a4 = a3*Scalar(1.1);
#endif
#if(UNROLL>4)
        Scalar a5 = a4*Scalar(1.3);
#endif
#if(UNROLL>5)
        Scalar a6 = a5*Scalar(1.1);
#endif
#if(UNROLL>6)
        Scalar a7 = a6*Scalar(1.1);
#endif
#if(UNROLL>7)
        Scalar a8 = a7*Scalar(1.1);
#endif


        for(int f = 0; f<F; f++) {
          a1 += b*a1;
#if(UNROLL>1)
          a2 += b*a2;
#endif
#if(UNROLL>2)
          a3 += b*a3;
#endif
#if(UNROLL>3)
          a4 += b*a4;
#endif
#if(UNROLL>4)
          a5 += b*a5;
#endif
#if(UNROLL>5)
          a6 += b*a6;
#endif
#if(UNROLL>6)
          a7 += b*a7;
#endif
#if(UNROLL>7)
          a8 += b*a8;
#endif


        }
#if(UNROLL==1)
        c += a1;
#endif
#if(UNROLL==2)
        c += a1+a2;
#endif
#if(UNROLL==3)
        c += a1+a2+a3;
#endif
#if(UNROLL==4)
        c += a1+a2+a3+a4;
#endif
#if(UNROLL==5)
        c += a1+a2+a3+a4+a5;
#endif
#if(UNROLL==6)
        c += a1+a2+a3+a4+a5+a6;
#endif
#if(UNROLL==7)
        c += a1+a2+a3+a4+a5+a6+a7;
#endif
#if(UNROLL==8)
        c += a1+a2+a3+a4+a5+a6+a7+a8;
#endif

      }
      C(i) = c ;
  });
  Kokkos::fence();
  }
  double seconds = timer.seconds();

  double bytes = 1.0*N*K*R*(2*sizeof(Scalar)+sizeof(int)) + 1.0*N*R*sizeof(Scalar);
  double flops = 1.0*N*K*R*(F*2*UNROLL + 2*(UNROLL-1));
  double gather_ops = 1.0*N*K*R*2;
  printf("SNKDRUF: %i %i %i %i %i %i %i Time: %lfs Bandwidth: %lfGiB/s GFlop/s: %lf GGather/s: %lf\n",sizeof(Scalar)/4,N,K,D,R,UNROLL,F,seconds,1.0*bytes/seconds/1024/1024/1024,1.e-9*flops/seconds,1.e-9*gather_ops/seconds);
}
};
