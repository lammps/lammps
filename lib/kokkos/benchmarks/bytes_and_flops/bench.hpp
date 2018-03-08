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

#include<Kokkos_Core.hpp>
#include<impl/Kokkos_Timer.hpp>

template<class Scalar, int Unroll,int Stride>
struct Run {
static void run(int N, int K, int R, int F, int T, int S);
};

template<class Scalar, int Stride>
struct RunStride {
static void run_1(int N, int K, int R, int F, int T, int S);
static void run_2(int N, int K, int R, int F, int T, int S);
static void run_3(int N, int K, int R, int F, int T, int S);
static void run_4(int N, int K, int R, int F, int T, int S);
static void run_5(int N, int K, int R, int F, int T, int S);
static void run_6(int N, int K, int R, int F, int T, int S);
static void run_7(int N, int K, int R, int F, int T, int S);
static void run_8(int N, int K, int R, int F, int T, int S);
static void run(int N, int K, int R, int U, int F, int T, int S);
};

#define STRIDE 1
#include<bench_stride.hpp>
#undef STRIDE
#define STRIDE 2
#include<bench_stride.hpp>
#undef STRIDE
#define STRIDE 4
#include<bench_stride.hpp>
#undef STRIDE
#define STRIDE 8
#include<bench_stride.hpp>
#undef STRIDE
#define STRIDE 16
#include<bench_stride.hpp>
#undef STRIDE
#define STRIDE 32
#include<bench_stride.hpp>
#undef STRIDE

template<class Scalar>
void run_stride_unroll(int N, int K, int R, int D, int U, int F, int T, int S) {
 if(D == 1)
   RunStride<Scalar,1>::run(N,K,R,U,F,T,S);
 if(D == 2)
   RunStride<Scalar,2>::run(N,K,R,U,F,T,S);
 if(D == 4)
   RunStride<Scalar,4>::run(N,K,R,U,F,T,S);
 if(D == 8)
   RunStride<Scalar,8>::run(N,K,R,U,F,T,S);
 if(D == 16)
   RunStride<Scalar,16>::run(N,K,R,U,F,T,S);
 if(D == 32)
   RunStride<Scalar,32>::run(N,K,R,U,F,T,S);
}

