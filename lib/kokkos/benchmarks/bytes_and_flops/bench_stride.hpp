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


#define UNROLL 1
#include<bench_unroll_stride.hpp>
#undef UNROLL
#define UNROLL 2
#include<bench_unroll_stride.hpp>
#undef UNROLL
#define UNROLL 3
#include<bench_unroll_stride.hpp>
#undef UNROLL
#define UNROLL 4
#include<bench_unroll_stride.hpp>
#undef UNROLL
#define UNROLL 5
#include<bench_unroll_stride.hpp>
#undef UNROLL
#define UNROLL 6
#include<bench_unroll_stride.hpp>
#undef UNROLL
#define UNROLL 7
#include<bench_unroll_stride.hpp>
#undef UNROLL
#define UNROLL 8
#include<bench_unroll_stride.hpp>
#undef UNROLL

template<class Scalar>
struct RunStride<Scalar,STRIDE> {
static void run_1(int N, int K, int R, int F, int T, int S) {
  Run<Scalar,1,STRIDE>::run(N,K,R,F,T,S);
}
static void run_2(int N, int K, int R, int F, int T, int S) {
  Run<Scalar,2,STRIDE>::run(N,K,R,F,T,S);
}
static void run_3(int N, int K, int R, int F, int T, int S) {
  Run<Scalar,3,STRIDE>::run(N,K,R,F,T,S);
}
static void run_4(int N, int K, int R, int F, int T, int S) {
  Run<Scalar,4,STRIDE>::run(N,K,R,F,T,S);
}
static void run_5(int N, int K, int R, int F, int T, int S) {
  Run<Scalar,5,STRIDE>::run(N,K,R,F,T,S);
}
static void run_6(int N, int K, int R, int F, int T, int S) {
  Run<Scalar,6,STRIDE>::run(N,K,R,F,T,S);
}
static void run_7(int N, int K, int R, int F, int T, int S) {
  Run<Scalar,7,STRIDE>::run(N,K,R,F,T,S);
}
static void run_8(int N, int K, int R, int F, int T, int S) {
  Run<Scalar,8,STRIDE>::run(N,K,R,F,T,S);
}

static void run(int N, int K, int R, int U, int F, int T, int S) {
  if(U==1) {
    run_1(N,K,R,F,T,S);
  }
  if(U==2) {
    run_2(N,K,R,F,T,S);
  }
  if(U==3) {
    run_3(N,K,R,F,T,S);
  }
  if(U==4) {
    run_4(N,K,R,F,T,S);
  }
  if(U==5) {
    run_5(N,K,R,F,T,S);
  }
  if(U==6) {
    run_6(N,K,R,F,T,S);
  }
  if(U==7) {
    run_7(N,K,R,F,T,S);
  }
  if(U==8) {
    run_8(N,K,R,F,T,S);
  } 
}
};

