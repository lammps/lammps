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

template<class Scalar, int UNROLL>
struct RunGather {
  static void run(int N, int K, int D, int R, int F);
};

#define UNROLL 1
#include<gather_unroll.hpp>
#undef UNROLL
#define UNROLL 2
#include<gather_unroll.hpp>
#undef UNROLL
#define UNROLL 3
#include<gather_unroll.hpp>
#undef UNROLL
#define UNROLL 4
#include<gather_unroll.hpp>
#undef UNROLL
#define UNROLL 5
#include<gather_unroll.hpp>
#undef UNROLL
#define UNROLL 6
#include<gather_unroll.hpp>
#undef UNROLL
#define UNROLL 7
#include<gather_unroll.hpp>
#undef UNROLL
#define UNROLL 8
#include<gather_unroll.hpp>
#undef UNROLL

template<class Scalar>
void run_gather_test(int N, int K, int D, int R, int U, int F) {
 if(U == 1)
   RunGather<Scalar,1>::run(N,K,D,R,F);
 if(U == 2)
   RunGather<Scalar,2>::run(N,K,D,R,F);
 if(U == 3)
   RunGather<Scalar,3>::run(N,K,D,R,F);
 if(U == 4)
   RunGather<Scalar,4>::run(N,K,D,R,F);
 if(U == 5)
   RunGather<Scalar,5>::run(N,K,D,R,F);
 if(U == 6)
   RunGather<Scalar,6>::run(N,K,D,R,F);
 if(U == 7)
   RunGather<Scalar,7>::run(N,K,D,R,F);
 if(U == 8)
   RunGather<Scalar,8>::run(N,K,D,R,F);
}
