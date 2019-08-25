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
#include <cuda/TestCuda_Category.hpp>

namespace Test {

__global__ void offset_streams(int* p) {
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if(idx<100) {
    p[idx]+=idx;
  }
}

namespace {
  struct FunctorRange {
    Kokkos::View<int*,Kokkos::CudaSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> a;
    FunctorRange(Kokkos::View<int*,Kokkos::CudaSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> a_):a(a_) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const int i) const {
      a(i)+=1;
    }
  };
  struct FunctorRangeReduce {
    Kokkos::View<int*,Kokkos::CudaSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> a;
    FunctorRangeReduce(Kokkos::View<int*,Kokkos::CudaSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> a_):a(a_) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const int i, int& lsum) const {
      lsum += a(i);
    }
  };
  struct FunctorMDRange {
    Kokkos::View<int*,Kokkos::CudaSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> a;
    FunctorMDRange(Kokkos::View<int*,Kokkos::CudaSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> a_):a(a_) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const int i, const int j) const {
      a(i*10+j)+=1;
    }
  };
  struct FunctorMDRangeReduce {
    Kokkos::View<int*,Kokkos::CudaSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> a;
    FunctorMDRangeReduce(Kokkos::View<int*,Kokkos::CudaSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> a_):a(a_) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const int i, const int j, int& lsum) const {
      lsum += a(i*10+j);
    }
  };
  struct FunctorTeam {
    Kokkos::View<int*,Kokkos::CudaSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> a;
    FunctorTeam(Kokkos::View<int*,Kokkos::CudaSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> a_):a(a_) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const Kokkos::TeamPolicy<Kokkos::Cuda>::member_type& team) const {
      int i = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,10),[&](const int j){
        a(i*10+j)+=1;
      });
    }
  };

  struct FunctorTeamReduce {
    Kokkos::View<int*,Kokkos::CudaSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> a;
    FunctorTeamReduce(Kokkos::View<int*,Kokkos::CudaSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> a_):a(a_) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const Kokkos::TeamPolicy<Kokkos::Cuda>::member_type& team, int& lsum) const {
      int i = team.league_rank();
      int team_sum;
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team,10),[&](const int j, int& tsum){
        tsum += a(i*10+j);
      },team_sum);
      Kokkos::single(Kokkos::PerTeam(team),[&]() {
        lsum += team_sum;
      });
    }
  };
}

// Test Interoperability with Cuda Streams
TEST_F( cuda, raw_cuda_streams )
{
  cudaStream_t stream;
  cudaStreamCreate(&stream);
  Kokkos::InitArguments arguments{-1,-1,-1, false};
  Kokkos::initialize(arguments);
  int* p;
  cudaMalloc(&p,sizeof(int)*100);

  {
  Kokkos::Cuda cuda0(stream);
  Kokkos::View<int*,Kokkos::CudaSpace>
    v(p,100);
  Kokkos::deep_copy(cuda0,v,5);
  int sum;

  Kokkos::parallel_for("Test::cuda::raw_cuda_stream::Range",
      Kokkos::RangePolicy<Kokkos::Cuda>(cuda0,0,100),FunctorRange(v));
  Kokkos::parallel_reduce("Test::cuda::raw_cuda_stream::RangeReduce",
      Kokkos::RangePolicy<Kokkos::Cuda,Kokkos::LaunchBounds<128,2>>(cuda0,0,100),FunctorRangeReduce(v),sum);
  cuda0.fence();
  ASSERT_EQ(600,sum);

  Kokkos::parallel_for("Test::cuda::raw_cuda_stream::MDRange",
      Kokkos::MDRangePolicy<Kokkos::Cuda,Kokkos::Rank<2>>(cuda0,{0,0},{10,10}),FunctorMDRange(v));
  Kokkos::parallel_reduce("Test::cuda::raw_cuda_stream::MDRangeReduce",
      Kokkos::MDRangePolicy<Kokkos::Cuda,Kokkos::Rank<2>,Kokkos::LaunchBounds<128,2>>(cuda0,{0,0},{10,10}),FunctorMDRangeReduce(v),sum);
  cuda0.fence();
  ASSERT_EQ(700,sum);

  Kokkos::parallel_for("Test::cuda::raw_cuda_stream::Team",
      Kokkos::TeamPolicy<Kokkos::Cuda>(cuda0,10,10),FunctorTeam(v));
  Kokkos::parallel_reduce("Test::cuda::raw_cuda_stream::Team",
      Kokkos::TeamPolicy<Kokkos::Cuda,Kokkos::LaunchBounds<128,2>>(cuda0,10,10),FunctorTeamReduce(v),sum);
  cuda0.fence();
  ASSERT_EQ(800,sum);

  }
  Kokkos::finalize();
  offset_streams<<<100,64,0,stream>>>(p);
  CUDA_SAFE_CALL( cudaDeviceSynchronize());
  cudaStreamDestroy(stream);

  int* h_p = new int[100];
  cudaMemcpy( h_p , p , sizeof(int)*100 , cudaMemcpyDefault );
  CUDA_SAFE_CALL( cudaDeviceSynchronize());
  int64_t sum = 0;
  int64_t sum_expect = 0;
  for(int i=0; i<100; i++) {
    sum += h_p[i];
    sum_expect += 8+i;
  }

  ASSERT_EQ(sum,sum_expect);
}
}
