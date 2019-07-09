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
#include <typeinfo>
#include <cmath>
#include <impl/Kokkos_Timer.hpp>

struct FillDevice {
  double value;
  Kokkos::View<double*,Kokkos::LayoutLeft,Kokkos::CudaSpace> a;
  FillDevice(const double& val, const Kokkos::View<double*,Kokkos::LayoutLeft,Kokkos::CudaSpace>& d_a):
     value(val),a(d_a){}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    a(i) = value;
  }
};

struct ComputeADevice {
  int iter;
  Kokkos::View<double*,Kokkos::LayoutLeft,Kokkos::CudaSpace> a;
  Kokkos::View<double*,Kokkos::LayoutLeft,Kokkos::CudaSpace> b;
  ComputeADevice(const int& iter_,
                 const Kokkos::View<double*,Kokkos::LayoutLeft,Kokkos::CudaSpace>& d_a,
                 const Kokkos::View<double*,Kokkos::LayoutLeft,Kokkos::CudaSpace>& d_b):
     iter(iter_),a(d_a),b(d_b){}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    for(int j=1;j<iter;j++) {
      a(i) += std::pow(b(i),1.0+1.0/iter);
    }
  }
};

struct ComputeAHost {
  Kokkos::View<double*,Kokkos::LayoutLeft,Kokkos::CudaHostPinnedSpace> a;
  Kokkos::View<double*,Kokkos::LayoutLeft,Kokkos::CudaHostPinnedSpace> b;
  ComputeAHost(  const Kokkos::View<double*,Kokkos::LayoutLeft,Kokkos::CudaHostPinnedSpace>& d_a,
                 const Kokkos::View<double*,Kokkos::LayoutLeft,Kokkos::CudaHostPinnedSpace>& d_b):
     a(d_a),b(d_b){}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    a(i) += b(i);
  }
};

struct MergeDevice {
  Kokkos::View<double*,Kokkos::LayoutLeft,Kokkos::CudaSpace> a;
  Kokkos::View<double*,Kokkos::LayoutLeft,Kokkos::CudaSpace> b;
  MergeDevice(
                 const Kokkos::View<double*,Kokkos::LayoutLeft,Kokkos::CudaSpace>& d_a,
                 const Kokkos::View<double*,Kokkos::LayoutLeft,Kokkos::CudaSpace>& d_b):
     a(d_a),b(d_b){}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    a(i) += b(i);
  }
};

int main(int argc, char * argv[]) {
  int size = 100000000;
  Kokkos::initialize();
  int synch = atoi(argv[1]);
  Kokkos::View<double*,Kokkos::LayoutLeft,Kokkos::CudaSpace> d_a("Device A",size);
  Kokkos::View<double*,Kokkos::LayoutLeft,Kokkos::CudaSpace> d_b("Device B",size);
  Kokkos::View<double*,Kokkos::LayoutLeft,Kokkos::CudaSpace> d_tmp("Device tmp",size);
  Kokkos::View<double*,Kokkos::LayoutLeft,Kokkos::CudaHostPinnedSpace> h_a("Host A",size);
  Kokkos::View<double*,Kokkos::LayoutLeft,Kokkos::CudaHostPinnedSpace> h_b("Host B",size);

  Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::Cuda>(0,size),FillDevice(0.0,d_a));
  Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::Cuda>(0,size),FillDevice(1.3513,d_b));
  Kokkos::fence();
  Kokkos::Timer timer;
  Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::Cuda>(0,size),ComputeADevice(20,d_a,d_b));

  if(synch==1)
    Kokkos::deep_copy(Kokkos::OpenMP(),h_b,d_b);
  if(synch==2)
    Kokkos::deep_copy(h_b,d_b);


  Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::OpenMP>(0,size),[=] (const int& i) {
    h_a(i) = 0.0;
  });
  Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::OpenMP>(0,size),ComputeAHost(h_a,h_b));
  Kokkos::OpenMP().fence();
  if(synch==1)
    Kokkos::deep_copy(Kokkos::OpenMP(), d_tmp,h_a);
  if(synch==2)
    Kokkos::deep_copy(d_tmp,h_a);
  Kokkos::fence();

  std::cout << "Time " << timer.seconds() << std::endl;
  Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::Cuda>(0,size),MergeDevice(d_a,d_tmp));

  Kokkos::deep_copy(h_a,d_a);
  std::cout << "h_a(0): " << h_a(0) << " ( Correct: 27.4154 )" << std::endl;
  Kokkos::finalize();
}



