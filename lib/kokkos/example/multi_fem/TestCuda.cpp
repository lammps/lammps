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

#include <TestBoxMeshFixture.hpp>
#include <Implicit.hpp>
#include <Nonlinear.hpp>
#include <Explicit.hpp>

#include <SparseLinearSystem.hpp>

#if defined( KOKKOS_ENABLE_CUDA )

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

CudaSparseSingleton & CudaSparseSingleton::singleton()
{ static CudaSparseSingleton s ; return s ; }

}
}

//----------------------------------------------------------------------------

void test_cuda_query( comm::Machine machine )
{
  const size_t comm_rank = comm::rank( machine );
  std::cout << "P" << comm_rank
            << ": Cuda device_count = "
            << Kokkos::Cuda::detect_device_count()
            << std::endl ;
}

//----------------------------------------------------------------------------

void test_cuda_fixture( comm::Machine machine ,
                        size_t nx , size_t ny , size_t nz )
{
  const size_t comm_rank = comm::rank( machine );
  const size_t comm_size = comm::size( machine );
  const size_t dev_count = Kokkos::Cuda::detect_device_count();
  const size_t dev_rank =
    dev_count && dev_count <= comm_size ? comm_rank % dev_count : 0 ;
  const size_t gang_count = 0 ;

  Kokkos::HostSpace::execution_space::initialize();
  Kokkos::Cuda::SelectDevice select_device( dev_rank );
  Kokkos::Cuda::initialize( select_device );
  test_box_fixture<Kokkos::Cuda>( machine , gang_count , nx , ny , nz );
  Kokkos::Cuda::finalize();
  Kokkos::HostSpace::execution_space::finalize();
}

//----------------------------------------------------------------------------

void test_cuda_implicit( comm::Machine machine , 
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t count_run )
{
  const size_t comm_rank = comm::rank( machine );
  const size_t comm_size = comm::size( machine );
  const size_t dev_count = Kokkos::Cuda::detect_device_count();
  const size_t dev_rank =
    dev_count && dev_count <= comm_size ? comm_rank % dev_count : 0 ;
  const size_t gang_count = 0 ;

  Kokkos::HostSpace::execution_space::initialize();
  Kokkos::Cuda::SelectDevice select_device( dev_rank );
  Kokkos::Cuda::initialize( select_device );
  HybridFEM::Implicit::driver<double,Kokkos::Cuda>( "Cuda" , machine , gang_count , elem_count_begin , elem_count_end , count_run );
  Kokkos::Cuda::finalize();
  Kokkos::HostSpace::execution_space::finalize();
}

//----------------------------------------------------------------------------

void test_cuda_explicit( comm::Machine machine , 
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t count_run )
{
  const size_t comm_rank = comm::rank( machine );
  const size_t comm_size = comm::size( machine );
  const size_t dev_count = Kokkos::Cuda::detect_device_count();
  const size_t dev_rank =
    dev_count && dev_count <= comm_size ? comm_rank % dev_count : 0 ;
  const size_t gang_count = 0 ;

  Kokkos::HostSpace::execution_space::initialize();
  Kokkos::Cuda::SelectDevice select_device( dev_rank );
  Kokkos::Cuda::initialize( select_device );
  Explicit::driver<double,Kokkos::Cuda>( "Cuda" , machine , gang_count , elem_count_begin , elem_count_end , count_run );
  Kokkos::Cuda::finalize();
  Kokkos::HostSpace::execution_space::finalize();
}

//----------------------------------------------------------------------------

void test_cuda_nonlinear( comm::Machine machine , 
                          size_t elem_count_begin ,
                          size_t elem_count_end ,
                          size_t count_run )
{
  const size_t comm_rank = comm::rank( machine );
  const size_t comm_size = comm::size( machine );
  const size_t dev_count = Kokkos::Cuda::detect_device_count();
  const size_t dev_rank =
    dev_count && dev_count <= comm_size ? comm_rank % dev_count : 0 ;
  const size_t gang_count = 0 ;

  Kokkos::HostSpace::execution_space::initialize();
  Kokkos::Cuda::SelectDevice select_device( dev_rank );
  Kokkos::Cuda::initialize( select_device );

  typedef Kokkos::Cuda device ;
  typedef FixtureElementHex8 hex8 ;
  HybridFEM::Nonlinear::driver<double,device,hex8>( "Cuda" , machine , gang_count , elem_count_begin , elem_count_end , count_run );
  Kokkos::Cuda::finalize();
  Kokkos::HostSpace::execution_space::finalize();
}

void test_cuda_nonlinear_quadratic( comm::Machine machine , 
                                    size_t elem_count_begin ,
                                    size_t elem_count_end ,
                                    size_t count_run )
{
  const size_t comm_rank = comm::rank( machine );
  const size_t comm_size = comm::size( machine );
  const size_t dev_count = Kokkos::Cuda::detect_device_count();
  const size_t dev_rank =
    dev_count && dev_count <= comm_size ? comm_rank % dev_count : 0 ;
  const size_t gang_count = 0 ;

  Kokkos::HostSpace::execution_space::initialize();
  Kokkos::Cuda::SelectDevice select_device( dev_rank );
  Kokkos::Cuda::initialize( select_device );

  typedef Kokkos::Cuda device ;
  typedef FixtureElementHex27 hex27 ;
  HybridFEM::Nonlinear::driver<double,device,hex27>( "Cuda" , machine , gang_count , elem_count_begin , elem_count_end , count_run );
  Kokkos::Cuda::finalize();
  Kokkos::HostSpace::execution_space::finalize();
}

//----------------------------------------------------------------------------

#endif  /* #if defined( KOKKOS_ENABLE_CUDA ) */

