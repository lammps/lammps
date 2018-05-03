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

// Must be included first on Intel-Phi systems due to
// redefinition of SEEK_SET in <mpi.h>.

#include <ParallelComm.hpp>

#include <iostream>
#include <stdexcept>
#include <limits>
#include <utility>

//----------------------------------------------------------------------------

#include <Kokkos_Core.hpp>

#include <BoxMeshFixture.hpp>
#include <TestBoxMeshFixture.hpp>
#include <Implicit.hpp>
#include <Nonlinear.hpp>
#include <Explicit.hpp>
#include <SparseLinearSystem.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void test_host_fixture( comm::Machine machine ,
                        size_t gang_count ,
                        size_t gang_worker_count ,
                        size_t nx , size_t ny , size_t nz )
{
  Kokkos::HostSpace::execution_space::initialize( gang_count * gang_worker_count );
  test_box_fixture<Kokkos::HostSpace::execution_space>( machine , gang_count , nx , ny , nz );
  Kokkos::HostSpace::execution_space::finalize();
}

//----------------------------------------------------------------------------

void test_host_implicit( comm::Machine machine ,
                         size_t gang_count ,
                         size_t gang_worker_count ,
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t count_run )
{
  Kokkos::HostSpace::execution_space::initialize( gang_count * gang_worker_count );
  HybridFEM::Implicit::driver<double,Kokkos::HostSpace::execution_space>( "Threads" , machine , gang_count , elem_count_begin , elem_count_end , count_run );
  Kokkos::HostSpace::execution_space::finalize();
}

//----------------------------------------------------------------------------

void test_host_explicit( comm::Machine machine ,
                         size_t gang_count ,
                         size_t gang_worker_count ,
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t count_run )
{
  Kokkos::HostSpace::execution_space::initialize( gang_count * gang_worker_count );
  Explicit::driver<double,Kokkos::HostSpace::execution_space>( "Threads" , machine , gang_count , elem_count_begin , elem_count_end , count_run );
  Kokkos::HostSpace::execution_space::finalize();
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void test_host_nonlinear( comm::Machine machine ,
                          size_t gang_count ,
                          size_t gang_worker_count ,
                          size_t elem_count_begin ,
                          size_t elem_count_end ,
                          size_t count_run )
{
  Kokkos::HostSpace::execution_space::initialize( gang_count * gang_worker_count );
  typedef FixtureElementHex8 hex8 ;
  typedef Kokkos::HostSpace::execution_space             device ;
  HybridFEM::Nonlinear::driver<double,device,hex8>( "Threads" , machine , gang_count , elem_count_begin , elem_count_end , count_run );
  Kokkos::HostSpace::execution_space::finalize();
}

void test_host_nonlinear_quadratic( comm::Machine machine ,
                                    size_t gang_count ,
                                    size_t gang_worker_count ,
                                    size_t elem_count_begin ,
                                    size_t elem_count_end ,
                                    size_t count_run )
{
  Kokkos::HostSpace::execution_space::initialize( gang_count * gang_worker_count );
  typedef FixtureElementHex27 hex27 ;
  typedef Kokkos::HostSpace::execution_space              device ;
  HybridFEM::Nonlinear::driver<double,device,hex27>( "Threads" , machine , gang_count , elem_count_begin , elem_count_end , count_run );
  Kokkos::HostSpace::execution_space::finalize();
}

//----------------------------------------------------------------------------


