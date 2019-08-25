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

#include <string>
#include <sstream>
#include <iostream>
#include <Kokkos_hwloc.hpp>

//----------------------------------------------------------------------------

void test_box_partition( bool print );

//----------------------------------------------------------------------------

void test_host_fixture( comm::Machine machine ,
                        size_t gang_count ,
                        size_t gang_worker_count ,
                        size_t nx , size_t ny , size_t nz );

void test_host_implicit( comm::Machine machine ,
                         size_t gang_count ,
                         size_t gang_worker_count ,
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t count_run );

void test_host_explicit( comm::Machine machine ,
                         size_t gang_count ,
                         size_t gang_worker_count ,
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t count_run );

void test_host_nonlinear( comm::Machine machine ,
                          size_t gang_count ,
                          size_t gang_worker_count ,
                          size_t elem_count_begin ,
                          size_t elem_count_end ,
                          size_t count_run );

void test_host_nonlinear_quadratic( comm::Machine machine ,
                                    size_t gang_count ,
                                    size_t gang_worker_count ,
                                    size_t elem_count_begin ,
                                    size_t elem_count_end ,
                                    size_t count_run );


//----------------------------------------------------------------------------

void test_cuda_query( comm::Machine );

void test_cuda_fixture( comm::Machine machine ,
                        size_t nx , size_t ny , size_t nz );

void test_cuda_implicit( comm::Machine machine ,
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t count_run );

void test_cuda_explicit( comm::Machine machine ,
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t count_run );

void test_cuda_nonlinear( comm:: Machine machine ,
                          size_t elem_count_begin ,
                          size_t elem_count_end ,
                          size_t count_run );

void test_cuda_nonlinear_quadratic( comm::Machine machine ,
                                    size_t elem_count_begin ,
                                    size_t elem_count_end ,
                                    size_t count_run );


//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace {

bool run_host( std::istream & input ,
               comm::Machine machine ,
               const size_t host_gang_count ,
               const size_t host_gang_worker_count )
{
  bool cmd_error = false ;

  std::string which ; input >> which ;

  if ( which == std::string("fixture") ) {

    size_t nx = 0 , ny = 0 , nz = 0 ;
    input >> nx >> ny >> nz ;
    test_host_fixture( machine , host_gang_count , host_gang_worker_count , nx , ny , nz );

  }
  else if ( which == std::string("explicit") ) {

    size_t mesh_node_begin = 100 ;
    size_t mesh_node_end   = 300 ;
    size_t run             =   1 ;
    input >> mesh_node_begin >> mesh_node_end >> run ;
    test_host_explicit( machine , host_gang_count , host_gang_worker_count , mesh_node_begin , mesh_node_end , run );

  }
  else if ( which == std::string("implicit") ) {

    size_t mesh_node_begin = 100 ;
    size_t mesh_node_end   = 300 ;
    size_t run             =   1 ;
    input >> mesh_node_begin >> mesh_node_end >> run ;
    test_host_implicit( machine , host_gang_count , host_gang_worker_count , mesh_node_begin , mesh_node_end , run );

  }
  else if ( which == std::string("nonlinear") ) {

    size_t mesh_node_begin = 100 ;
    size_t mesh_node_end   = 300 ;
    size_t run             =   1 ;
    input >> mesh_node_begin >> mesh_node_end >> run ;
    test_host_nonlinear( machine , host_gang_count , host_gang_worker_count , mesh_node_begin , mesh_node_end , run );

  }
  else if ( which == std::string("nonlinear_quadratic") ) {

    size_t mesh_node_begin = 100 ;
    size_t mesh_node_end   = 300 ;
    size_t run             =   1 ;
    input >> mesh_node_begin >> mesh_node_end >> run ;
    test_host_nonlinear_quadratic( machine , host_gang_count , host_gang_worker_count , mesh_node_begin , mesh_node_end , run );

  }
  else {
    cmd_error = true ;
  }

  return cmd_error ;
}

#if defined( KOKKOS_ENABLE_CUDA )
bool run_cuda( std::istream & input , comm::Machine machine )
{
  bool cmd_error = false ;

  std::string which ; input >> which ;

  if ( which == std::string("fixture") ) {

    size_t nx = 0 , ny = 0 , nz = 0 ;
    input >> nx >> ny >> nz ;
    test_cuda_fixture( machine , nx , ny , nz );

  }
  else if ( which == std::string("explicit") ) {

    size_t mesh_node_begin = 100 ;
    size_t mesh_node_end   = 300 ;
    size_t run             =   1 ;
    input >> mesh_node_begin >> mesh_node_end >> run ;
    test_cuda_explicit( machine , mesh_node_begin , mesh_node_end , run );

  }
  else if ( which == std::string("implicit") ) {

    size_t mesh_node_begin = 100 ;
    size_t mesh_node_end   = 300 ;
    size_t run             =   1 ;
    input >> mesh_node_begin >> mesh_node_end >> run ;
    test_cuda_implicit( machine , mesh_node_begin , mesh_node_end , run );

  }
  else if ( which == std::string("nonlinear") ) {

    size_t mesh_node_begin = 100 ;
    size_t mesh_node_end   = 300 ;
    size_t run             =   1 ;
    input >> mesh_node_begin >> mesh_node_end >> run ;
    test_cuda_nonlinear( machine , mesh_node_begin , mesh_node_end , run );

  }
  else if ( which == std::string("nonlinear_quadratic") ) {

    size_t mesh_node_begin = 100 ;
    size_t mesh_node_end   = 300 ;
    size_t run             =   1 ;
    input >> mesh_node_begin >> mesh_node_end >> run ;
    test_cuda_nonlinear_quadratic( machine , mesh_node_begin , mesh_node_end , run );

  }
  else {
    cmd_error = true ;
  }

  return cmd_error ;
}
#endif

void run( const std::string & argline , comm::Machine machine )
{
  const unsigned numa_count       = Kokkos::hwloc::get_available_numa_count();
  const unsigned cores_per_numa   = Kokkos::hwloc::get_available_cores_per_numa();
  const unsigned threads_per_core = Kokkos::hwloc::get_available_threads_per_core();

  std::istringstream input( argline );

  bool cmd_error = false ;

  std::string which ; input >> which ;

  if ( which == std::string("query") ) {
    std::cout << "P" << comm::rank( machine )
              << ": hwloc { NUMA[" << numa_count << "]"
              << " CORE[" << cores_per_numa << "]"
              << " PU[" << threads_per_core << "] }"
              << std::endl ;
#if defined( KOKKOS_ENABLE_CUDA )
    test_cuda_query( machine );
#endif
  }
  else if ( which == std::string("partition") ) {
    if ( 0 == comm::rank( machine ) ) {
      test_box_partition( false /* print flag */ );
    }
  }
  else {
    if ( which == std::string("host") ) {
      size_t host_gang_count = 0 ;
      size_t host_gang_worker_count = 1 ;

      input >> host_gang_count ;
      input >> host_gang_worker_count ;

      cmd_error = run_host( input , machine , host_gang_count , host_gang_worker_count );
    }
    else if ( which == std::string("host-all") ) {
      size_t host_gang_count        = numa_count ;
      size_t host_gang_worker_count = cores_per_numa * threads_per_core ;

      cmd_error = run_host( input , machine , host_gang_count , host_gang_worker_count );
    }
    else if ( which == std::string("host-most") ) {
      size_t host_gang_count        = numa_count ;
      size_t host_gang_worker_count = ( cores_per_numa - 1 ) * threads_per_core ;

      cmd_error = run_host( input , machine , host_gang_count , host_gang_worker_count );
    }
#if defined( KOKKOS_ENABLE_CUDA )
    else if ( which == std::string("cuda") ) {
      cmd_error = run_cuda( input , machine );
    }
#endif
    else {
      cmd_error = true ;
    }
  }

  if ( cmd_error && 0 == comm::rank( machine ) ) {
    std::cout << "Expecting command line with" << std::endl
              << "    query" << std::endl
              << "    partition" << std::endl
              << "    host NumNumaNode NumThreadPerNode <test>" << std::endl
              << "    host-all <test>" << std::endl
              << "    host-most <test>" << std::endl
              << "    cuda <test>" << std::endl
              << "where <test> is" << std::endl
              << "    fixture   NumElemX NumElemY NumElemZ" << std::endl
              << "    implicit  NumElemBegin NumElemEnd NumRun" << std::endl
              << "    explicit  NumElemBegin NumElemEnd NumRun" << std::endl
              << "    nonlinear NumElemBegin NumElemEnd NumRun" << std::endl
              << "    nonlinear_quadratic NumElemBegin NumElemEnd NumRun" << std::endl ;

  }
}

} // namespace

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

int main( int argc , char ** argv )
{
  comm::Machine machine = comm::Machine::init( & argc , & argv );

  const unsigned comm_rank = comm::rank( machine );

  const std::string argline = comm::command_line( machine , argc , argv );

  try {
    run( argline , machine );
  }
  catch( const std::exception & x ) {
    std::cerr << "P" << comm_rank << " throw: " << x.what() << std::endl ;
  }
  catch( ... ) {
    std::cerr << "P" << comm_rank << " throw: unknown exception" << std::endl ;
  }

  comm::Machine::finalize();

  return 0 ;
}

