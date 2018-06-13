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

#include <cmath>

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <utility>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <Kokkos_Core.hpp>

#include <WrapMPI.hpp>
#include <fenl.hpp>

// For vtune
#include <sys/types.h>
#include <unistd.h>

//----------------------------------------------------------------------------

enum { CMD_USE_THREADS = 0
     , CMD_USE_NUMA
     , CMD_USE_CORE_PER_NUMA
     , CMD_USE_CUDA
     , CMD_USE_ROCM
     , CMD_USE_OPENMP
     , CMD_USE_CUDA_DEV
     , CMD_USE_FIXTURE_X
     , CMD_USE_FIXTURE_Y
     , CMD_USE_FIXTURE_Z
     , CMD_USE_FIXTURE_BEGIN
     , CMD_USE_FIXTURE_END
     , CMD_USE_FIXTURE_QUADRATIC
     , CMD_USE_ATOMIC
     , CMD_USE_TRIALS
     , CMD_VTUNE
     , CMD_PRINT
     , CMD_ECHO
     , CMD_ERROR
     , CMD_COUNT };

void print_cmdline( std::ostream & s , const int cmd[] )
{
  if ( cmd[ CMD_USE_THREADS ] ) {
    s << " Threads(" << cmd[ CMD_USE_THREADS ]
      << ") NUMA(" << cmd[ CMD_USE_NUMA ]
      << ") CORE_PER_NUMA(" << cmd[ CMD_USE_CORE_PER_NUMA ]
      << ")" ;
  }
  if ( cmd[ CMD_USE_OPENMP ] ) {
    s << " OpenMP(" << cmd[ CMD_USE_OPENMP ]
      << ") NUMA(" << cmd[ CMD_USE_NUMA ]
      << ") CORE_PER_NUMA(" << cmd[ CMD_USE_CORE_PER_NUMA ]
      << ")" ;
  }
  if ( cmd[ CMD_USE_FIXTURE_X ] ) {
    s << " Fixture(" << cmd[ CMD_USE_FIXTURE_X ]
      << "x" << cmd[ CMD_USE_FIXTURE_Y ]
      << "x" << cmd[ CMD_USE_FIXTURE_Z ]
      << ")" ;
  }
  if ( cmd[ CMD_USE_FIXTURE_BEGIN ] ) {
    s << " Fixture( " << cmd[ CMD_USE_FIXTURE_BEGIN ]
      << " .. " << cmd[ CMD_USE_FIXTURE_END ]
      << " )" ;
  }
  if ( cmd[ CMD_USE_FIXTURE_QUADRATIC ] ) {
    s << " Quadratic-Element" ;
  }
  if ( cmd[ CMD_USE_CUDA ] ) {
    s << " CUDA(" << cmd[ CMD_USE_CUDA_DEV ] << ")" ;
  }
  if ( cmd[ CMD_USE_ROCM ] ) {
    s << " ROCM" ;
  }
  if ( cmd[ CMD_USE_ATOMIC ] ) {
    s << " ATOMIC" ;
  }
  if ( cmd[ CMD_USE_TRIALS ] ) {
    s << " TRIALS(" << cmd[ CMD_USE_TRIALS ] << ")" ;
  }
  if ( cmd[ CMD_VTUNE ] ) {
    s << " VTUNE" ;
  }
  if ( cmd[ CMD_PRINT ] ) {
    s << " PRINT" ;
  }
  s << std::endl ;
}

void print_perf_value( std::ostream & s , const std::vector<size_t> & widths,  const Kokkos::Example::FENL::Perf & perf )
{
  int i=0;
  s << std::setw(widths[i++]) << perf.global_elem_count << " ,";
  s << std::setw(widths[i++]) << perf.global_node_count << " ,";
  s << std::setw(widths[i++]) << perf.newton_iter_count << " ,";
  s << std::setw(widths[i++]) << perf.cg_iter_count << " ,";
  s << std::setw(widths[i++]) << perf.map_ratio << " ,";
  s << std::setw(widths[i++]) << ( perf.fill_node_set * 1000.0 ) / perf.global_node_count << " ,";
  s << std::setw(widths[i++]) << ( perf.scan_node_count * 1000.0 ) / perf.global_node_count << " ,";
  s << std::setw(widths[i++]) << ( perf.fill_graph_entries * 1000.0 ) / perf.global_node_count << " ,";
  s << std::setw(widths[i++]) << ( perf.sort_graph_entries * 1000.0 ) / perf.global_node_count << " ,";
  s << std::setw(widths[i++]) << ( perf.fill_element_graph * 1000.0 ) / perf.global_node_count << " ,";
  s << std::setw(widths[i++]) << ( perf.create_sparse_matrix * 1000.0 ) / perf.global_node_count << " ,";
  s << std::setw(widths[i++]) << ( perf.fill_time * 1000.0 ) / perf.global_node_count << " ,";
  s << std::setw(widths[i++]) << ( perf.bc_time * 1000.0 ) / perf.global_node_count << " ,";
  s << std::setw(widths[i++]) << ( ( perf.matvec_time * 1000.0 ) / perf.cg_iter_count ) / perf.global_node_count << " ,";
  s << std::setw(widths[i++]) << ( ( perf.cg_time * 1000.0 ) / perf.cg_iter_count ) / perf.global_node_count << " ,";
  s << std::setw(widths[i])   << perf.error_max;
  s << std::endl ;
}

template< class Device , Kokkos::Example::BoxElemPart::ElemOrder ElemOrder >
void run( MPI_Comm comm , const int cmd[] )
{
  int comm_rank = 0 ;

#if defined( KOKKOS_ENABLE_MPI )
  MPI_Comm_rank( comm , & comm_rank );
#else
  comm = 0 ;
#endif


  if ( 0 == comm_rank ) {
    if ( cmd[ CMD_USE_THREADS ] ) { std::cout << "THREADS , " << cmd[ CMD_USE_THREADS ] ; }
    else if ( cmd[ CMD_USE_OPENMP ] ) { std::cout << "OPENMP , " << cmd[ CMD_USE_OPENMP ] ; }
    else if ( cmd[ CMD_USE_CUDA ] ) { std::cout << "CUDA" ; }
    else if ( cmd[ CMD_USE_ROCM ] ) { std::cout << "ROCM" ; }

    if ( cmd[ CMD_USE_FIXTURE_QUADRATIC ] ) { std::cout << " , QUADRATIC-ELEMENT" ; }
    else { std::cout << " , LINEAR-ELEMENT" ; }

    if ( cmd[ CMD_USE_ATOMIC ] ) { std::cout << " , USING ATOMICS" ; }
  }

  std::vector< std::pair<std::string,std::string> > headers;


  headers.push_back(std::make_pair("ELEMS","count"));
  headers.push_back(std::make_pair("NODES","count"));
  headers.push_back(std::make_pair("NEWTON","iter"));
  headers.push_back(std::make_pair("CG","iter"));
  headers.push_back(std::make_pair("MAP_RATIO","ratio"));
  headers.push_back(std::make_pair("SET_FILL/NODE","millisec"));
  headers.push_back(std::make_pair("SCAN/NODE","millisec"));
  headers.push_back(std::make_pair("GRAPH_FILL/NODE","millisec"));
  headers.push_back(std::make_pair("SORT/NODE","millisec"));
  headers.push_back(std::make_pair("ELEM_GRAPH_FILL/NODE","millisec"));
  headers.push_back(std::make_pair("MATRIX_CREATE/NODE","millisec"));
  headers.push_back(std::make_pair("MATRIX_FILL/NODE","millisec"));
  headers.push_back(std::make_pair("BOUNDARY/NODE","millisec"));
  headers.push_back(std::make_pair("MAT_VEC/ITER/ROW","millisec"));
  headers.push_back(std::make_pair("CG/ITER/ROW","millisec"));
  headers.push_back(std::make_pair("ERROR","ratio"));

  // find print widths
  size_t min_width = 10;
  std::vector< size_t > widths(headers.size());
  for (size_t i=0, ie=headers.size(); i<ie; ++i)
    widths[i] = std::max(min_width, headers[i].first.size()+1);

  // print column headers
  if ( 0 == comm_rank ) {
    std::cout << std::endl ;
    for (size_t i=0; i<headers.size(); ++i)
      std::cout << std::setw(widths[i]) << headers[i].first << " ,";
    std::cout << "\b\b  " << std::endl;
    for (size_t i=0; i<headers.size(); ++i)
      std::cout << std::setw(widths[i]) << headers[i].second << " ,";
    std::cout << "\b\b  " << std::endl;

    std::cout << std::scientific;
    std::cout.precision(3);
  }

  if ( cmd[ CMD_USE_FIXTURE_BEGIN ] ) {
    for ( int i = cmd[CMD_USE_FIXTURE_BEGIN] ; i < cmd[CMD_USE_FIXTURE_END] * 2 ; i *= 2 ) {
      int nelem[3] ;
      nelem[0] = std::max( 1 , (int) cbrt( ((double) i) / 2.0 ) );
      nelem[1] = 1 + nelem[0] ;
      nelem[2] = 2 * nelem[0] ;

      const Kokkos::Example::FENL::Perf perf =
        cmd[ CMD_USE_FIXTURE_QUADRATIC ]
        ? Kokkos::Example::FENL::fenl< Device , Kokkos::Example::BoxElemPart::ElemQuadratic >
            ( comm , cmd[CMD_PRINT], cmd[CMD_USE_TRIALS], cmd[CMD_USE_ATOMIC], nelem )
        : Kokkos::Example::FENL::fenl< Device , Kokkos::Example::BoxElemPart::ElemLinear >
            ( comm , cmd[CMD_PRINT], cmd[CMD_USE_TRIALS], cmd[CMD_USE_ATOMIC], nelem )
        ;

      if ( 0 == comm_rank ) print_perf_value( std::cout , widths, perf );
    }
  }
  else {
    int nelem[3] = { cmd[ CMD_USE_FIXTURE_X ] ,
                     cmd[ CMD_USE_FIXTURE_Y ] ,
                     cmd[ CMD_USE_FIXTURE_Z ] };

    const Kokkos::Example::FENL::Perf perf =
      cmd[ CMD_USE_FIXTURE_QUADRATIC ]
      ? Kokkos::Example::FENL::fenl< Device , Kokkos::Example::BoxElemPart::ElemQuadratic >
          ( comm , cmd[CMD_PRINT], cmd[CMD_USE_TRIALS], cmd[CMD_USE_ATOMIC], nelem )
      : Kokkos::Example::FENL::fenl< Device , Kokkos::Example::BoxElemPart::ElemLinear >
          ( comm , cmd[CMD_PRINT], cmd[CMD_USE_TRIALS], cmd[CMD_USE_ATOMIC], nelem )
      ;

    if ( 0 == comm_rank ) print_perf_value( std::cout , widths, perf );
  }
}

//----------------------------------------------------------------------------

int main( int argc , char ** argv )
{
  int comm_rank = 0 ;

#if defined( KOKKOS_ENABLE_MPI )
  MPI_Init( & argc , & argv );
  MPI_Comm comm = MPI_COMM_WORLD ;
  MPI_Comm_rank( comm , & comm_rank );
#else
  MPI_Comm comm = 0 ;
  (void) comm ; // suppress warning
#endif

  Kokkos::initialize(argc,argv);
  int cmdline[ CMD_COUNT ] ;

  for ( int i = 0 ; i < CMD_COUNT ; ++i ) cmdline[i] = 0 ;

  if ( 0 == comm_rank ) {
    for ( int i = 1 ; i < argc ; ++i ) {
      if ( 0 == strcasecmp( argv[i] , "fixture" ) ) {
        sscanf( argv[++i] , "%dx%dx%d" ,
                cmdline + CMD_USE_FIXTURE_X ,
                cmdline + CMD_USE_FIXTURE_Y ,
                cmdline + CMD_USE_FIXTURE_Z );
      }
      else if ( 0 == strcasecmp( argv[i] , "fixture-range" ) ) {
        sscanf( argv[++i] , "%d..%d" ,
                cmdline + CMD_USE_FIXTURE_BEGIN ,
                cmdline + CMD_USE_FIXTURE_END );
      }
      else if ( 0 == strcasecmp( argv[i] , "fixture-quadratic" ) ) {
        cmdline[ CMD_USE_FIXTURE_QUADRATIC ] = 1 ;
      }
      else if ( 0 == strcasecmp( argv[i] , "atomic" ) ) {
        cmdline[ CMD_USE_ATOMIC ] = 1 ;
      }
      else if ( 0 == strcasecmp( argv[i] , "trials" ) ) {
        cmdline[ CMD_USE_TRIALS ] = atoi( argv[++i] ) ;
      }
      else if ( 0 == strcasecmp( argv[i] , "vtune" ) ) {
        cmdline[ CMD_VTUNE ] = 1 ;
      }
      else if ( 0 == strcasecmp( argv[i] , "print" ) ) {
        cmdline[ CMD_PRINT ] = 1 ;
      }
      else if ( 0 == strcasecmp( argv[i] , "echo" ) ) {
        cmdline[ CMD_ECHO ] = 1 ;
      }
      else {
        cmdline[ CMD_ERROR ] = 1 ;

        std::cerr << "Unrecognized command line argument #" << i << ": " << argv[i] << std::endl ;
      }
    }

    if ( cmdline[ CMD_ECHO ] && 0 == comm_rank ) { print_cmdline( std::cout , cmdline ); }
  }

#if defined( KOKKOS_ENABLE_MPI )
  MPI_Bcast( cmdline , CMD_COUNT , MPI_INT , 0 , comm );
#endif

  if ( cmdline[ CMD_VTUNE ] ) {
    std::stringstream cmd;
    pid_t my_os_pid=getpid();
    const std::string vtune_loc =
      "/usr/local/intel/vtune_amplifier_xe_2013/bin64/amplxe-cl";
    const std::string output_dir = "./vtune/vtune.";
    const int p_rank = comm_rank;
    cmd << vtune_loc
        << " -collect hotspots -result-dir " << output_dir << p_rank
        << " -target-pid " << my_os_pid << " &";
    if (p_rank == 0)
      std::cout << cmd.str() << std::endl;
    system(cmd.str().c_str());
    system("sleep 10");
  }

  if ( ! cmdline[ CMD_ERROR ] && ! cmdline[ CMD_ECHO ] ) {

    if ( ! cmdline[ CMD_USE_TRIALS ] ) { cmdline[ CMD_USE_TRIALS ] = 1 ; }

    if ( ! cmdline[ CMD_USE_FIXTURE_X ] && ! cmdline[ CMD_USE_FIXTURE_BEGIN ] ) {
      cmdline[ CMD_USE_FIXTURE_X ] = 2 ;
      cmdline[ CMD_USE_FIXTURE_Y ] = 2 ;
      cmdline[ CMD_USE_FIXTURE_Z ] = 2 ;
    }

    run< Kokkos::DefaultExecutionSpace , Kokkos::Example::BoxElemPart::ElemLinear >( comm , cmdline );

  }

  Kokkos::finalize();
#if defined( KOKKOS_ENABLE_MPI )
  MPI_Finalize();
#endif

  return cmdline[ CMD_ERROR ] ? -1 : 0 ;
}

