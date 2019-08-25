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

#ifndef PARALLELCOMM_HPP
#define PARALLELCOMM_HPP

//------------------------------------------------------------------------

#include <Kokkos_Macros.hpp>

//------------------------------------------------------------------------

#if defined( KOKKOS_ENABLE_MPI )

#include <mpi.h>
#include <string>

namespace comm {

struct Machine {
  MPI_Comm mpi_comm ;

  Machine() : mpi_comm( MPI_COMM_NULL ) {}

  Machine( const Machine & rhs )
    : mpi_comm( rhs.mpi_comm ) {}

  Machine( MPI_Comm c ) : mpi_comm( c ) {}

  static Machine init( int * argc , char *** argv )
  {
    MPI_Init( argc , argv );
    return Machine( MPI_COMM_WORLD );
  }

  static void finalize() { MPI_Finalize(); }
};

inline
unsigned  size( Machine machine )
{
  int np ; MPI_Comm_size( machine.mpi_comm , & np ); return np ;
}

inline
unsigned  rank( Machine machine )
{
  int ip ; MPI_Comm_rank( machine.mpi_comm , & ip ); return ip ;
}

inline
double max( Machine machine , double local )
{
  double global = 0;
  MPI_Allreduce( & local , & global , 1 , MPI_DOUBLE , MPI_MAX , machine.mpi_comm );
  return global ;
}

inline
std::string command_line( Machine machine , const int argc , const char * const * const argv )
{
  std::string argline ;

  if ( 0 == rank( machine ) ) {
    for ( int i = 1 ; i < argc ; ++i ) {
      argline.append(" ").append( argv[i] );
    }
  }

  int length = argline.length();
  MPI_Bcast( & length , 1 , MPI_INT , 0 , machine.mpi_comm );
  argline.resize( length , ' ' );
  MPI_Bcast( (void*) argline.data() , length , MPI_CHAR , 0 , machine.mpi_comm );

  return argline ;
}

}

#else /* ! defined( KOKKOS_ENABLE_MPI ) */

#include <string>

namespace comm {

// Stub for non-parallel

struct Machine {
  static Machine init( int * , char *** )
  { return Machine(); }

  static void finalize() {}
};

inline
unsigned  size( Machine ) { return 1 ; }

inline
unsigned  rank( Machine ) { return 0 ; }

inline
double max( Machine , double local )
{ return local ; }

inline
std::string command_line( Machine machine , const int argc , const char * const * const argv )
{
  std::string argline ;

  if ( 0 == rank( machine ) ) {
    for ( int i = 1 ; i < argc ; ++i ) {
      argline.append(" ").append( argv[i] );
    }
  }

  return argline ;
}

}

#endif /* ! defined( KOKKOS_ENABLE_MPI ) */

//------------------------------------------------------------------------

#endif /* #ifndef PARALLELCOMM_HPP */


