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

#include <cstring>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include <Kokkos_Core.hpp>

#include <sort_array.hpp>


int main( int argc , char ** argv )
{
#if defined( KOKKOS_ENABLE_CUDA ) || defined( KOKKOS_ENABLE_THREADS ) || defined( KOKKOS_ENABLE_OPENMP )
  Kokkos::initialize( argc , argv );

  int length_array = 100000 ;

  for ( int i = 0 ; i < argc ; ++i ) {
    if ( 0 == strcmp( argv[i] , "length_array" ) ) {
      length_array = atoi( argv[i+1] );
    }
  }

  int length_total_array  = length_array * 100;

#if defined( KOKKOS_ENABLE_CUDA )
  if ( Kokkos::Cuda::is_initialized() ) {
    std::cout << "Kokkos::Cuda" << std::endl ;
    Example::sort_array< Kokkos::Cuda >( length_array , length_total_array );
  }
#endif

#if defined( KOKKOS_ENABLE_THREADS )
  if ( Kokkos::Threads::is_initialized() ) {
    std::cout << "Kokkos::Threads" << std::endl ;
    Example::sort_array< Kokkos::Threads >( length_array , length_total_array );
  }
#endif

#if defined( KOKKOS_ENABLE_OPENMP )
  if ( Kokkos::OpenMP::is_initialized() ) {
    std::cout << "Kokkos::OpenMP" << std::endl ;
    Example::sort_array< Kokkos::OpenMP >( length_array , length_total_array );
  }
#endif

  Kokkos::finalize();
#endif

  return 0 ;
}

