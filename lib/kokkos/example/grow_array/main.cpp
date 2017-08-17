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

#include <iostream>
#include <sstream>

#include <Kokkos_Core.hpp>

#include <grow_array.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

int main( int argc , char ** argv )
{
  int num_threads = 4 ;
  int use_numa = 1 ;
  int use_core = 1 ;
  int length_array  = 1000000 ;
  int span_values = 100000000 ;


  if ( Kokkos::hwloc::available() ) {
    use_numa = Kokkos::hwloc::get_available_numa_count();
    use_core = Kokkos::hwloc::get_available_cores_per_numa() - 1 ;
    num_threads = use_numa * use_core * Kokkos::hwloc::get_available_threads_per_core();
  }

#if defined( KOKKOS_ENABLE_SERIAL )
  {
    std::cout << "Kokkos::Serial" << std::endl ;
    // The Serial device accepts these arguments, though it may ignore them.
    Kokkos::Serial::initialize( num_threads , use_numa , use_core );
    Example::grow_array< Kokkos::Serial >( length_array , span_values );
    Kokkos::Serial::finalize ();
  }
#endif // defined( KOKKOS_ENABLE_SERIAL )

#if defined( KOKKOS_ENABLE_THREADS )
  {
    std::cout << "Kokkos::Threads" << std::endl ;
    Kokkos::Threads::initialize( num_threads , use_numa , use_core );
    Example::grow_array< Kokkos::Threads >( length_array , span_values );
    Kokkos::Threads::finalize();
  }
#endif

#if defined( KOKKOS_ENABLE_OPENMP )
  {
    std::cout << "Kokkos::OpenMP" << std::endl ;
    Kokkos::OpenMP::initialize();
    Example::grow_array< Kokkos::OpenMP >( length_array , span_values );
    Kokkos::OpenMP::finalize();
  }
#endif

#if defined( KOKKOS_ENABLE_CUDA )
  {
    std::cout << "Kokkos::Cuda" << std::endl ;
    Kokkos::HostSpace::execution_space::initialize(1);
    Kokkos::Cuda::initialize();
    Example::grow_array< Kokkos::Cuda >( length_array , span_values );
    Kokkos::Cuda::finalize();
    Kokkos::HostSpace::execution_space::finalize();
  }
#endif

  return 0 ;
}

