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

#include <utility>
#include <iostream>

#include <Kokkos_Core.hpp>

#include <feint_fwd.hpp>

int main()
{
#if defined( KOKKOS_ENABLE_THREADS )
  {
    // Use 4 cores per NUMA region, unless fewer available

    const unsigned use_numa_count     = Kokkos::hwloc::get_available_numa_count();
    const unsigned use_cores_per_numa = std::min( 4u , Kokkos::hwloc::get_available_cores_per_numa() );

    Kokkos::Threads::initialize( use_numa_count * use_cores_per_numa );

    std::cout << "feint< Threads , NotUsingAtomic >" << std::endl ;
    Kokkos::Example::feint< Kokkos::Threads , false >();

    std::cout << "feint< Threads , Usingtomic >" << std::endl ;
    Kokkos::Example::feint< Kokkos::Threads , true  >();

    Kokkos::Threads::finalize();
  }
#endif

#if defined( KOKKOS_ENABLE_OPENMP )
  {

    int num_threads  = 0;
    if ( Kokkos::hwloc::available() ) {
      // Use 4 cores per NUMA region, unless fewer available
      const unsigned use_numa_count     = Kokkos::hwloc::get_available_numa_count();
      const unsigned use_cores_per_numa = std::min( 4u , Kokkos::hwloc::get_available_cores_per_numa() );
      num_threads = use_numa_count * use_cores_per_numa;

    }
    else {
      #pragma omp parallel
      {
        #pragma omp atomic
        ++num_threads;
      }
      num_threads = std::max(4, num_threads/4);
    }


    Kokkos::OpenMP::initialize( num_threads );

    std::cout << "feint< OpenMP , NotUsingAtomic >" << std::endl ;
    Kokkos::Example::feint< Kokkos::OpenMP , false >();

    std::cout << "feint< OpenMP , Usingtomic >" << std::endl ;
    Kokkos::Example::feint< Kokkos::OpenMP , true  >();

    Kokkos::OpenMP::finalize();
  }
#endif

#if defined( KOKKOS_ENABLE_CUDA )
  {
    // Initialize Host mirror device
    Kokkos::HostSpace::execution_space::initialize(1);
    const unsigned device_count = Kokkos::Cuda::detect_device_count();

    // Use the last device:
    Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice(device_count-1) );

    std::cout << "feint< Cuda , NotUsingAtomic >" << std::endl ;
    Kokkos::Example::feint< Kokkos::Cuda , false >();

    std::cout << "feint< Cuda , UsingAtomic >" << std::endl ;
    Kokkos::Example::feint< Kokkos::Cuda , true  >();

    Kokkos::Cuda::finalize();
    Kokkos::HostSpace::execution_space::finalize();

  }
#endif
}

