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

#if 0

#include <cstdlib>
#include <cstring>

#include <ParallelMachine.hpp>

#include <Kokkos_Core.hpp>

#if ! defined( KOKKOS_ENABLE_MPI )
#define MPI_COMM_NULL 0
#endif

//------------------------------------------------------------------------

namespace Parallel {

Machine::Machine( int * argc , char *** argv )
  : m_mpi_comm( MPI_COMM_NULL )
  , m_mpi_size(0)
  , m_mpi_rank(0)
  , m_mpi_gpu(0)
{

#if defined( KOKKOS_ENABLE_CUDA )
  //------------------------------------
  // Might be using a Cuda aware version of MPI.
  // Must select Cuda device before initializing MPI.
  {
    int i = 1 ;
    for ( ; i < *argc && strcmp((*argv)[i],"mpi_cuda") ; ++i );

    if ( i < *argc ) {
      // Determine, if possible, what will be the node-local
      // rank of the MPI process once MPI has been initialized.
      // This rank is needed to set the Cuda device before 'mvapich'
      // is initialized.

      const char * const mvapich_local_rank = getenv("MV2_COMM_WORLD_LOCAL_RANK");
      const char * const slurm_local_rank   = getenv("SLURM_LOCALID");

      const int pre_mpi_local_rank =
        0 != mvapich_local_rank ? atoi( mvapich_local_rank ) : (
        0 != slurm_local_rank   ? atoi( slurm_local_rank ) : (
        -1 ) );

      if ( 0 <= pre_mpi_local_rank ) {

        const int ngpu = Kokkos::Cuda::detect_device_count();

        const int cuda_device_rank = pre_mpi_local_rank % ngpu ;

        Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice( cuda_device_rank ) );

        m_mpi_gpu = 1 ;
      }
    }
  }
#endif

  //------------------------------------

#if defined( KOKKOS_ENABLE_MPI )
  MPI_Init( argc , argv );
  m_mpi_comm = MPI_COMM_WORLD ;
  MPI_Comm_size( m_mpi_comm , & m_mpi_size );
  MPI_Comm_rank( m_mpi_comm , & m_mpi_rank );
#endif

  // Query hwloc after MPI initialization to allow MPI binding:
  //------------------------------------
  // Request to use host device:
  {
    int i = 1 ;
    for ( ; i < *argc && strcmp((*argv)[i],"host") ; ++i );

    if ( i < *argc ) {

      unsigned team_count       = Kokkos::hwloc::get_available_numa_count();
      unsigned threads_per_team = Kokkos::hwloc::get_available_cores_per_numa() *
                                  Kokkos::hwloc::get_available_threads_per_core();

      if ( i + 2 < *argc ) {
        team_count       = atoi( (*argv)[i+1] );
        threads_per_team = atoi( (*argv)[i+2] );
      }

      Kokkos::Threads::initialize( team_count * threads_per_team );
    }
  }

#if defined( KOKKOS_ENABLE_CUDA )
  //------------------------------------
  // Request to use Cuda device and not already initialized.
  if ( ! m_mpi_gpu ) {
    int i = 1 ;
    for ( ; i < *argc && strcmp((*argv)[i],"mpi_cuda") && strcmp((*argv)[i],"cuda") ; ++i );

    if ( i < *argc ) {

      const int ngpu = Kokkos::Cuda::detect_device_count();

      const int cuda_device_rank = m_mpi_rank % ngpu ;

      Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice( cuda_device_rank ) );
    }
  }
#endif

}

Machine::~Machine()
{
  Kokkos::Threads::finalize();
#if defined( KOKKOS_ENABLE_CUDA )
  Kokkos::Cuda::finalize();
#endif
#if defined( KOKKOS_ENABLE_MPI )
  MPI_Finalize();
#endif
}

void Machine::print_configuration( std::ostream & msg ) const
{
  msg << "MPI [ " << m_mpi_rank << " / " << m_mpi_size << " ]" << std::endl ;
  Kokkos::Threads::print_configuration( msg );
#if defined( KOKKOS_ENABLE_CUDA )
  Kokkos::Cuda::print_configuration( msg );
#endif
}

}

#endif /* #if 0 */

