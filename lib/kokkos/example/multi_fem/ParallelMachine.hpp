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

#error "ParallelMachine"

#ifndef PARALLELMACHINE_HPP
#define PARALLELMACHINE_HPP

//------------------------------------------------------------------------

#include <iosfwd>

#include <Kokkos_Core.hpp>

//------------------------------------------------------------------------

#if defined( KOKKOS_ENABLE_MPI )
#include <mpi.h>
#else
  typedef int MPI_Comm ;
#endif

//------------------------------------------------------------------------
//------------------------------------------------------------------------

namespace Parallel {

/** \brief  Hybrid parallel machine with MPI+Kokkos::Threads or MPI+Kokkos::Cuda.
 *
 *  Initialization of MPI and Kokkos device has interdependencies which this
 *  class manages.  The command line and environment variables are queried to initialize
 *  the Threads or Cuda device:
 *
 *    1)  cuda               : initializes Cuda device
 *    2)  host               : initializes Threads device with all hwloc detected cores.
 *    3)  host #gang #worker : initializes Threads with specified
 */
class Machine {
private:

  MPI_Comm m_mpi_comm ;
  int      m_mpi_size ;
  int      m_mpi_rank ;
  unsigned m_mpi_gpu ;
  unsigned m_gpu_arch ;

  Machine();
  Machine( const Machine & );
  Machine & operator = ( const Machine & );

public:

  /** \brief  Coordinated initialize MPI, Cuda, or Threads devices from 'main'.  */
  Machine( int * argc , char *** argv );

  ~Machine();

  MPI_Comm mpi_comm() const { return m_mpi_comm ; }

  int mpi_size() const { return m_mpi_size ; }
  int mpi_rank() const { return m_mpi_rank ; }

  /** \brief  If using MPI that can directly operate on GPU memory */
  bool mpi_gpu() const { return m_mpi_gpu ; }

  /** \brief  If using GPU then what architecture */
  unsigned gpu_arch() const { return m_gpu_arch ; }

  void print_configuration( std::ostream & ) const ;
};

}

//------------------------------------------------------------------------

#endif /* #ifndef PARALLELMACHINE_HPP */


