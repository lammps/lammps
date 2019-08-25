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

#ifndef KOKKOS_EXAMPLE_WRAP_MPI
#define KOKKOS_EXAMPLE_WRAP_MPI

#include <Kokkos_Macros.hpp>
#include <string>

#if defined( KOKKOS_ENABLE_MPI )

#include <mpi.h>

namespace Kokkos {
namespace Example {

inline
double all_reduce( double value , MPI_Comm comm )
{
  double local = value ;
  MPI_Allreduce( & local , & value , 1 , MPI_DOUBLE , MPI_SUM , comm );
  return value ;
}

inline
double all_reduce_max( double value , MPI_Comm comm )
{
  double local = value ;
  MPI_Allreduce( & local , & value , 1 , MPI_DOUBLE , MPI_MAX , comm );
  return value ;
}

} // namespace Example
} // namespace Kokkos

#elif ! defined( KOKKOS_ENABLE_MPI )

/* Wrap the the MPI_Comm type and heavily used MPI functions
 * to reduce the number of '#if defined( KOKKOS_ENABLE_MPI )'
 * blocks which have to be sprinkled throughout the examples.
 */

typedef int MPI_Comm ;

inline int MPI_Comm_size( MPI_Comm , int * size ) { *size = 1 ; return 0 ; }
inline int MPI_Comm_rank( MPI_Comm , int * rank ) { *rank = 0 ; return 0 ; }
inline int MPI_Barrier( MPI_Comm ) { return 0; }

namespace Kokkos {
namespace Example {

inline
double all_reduce( double value , MPI_Comm ) { return value ; }

inline
double all_reduce_max( double value , MPI_Comm ) { return value ; }

} // namespace Example
} // namespace Kokkos

#endif /* ! defined( KOKKOS_ENABLE_MPI ) */
#endif /* #ifndef KOKKOS_EXAMPLE_WRAP_MPI */

