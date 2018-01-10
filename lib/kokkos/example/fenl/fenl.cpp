/*
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
*/

#include <HexElement.hpp>
#include <fenl_impl.hpp>

namespace Kokkos {
namespace Example {
namespace FENL {

#if defined( KOKKOS_ENABLE_THREADS )

template
Perf fenl< Kokkos::Threads , Kokkos::Example::BoxElemPart::ElemLinear >(
  MPI_Comm comm ,
  const int use_print ,
  const int use_trials ,
  const int use_atomic ,
  const int global_elems[] );


template
Perf fenl< Kokkos::Threads , Kokkos::Example::BoxElemPart::ElemQuadratic >(
  MPI_Comm comm ,
  const int use_print ,
  const int use_trials ,
  const int use_atomic ,
  const int global_elems[] );

#endif


#if defined (KOKKOS_ENABLE_OPENMP)

template
Perf fenl< Kokkos::OpenMP , Kokkos::Example::BoxElemPart::ElemLinear >(
  MPI_Comm comm ,
  const int use_print ,
  const int use_trials ,
  const int use_atomic ,
  const int global_elems[] );


template
Perf fenl< Kokkos::OpenMP , Kokkos::Example::BoxElemPart::ElemQuadratic >(
  MPI_Comm comm ,
  const int use_print ,
  const int use_trials ,
  const int use_atomic ,
  const int global_elems[] );

#endif

#if defined( KOKKOS_ENABLE_CUDA )

template
Perf fenl< Kokkos::Cuda , Kokkos::Example::BoxElemPart::ElemLinear >(
  MPI_Comm comm ,
  const int use_print ,
  const int use_trials ,
  const int use_atomic ,
  const int global_elems[] );


template
Perf fenl< Kokkos::Cuda , Kokkos::Example::BoxElemPart::ElemQuadratic >(
  MPI_Comm comm ,
  const int use_print ,
  const int use_trials ,
  const int use_atomic ,
  const int global_elems[] );

#endif

#if defined( KOKKOS_ENABLE_ROCM )

template
Perf fenl< Kokkos::Experimental::ROCm , Kokkos::Example::BoxElemPart::ElemLinear >(
  MPI_Comm comm ,
  const int use_print ,
  const int use_trials ,
  const int use_atomic ,
  const int global_elems[] );


template
Perf fenl< Kokkos::Experimental::ROCm , Kokkos::Example::BoxElemPart::ElemQuadratic >(
  MPI_Comm comm ,
  const int use_print ,
  const int use_trials ,
  const int use_atomic ,
  const int global_elems[] );

#endif


} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos */

