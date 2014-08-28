/*
//@HEADER
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

/// \file Kokkos_Vectorization.hpp
/// \brief Declaration and definition of Kokkos::Vectorization interface.
#ifndef KOKKOS_VECTORIZATION_HPP
#define KOKKOS_VECTORIZATION_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_ExecPolicy.hpp>

namespace Kokkos {
template<class Device, int N>
struct Vectorization {
  typedef Kokkos::TeamPolicy< Device >       team_policy ;
  typedef typename team_policy::member_type  team_member ;

  enum {increment = 1};

  KOKKOS_FORCEINLINE_FUNCTION
  static int begin() { return 0;}

  KOKKOS_FORCEINLINE_FUNCTION
  static int thread_rank(const team_member &dev) {
    return dev.team_rank();
  }

  KOKKOS_FORCEINLINE_FUNCTION
  static int team_rank(const team_member &dev) {
    return dev.team_rank()/increment;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  static int team_size(const team_member &dev) {
    return dev.team_size()/increment;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  static int global_thread_rank(const team_member &dev) {
    return (dev.league_rank()*dev.team_size()+dev.team_rank());
  }

  KOKKOS_FORCEINLINE_FUNCTION
  static bool is_lane_0(const team_member &dev) {
    return true;
  }

  template<class Scalar>
  KOKKOS_FORCEINLINE_FUNCTION
  static Scalar reduce(const Scalar& val) {
    return val;
  }
};
}

#if defined( KOKKOS_HAVE_CUDA )
#include <Cuda/Kokkos_Cuda_Vectorization.hpp>
#endif

#endif
