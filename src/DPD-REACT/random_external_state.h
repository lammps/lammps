/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/*
   This random_external_state.h file was derrived from the Kokkos
   file algorithms/src/Kokkos_Random.hpp and adapted to work
   without Kokkos installed, as well as being converted to a form
   that has no internal state.  All RNG state information is kept
   outside this "class", and is passed in by reference by the caller.
 */
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

#ifndef LMP_RANDOM_EXTERNALSTATE_H
#define LMP_RANDOM_EXTERNALSTATE_H

#include "accelerator_kokkos.h"
#include <cmath>

/// \file random_external_state.h
/// \brief Pseudorandom number generators
///
/// These generators are based on Vigna, Sebastiano (2014). "An
/// experimental exploration of Marsaglia's xorshift generators,
/// scrambled."  See: http://arxiv.org/abs/1402.6246

// A replacement for the Kokkos Random_XorShift64 class that uses
// an external state variable, instead of a class member variable.
namespace random_external_state {
typedef uint64_t es_RNG_t;

constexpr uint32_t MAX_URAND = 0xffffffffU;
constexpr uint64_t MAX_URAND64 = 0xffffffffffffffffULL - 1;

LAMMPS_INLINE
uint32_t es_urand(es_RNG_t &state_)
{
  state_ ^= state_ >> 12;
  state_ ^= state_ << 25;
  state_ ^= state_ >> 27;

  es_RNG_t tmp = state_ * 2685821657736338717ULL;
  tmp = tmp >> 16;
  return static_cast<uint32_t>(tmp & MAX_URAND);
}

LAMMPS_INLINE
uint64_t es_urand64(es_RNG_t &state_)
{
  state_ ^= state_ >> 12;
  state_ ^= state_ << 25;
  state_ ^= state_ >> 27;
  return (state_ * 2685821657736338717ULL) - 1;
}

LAMMPS_INLINE
int es_rand(es_RNG_t &state_)
{
  return static_cast<int>(es_urand(state_) / 2);
}

LAMMPS_INLINE
double es_drand(es_RNG_t &state_)
{
  return static_cast<double>(es_urand64(state_)) / static_cast<double>(MAX_URAND64);
}

//Marsaglia polar method for drawing a standard normal distributed random number
LAMMPS_INLINE
double es_normal(es_RNG_t &state_)
{
  double S, U;
  do {
    U = 2.0 * es_drand(state_) - 1.0;
    const double V = 2.0 * es_drand(state_) - 1.0;
    S = U * U + V * V;
  } while ((S >= 1.0) || (S == 0.0));
  return U * sqrt(-2.0 * log(S) / S);
}

LAMMPS_INLINE
double es_normalPair(es_RNG_t &state_, double &second)
{
  double S, U, V;
  do {
    U = 2.0 * es_drand(state_) - 1.0;
    V = 2.0 * es_drand(state_) - 1.0;
    S = U * U + V * V;
  } while ((S >= 1.0) || (S == 0.0));
  const double fac = sqrt(-2.0 * log(S) / S);
  second = V * fac;
  return U * fac;
}

// Use es_init() to init a serial RNG, that is then
// used to generate the initial state of your k parallel
// RNGs with k calls to genNextParallelState()
LAMMPS_INLINE
void es_init(es_RNG_t &serial_state, uint64_t seed)
{
  if (seed == 0) seed = uint64_t(1318319);
  serial_state = seed;
  for (int i = 0; i < 17; i++) es_rand(serial_state);
}

// Call genNextParallelState() once for each RNG to generate
// the initial state for that RNG.
LAMMPS_INLINE
void es_genNextParallelState(es_RNG_t &serial_state, es_RNG_t &new_state)
{
  int n1 = es_rand(serial_state);
  int n2 = es_rand(serial_state);
  int n3 = es_rand(serial_state);
  int n4 = es_rand(serial_state);
  new_state = ((((static_cast<es_RNG_t>(n1)) & 0xffff) << 00) |
               (((static_cast<es_RNG_t>(n2)) & 0xffff) << 16) |
               (((static_cast<es_RNG_t>(n3)) & 0xffff) << 32) |
               (((static_cast<es_RNG_t>(n4)) & 0xffff) << 48));
}
}    // namespace random_external_state

#endif
